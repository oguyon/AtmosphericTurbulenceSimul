// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cuda_psf.cu
 * @brief   CUDA/cuFFT open-loop and closed-loop AO science PSF simulation engine
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <cuda_runtime.h>
#include <cufft.h>

#include "atmturb_cuda.h"
#include "atmturb_cuda_psf.h"

#define ATMTURB_AO_CUDA_WF_CHUNK 256

typedef struct
{
    const float  *d_wf;
    const float  *d_amp;
    const float  *d_mask;
    cufftComplex *d_in_grid;
    int           pup_size;
    int           n_fft;
    int           off_x;
    int           off_y;
    int           count;
    int           wf_bcast;
    int           amp_bcast;
    float         phase_scale;
} atmturb_ao_open_embed_kargs_t;

typedef struct
{
    const cufftComplex *d_out_grid;
    float              *d_psf;
    float              *d_peaks;
    int                 n_fft;
    int                 n_psf;
    int                 count;
    int                 do_accum;
    int                 write_psf;
    int                 idx_peak0;
    int                 slot0;
    int                 idx_peak1;
    int                 slot1;
    float               norm;
} atmturb_ao_accum_kargs_t;

typedef struct
{
    const float  *d_wf_slice;
    const float  *d_swf_slice;
    const float  *d_amp_slice;
    const float  *d_mask;
    const float  *d_read_slot;
    float        *d_write_slot;
    float        *d_dm;
    float        *d_pint;
    float        *d_perr;
    cufftComplex *d_in_grid;
    int           pup_size;
    int           n_fft;
    int           off_x;
    int           off_y;
    int           loop_mode;
    int           delay;
    float         opd_scale;
    float         k_sci;
    float         keep;
    float         gain;
    float         kp;
    float         ki;
    float         kd_dt;
    float         dt;
} atmturb_ao_ctrl_kargs_t;

static cudaStream_t s_psf_stream       = 0;
static cufftHandle  s_psf_plan_1       = 0;
static int          s_psf_plan_1_nfft  = 0;
static cufftHandle  s_psf_plan_b       = 0;
static int          s_psf_plan_b_nfft  = 0;
static int          s_psf_plan_b_batch = 0;

/**
 * atmturb_cuda_atomic_max_pos - Reduce non-negative IEEE-754 float across warp into global max
 */
__device__ __forceinline__ void atmturb_cuda_atomic_max_pos(
    float *dst,
    float  val)
{
    for (int offset = 16; offset > 0; offset >>= 1)
    {
        float other = __shfl_down_sync(0xffffffffU, val, offset);
        if (other > val)
        {
            val = other;
        }
    }
    if ((threadIdx.x & 31) == 0 && val > 0.0f)
    {
        atomicMax((int *) dst, __float_as_int(val));
    }
}

/**
 * atmturb_cuda_ao_open_embed_kernel - Embed flat or open-loop masked pupil phasors into FFT grid
 */
__global__ static void atmturb_cuda_ao_open_embed_kernel(
    atmturb_ao_open_embed_kargs_t a)
{
    int x = (int) (blockIdx.x * blockDim.x + threadIdx.x);
    int y = (int) (blockIdx.y * blockDim.y + threadIdx.y);
    int b = (int) blockIdx.z;
    if (x >= a.pup_size || y >= a.pup_size || b >= a.count)
    {
        return;
    }

    size_t idx  = (size_t) y * (size_t) a.pup_size + (size_t) x;
    size_t ntot = (size_t) a.pup_size * (size_t) a.pup_size;
    size_t nf2  = (size_t) a.n_fft * (size_t) a.n_fft;
    float  m    = a.d_mask[idx];
    if (a.d_amp != NULL)
    {
        size_t a_off = a.amp_bcast ? idx : ((size_t) b * ntot + idx);
        m *= a.d_amp[a_off];
    }
    float re = m;
    float im = 0.0f;

    if (a.d_wf != NULL)
    {
        size_t p_off = a.wf_bcast ? idx : ((size_t) b * ntot + idx);
        float  ph    = a.d_wf[p_off] * a.phase_scale;
        sincosf(ph, &im, &re);
        re *= m;
        im *= m;
    }

    size_t dst_idx = (size_t) b * nf2 +
                     (size_t) (a.off_y + y) * (size_t) a.n_fft + (size_t) (a.off_x + x);
    a.d_in_grid[dst_idx] = make_cuFloatComplex(re, im);
}

/**
 * atmturb_cuda_ao_ctrl_embed_kernel - Fused closed-loop DM update and science phasor embedding
 */
__global__ static void atmturb_cuda_ao_ctrl_embed_kernel(
    atmturb_ao_ctrl_kargs_t a)
{
    int x = (int) (blockIdx.x * blockDim.x + threadIdx.x);
    int y = (int) (blockIdx.y * blockDim.y + threadIdx.y);
    if (x >= a.pup_size || y >= a.pup_size)
    {
        return;
    }

    size_t idx     = (size_t) y * (size_t) a.pup_size + (size_t) x;
    float  m       = a.d_mask[idx];
    float  wf_val  = (a.d_wf_slice != NULL) ? a.d_wf_slice[idx] : 0.0f;
    float  dm_cur  = a.d_dm[idx];
    float  res     = wf_val * a.opd_scale - dm_cur;
    float  cur_err = res * m;
    float  e       = (a.delay == 0) ? cur_err : a.d_read_slot[idx];

    if (a.delay > 0)
    {
        a.d_write_slot[idx] = cur_err;
    }

    if (a.loop_mode == 1)
    {
        a.d_dm[idx] = a.keep * dm_cur + a.gain * e;
    }
    else if (a.loop_mode == 2)
    {
        float pint    = a.d_pint[idx];
        float perr    = a.d_perr[idx];
        float der     = (e - perr) * a.kd_dt;
        float pi      = a.keep * pint + e * a.dt;
        a.d_pint[idx] = pi;
        a.d_dm[idx]   = a.kp * e + a.ki * pi + der;
        a.d_perr[idx] = e;
    }

    float amp    = (a.d_amp_slice != NULL) ? (m * a.d_amp_slice[idx]) : m;
    float ph_sci = (a.d_swf_slice != NULL) ? (a.d_swf_slice[idx] - dm_cur * a.k_sci)
                                           : (res * a.k_sci);
    float re     = 0.0f;
    float im     = 0.0f;
    sincosf(ph_sci, &im, &re);
    size_t dst_idx       = (size_t) (a.off_y + y) * (size_t) a.n_fft + (size_t) (a.off_x + x);
    a.d_in_grid[dst_idx] = make_cuFloatComplex(re * amp, im * amp);
}

/**
 * atmturb_cuda_ao_accum_kernel - Extract FFT-shifted Nyquist PSF intensity and endpoint peaks
 */
__global__ static void atmturb_cuda_ao_accum_kernel(
    atmturb_ao_accum_kargs_t a)
{
    int  lane_idx = (int) (blockIdx.x * blockDim.x + threadIdx.x);
    int  vj       = (int) (blockIdx.y * blockDim.y + threadIdx.y);
    bool valid    = (lane_idx < a.n_psf && vj < a.n_psf);

    int    half_fft = a.n_fft >> 1;
    int    u_start  = (a.n_fft - a.n_psf) >> 1;
    int    v_start  = u_start;
    int    u        = u_start + lane_idx;
    int    v        = v_start + vj;
    int    fx       = (u < half_fft) ? (u + half_fft) : (u - half_fft);
    int    fy       = (v < half_fft) ? (v + half_fft) : (v - half_fft);
    size_t nf2      = (size_t) a.n_fft * (size_t) a.n_fft;
    size_t fft_idx  = (size_t) fy * (size_t) a.n_fft + (size_t) fx;
    float  sum      = 0.0f;
    float  last_val = 0.0f;

    for (int b = 0; b < a.count; b++)
    {
        float val = 0.0f;
        if (valid)
        {
            cufftComplex z = a.d_out_grid[(size_t) b * nf2 + fft_idx];
            val            = (z.x * z.x + z.y * z.y) * a.norm;
            sum           += val;
            last_val       = val;
        }
        if (b == a.idx_peak0)
        {
            atmturb_cuda_atomic_max_pos(&a.d_peaks[a.slot0], val);
        }
        if (b == a.idx_peak1)
        {
            atmturb_cuda_atomic_max_pos(&a.d_peaks[a.slot1], val);
        }
    }

    if (valid && a.write_psf)
    {
        size_t psf_idx = (size_t) vj * (size_t) a.n_psf + (size_t) lane_idx;
        if (a.do_accum)
        {
            a.d_psf[psf_idx] += sum;
        }
        else
        {
            a.d_psf[psf_idx] = last_val;
        }
    }
}

/**
 * atmturb_cuda_psf_cleanup - Release cached cuFFT plans and stream for GPU AO PSF engine
 */
void atmturb_cuda_psf_cleanup(void)
{
    if (s_psf_plan_1 != 0)
    {
        cufftDestroy(s_psf_plan_1);
        s_psf_plan_1      = 0;
        s_psf_plan_1_nfft = 0;
    }
    if (s_psf_plan_b != 0)
    {
        cufftDestroy(s_psf_plan_b);
        s_psf_plan_b       = 0;
        s_psf_plan_b_nfft  = 0;
        s_psf_plan_b_batch = 0;
    }
    if (s_psf_stream != 0)
    {
        cudaStreamSynchronize(s_psf_stream);
        cudaStreamDestroy(s_psf_stream);
        s_psf_stream = 0;
    }
}

/**
 * atmturb_cuda_psf_ensure_plans - Create or reuse persistent stream and cuFFT 2D plans
 */
static int atmturb_cuda_psf_ensure_plans(
    int n_fft,
    int batch)
{
    if (s_psf_stream == 0 && cudaStreamCreate(&s_psf_stream) != cudaSuccess)
    {
        return -1;
    }

    if (s_psf_plan_1 == 0 || s_psf_plan_1_nfft != n_fft)
    {
        if (s_psf_plan_1 != 0)
        {
            cufftDestroy(s_psf_plan_1);
            s_psf_plan_1 = 0;
        }
        if (cufftPlan2d(&s_psf_plan_1, n_fft, n_fft, CUFFT_C2C) != CUFFT_SUCCESS ||
            cufftSetStream(s_psf_plan_1, s_psf_stream) != CUFFT_SUCCESS)
        {
            return -1;
        }
        s_psf_plan_1_nfft = n_fft;
    }

    if (batch > 1 &&
        (s_psf_plan_b == 0 || s_psf_plan_b_nfft != n_fft || s_psf_plan_b_batch != batch))
    {
        if (s_psf_plan_b != 0)
        {
            cufftDestroy(s_psf_plan_b);
            s_psf_plan_b = 0;
        }
        int dims[2] = {n_fft, n_fft};
        int n2      = n_fft * n_fft;
        if (cufftPlanMany(&s_psf_plan_b, 2, dims, NULL, 1, n2, NULL, 1, n2,
                          CUFFT_C2C, batch) != CUFFT_SUCCESS ||
            cufftSetStream(s_psf_plan_b, s_psf_stream) != CUFFT_SUCCESS)
        {
            return -1;
        }
        s_psf_plan_b_nfft  = n_fft;
        s_psf_plan_b_batch = batch;
    }

    return 0;
}

/**
 * atmturb_cuda_ao_alloc_opt - Allocate optional chunked or broadcast 2D device buffer
 */
static int atmturb_cuda_ao_alloc_opt(
    const float *h_src,
    long         stride,
    size_t       chunk_len,
    size_t       ntot,
    float      **d_out)
{
    *d_out = NULL;
    if (h_src == NULL)
    {
        return 0;
    }
    size_t n_frames = (stride == 0) ? 1U : chunk_len;
    if (cudaMalloc((void **) d_out, sizeof(float) * n_frames * ntot) != cudaSuccess)
    {
        return -1;
    }
    if (stride == 0 &&
        cudaMemcpyAsync(*d_out, h_src, sizeof(float) * ntot,
                        cudaMemcpyHostToDevice, s_psf_stream) != cudaSuccess)
    {
        cudaFree(*d_out);
        *d_out = NULL;
        return -1;
    }
    return 0;
}

/**
 * atmturb_cuda_ao_upload_chunk - Copy frame chunk when optional buffer is 3D (`stride > 0`)
 */
static int atmturb_cuda_ao_upload_chunk(
    float       *d_dst,
    const float *h_src,
    long         stride,
    long         t0,
    long         clen,
    size_t       ntot)
{
    if (d_dst == NULL || h_src == NULL || stride == 0)
    {
        return 0;
    }
    return (cudaMemcpyAsync(d_dst, h_src + (size_t) t0 * ntot,
                            sizeof(float) * (size_t) clen * ntot,
                            cudaMemcpyHostToDevice, s_psf_stream) == cudaSuccess) ? 0 : -1;
}

/**
 * atmturb_cuda_ao_run_open - Execute batched open-loop AO PSF accumulation on GPU
 */
static int atmturb_cuda_ao_run_open(
    const atmturb_cuda_ao_req_t *req,
    const float                 *d_mask,
    cufftComplex                *d_in_grid,
    cufftComplex                *d_out_grid,
    float                       *d_psf,
    float                       *d_peaks,
    int                          batch)
{
    size_t       ntot       = (size_t) req->pup_size * (size_t) req->pup_size;
    const float *h_pha      = (req->swf_data != NULL) ? req->swf_data : req->wf_data;
    long         pha_stride = (req->swf_data != NULL) ? req->swf_stride : (long) ntot;
    float       *d_wf_batch = NULL;
    float       *d_amp_b    = NULL;
    int          ret        = -1;

    if (atmturb_cuda_ao_alloc_opt(h_pha, pha_stride, (size_t) batch, ntot, &d_wf_batch) != 0 ||
        atmturb_cuda_ao_alloc_opt(req->amp_data, req->amp_stride, (size_t) batch, ntot,
                                  &d_amp_b) != 0)
    {
        goto cleanup;
    }

    {
        float open_scale = (req->swf_data != NULL)
                               ? 1.0f
                               : (float) (req->params->lambda_ref_m / req->params->lambda_sci_m);
        float norm       = 1.0f / ((float) req->n_fft * (float) req->n_fft);
        int   off        = (int) ((req->n_fft - req->pup_size) / 2);
        dim3  blk(32, 8);
        dim3  grd_psf(((unsigned int) req->n_psf + 31U) / 32U,
                      ((unsigned int) req->n_psf + 7U) / 8U);

        for (long t0 = 0; t0 < req->nbframes; t0 += batch)
        {
            int count = (t0 + batch <= req->nbframes) ? batch : (int) (req->nbframes - t0);
            if (atmturb_cuda_ao_upload_chunk(d_wf_batch, h_pha, pha_stride,
                                             t0, count, ntot) != 0 ||
                atmturb_cuda_ao_upload_chunk(d_amp_b, req->amp_data, req->amp_stride,
                                             t0, count, ntot) != 0)
            {
                goto cleanup;
            }

            dim3 grd_pup(((unsigned int) req->pup_size + 31U) / 32U,
                         ((unsigned int) req->pup_size + 7U) / 8U, (unsigned int) count);
            atmturb_ao_open_embed_kargs_t ek = {
                d_wf_batch, d_amp_b, d_mask, d_in_grid, (int) req->pup_size, (int) req->n_fft,
                off, off, count, (pha_stride == 0) ? 1 : 0, (req->amp_stride == 0) ? 1 : 0,
                open_scale
            };
            atmturb_cuda_ao_open_embed_kernel<<<grd_pup, blk, 0, s_psf_stream>>>(ek);

            cufftHandle plan = (batch > 1) ? s_psf_plan_b : s_psf_plan_1;
            if (cufftExecC2C(plan, d_in_grid, d_out_grid, CUFFT_FORWARD) != CUFFT_SUCCESS)
            {
                goto cleanup;
            }

            int idx0 = (t0 == 0) ? 0 : -1;
            int idx1 = (req->nbframes > 1 && (req->nbframes - 1) >= t0 &&
                        (req->nbframes - 1) < t0 + count) ? (int) (req->nbframes - 1 - t0) : -1;
            atmturb_ao_accum_kargs_t ak = {
                d_out_grid, d_psf, d_peaks, (int) req->n_fft, (int) req->n_psf, count,
                req->do_accum, 1, idx0, 1, idx1, 2, norm
            };
            atmturb_cuda_ao_accum_kernel<<<grd_psf, blk, 0, s_psf_stream>>>(ak);
        }
    }
    ret = 0;

cleanup:
    cudaFree(d_wf_batch);
    cudaFree(d_amp_b);
    return ret;
}

/**
 * atmturb_cuda_ao_run_closed - Execute closed-loop AO controller and PSF synthesis on GPU
 */
static int atmturb_cuda_ao_run_closed(
    const atmturb_cuda_ao_req_t *req,
    const float                 *d_mask,
    cufftComplex                *d_in_grid,
    cufftComplex                *d_out_grid,
    float                       *d_psf,
    float                       *d_peaks)
{
    const atmturb_ao_params_t *p     = req->params;
    size_t                     ntot  = (size_t) req->pup_size * (size_t) req->pup_size;
    int                        delay = (p->loop_delay < 0) ? 0
                                                           : ((p->loop_delay > 3) ? 3
                                                                                  : p->loop_delay);
    int                        depth = delay + 1;
    long                       wchk  = (req->nbframes < ATMTURB_AO_CUDA_WF_CHUNK)
                                           ? req->nbframes : ATMTURB_AO_CUDA_WF_CHUNK;

    size_t state_floats = (size_t) (3 + depth) * ntot;
    float *d_state      = NULL;
    float *d_wf_chunk   = NULL;
    float *d_swf_chunk  = NULL;
    float *d_amp_chunk  = NULL;
    int    ret          = -1;

    if (cudaMalloc((void **) &d_state, sizeof(float) * state_floats) != cudaSuccess ||
        cudaMemsetAsync(d_state, 0, sizeof(float) * state_floats, s_psf_stream) != cudaSuccess ||
        atmturb_cuda_ao_alloc_opt(req->wf_data, (long) ntot, (size_t) wchk, ntot,
                                  &d_wf_chunk) != 0 ||
        atmturb_cuda_ao_alloc_opt(req->swf_data, req->swf_stride, (size_t) wchk, ntot,
                                  &d_swf_chunk) != 0 ||
        atmturb_cuda_ao_alloc_opt(req->amp_data, req->amp_stride, (size_t) wchk, ntot,
                                  &d_amp_chunk) != 0)
    {
        goto cleanup;
    }

    {
        float dt   = (p->time_step_s > 0.0) ? (float) p->time_step_s : 0.001f;
        float norm = 1.0f / ((float) req->n_fft * (float) req->n_fft);
        int   off  = (int) ((req->n_fft - req->pup_size) / 2);
        int   head = 0;
        dim3  blk(32, 8);
        dim3  grd_pup(((unsigned int) req->pup_size + 31U) / 32U,
                      ((unsigned int) req->pup_size + 7U) / 8U);
        dim3  grd_psf(((unsigned int) req->n_psf + 31U) / 32U,
                      ((unsigned int) req->n_psf + 7U) / 8U);

        atmturb_ao_ctrl_kargs_t ck = {
            NULL, NULL, NULL, d_mask, NULL, NULL, d_state, d_state + ntot, d_state + 2 * ntot,
            d_in_grid, (int) req->pup_size, (int) req->n_fft, off, off, p->loop_mode, delay,
            (float) (p->lambda_ref_m / (2.0 * M_PI)), (float) (2.0 * M_PI / p->lambda_sci_m),
            1.0f - (float) p->leak, (float) p->gain, (float) p->Kp, (float) p->Ki,
            (float) p->Kd / dt, dt
        };

        for (long t0 = 0; t0 < req->nbframes; t0 += wchk)
        {
            long clen = (t0 + wchk <= req->nbframes) ? wchk : (req->nbframes - t0);
            if (atmturb_cuda_ao_upload_chunk(d_wf_chunk, req->wf_data, (long) ntot,
                                             t0, clen, ntot) != 0 ||
                atmturb_cuda_ao_upload_chunk(d_swf_chunk, req->swf_data, req->swf_stride,
                                             t0, clen, ntot) != 0 ||
                atmturb_cuda_ao_upload_chunk(d_amp_chunk, req->amp_data, req->amp_stride,
                                             t0, clen, ntot) != 0)
            {
                goto cleanup;
            }

            for (long c = 0; c < clen; c++)
            {
                long   t       = t0 + c;
                int    r_idx   = (head - delay + depth) % depth;
                size_t swf_off = (req->swf_stride == 0) ? 0U : ((size_t) c * ntot);
                size_t amp_off = (req->amp_stride == 0) ? 0U : ((size_t) c * ntot);
                ck.d_wf_slice   = (d_wf_chunk != NULL) ? (d_wf_chunk + (size_t) c * ntot) : NULL;
                ck.d_swf_slice  = (d_swf_chunk != NULL) ? (d_swf_chunk + swf_off) : NULL;
                ck.d_amp_slice  = (d_amp_chunk != NULL) ? (d_amp_chunk + amp_off) : NULL;
                ck.d_read_slot  = d_state + (size_t) (3 + r_idx) * ntot;
                ck.d_write_slot = d_state + (size_t) (3 + head) * ntot;
                head            = (head + 1) % depth;

                atmturb_cuda_ao_ctrl_embed_kernel<<<grd_pup, blk, 0, s_psf_stream>>>(ck);
                if (cufftExecC2C(s_psf_plan_1, d_in_grid, d_out_grid,
                                 CUFFT_FORWARD) != CUFFT_SUCCESS)
                {
                    goto cleanup;
                }

                int idx0 = (t == 0) ? 0 : -1;
                int idx1 = (req->nbframes > 1 && t == req->nbframes - 1) ? 0 : -1;
                atmturb_ao_accum_kargs_t ak = {
                    d_out_grid, d_psf, d_peaks, (int) req->n_fft, (int) req->n_psf, 1,
                    req->do_accum, 1, idx0, 1, idx1, 2, norm
                };
                atmturb_cuda_ao_accum_kernel<<<grd_psf, blk, 0, s_psf_stream>>>(ak);
            }
        }
    }
    ret = 0;

cleanup:
    cudaFree(d_state);
    cudaFree(d_wf_chunk);
    cudaFree(d_swf_chunk);
    cudaFree(d_amp_chunk);
    return ret;
}

/**
 * atmturb_cuda_ao_sim_run - Execute open-loop or closed-loop AO PSF synthesis on CUDA GPU
 */
int atmturb_cuda_ao_sim_run(
    const atmturb_cuda_ao_req_t *req)
{
    if (req == NULL || req->params == NULL || req->pupil_mask == NULL || req->psf_out == NULL ||
        req->pup_size <= 0 || req->n_fft < 2 * req->pup_size || req->n_psf <= 0 ||
        req->nbframes <= 0 || req->params->loop_mode < 0 || req->params->loop_mode > 2 ||
        !atmturb_cuda_device_available())
    {
        return -1;
    }

    int batch = 1;
    if (req->params->loop_mode == 0)
    {
        batch = (req->n_fft <= 512) ? 16 : ((req->n_fft <= 1024) ? 8 : 4);
        if ((long) batch > req->nbframes)
        {
            batch = (int) req->nbframes;
        }
    }
    if (atmturb_cuda_psf_ensure_plans((int) req->n_fft, batch) != 0)
    {
        return -1;
    }

    size_t        ntot       = (size_t) req->pup_size * (size_t) req->pup_size;
    size_t        nf2        = (size_t) req->n_fft * (size_t) req->n_fft;
    size_t        npsf2      = (size_t) req->n_psf * (size_t) req->n_psf;
    float        *d_mask     = NULL;
    float        *d_psf      = NULL;
    float        *d_peaks    = NULL;
    cufftComplex *d_in_grid  = NULL;
    cufftComplex *d_out_grid = NULL;
    int           ret        = -1;

    if (cudaMalloc((void **) &d_mask, sizeof(float) * ntot) != cudaSuccess ||
        cudaMalloc((void **) &d_psf, sizeof(float) * npsf2) != cudaSuccess ||
        cudaMalloc((void **) &d_peaks, sizeof(float) * 3U) != cudaSuccess ||
        cudaMalloc((void **) &d_in_grid, sizeof(cufftComplex) * (size_t) batch * nf2)
            != cudaSuccess ||
        cudaMalloc((void **) &d_out_grid, sizeof(cufftComplex) * (size_t) batch * nf2)
            != cudaSuccess)
    {
        goto cleanup;
    }

    if (cudaMemcpyAsync(d_mask, req->pupil_mask, sizeof(float) * ntot,
                        cudaMemcpyHostToDevice, s_psf_stream) != cudaSuccess ||
        cudaMemsetAsync(d_psf, 0, sizeof(float) * npsf2, s_psf_stream) != cudaSuccess ||
        cudaMemsetAsync(d_peaks, 0, sizeof(float) * 3U, s_psf_stream) != cudaSuccess ||
        cudaMemsetAsync(d_in_grid, 0, sizeof(cufftComplex) * (size_t) batch * nf2,
                        s_psf_stream) != cudaSuccess)
    {
        goto cleanup;
    }

    {
        int  off = (int) ((req->n_fft - req->pup_size) / 2);
        dim3 blk(32, 8);
        dim3 grd_pup(((unsigned int) req->pup_size + 31U) / 32U,
                     ((unsigned int) req->pup_size + 7U) / 8U, 1U);
        dim3 grd_psf(((unsigned int) req->n_psf + 31U) / 32U,
                     ((unsigned int) req->n_psf + 7U) / 8U);
        atmturb_ao_open_embed_kargs_t dl_ek = {
            NULL, NULL, d_mask, d_in_grid, (int) req->pup_size, (int) req->n_fft,
            off, off, 1, 0, 0, 0.0f
        };
        atmturb_ao_accum_kargs_t dl_ak = {
            d_out_grid, d_psf, d_peaks, (int) req->n_fft, (int) req->n_psf, 1, 0, 0, 0, 0, -1, 0,
            1.0f / ((float) req->n_fft * (float) req->n_fft)
        };
        atmturb_cuda_ao_open_embed_kernel<<<grd_pup, blk, 0, s_psf_stream>>>(dl_ek);
        if (cufftExecC2C(s_psf_plan_1, d_in_grid, d_out_grid, CUFFT_FORWARD) != CUFFT_SUCCESS)
        {
            goto cleanup;
        }
        atmturb_cuda_ao_accum_kernel<<<grd_psf, blk, 0, s_psf_stream>>>(dl_ak);
    }

    if (req->params->loop_mode == 0)
    {
        if (atmturb_cuda_ao_run_open(req, d_mask, d_in_grid, d_out_grid, d_psf,
                                     d_peaks, batch) != 0)
        {
            goto cleanup;
        }
    }
    else if (atmturb_cuda_ao_run_closed(req, d_mask, d_in_grid, d_out_grid, d_psf, d_peaks) != 0)
    {
        goto cleanup;
    }

    {
        float h_peaks[3] = {0.0f, 0.0f, 0.0f};
        if (cudaMemcpyAsync(req->psf_out, d_psf, sizeof(float) * npsf2,
                            cudaMemcpyDeviceToHost, s_psf_stream) != cudaSuccess ||
            cudaMemcpyAsync(h_peaks, d_peaks, sizeof(float) * 3U,
                            cudaMemcpyDeviceToHost, s_psf_stream) != cudaSuccess ||
            cudaStreamSynchronize(s_psf_stream) != cudaSuccess)
        {
            goto cleanup;
        }

        double dl_peak = (double) h_peaks[0];
        if (req->dl_peak_out != NULL)
        {
            *req->dl_peak_out = dl_peak;
        }
        if (req->s_endpoints != NULL)
        {
            double s0           = (dl_peak > 0.0) ? ((double) h_peaks[1] / dl_peak) : 0.0;
            double s1           = (req->nbframes == 1)
                                      ? s0
                                      : ((dl_peak > 0.0) ? ((double) h_peaks[2] / dl_peak) : 0.0);
            req->s_endpoints[0] = s0;
            req->s_endpoints[1] = s1;
        }
    }
    ret = 0;

cleanup:
    cudaFree(d_mask);
    cudaFree(d_psf);
    cudaFree(d_peaks);
    cudaFree(d_in_grid);
    cudaFree(d_out_grid);
    return ret;
}
