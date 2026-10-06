# WFpropagate

**Architectural Layer**: Level 1 (Diffractive Optical Propagation)

## Purpose

The `WFpropagate` module implements Fresnel optical wavefront propagation using angular spectrum and Fresnel kernel methods. It provides single-plane propagation, multi-plane propagation cubes, and optical train simulations (such as Lyot coronagraphs) with diffraction effects.

## Components

- `WFpropagate.c` / `WFpropagate.h`: Public API module initialization and milk CLI command registration.
- `wfprop_fresnel.c`: Fresnel wavefront propagation kernels, transfer function initialization, and plane-to-plane propagation.
- `wfprop_cube.c`: Multi-distance propagation creating amplitude and phase data cubes.
- `wfprop_lyot.c`: Optical train simulation including focal plane masks, Lyot stops, and multi-plane Fresnel diffraction.

## Public Headers

- `WFpropagate.h`:
  - `int init_WFpropagate(void)`
  - `int Fresnel_propagate_wavefront(char *in, char *out, double PUPIL_SCALE, double z, double lambda)`
  - `int Init_Fresnel_propagate_wavefront(char *Cim, long size, double PUPIL_SCALE, double z, double lambda, double FPMASKRAD, int Precision)`
  - `int Fresnel_propagate_wavefront1(char *in, char *out, char *Cin)`
  - `long Fresnel_propagate_cube(char *IDcin_name, char *IDout_name_amp, char *IDout_name_pha, double PUPIL_SCALE, double zstart, double zend, long NBzpts, double lambda)`
  - `long WFpropagate_run(void)`

## Dependencies

- `CLIcore` (image memory management, command registration)
- `milkCOREMODarith`, `milkCOREMODmemory`, `milkfft`
- `OpticsMaterials`
