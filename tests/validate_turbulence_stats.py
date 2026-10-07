#!/usr/bin/env python3
"""Statistical validator for milkatmturb turbulence outputs.

Each subcommand loads one or more FITS files, computes a statistic, compares it with a
theoretical expectation, prints a one-line verdict and exits with status 0 (pass) or 1 (fail).
Exit status 77 means a required Python module is missing (treated as SKIP by the runner).

Subcommands:
  finite  - data are finite and not flat
  sf      - phase structure function vs discrete-PSD and/or continuous von Karman theory
  ratio   - RMS ratio between two cubes
  same    - two cubes agree to a relative tolerance
  shift   - mean frame-to-frame displacement from cross-correlation peak
  tilt    - Zernike tip/tilt variance over a circular pupil vs Noll / von Karman theory
  repeat  - correlation between frames separated by a fixed lag
  scint   - scintillation index and mean intensity of an amplitude cube
  corr    - correlation coefficient between planes of a cube
"""
import argparse
import math
import sys

try:
    import numpy as np
    from astropy.io import fits
    from scipy import special
except ImportError as exc:  # pragma: no cover
    print(f"SKIP: missing python module ({exc})")
    sys.exit(77)

ARCSEC = math.pi / 180.0 / 3600.0


# ----------------------------------------------------------------------------------------------
# I/O helpers
# ----------------------------------------------------------------------------------------------
def load_cube(fname):
    """Load a FITS file as a float64 3D array (nframes, ny, nx)."""
    with fits.open(fname) as hdul:
        data = np.asarray(hdul[0].data, dtype=np.float64)
    if data.ndim == 2:
        data = data[np.newaxis, :, :]
    return data


def verdict(ok, msg):
    """Print verdict line and return process exit status."""
    print(("PASS: " if ok else "FAIL: ") + msg)
    return 0 if ok else 1


# ----------------------------------------------------------------------------------------------
# Theory
# ----------------------------------------------------------------------------------------------
def r0_from_seeing(seeing_arcsec, lam_m):
    """Fried parameter [m] from Kolmogorov seeing FWHM [arcsec] at wavelength lam_m."""
    return 0.98 * lam_m / (seeing_arcsec * ARCSEC)


def sf_vonkarman(r, r0, L0):
    """Continuous phase structure function [rad^2]; L0 <= 0 means Kolmogorov."""
    r = np.asarray(r, dtype=np.float64)
    if L0 is None or L0 <= 0:
        return 6.88 * (r / r0) ** (5.0 / 3.0)
    x = 2.0 * math.pi * r / L0
    bracket = 1.0 - (2.0 ** (1.0 / 6.0) / special.gamma(5.0 / 6.0)) * x ** (5.0 / 6.0) \
        * special.kv(5.0 / 6.0, x)
    return 0.17253 * (L0 / r0) ** (5.0 / 3.0) * bracket


def sf_discrete(lags_pix, nmaster, dx_m, r0, L0, l0):
    """Ensemble structure function of an FFT screen (master grid nmaster, pixel dx_m).

    Uses the analytic per-mode variance of the generator:
    sigma_k^2 = 0.023 (r0/dx)^(-5/3) N^(5/3) (k^2 + k0^2)^(-11/6) exp(-k^2/km^2), DC = 0.
    lags_pix are in master pixels (along x).
    """
    k = np.fft.fftfreq(nmaster) * nmaster
    kx, ky = np.meshgrid(k, k)
    k2 = kx ** 2 + ky ** 2
    k0 = nmaster * dx_m / L0 if (L0 is not None and L0 > 0) else 0.0
    with np.errstate(divide="ignore"):
        psd = 0.023 * (r0 / dx_m) ** (-5.0 / 3.0) * nmaster ** (5.0 / 3.0) \
            * (k2 + k0 ** 2) ** (-11.0 / 6.0)
    if l0 is not None and l0 > 0:
        km = (5.92 / (2.0 * math.pi)) * nmaster * dx_m / l0
        psd *= np.exp(-k2 / km ** 2)
    psd[0, 0] = 0.0
    cov = np.real(np.fft.ifft2(psd)) * nmaster * nmaster
    lags = np.asarray(lags_pix, dtype=np.int64)
    return 2.0 * (cov[0, 0] - cov[0, lags % nmaster])


def measure_sf(cube, lags):
    """Empirical structure function along x and y, averaged over frames."""
    out = []
    for lag in lags:
        dx = cube[:, :, lag:] - cube[:, :, :-lag]
        dy = cube[:, lag:, :] - cube[:, :-lag, :]
        out.append(0.5 * (np.mean(dx ** 2) + np.mean(dy ** 2)))
    return np.asarray(out)


# ----------------------------------------------------------------------------------------------
# Subcommands
# ----------------------------------------------------------------------------------------------
def cmd_finite(args):
    cube = load_cube(args.cube)
    ok_finite = bool(np.all(np.isfinite(cube)))
    std = float(np.std(cube))
    return verdict(ok_finite and std > args.min_std,
                   f"{args.cube}: finite={ok_finite} std={std:.4g} (min {args.min_std})")


def resolve_r0(args):
    """r0 [m] from --r0 or from --seeing/--lam, scaled by airmass (cos z)^(3/5)."""
    r0 = args.r0 if args.r0 else r0_from_seeing(args.seeing, args.lam)
    return r0 * math.cos(args.zenith) ** 0.6


def cmd_sf(args):
    cube = load_cube(args.cube)
    if args.subtract_piston:
        cube = cube - cube.mean(axis=(1, 2), keepdims=True)
    npix = cube.shape[-1]
    lag_max = args.lag_max if args.lag_max else max(args.lag_min + 1, npix // 16)
    lags = np.arange(args.lag_min, lag_max + 1)
    dmeas = measure_sf(cube, lags)
    r0 = resolve_r0(args)
    r_m = lags * args.pixscale
    if args.reference == "discrete":
        os_ = args.oversample
        dref = sf_discrete(lags * os_, args.master_size, args.pixscale / os_, r0, args.L0, args.l0)
    else:
        dref = sf_vonkarman(r_m, r0, args.L0)
    ratio = dmeas / dref
    r0_fit = r0 * float(np.mean(ratio)) ** (-0.6)
    worst = float(np.max(np.abs(ratio - 1.0)))
    err = abs(r0_fit / r0 - 1.0)
    ok = err <= args.tol and (args.lag_tol <= 0 or worst <= args.lag_tol)
    return verdict(ok, f"sf[{args.reference}] r0 expected {r0:.5g} m, fitted {r0_fit:.5g} m "
                       f"(err {100 * err:.2f}%, tol {100 * args.tol:.1f}%), "
                       f"worst lag dev {100 * worst:.1f}% over lags {lags[0]}..{lags[-1]}")


def cmd_ratio(args):
    a = load_cube(args.cube_a)
    b = load_cube(args.cube_b)
    if args.subtract_piston:
        a = a - a.mean(axis=(1, 2), keepdims=True)
        b = b - b.mean(axis=(1, 2), keepdims=True)
    meas = float(np.std(a) / np.std(b))
    if args.variance:
        meas = meas ** 2
    err = abs(meas / args.expect - 1.0)
    kind = "variance" if args.variance else "rms"
    return verdict(err <= args.tol, f"{kind} ratio {meas:.5g}, expected {args.expect:.5g} "
                                    f"(err {100 * err:.2f}%, tol {100 * args.tol:.1f}%)")


def cmd_same(args):
    a = load_cube(args.cube_a)
    b = load_cube(args.cube_b)
    if a.shape != b.shape:
        return verdict(False, f"shape mismatch {a.shape} vs {b.shape}")
    rel = float(np.max(np.abs(a - b)) / max(np.std(a), 1e-30))
    return verdict(rel <= args.tol, f"max |a-b| / std(a) = {rel:.3g} (tol {args.tol:.3g})")


def xcorr_peak(f0, f1):
    """Sub-pixel (dx, dy) such that f1(x) ~ f0(x - d), via FFT cross-correlation."""
    f0 = f0 - f0.mean()
    f1 = f1 - f1.mean()
    cc = np.real(np.fft.ifft2(np.conj(np.fft.fft2(f0)) * np.fft.fft2(f1)))
    ny, nx = cc.shape
    iy, ix = np.unravel_index(np.argmax(cc), cc.shape)

    def parabolic(cm, c0, cp):
        den = cm - 2.0 * c0 + cp
        return 0.0 if den == 0 else 0.5 * (cm - cp) / den

    dx = ix + parabolic(cc[iy, (ix - 1) % nx], cc[iy, ix], cc[iy, (ix + 1) % nx])
    dy = iy + parabolic(cc[(iy - 1) % ny, ix], cc[iy, ix], cc[(iy + 1) % ny, ix])
    dx = dx - nx if dx > nx / 2 else dx
    dy = dy - ny if dy > ny / 2 else dy
    return dx, dy


def cmd_shift(args):
    if args.cube_b:
        a = load_cube(args.cube)
        b = load_cube(args.cube_b)
        pairs = [(a[i], b[i]) for i in range(a.shape[0])]
    else:
        cube = load_cube(args.cube)
        pairs = [(cube[i], cube[i + args.step]) for i in range(cube.shape[0] - args.step)]
    d = np.array([xcorr_peak(p0, p1) for p0, p1 in pairs])
    mx, my = float(d[:, 0].mean()), float(d[:, 1].mean())
    err = math.hypot(mx - args.dx, my - args.dy)
    return verdict(err <= args.tol, f"mean shift ({mx:.3f}, {my:.3f}) px, expected "
                                    f"({args.dx:.3f}, {args.dy:.3f}) (err {err:.3f}, "
                                    f"tol {args.tol:.3f})")


def tilt_variance_theory(D, r0, L0):
    """Two-axis Zernike tip+tilt phase variance [rad^2] over a circular pupil of diameter D."""
    kolmo = 0.896 * (D / r0) ** (5.0 / 3.0)
    if L0 is None or L0 <= 0:
        return kolmo
    # Integrate von Karman PSD against the Zernike tilt filter (Noll 1976 / Sasiela).
    f = np.logspace(-6, 3, 20000) / D
    psd = 0.023 * r0 ** (-5.0 / 3.0) * (f ** 2 + 1.0 / L0 ** 2) ** (-11.0 / 6.0)
    x = math.pi * D * f
    filt = 4.0 * (2.0 * special.jv(2, x) / x) ** 2  # (n+1)=2 per mode, two modes
    integrate = getattr(np, "trapezoid", None) or np.trapz
    return float(integrate(psd * filt * 2.0 * math.pi * f, f))


def cmd_tilt(args):
    cube = load_cube(args.cube)
    _, ny, nx = cube.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    rad = (min(nx, ny) / 2.0)
    xn = (xx - (nx - 1) / 2.0) / rad
    yn = (yy - (ny - 1) / 2.0) / rad
    mask = (xn ** 2 + yn ** 2) <= 1.0
    z2, z3 = 2.0 * xn[mask], 2.0 * yn[mask]
    var = 0.0
    for frame in cube:
        v = frame[mask] - frame[mask].mean()
        a2 = np.sum(v * z2) / np.sum(z2 * z2)
        a3 = np.sum(v * z3) / np.sum(z3 * z3)
        var += a2 * a2 + a3 * a3
    var /= cube.shape[0]
    D = 2.0 * rad * args.pixscale
    expect = tilt_variance_theory(D, resolve_r0(args), args.L0)
    err = abs(var / expect - 1.0)
    return verdict(err <= args.tol, f"tip/tilt variance {var:.4g} rad^2, expected {expect:.4g} "
                                    f"(err {100 * err:.1f}%, tol {100 * args.tol:.0f}%)")


def cmd_repeat(args):
    cube = load_cube(args.cube)
    cube = cube - cube.mean(axis=(1, 2), keepdims=True)
    n = cube.shape[0] - args.lag
    if n <= 0:
        return verdict(False, f"cube too short ({cube.shape[0]}) for lag {args.lag}")
    cc = [np.sum(cube[i] * cube[i + args.lag]) /
          math.sqrt(np.sum(cube[i] ** 2) * np.sum(cube[i + args.lag] ** 2)) for i in range(n)]
    worst = float(np.max(np.abs(cc)))
    return verdict(worst <= args.max_corr, f"max |corr| at lag {args.lag} frames = {worst:.3f} "
                                           f"(max {args.max_corr:.3f})")


def cmd_scint(args):
    amp = load_cube(args.cube)
    inten = amp ** 2
    mean_i = float(inten.mean())
    sig2 = float(inten.var() / mean_i ** 2)
    ok = abs(mean_i - 1.0) <= args.mean_tol
    msg = f"mean intensity {mean_i:.4f} (tol {args.mean_tol}), sigma_I^2 = {sig2:.4g}"
    if args.expect is not None:
        err = abs(sig2 / args.expect - 1.0)
        ok = ok and err <= args.tol
        msg += f", expected {args.expect:.4g} (err {100 * err:.1f}%, tol {100 * args.tol:.0f}%)"
    return verdict(ok, msg)


def cmd_corr(args):
    data = load_cube(args.cube)
    planes = data.reshape(data.shape[0], -1)
    a = planes[args.plane_a] - planes[args.plane_a].mean()
    b = planes[args.plane_b] - planes[args.plane_b].mean()
    rho = float(np.sum(a * b) / math.sqrt(np.sum(a * a) * np.sum(b * b)))
    return verdict(abs(rho) <= args.max_corr, f"corr(plane {args.plane_a}, plane {args.plane_b})"
                                              f" = {rho:.4f} (max |corr| {args.max_corr})")


# ----------------------------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------------------------
def add_r0_args(p):
    p.add_argument("--r0", type=float, default=0.0, help="Fried parameter [m] at zenith")
    p.add_argument("--seeing", type=float, default=0.0, help="seeing [arcsec] (if no --r0)")
    p.add_argument("--lam", type=float, default=0.5e-6, help="reference wavelength [m]")
    p.add_argument("--zenith", type=float, default=0.0, help="zenith angle [rad]")
    p.add_argument("--L0", type=float, default=0.0, help="outer scale [m], <=0 infinite")
    p.add_argument("--pixscale", type=float, required=True, help="pupil pixel [m]")


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("finite")
    p.add_argument("cube")
    p.add_argument("--min-std", type=float, default=0.0)
    p.set_defaults(func=cmd_finite)

    p = sub.add_parser("sf")
    p.add_argument("cube")
    add_r0_args(p)
    p.add_argument("--l0", type=float, default=0.0, help="inner scale [m]")
    p.add_argument("--reference", choices=("continuous", "discrete"), default="continuous")
    p.add_argument("--master-size", type=int, default=0, help="master grid size (discrete)")
    p.add_argument("--oversample", type=int, default=1, help="master px per pupil px")
    p.add_argument("--lag-min", type=int, default=2)
    p.add_argument("--lag-max", type=int, default=0)
    p.add_argument("--tol", type=float, default=0.05, help="relative tolerance on r0")
    p.add_argument("--lag-tol", type=float, default=0.0, help="max per-lag D deviation")
    p.add_argument("--subtract-piston", action="store_true")
    p.set_defaults(func=cmd_sf)

    p = sub.add_parser("ratio")
    p.add_argument("cube_a")
    p.add_argument("cube_b")
    p.add_argument("--expect", type=float, required=True)
    p.add_argument("--tol", type=float, default=0.03)
    p.add_argument("--variance", action="store_true")
    p.add_argument("--subtract-piston", action="store_true")
    p.set_defaults(func=cmd_ratio)

    p = sub.add_parser("same")
    p.add_argument("cube_a")
    p.add_argument("cube_b")
    p.add_argument("--tol", type=float, default=1e-4)
    p.set_defaults(func=cmd_same)

    p = sub.add_parser("shift")
    p.add_argument("cube")
    p.add_argument("--cube-b", default="", help="compare frame i of cube with frame i of B")
    p.add_argument("--step", type=int, default=1)
    p.add_argument("--dx", type=float, default=0.0)
    p.add_argument("--dy", type=float, default=0.0)
    p.add_argument("--tol", type=float, default=0.2)
    p.set_defaults(func=cmd_shift)

    p = sub.add_parser("tilt")
    p.add_argument("cube")
    add_r0_args(p)
    p.add_argument("--tol", type=float, default=0.15)
    p.set_defaults(func=cmd_tilt)

    p = sub.add_parser("repeat")
    p.add_argument("cube")
    p.add_argument("--lag", type=int, required=True)
    p.add_argument("--max-corr", type=float, default=0.1)
    p.set_defaults(func=cmd_repeat)

    p = sub.add_parser("scint")
    p.add_argument("cube")
    p.add_argument("--expect", type=float, default=None)
    p.add_argument("--tol", type=float, default=0.2)
    p.add_argument("--mean-tol", type=float, default=0.01)
    p.set_defaults(func=cmd_scint)

    p = sub.add_parser("corr")
    p.add_argument("cube")
    p.add_argument("--plane-a", type=int, default=0)
    p.add_argument("--plane-b", type=int, default=1)
    p.add_argument("--max-corr", type=float, default=0.05)
    p.set_defaults(func=cmd_corr)
    return ap


def main():
    args = build_parser().parse_args()
    if getattr(args, "r0", None) == 0.0 and getattr(args, "seeing", 0.0) == 0.0 \
            and args.cmd in ("sf", "tilt"):
        print("FAIL: --r0 or --seeing is required")
        return 1
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
