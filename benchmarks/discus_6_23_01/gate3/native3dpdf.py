"""Line-by-line numpy emulation of DISCUS_SUITE 6.23.01 `output -> value 3DPDF` for a 3D cubic-index grid.

Source (tag v.6.23.01, commit 19a214a):
  graphic.f90  'run':  dsi = dsi3d  (for set aver != 0, dsi3d = dsi; four_weight in fourier_sup.f90)
  graphic.f90  out_prep_3dpdf / out_prep_3dpdf_3d:
      maptofftfd_33D_R_C_DP   : forward map, centre index  i0 = INT(n/2) + 1 - MOD(n-1,2)  (1-based) -> FFT index 1
      fftw_plan_dft_3d(FFTW_FORWARD)   (sign exp(-2 pi i m k / n), unnormalized)
      out /= sqrt(n1 n2 n3)
      fftfdtoline_33D_C_R_DP  : reverse map with axis reversal; real part only
      out_vi(j,i)  = 1 / (vi(j,i) (n_i - 1))       written as data_increment_vector
      out_eck(i,1) = - out_vi * INT((n-1)/2)        written as data_corner
  qval.f90: value 3DPDF -> REAL(rpdf);  write_unified_data: data_values * out_scale (out_scale = 1 unless set).
Only identity dsort is supported (all three dimensions equal, or already in decreasing order).
"""
import numpy as np


def _fwd_index(n):
    l = (n - 1) % 2
    i = np.arange(1, n + 1)
    return (i + n // 2 - l) % n          # 0-based FFT index of input index i (1-based)


def _rev_index(n):
    l = (n - 1) % 2
    ii = np.arange(1, n + 1)
    return n + 1 - (((ii + n // 2 - l - 1) % n) + 1) - 1   # 0-based output index of FFT index ii


def bragg_reweight(I, corner, vi_diag, ncells, tol=1e-4):
    """four_weight (fourier_sup.f90) for set aver = 0: Bragg voxels (all |hkl - nint| < 1e-4) of dsi3d are divided by
    scalef = ncells * voxel / V*_cell = ncells * |det(vi)| (fractional voxel volume), if |scalef - 1| > 1e-5."""
    scalef = ncells * float(np.prod(vi_diag))
    if abs(scalef - 1.0) <= 1e-5:
        return I.copy(), scalef
    ax = [c + v * np.arange(k) for c, v, k in zip(corner, vi_diag, I.shape)]
    isint = [np.abs(a - np.rint(a)) < tol for a in ax]
    J = I.copy()
    m = np.ix_(*isint)
    J[m] = J[m] / scalef
    return J, scalef


def native_3dpdf(I):
    n = I.shape
    assert list(n) == sorted(n, reverse=True), "only identity dsort emulated"
    a = np.zeros(n, complex)
    f = [_fwd_index(k) for k in n]
    a[np.ix_(*f)] = I
    out = np.fft.fftn(a) / np.sqrt(np.prod(n))
    res = np.empty(n)
    r = [_rev_index(k) for k in n]
    res[np.ix_(*r)] = out.real
    return res


def native_axes(n, vi_diag):
    """(corner, increment) as DISCUS writes them for value 3DPDF."""
    inc = np.array([1.0 / (v * (k - 1)) for v, k in zip(vi_diag, n)])
    return -inc * np.array([(k - 1) // 2 for k in n]), inc


def true_axes(n, vi_diag):
    """(corner, increment) of the same values on the discrete-Fourier grid: step 1/(n dH); u = 0 at the output index
    that receives FFT index 0."""
    inc = np.array([1.0 / (v * k) for v, k in zip(vi_diag, n)])
    p0 = np.array([int(_rev_index(k)[0]) for k in n])
    return -inc * p0, inc
