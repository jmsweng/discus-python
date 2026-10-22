import numpy as np
import pytest

from discuspdf import ReciprocalVolume, delta_pdf_from_diffuse


def axes(n, dq=.2):
    return (np.arange(n) - n // 2) * dq


@pytest.mark.parametrize("n", [15, 16])
def test_constant_intensity_fft(n):
    q = axes(n)
    result = delta_pdf_from_diffuse(np.ones((n, n, n)), (q, q, q))
    expected = (n * .2)**3 / (2 * np.pi)**3
    assert result.complex_volume.shape == (n, n, n)
    assert np.iscomplexobj(result.complex_volume)
    assert np.isclose(result.complex_volume[n//2, n//2, n//2].real, expected)
    assert result.diagnostics.imaginary_ratio < 1e-12
    assert result.diagnostics.friedel_relative_error < 1e-12


@pytest.mark.parametrize("n", [15, 16])
def test_analytic_cosine_pair(n):
    q = axes(n)
    dq = .2
    r0 = 2*np.pi / (n*dq) * 2
    intensity = 2 + np.cos(q[:, None, None] * r0) * np.ones((1,n,n))
    result = delta_pdf_from_diffuse(intensity, (q,q,q))
    line = np.abs(result.real_volume[:, n//2, n//2])
    assert line[n//2 + 2] > .1 * line[n//2]
    assert line[n//2 - 2] > .1 * line[n//2]
    assert result.diagnostics.imaginary_ratio < 1e-10


@pytest.mark.parametrize("n", [15, 16])
def test_unequal_axis_steps(n):
    qs = (axes(n,.1), axes(n,.15), axes(n,.25))
    result = delta_pdf_from_diffuse(np.ones((n,n,n)), qs)
    assert np.allclose(result.diagnostics.dq, (.1, .15, .25))
    assert np.isclose(result.r_axes[0][1] - result.r_axes[0][0], 2*np.pi/(n*.1))


def test_complex_and_threshold():
    n=16
    q=axes(n)
    a=np.zeros((n,n,n))
    a[n//2+1,n//2,n//2]=1.
    with pytest.raises(ValueError, match="exceeds"):
        delta_pdf_from_diffuse(a, (q,q,q), imaginary_tolerance=1e-8)
    with pytest.warns(RuntimeWarning):
        result=delta_pdf_from_diffuse(a, (q,q,q), imaginary_tolerance=1e-8,
                                      raise_on_imaginary=False)
    assert result.diagnostics.imaginary_ratio > 1e-8
    assert np.max(np.abs(result.complex_volume.imag)) > 0


def test_metadata_gate():
    q=axes(16)
    vol=ReciprocalVolume(np.ones((16,)*3), (q,q,q), "cartesian_q",
                         "1/angstrom", preprocessing=("bragg_excluded",))
    assert delta_pdf_from_diffuse(vol).diagnostics.imaginary_ratio < 1e-12
    missing=ReciprocalVolume(np.ones((16,)*3), (q,q,q), "cartesian_q", "1/angstrom")
    with pytest.raises(ValueError, match="Bragg"):
        delta_pdf_from_diffuse(missing)
