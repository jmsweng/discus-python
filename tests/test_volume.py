import numpy as np
import pytest
from discuspdf import delta_pdf_from_diffuse


def test_fft_peak_for_constant_intensity():
    n = 16
    q = (np.arange(n) - n//2) * 0.2
    intensity = np.ones((n, n, n))
    volume, r = delta_pdf_from_diffuse(intensity, (q, q, q))
    assert len(r) == 3
    center = n//2
    assert np.argmax(np.abs(volume)) == np.ravel_multi_index((center,)*3, (n,)*3)
    assert np.isfinite(volume).all()


def test_reject_rlu_axis():
    n = 4
    q = np.arange(n) * 0.2
    with pytest.raises(ValueError, match="centered"):
        delta_pdf_from_diffuse(np.ones((n,n,n)), (q,q,q))
