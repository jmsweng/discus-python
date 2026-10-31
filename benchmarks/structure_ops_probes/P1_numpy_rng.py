"""P1: Record NumPy RNG values (record only; tests never pin these)."""
import numpy as np

print(f"numpy.__version__ = {np.__version__}")
print(f"default_rng(1234).choice(np.arange(10), size=4, replace=False) = {np.random.default_rng(1234).choice(np.arange(10), size=4, replace=False)}")
print(f"default_rng(1234).normal(0,1,(2,3)) = {np.random.default_rng(1234).normal(0, 1, (2, 3))}")
print(f"default_rng(1234).random(3) = {np.random.default_rng(1234).random(3)}")
