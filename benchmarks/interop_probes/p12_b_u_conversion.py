"""P12: B/U conversion (gemmi.u_to_b vs 8*pi^2*U)."""
import math

def main():
    print("=== P12: B/U conversion ===")
    import gemmi

    # gemmi does not have a u_to_b function
    print(f"gemmi has u_to_b: {hasattr(gemmi, 'u_to_b')}")

    u = 0.01
    b_manual = 8 * math.pi**2 * u
    print(f"8*pi^2*{u} = {b_manual}")

    # Check if there's any B-related function in gemmi
    b_funcs = [a for a in dir(gemmi) if 'b' in a.lower() and ('factor' in a.lower() or 'conv' in a.lower())]
    print(f"gemmi B-related functions: {b_funcs}")

    # The conversion is simply B = 8*pi^2*U (in Angstrom^2 when U is in Angstrom^2)
    # Verify with a known value
    u2 = 0.008
    b2 = 8 * math.pi**2 * u2
    print(f"8*pi^2*{u2} = {b2}")

    print("=== P12 done ===")

if __name__ == "__main__":
    main()
