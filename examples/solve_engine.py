"""Run from project root: PYTHONPATH=src python examples/solve_engine.py."""

import numpy as np

from core.engine_solver import OpposedPistonEngineCycle


def main() -> None:
    """Save V [m³], p [Pa] for an illustrative 720-degree transient."""
    engine = OpposedPistonEngineCycle(
        bore=0.08,
        stroke=0.08,
        con_rod=0.14,
        rpm=2000.0,
        fuel_mass=2e-5,
        air_fuel_ratio=30.0,
        theta_s=np.deg2rad(-10),
        delta_theta_b=np.deg2rad(50),
        m=2.0,
    )
    result = engine.solve_cycle()
    np.savetxt("engine_pv.csv", result.pv, delimiter=",",
               header="volume_m3,pressure_Pa", comments="")
    print(f"EOS relative residual: {result.eos_relative_error:.3e}")
    print(f"Peak pressure: {result.pressure.max() / 1e5:.2f} bar")


if __name__ == "__main__":
    main()
