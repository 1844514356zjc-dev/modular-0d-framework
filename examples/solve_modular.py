"""Run three compositions with PYTHONPATH=src python examples/solve_modular.py.

Core configurations use SI magnitudes and radians, without UI imports.
"""

from __future__ import annotations

import numpy as np

from core.engine_solver import EngineConfig, EngineSimulator
from core.heat_models import OttoHeatRelease, WiebeCombustion
from core.mechanics import CrankMechanism, OttoCycleMechanism
from core.valves_ports import PoppetValve, Reservoir
from core.working_fluids import IdealGas


def main() -> None:
    """Solve dictionary ports, typed valves, and an exact closed Otto cycle."""
    gas = IdealGas(gamma=1.4)
    crank = CrankMechanism(compression_ratio=10)
    mass = 0.0003
    temperature = 800.0
    initial_pressure = gas.evaluate(
        temperature, mass / crank.get_volume(0)
    ).pressure
    valves = tuple(
        PoppetValve(
            name=name,
            reservoir=reservoir,
            opening_theta=np.deg2rad(opening),
            closing_theta=np.deg2rad(closing),
        )
        for name, reservoir, opening, closing in (
            ("intake", Reservoir(1.05e5, 320), 350, 580),
            ("exhaust", Reservoir(1e5, 700), 140, 370),
        )
    )
    four_stroke = EngineConfig(
        mechanism=crank,
        boundaries=valves,
        fluid=gas,
        combustion=WiebeCombustion(416.5, period=4 * np.pi),
        initial_state=(initial_pressure, temperature, mass),
    )
    # Serializable dicts can also mix constructed objects and nested specs.
    opposed = {
        "mechanism": {"type": "opposed"},
        "fluid": {"type": "ideal"},
        "combustion": {"type": "wiebe", "heat_per_event": 416.5},
        "boundaries": [
            {
                "type": "port",
                "name": "intake",
                "height": 0.0096,
                "width_total": 0.08,
                "edge_position": 0.0656,
                "reservoir": {"pressure": 1.3e5, "temperature": 320},
            },
            {
                "type": "port",
                "name": "exhaust",
                "height": 0.012,
                "width_total": 0.08,
                "edge_position": 0.06,
                "reservoir": {"pressure": 1e5, "temperature": 700},
            },
        ],
    }
    opposed_engine = EngineSimulator(opposed)
    opposed_state = opposed_engine.fluid.evaluate(
        temperature, mass / opposed_engine.mechanism.get_volume(0)
    )
    opposed_result = opposed_engine.solve_cycle(
        [opposed_state.pressure, temperature, mass], samples=361
    )
    print("Opposed ports:", opposed_result.work[-1], "J")
    valve_result = EngineSimulator(four_stroke).solve_cycle(samples=361)
    print("Four-stroke valves:", valve_result.work[-1], "J")
    volume_law = OttoCycleMechanism(crank)
    charge = gas.from_pressure_temperature(1e5, 300)
    otto = EngineConfig(
        mechanism=volume_law,
        fluid=gas,
        combustion=OttoHeatRelease(500, 10, gas.gamma),
        initial_state=(1e5, 300, charge.density * volume_law.get_volume(0)),
    )
    cycle = EngineSimulator(otto).solve_cycle(samples=361)
    print("Otto efficiency:", cycle.work[-1] / 500)
    print("Analytic efficiency:", 1 - 10 ** (1 - gas.gamma))
    print("Relative energy residual:", cycle.energy_relative_error)


if __name__ == "__main__":
    main()
