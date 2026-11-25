"""Invariants for interchangeable mechanisms, boundaries and EOS."""

import json

import numpy as np
import pytest

from core.engine_solver import EngineConfig, EngineSimulator
from core.heat_models import OttoHeatRelease
from core.mechanics import (
    BaseMechanism,
    CrankMechanism,
    OpposedPistonMechanism,
    OttoCycleMechanism,
)
from core.valves_ports import (
    PistonMaskedPort,
    PoppetValve,
    Reservoir,
    nozzle_mass_flow,
)
from core.working_fluids import CoolPropGas, IdealGas
from ui.engine_dashboard import EngineSettings, json_export, simulate


@pytest.mark.parametrize("kind", [CrankMechanism, OpposedPistonMechanism])
def test_mechanism_derivative_and_displacement(
    kind: type[CrankMechanism],
) -> None:
    """Check analytic dV/dtheta [m³/rad], ratio and total piston travel."""
    mechanism = kind(compression_ratio=12)
    assert mechanism.get_volume(np.pi) / mechanism.get_volume(
        0
    ) == pytest.approx(12)
    assert mechanism.swept_volume == pytest.approx(
        mechanism.piston_count * mechanism.piston_area * mechanism.stroke
    )
    for theta in np.linspace(0.1, 6.1, 19):
        step = 1e-6
        numeric = (
            mechanism.get_volume(theta + step)
            - mechanism.get_volume(theta - step)
        ) / (2 * step)
        assert mechanism.get_dV_dtheta(theta) == pytest.approx(
            numeric, rel=1e-7
        )
    with pytest.raises(ValueError):
        kind(con_rod=0.01)


def test_cam_lift_timing_and_area() -> None:
    """Check curtain/throat limitation and a 720°, not 360°, cam period."""
    valve = PoppetValve(
        name="valve",
        reservoir=Reservoir(),
        diameter=0.03,
        max_lift=0.02,
        opening_theta=0,
        closing_theta=np.pi,
    )
    assert valve.get_lift(np.pi / 2) == pytest.approx(0.02)
    assert valve.get_geometric_area(np.pi / 2) == pytest.approx(
        np.pi * 0.03**2 / 4
    )
    assert valve.get_lift(np.pi) == 0
    assert valve.get_lift(2.5 * np.pi) == 0
    assert valve.get_lift(4.5 * np.pi) == pytest.approx(0.02)
    wrapped = PoppetValve(
        name="wrapped",
        reservoir=Reservoir(),
        opening_theta=np.deg2rad(700),
        closing_theta=np.deg2rad(20),
    )
    assert wrapped.get_lift(0) == pytest.approx(wrapped.max_lift)


def test_nozzle_inward_sign_choking_and_cd() -> None:
    """Signed nozzle [kg/s] is linear in Cd-corrected area; sonic flow caps."""
    fluid = IdealGas()
    reservoir = fluid.from_pressure_temperature(2e5, 320)
    low = fluid.from_pressure_temperature(1e4, 600)
    medium = fluid.from_pressure_temperature(2e4, 600)
    high = fluid.from_pressure_temperature(4e5, 600)
    first = nozzle_mass_flow(1e-4, low, reservoir)
    assert first > 0
    assert nozzle_mass_flow(1e-4, medium, reservoir) == pytest.approx(first)
    assert nozzle_mass_flow(2e-4, low, reservoir) == pytest.approx(2 * first)
    assert nozzle_mass_flow(1e-4, high, reservoir) < 0


def test_custom_mechanism_closed_invariants() -> None:
    """A volume-only extension satisfies closed adiabatic invariants."""

    class HarmonicVolume(BaseMechanism):
        """Independent sinusoidal test volume, in m³ and radians."""

        def get_volume(self, theta: float) -> float:
            """Return positive volume [m³]."""
            return float(1e-3 * (1.5 - 0.5 * np.cos(theta)))

        def get_dV_dtheta(self, theta: float) -> float:
            """Return analytic derivative [m³/rad]."""
            return float(0.5e-3 * np.sin(theta))

    gas = IdealGas(gamma=1.4)
    engine = EngineSimulator(
        EngineConfig(
            mechanism=HarmonicVolume(),
            fluid=gas,
            theta_start=0.3,
            theta_end=6.1,
        )
    )
    result = engine.solve_cycle(samples=121)
    invariant = result.pressure * result.volume**gas.gamma
    np.testing.assert_allclose(invariant, invariant[0], rtol=3e-7)
    np.testing.assert_allclose(result.mass, result.mass[0])
    assert result.boundary_flows == {}
    assert result.energy_relative_error < 1e-7


def test_dict_configuration_and_mixed_boundaries() -> None:
    """A JSON configuration combines ports and valves in one cylinder."""
    specification = {
        "mechanism": {"type": "crank"},
        "fluid": {"type": "ideal"},
        "combustion": {"type": "none"},
        "boundaries": [
            {
                "type": "port",
                "name": "side_port",
                "height": 0.005,
                "width_total": 0.04,
                "edge_position": 0.06,
                "reservoir": {"pressure": 1.3e5, "temperature": 320},
            },
            {
                "type": "poppet",
                "name": "head_valve",
                "opening_theta": 1,
                "closing_theta": 5,
                "reservoir": {"pressure": 1e5, "temperature": 700},
            },
        ],
    }
    snapshot = json.dumps(specification, sort_keys=True)
    engine = EngineSimulator(specification)
    assert isinstance(engine.config.boundaries[0], PistonMaskedPort)
    assert isinstance(engine.config.boundaries[1], PoppetValve)
    result = engine.solve_cycle(samples=181)
    assert set(result.boundary_flows) == {"side_port", "head_valve"}
    assert np.ptp(result.mass) > 0
    assert result.energy_relative_error < 2e-6
    assert json.dumps(specification, sort_keys=True) == snapshot
    with pytest.raises(ValueError, match="Unknown mechanism"):
        EngineSimulator({"mechanism": {"type": "unknown"}})


def test_exact_air_standard_otto() -> None:
    """Compare a complete closed heat-engine cycle with eta=1-r^(1-gamma)."""
    ratio, gamma, heat = 10.0, 1.4, 500.0
    mechanism = OttoCycleMechanism(CrankMechanism(compression_ratio=ratio))
    fluid = IdealGas(gamma=gamma)
    initial = fluid.from_pressure_temperature(1e5, 300)
    config = EngineConfig(
        mechanism=mechanism,
        fluid=fluid,
        combustion=OttoHeatRelease(heat, ratio, gamma),
        initial_state=(1e5, 300, initial.density * mechanism.get_volume(0)),
    )
    result = EngineSimulator(config).solve_cycle(samples=361)
    assert result.work[-1] / heat == pytest.approx(
        1 - ratio ** (1 - gamma), rel=2e-6
    )
    assert result.temperature[-1] == pytest.approx(300, rel=2e-6)
    assert result.pressure[-1] == pytest.approx(1e5, rel=2e-6)
    assert result.energy_relative_error < 1e-7
    with pytest.raises(ValueError, match="Otto benchmark"):
        EngineConfig(mechanism=CrankMechanism(), combustion=OttoHeatRelease())


def test_coolprop_derivatives_and_adiabatic_entropy() -> None:
    """Check real-EOS derivatives and closed isentropic motion."""
    fluid = CoolPropGas("Nitrogen")
    temperature, density, step = 500.0, 20.0, 1e-3
    state = fluid.evaluate(temperature, density)
    dpt = (
        fluid.evaluate(temperature + step, density).pressure
        - fluid.evaluate(temperature - step, density).pressure
    ) / (2 * step)
    durho = (
        fluid.evaluate(temperature, density + step).internal_energy
        - fluid.evaluate(temperature, density - step).internal_energy
    ) / (2 * step)
    assert state.dp_dT == pytest.approx(dpt, rel=1e-7)
    assert state.du_drho == pytest.approx(durho, rel=1e-6)
    engine = EngineSimulator(
        EngineConfig(mechanism=CrankMechanism(), fluid=fluid)
    )
    result = engine.solve_cycle(samples=61)
    entropy = [
        fluid.evaluate(t, m / v).entropy
        for t, m, v in zip(result.temperature, result.mass, result.volume)
    ]
    np.testing.assert_allclose(entropy, entropy[0], atol=2e-4, rtol=0)
    assert result.eos_relative_error < 2e-6
    assert result.energy_relative_error < 2e-6
    with pytest.raises(ValueError, match="gas phase"):
        CoolPropGas("Water").from_pressure_temperature(1e5, 300)


@pytest.mark.parametrize(
    "settings",
    [
        EngineSettings(topology="four_stroke", compression_ratio=10),
        EngineSettings(topology="closed"),
        EngineSettings(
            topology="closed",
            closed_model="Adiabatic gas spring",
            fluid="Nitrogen",
        ),
    ],
)
def test_topology_reports(settings: EngineSettings) -> None:
    """UI presets run and export the selected physical topology."""
    payload = simulate(settings)
    report = json.loads(json_export(payload))
    assert payload["energy_error"] < 2e-6
    assert report["settings"]["topology"] == settings.topology
    frame = payload["frame"]
    if settings.topology == "closed":
        assert np.ptp(frame.mass_kg) == 0
        assert np.all(frame.intake_flow_kg_s == 0)
        assert report["metrics"]["eta_tr"] is None
        if settings.closed_model == "Adiabatic gas spring":
            assert report["thermal_efficiency"] is None
    else:
        assert frame.intake_lift_m.max() > 0
        assert frame.intake_lift_m.iloc[0] == 0
        assert payload["fuel_energy_J"] == pytest.approx(20e-6 * 42.5e6)
