"""Physical invariants for the SI opposed-piston single-zone model."""

import numpy as np
import pytest
from scipy.integrate import quad

from core.engine_solver import OpposedPistonEngineCycle, Port


def test_geometry() -> None:
    """Check combined displacement and analytic derivative [m³/rad]."""
    engine = OpposedPistonEngineCycle()
    v_min, _ = engine.volume(0)
    v_max, _ = engine.volume(np.pi)
    assert v_max / v_min == pytest.approx(engine.compression_ratio)
    for angle in np.linspace(0.1, 6.1, 21):
        step = 1e-6
        numerical = (engine.volume(angle + step)[0]
                     - engine.volume(angle - step)[0]) / (2 * step)
        assert engine.volume(angle)[1] == pytest.approx(numerical, rel=1e-7)


def test_wiebe_energy() -> None:
    """Normalized Wiebe integral equals fuel energy per firing [J]."""
    engine = OpposedPistonEngineCycle()
    integral, _ = quad(engine.combustion_heat_release, engine.theta_s,
                       engine.theta_s + engine.delta_theta_b)
    assert integral == pytest.approx(engine.heat_per_event)
    assert engine.combustion_heat_release(2.0) == 0


def test_nozzle_direction_and_choking() -> None:
    """Check pressure reversal and downstream-independent sonic flow [kg/s]."""
    engine = OpposedPistonEngineCycle()
    port = engine.intake
    assert engine.mass_flow(0, 1e5, 600, port) == 0
    assert engine.mass_flow(np.pi, port.pressure, 600, port) == 0
    flow = engine.mass_flow(np.pi, 1e4, 600, port)
    assert flow > 0
    assert engine.mass_flow(np.pi, 2e4, 600, port) == pytest.approx(flow)
    assert engine.mass_flow(np.pi, 2e5, 600, port) < 0


def test_closed_adiabatic_invariants() -> None:
    """Verify pV**gamma, mass and EOS in a closed motored cylinder [SI]."""
    closed = Port(0.09, 0.01, 0.01, 1e5, 300)
    engine = OpposedPistonEngineCycle(
        fuel_mass=0, heat_transfer_multiplier=0, intake=closed, exhaust=closed
    )
    result = engine.solve_cycle(samples=181)
    invariant = result.pressure * result.volume**engine.gamma
    np.testing.assert_allclose(invariant, invariant[0], rtol=3e-7)
    np.testing.assert_allclose(result.mass, result.mass[0], rtol=1e-12)
    assert result.eos_relative_error < 3e-7


def test_fired_cycle() -> None:
    """Integrate full 720-degree fired/open trajectory and verify EOS."""
    result = OpposedPistonEngineCycle().solve_cycle(samples=361)
    assert result.pv.shape == (361, 2)
    assert result.theta[-1] == pytest.approx(4 * np.pi)
    assert result.eos_relative_error < 2e-6
    assert np.ptp(result.mass) > 0
    assert np.all(result.temperature > 0)


def test_reject_inconsistent_state() -> None:
    """Reject initial mass inconsistent with ideal-gas closure [kg]."""
    with pytest.raises(ValueError, match='pV=mRT'):
        OpposedPistonEngineCycle().solve_cycle([1e5, 300, 1])
