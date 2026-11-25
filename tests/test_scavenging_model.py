"""Analytical geometric and mass-integration checks in SI units."""

import numpy as np
import pytest

from core.engine_solver import OpposedPistonEngineCycle
from core.scavenging_model import PortGeometry, ScavengingEvaluator


def test_rectangular_overlap() -> None:
    """Check covered, partial and saturated areas [m²], including Cd."""
    port = PortGeometry(0.01, 0.05, 1.0, 5.0)
    assert port.get_effective_area(2, 0.02, 0.03) == 0
    assert port.get_effective_area(2, 0.035, 0.03) == pytest.approx(
        0.65 * 0.05 * 0.005
    )
    assert port.get_effective_area(2, 0.08, 0.03) == pytest.approx(
        0.65 * 0.05 * 0.01
    )
    assert port.get_effective_area(0, 0.08, 0.03) == 0
    assert port.get_effective_area(5, 0.08, 0.03) == 0


def test_wrapped_timing() -> None:
    """A timing gate crossing TDC repeats each revolution [rad]."""
    port = PortGeometry(0.01, 0.05, np.deg2rad(350), np.deg2rad(20))
    for angle in (0, 2 * np.pi, -2 * np.pi, np.deg2rad(355)):
        assert port.get_effective_area(angle, 0.05, 0.03) > 0
    assert port.get_effective_area(np.pi, 0.05, 0.03) == 0


def test_kinematics_compatibility() -> None:
    """Geometry matches engine_solver's existing port area times Cd [m²]."""
    engine = OpposedPistonEngineCycle()
    original = engine.intake
    port = PortGeometry(original.height, original.width, 0, 2 * np.pi,
                        original.discharge_coefficient)
    for angle in np.linspace(0, 4 * np.pi, 81):
        displacement, _ = engine.piston_motion(angle)
        assert port.get_effective_area(
            angle, displacement, original.top
        ) == pytest.approx(
            original.discharge_coefficient * engine.port_area(angle, original)
        )


def test_constant_flow_metrics() -> None:
    """One revolution at 60 RPM lasts 1 s; verify all metric definitions."""
    result = ScavengingEvaluator.evaluate_cycle_metrics(
        [0, 1, 2 * np.pi], [0.002] * 3, [0.001] * 3,
        0.0015, 0.0005, 0.001, 1.2, rpm=60,
    )
    assert result == pytest.approx({
        'DR': 0.002 / 0.0012, 'eta_tr': 0.75, 'eta_sc': 0.75,
        'm_del': 0.002, 'm_exhaust': 0.001,
    })


def test_reverse_flow_and_interpolated_boundaries() -> None:
    """Integrate only the positive triangle, not clipped sample trapezoids."""
    result = ScavengingEvaluator.evaluate_cycle_metrics(
        [0, 2], [-1, 1], [1, -1], 0, 1, 1, 1,
        rpm=60 / (2 * np.pi), scavenging_interval=(0.5, 1.5),
    )
    assert result['m_del'] == pytest.approx(0.125)
    assert result['m_exhaust'] == pytest.approx(0.125)


def test_zero_denominators() -> None:
    """Undefined efficiencies are NaN, while zero delivery has DR=0."""
    result = ScavengingEvaluator.evaluate_cycle_metrics(
        [0, 1], [0, 0], [0, 0], 0, 0, 1, 1, rpm=60
    )
    assert result['DR'] == 0
    assert np.isnan(result['eta_tr'])
    assert np.isnan(result['eta_sc'])


@pytest.mark.parametrize('angles, flows, speed, interval', [
    ([0, 0], [1, 1], 60, None),
    ([0, 1], [1], 60, None),
    ([0, np.nan], [1, 1], 60, None),
    ([0, 1], [1, 1], 0, None),
    ([0, 1], [1, 1], 60, (-1, 1)),
])
def test_invalid_metrics(angles: list[float], flows: list[float],
                         speed: float,
                         interval: tuple[float, float] | None) -> None:
    """Reject invalid sample coverage, shapes, speed and nonfinite data."""
    with pytest.raises(ValueError):
        ScavengingEvaluator.evaluate_cycle_metrics(
            angles, flows, [1, 1], 0, 0, 1, 1, rpm=speed,
            scavenging_interval=interval,
        )


@pytest.mark.parametrize('height, opening, closing, cd', [
    (0, 0, 1, 0.65), (0.01, 0, 0, 0.65),
    (0.01, 0, 7, 0.65), (0.01, 0, 1, 1.1),
])
def test_invalid_geometry(height: float, opening: float, closing: float,
                          cd: float) -> None:
    """Reject nonphysical port geometry and ambiguous timing [SI]."""
    with pytest.raises(ValueError):
        PortGeometry(height, 0.05, opening, closing, cd)
