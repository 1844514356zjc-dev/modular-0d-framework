"""Piston-controlled rectangular ports and cycle scavenging metrics.

All angles are radians; all other dimensional quantities use SI. Definitions
follow MIT 2.61, Lecture 8, Gas Exchange Processes:
https://web.mit.edu/2.61/www/Lecture%20notes/Lec.%2008%20Gas%20exchange%20processes.pdf
Port area follows rectangular overlap geometry, not an empirical flow law.
No species transport is solved here: retained fresh and residual masses must
come from a species model or measurements, not total cylinder mass alone.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.integrate import trapezoid


@dataclass(frozen=True)
class PortGeometry:
    """Rectangular window with a repeating angular timing gate.

    Args:
        height: Axial port height [m], positive.
        width_total: Sum of circumferential window widths [m], positive.
        opening_theta: Opening angle [rad], measured from piston TDC.
        closing_theta: Closing angle [rad]. May be below opening_theta for
            a gate crossing TDC, e.g. 350 to 20 degrees. The gate repeats
            every 2*pi; a difference of exactly 2*pi denotes always open.
        discharge_coeff: Area correction Cd [-], in (0, 1].

    Notes:
        Coordinates increase from TDC toward BDC for each piston separately.
        Both piston displacement and port edge use that same coordinate.
        Unequal piston phasing belongs in the supplied piston kinematics.
        Equal opening/closing angles are ambiguous and rejected.
    """

    height: float
    width_total: float
    opening_theta: float
    closing_theta: float
    discharge_coeff: float = 0.65

    def __post_init__(self) -> None:
        """Validate SI geometry, finite timing and dimensionless Cd."""
        values = (self.height, self.width_total, self.opening_theta,
                  self.closing_theta, self.discharge_coeff)
        if not np.all(np.isfinite(values)):
            raise ValueError("Port parameters must be finite.")
        if self.height <= 0 or self.width_total <= 0:
            raise ValueError("Port height and width must be positive.")
        if not 0 < self.discharge_coeff <= 1:
            raise ValueError("discharge_coeff must be in (0, 1].")
        span = self.closing_theta - self.opening_theta
        if span == 0 or abs(span) > 2 * np.pi:
            raise ValueError("Timing span must be nonzero and at most 2*pi.")

    def get_effective_area(
        self, theta: float, piston_displacement: float, edge_position: float,
    ) -> float:
        """Compute discharge-corrected uncovered area [m²].

        Rectangular overlap gives h_open=clip(x-edge, 0, height),
        A_geo=width_total*h_open and A_eff=Cd*A_geo inside the timing gate.
        This method returns A_eff, including Cd as requested. Do not apply
        another discharge coefficient in a downstream nozzle calculation.

        Args:
            theta: Current crank angle [rad], may be unwrapped.
            piston_displacement: Piston crown distance from its TDC [m].
            edge_position: Distance from TDC to port's upper edge [m].

        Returns:
            Effective area [m²], zero outside the open/close timing gate.
            Gate convention is opening inclusive, closing exclusive.

        Raises:
            ValueError: Nonfinite coordinates or negative distances.
        """
        values = (theta, piston_displacement, edge_position)
        if not np.all(np.isfinite(values)):
            raise ValueError("Angle and coordinates must be finite.")
        if piston_displacement < 0 or edge_position < 0:
            raise ValueError("Distances from TDC must be nonnegative.")
        span = self.closing_theta - self.opening_theta
        if abs(span) != 2 * np.pi:
            duration = span % (2 * np.pi)
            phase = (theta - self.opening_theta) % (2 * np.pi)
            if phase >= duration:
                return 0.0
        axial_opening = np.clip(
            piston_displacement - edge_position, 0.0, self.height
        )
        return float(self.discharge_coeff * self.width_total * axial_opening)


class ScavengingEvaluator:
    """Evaluate supplied species inventories using MIT 2.61 definitions.

    Intake flow is positive INTO the cylinder; exhaust flow is positive OUT
    of the cylinder. Negative values represent reverse flow. Delivered mass
    counts gross positive fresh intake flow, not net flow. Exhaust flow alone
    cannot identify fresh-air short circuiting. The intake plenum is assumed
    to contain only fresh air; exhaust backflow is not counted as delivery.
    """

    @staticmethod
    def evaluate_cycle_metrics(
        theta: Sequence[float], m_dot_in: Sequence[float],
        m_dot_out: Sequence[float], m_trap: float, m_res: float,
        V_d: float, rho_ref: float, *, rpm: float,
        scavenging_interval: tuple[float, float] | None = None,
    ) -> dict[str, float]:
        """Integrate flow in time and calculate dimensionless efficiencies.

        Formulas (MIT 2.61 Lecture 8):
        m_del=integral(max(m_dot_in, 0)*dtheta/omega), omega=2*pi*rpm/60;
        DR=m_del/(rho_ref*V_d), eta_tr=m_trap/m_del,
        eta_sc=m_trap/(m_trap+m_res). Piecewise-linear flow interpolation
        with inserted zero crossings makes positive-flow integration exact
        for the supplied interpolant, including partial boundary intervals.

        Args:
            theta: Strictly increasing unwrapped angles [rad], length >=2.
                Supply one cycle or one sweep matching the final inventories.
                A 720-degree two-stroke trace contains two sweeps: select the
                relevant sweep rather than mixing two deliveries with one
                sweep's retained mass. Unwrap angles before calling if needed.
            m_dot_in: Signed fresh intake flow samples [kg/s].
            m_dot_out: Signed exhaust flow samples [kg/s], positive outward.
                Negate engine_solver.mass_flow for its exhaust Port, since
                that solver uses positive inward for both ports.
            m_trap: Retained fresh air from this delivery at sweep end [kg].
            m_res: Residual burned gas at the same instant [kg].
            V_d: Combined cylinder swept volume [m³].
            rho_ref: Reference fresh intake density [kg/m³].
            rpm: Constant crankshaft speed [rev/min], positive; required for
                dimensional conversion of kg/s versus crank angle.
            scavenging_interval: Optional (start, end) unwrapped angles [rad]
                within theta coverage. Defaults to the entire input trace.
                Boundaries are interpolated, never silently clipped.

        Returns:
            Keys: DR, eta_tr, eta_sc [-], m_del [kg], m_exhaust [kg].
            m_exhaust is gross positive exhaust flow integrated over the same
            interval and is diagnostic only. Zero-denominator efficiencies
            are NaN; zero delivery gives DR=0. Ratios are never clipped;
            eta_tr>1 indicates mismatched inventories or cycle boundaries.

        Raises:
            ValueError: Invalid shapes, nonfinite data, invalid masses, speed,
                reference state, angle ordering or interval coverage.
        """
        angles = np.asarray(theta, dtype=float)
        intake = np.asarray(m_dot_in, dtype=float)
        exhaust = np.asarray(m_dot_out, dtype=float)
        if (angles.ndim != 1 or angles.size < 2
                or intake.shape != angles.shape
                or exhaust.shape != angles.shape):
            raise ValueError("Angles and flows require matching 1D shapes.")
        if not all(np.all(np.isfinite(a)) for a in (angles, intake, exhaust)):
            raise ValueError("Angles and flows must be finite.")
        if np.any(np.diff(angles) <= 0):
            raise ValueError("Angles must increase strictly without wrapping.")
        if not np.all(np.isfinite((m_trap, m_res, V_d, rho_ref, rpm))):
            raise ValueError("Masses, reference state and RPM must be finite.")
        if min(m_trap, m_res) < 0 or min(V_d, rho_ref, rpm) <= 0:
            raise ValueError("Invalid mass, reference state or RPM.")
        bounds = np.asarray(
            (angles[0], angles[-1]) if scavenging_interval is None
            else scavenging_interval, dtype=float,
        )
        if (bounds.shape != (2,) or not np.all(np.isfinite(bounds))
                or bounds[0] >= bounds[1] or bounds[0] < angles[0]
                or bounds[1] > angles[-1]):
            raise ValueError("Interval must lie inside the sampled range.")
        start, end = bounds
        selected = np.concatenate((
            [start], angles[(angles > start) & (angles < end)], [end]
        ))
        omega = 2 * np.pi * rpm / 60
        totals: list[float] = []
        for flow in (intake, exhaust):
            values = np.interp(selected, angles, flow)
            left, right = values[:-1], values[1:]
            crosses = ((left < 0) & (right > 0)
                       | (left > 0) & (right < 0))
            roots = (selected[:-1][crosses]
                     - left[crosses] * np.diff(selected)[crosses]
                     / (right[crosses] - left[crosses]))
            grid = np.sort(np.concatenate((selected, roots)))
            positive = np.maximum(np.interp(grid, angles, flow), 0.0)
            totals.append(float(trapezoid(positive, grid) / omega))
        delivered, expelled = totals
        total_trapped = m_trap + m_res
        return {
            "DR": delivered / (rho_ref * V_d),
            "eta_tr": m_trap / delivered if delivered > 0 else float("nan"),
            "eta_sc": (m_trap / total_trapped if total_trapped > 0
                       else float("nan")),
            "m_del": delivered,
            "m_exhaust": expelled,
        }
