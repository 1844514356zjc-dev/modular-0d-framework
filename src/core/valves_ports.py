"""Composable gas boundaries: cam valves, masked ports and signed nozzle flow.

Geometry: Heywood, Internal Combustion Engine Fundamentals (1988).
Flow: https://www.grc.nasa.gov/www/k-12/airplane/mflchk.html.
Areas are geometric until Cd is applied once in get_effective_area.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Callable

import numpy as np

from core.working_fluids import GasState


@dataclass(frozen=True)
class Reservoir:
    """Stationary, infinite plenum: pressure [Pa] and temperature [K]."""

    pressure: float = 1e5
    temperature: float = 320.0

    def __post_init__(self) -> None:
        """Reject invalid SI stagnation states."""
        if (
            not np.all(np.isfinite((self.pressure, self.temperature)))
            or min(self.pressure, self.temperature) <= 0
        ):
            raise ValueError(
                "Reservoir pressure and temperature must be positive."
            )


@dataclass(frozen=True, kw_only=True)
class GasExchangeBoundary(ABC):
    """Area/pressure boundary contract; all angles in radians.

    Args:
        name: Unique flow channel identifier, e.g. 'intake' or 'exhaust'.
        reservoir: Infinite plenum supplying upstream stagnation conditions.
        discharge_coeff: Cd [-], in (0, 1].
    """

    name: str
    reservoir: Reservoir
    discharge_coeff: float = 0.65

    def __post_init__(self) -> None:
        """Validate channel identity and discharge coefficient [-]."""
        if (
            not self.name
            or not np.isfinite(self.discharge_coeff)
            or not 0 < self.discharge_coeff <= 1
        ):
            raise ValueError("Boundary needs a name and Cd in (0, 1].")

    @abstractmethod
    def get_geometric_area(self, theta: float) -> float:
        """Return uncovered area [m²] at theta [rad], without Cd."""

    def get_effective_area(self, theta: float) -> float:
        """Return Cd*A_geo [m²], applying discharge correction once."""
        area = self.get_geometric_area(theta)
        if not np.isfinite(area) or area < 0:
            raise ValueError("Boundary returned an invalid geometric area.")
        return self.discharge_coeff * area


@dataclass(frozen=True, kw_only=True)
class PistonMaskedPort(GasExchangeBoundary):
    """Rectangular port uncovered by an injected piston displacement law.

    Args:
        height: Port axial height [m].
        width_total: Total circumferential width [m].
        edge_position: Upper edge distance from the piston TDC [m].
        piston_displacement: Callable(theta [rad]) -> displacement [m].

    Notes:
        The callback makes the port independent of any particular mechanism.
        A_geo=width_total*clip(x-edge_position, 0, height).
    """

    height: float
    width_total: float
    edge_position: float
    piston_displacement: Callable[[float], float]

    def __post_init__(self) -> None:
        """Validate SI window dimensions and the kinematic callback."""
        super().__post_init__()
        values = (self.height, self.width_total, self.edge_position)
        if (
            not np.all(np.isfinite(values))
            or min(values[:2]) <= 0
            or self.edge_position < 0
            or not callable(self.piston_displacement)
        ):
            raise ValueError(
                "Invalid masked-port geometry or displacement law."
            )

    def get_geometric_area(self, theta: float) -> float:
        """Return rectangular overlap area [m²] for theta [rad]."""
        displacement = self.piston_displacement(theta)
        if not np.isfinite(displacement):
            raise ValueError("Piston displacement must be finite.")
        return float(
            self.width_total
            * np.clip(displacement - self.edge_position, 0, self.height)
        )


@dataclass(frozen=True, kw_only=True)
class PoppetValve(GasExchangeBoundary):
    """Cam-driven valve: min(pi*diameter*lift, pi*diameter²/4).

    Args:
        diameter: Valve head/throat diameter [m].
        max_lift: Maximum cam lift [m].
        opening_theta: Opening crank angle [rad].
        closing_theta: Closing crank angle [rad]; may wrap across period.
        period: Crank period [rad], default 4*pi for a four-stroke engine.
        cam_lift: Optional custom callable(theta [rad]) -> lift [m].
            Lift is gated by timing and limited to [0, max_lift].

    Notes:
        Default lift=max_lift*sin(pi*phase/duration)². This smooth kinematic
        cam prescription has zero lift and velocity at both timing edges.
        Valve inertia, seat leakage and variable Cd are not modeled.
    """

    diameter: float = 0.03
    max_lift: float = 0.008
    opening_theta: float = 0.0
    closing_theta: float = np.pi
    period: float = 4 * np.pi
    cam_lift: Callable[[float], float] | None = None

    def __post_init__(self) -> None:
        """Validate a nonempty cam interval and SI valve dimensions."""
        super().__post_init__()
        values = (
            self.diameter,
            self.max_lift,
            self.opening_theta,
            self.closing_theta,
            self.period,
        )
        span = self.closing_theta - self.opening_theta
        if (
            not np.all(np.isfinite(values))
            or self.diameter <= 0
            or self.max_lift <= 0
            or self.period <= 0
            or span == 0
            or abs(span) >= self.period
            or (self.cam_lift is not None and not callable(self.cam_lift))
        ):
            raise ValueError("Invalid poppet-valve geometry or cam timing.")

    def get_lift(self, theta: float) -> float:
        """Return gated cam lift [m], periodic in crank angle [rad]."""
        if not np.isfinite(theta):
            raise ValueError("Valve angle must be finite.")
        duration = (self.closing_theta - self.opening_theta) % self.period
        phase = (theta - self.opening_theta) % self.period
        if phase >= duration:
            return 0.0
        lift = (
            self.max_lift * np.sin(np.pi * phase / duration) ** 2
            if self.cam_lift is None
            else self.cam_lift(theta)
        )
        if not np.isfinite(lift):
            raise ValueError("Cam lift must be finite.")
        return float(np.clip(lift, 0, self.max_lift))

    def get_geometric_area(self, theta: float) -> float:
        """Return the smaller curtain or throat area [m²] (Heywood)."""
        return float(
            min(
                np.pi * self.diameter * self.get_lift(theta),
                np.pi * self.diameter**2 / 4,
            )
        )


def nozzle_mass_flow(
    effective_area: float, cylinder: GasState, reservoir: GasState
) -> float:
    """Return inward-positive quasi-steady nozzle flow [kg/s].

    Args:
        effective_area: Cd-corrected opening [m²]; Cd is not reapplied.
        cylinder: Cylinder SI thermodynamic snapshot.
        reservoir: External SI stagnation snapshot of the same working fluid.

    Returns:
        Signed ideal-gas isentropic nozzle rate, including sonic choking and
        reversal. For CoolProp gases this is explicitly a locally idealized
        nozzle approximation using upstream cp/cv and R_eff=p/(rho*T).
        Cylinder energy/EOS still use full real-gas caloric derivatives.
    """
    if not np.isfinite(effective_area) or effective_area < 0:
        raise ValueError("Effective nozzle area must be nonnegative.")
    if effective_area == 0 or cylinder.pressure == reservoir.pressure:
        return 0.0
    inward = reservoir.pressure > cylinder.pressure
    upstream, downstream = (
        (reservoir, cylinder) if inward else (cylinder, reservoir)
    )
    gamma = upstream.cp / upstream.cv
    if gamma <= 1 or not np.isfinite(gamma):
        raise ValueError("Nozzle approximation requires cp/cv > 1.")
    critical = (2 / (gamma + 1)) ** (gamma / (gamma - 1))
    ratio = max(downstream.pressure / upstream.pressure, critical)
    factor = (
        2
        * gamma
        / (gamma - 1)
        * (ratio ** (2 / gamma) - ratio ** ((gamma + 1) / gamma))
    )
    r_effective = upstream.pressure / (upstream.density * upstream.temperature)
    magnitude = (
        effective_area
        * upstream.pressure
        / np.sqrt(r_effective * upstream.temperature)
        * np.sqrt(max(factor, 0.0))
    )
    return float(magnitude if inward else -magnitude)
