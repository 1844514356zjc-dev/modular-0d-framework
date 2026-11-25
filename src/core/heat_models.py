"""Injectable heat sources [J/rad] and wall losses [W].

Sources: normalized Wiebe (Heywood, 1988); simplified Woschni, SAE 670931.
Heat sources exclude injection mass and chemistry. Angle units are radians.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np


class HeatReleaseModel(Protocol):
    """Protocol for any angle-prescribed heat input or removal law."""

    def heat_release(self, theta: float) -> float:
        """Return dQ/dtheta [J/rad], positive into the control volume."""
        ...


@dataclass(frozen=True)
class NoHeatRelease:
    """Adiabatic/motored heat source: zero J/rad at every angle."""

    def heat_release(self, theta: float) -> float:
        """Return zero heat release [J/rad] for any theta [rad]."""
        return 0.0


@dataclass(frozen=True)
class WiebeCombustion:
    """Normalized periodic Wiebe law, integral equal to heat_per_event.

    Args:
        heat_per_event: Released heat per firing [J], nonnegative.
        theta_s: Burn start [rad], may precede zero.
        duration: Burn duration [rad], in (0, period).
        m: Shape exponent [-], nonnegative.
        a: Completeness coefficient [-], positive.
        period: Firing interval [rad], e.g. 2*pi or 4*pi.
    """

    heat_per_event: float = 833.0
    theta_s: float = np.deg2rad(-10)
    duration: float = np.deg2rad(50)
    m: float = 2.0
    a: float = 6.908
    period: float = 2 * np.pi

    def __post_init__(self) -> None:
        """Validate finite heat, exponent and angular support [SI]."""
        values = (
            self.heat_per_event,
            self.theta_s,
            self.duration,
            self.m,
            self.a,
            self.period,
        )
        if (
            not np.all(np.isfinite(values))
            or self.heat_per_event < 0
            or self.m < 0
            or self.a <= 0
            or not 0 < self.duration < self.period
        ):
            raise ValueError("Invalid Wiebe heat-release configuration.")

    def heat_release(self, theta: float) -> float:
        """Return normalized Q*a*(m+1)/duration*z^m*exp(-a*z^(m+1))."""
        phase = (theta - self.theta_s) % self.period
        if phase >= self.duration:
            return 0.0
        z = phase / self.duration
        return float(
            self.heat_per_event
            * self.a
            * (self.m + 1)
            / self.duration
            * z**self.m
            * np.exp(-self.a * z ** (self.m + 1))
            / (-np.expm1(-self.a))
        )


@dataclass(frozen=True)
class OttoHeatRelease:
    """Isochoric heat pulses on OttoCycleMechanism's constant-volume dwells.

    Args:
        heat_input: Heat addition [J] at pi..2*pi, positive.
        compression_ratio: Vmax/Vmin [-], >1.
        gamma: Constant ideal-gas cp/cv [-], >1.

    Notes:
        Heat rejection is Q_in/r^(gamma-1), giving eta=1-r^(1-gamma).
        This analytic benchmark requires the matching calorically perfect
        IdealGas and OttoCycleMechanism, with wall heat transfer disabled.
    """

    heat_input: float = 500.0
    compression_ratio: float = 10.0
    gamma: float = 1.35

    def __post_init__(self) -> None:
        """Validate finite benchmark heat [J] and ideal-gas constants [-]."""
        values = (self.heat_input, self.compression_ratio, self.gamma)
        if (
            not np.all(np.isfinite(values))
            or self.heat_input <= 0
            or self.compression_ratio <= 1
            or self.gamma <= 1
        ):
            raise ValueError("Invalid air-standard Otto heat configuration.")

    def heat_release(self, theta: float) -> float:
        """Return smooth isochoric addition/removal [J/rad], period 4*pi."""
        angle = theta % (4 * np.pi)
        if np.pi <= angle < 2 * np.pi:
            return float(2 * self.heat_input / np.pi * np.sin(angle) ** 2)
        if 3 * np.pi <= angle < 4 * np.pi:
            rejection = self.heat_input / self.compression_ratio ** (
                self.gamma - 1
            )
            return float(-2 * rejection / np.pi * np.sin(angle) ** 2)
        return 0.0


@dataclass(frozen=True)
class WoschniHeatTransfer:
    """Simplified Woschni convective wall loss for piston engines.

    Args:
        bore: Cylinder diameter [m].
        stroke: Single-piston travel [m].
        rpm: Rotational speed [rev/min].
        wall_temperature: Fixed wall temperature [K].
        multiplier: Empirical calibration factor [-], nonnegative.
    """

    bore: float
    stroke: float
    rpm: float
    wall_temperature: float = 450.0
    multiplier: float = 1.0

    def __post_init__(self) -> None:
        """Validate finite geometric and thermal inputs [SI]."""
        values = (
            self.bore,
            self.stroke,
            self.rpm,
            self.wall_temperature,
            self.multiplier,
        )
        if (
            not np.all(np.isfinite(values))
            or min(values[:4]) <= 0
            or self.multiplier < 0
        ):
            raise ValueError("Invalid Woschni heat-transfer configuration.")

    def __call__(
        self, theta: float, pressure: float, temperature: float, volume: float
    ) -> float:
        """Return h*Aw*(T-Twall) [W], positive outward.

        h=3.26*B^(-0.2)*p_kPa^0.8*T^(-0.55)*w^0.8 [W/(m² K)],
        w=2.28*2*S*rpm/60 [m/s], Aw=2*Ap+4*V/B [m²].
        """
        speed = 2.28 * 2 * self.stroke * self.rpm / 60
        coefficient = (
            self.multiplier
            * 3.26
            * self.bore ** (-0.2)
            * (pressure / 1000) ** 0.8
            * temperature ** (-0.55)
            * speed**0.8
        )
        area = np.pi * self.bore**2 / 2 + 4 * volume / self.bore
        return float(
            coefficient * area * (temperature - self.wall_temperature)
        )
