"""Volume laws in SI and radians; slider-crank geometry follows Heywood.

Reference: Internal Combustion Engine Fundamentals (1988), engine geometry.
Positive displacement points from TDC toward BDC. Stroke is per piston.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import ClassVar

import numpy as np


class BaseMechanism(ABC):
    """Minimum contract for any prescribed moving control volume."""

    @abstractmethod
    def get_volume(self, theta: float) -> float:
        """Return positive control-volume size [m³] at angle theta [rad]."""

    @abstractmethod
    def get_dV_dtheta(self, theta: float) -> float:
        """Return the analytic volume derivative [m³/rad] at theta [rad]."""


@dataclass(frozen=True)
class CrankMechanism(BaseMechanism):
    """Single-piston slider-crank geometry, minimum volume at theta=0.

    Args:
        bore: Cylinder diameter [m].
        stroke: Travel of one piston [m].
        con_rod: Rod length [m], strictly greater than stroke/2.
        compression_ratio: Maximum/minimum cylinder volume [-], >1.
    """

    bore: float = 0.08
    stroke: float = 0.08
    con_rod: float = 0.14
    compression_ratio: float = 16.0
    piston_count: ClassVar[int] = 1

    def __post_init__(self) -> None:
        """Validate finite SI dimensions and a nonsingular crank linkage."""
        values = (self.bore, self.stroke, self.con_rod, self.compression_ratio)
        if (
            not np.all(np.isfinite(values))
            or min(values[:3]) <= 0
            or self.con_rod <= self.stroke / 2
            or self.compression_ratio <= 1
        ):
            raise ValueError(
                "Invalid slider-crank geometry or compression ratio."
            )

    @property
    def piston_area(self) -> float:
        """Return circular crown area pi*B²/4 [m²]."""
        return float(np.pi * self.bore**2 / 4)

    @property
    def swept_volume(self) -> float:
        """Return total swept volume n*Ap*S [m³]."""
        return self.piston_count * self.piston_area * self.stroke

    @property
    def clearance_volume(self) -> float:
        """Return minimum volume Vd/(compression_ratio-1) [m³]."""
        return self.swept_volume / (self.compression_ratio - 1)

    def piston_motion(self, theta: float) -> tuple[float, float]:
        """Return one-piston x [m] and dx/dtheta [m/rad] (Heywood).

        x=r*(1-cos(theta))+l-sqrt(l²-r²*sin²(theta)). Theta is [rad].
        """
        if not np.isfinite(theta):
            raise ValueError("Crank angle must be finite.")
        r = self.stroke / 2
        sine, cosine = np.sin(theta), np.cos(theta)
        root = np.sqrt(self.con_rod**2 - r**2 * sine**2)
        return (
            float(r * (1 - cosine) + self.con_rod - root),
            float(r * sine + r**2 * sine * cosine / root),
        )

    def get_piston_displacement(self, theta: float) -> float:
        """Return crown displacement from TDC [m] at angle theta [rad]."""
        return self.piston_motion(theta)[0]

    def get_volume(self, theta: float) -> float:
        """Return Vc+n*Ap*x [m³], with theta in radians."""
        return (
            self.clearance_volume
            + self.piston_count
            * self.piston_area
            * self.get_piston_displacement(theta)
        )

    def get_dV_dtheta(self, theta: float) -> float:
        """Return n*Ap*dx/dtheta [m³/rad]."""
        return (
            self.piston_count * self.piston_area * self.piston_motion(theta)[1]
        )


class OpposedPistonMechanism(CrankMechanism):
    """Two identical synchronous opposed pistons; total stroke is 2*S.

    Constructor parameters and units match CrankMechanism. No crank phasing
    or proprietary non-crank single-stroke motion law is assumed.
    """

    piston_count: ClassVar[int] = 2


@dataclass(frozen=True)
class OttoCycleMechanism(BaseMechanism):
    """Four-segment volume schedule for an exact air-standard Otto benchmark.

    Args:
        crank: Single-piston geometry supplying Vmin and Vmax [m³].

    Notes:
        0..pi: compression, pi..2pi: constant Vmin, 2pi..3pi: expansion,
        3pi..4pi: constant Vmax. Theta [rad] is a process coordinate, not a
        physical crankshaft during the constant-volume dwells (Heywood).
    """

    crank: CrankMechanism

    def get_volume(self, theta: float) -> float:
        """Return the piecewise isochoric/slider-crank volume [m³]."""
        angle = theta % (4 * np.pi)
        if angle < np.pi:
            return self.crank.get_volume(np.pi + angle)
        if angle < 2 * np.pi:
            return self.crank.clearance_volume
        if angle < 3 * np.pi:
            return self.crank.get_volume(angle - 2 * np.pi)
        return self.crank.clearance_volume + self.crank.swept_volume

    def get_dV_dtheta(self, theta: float) -> float:
        """Return dV/dtheta [m³/rad], zero on the isochoric segments."""
        angle = theta % (4 * np.pi)
        if angle < np.pi:
            return self.crank.get_dV_dtheta(np.pi + angle)
        if 2 * np.pi <= angle < 3 * np.pi:
            return self.crank.get_dV_dtheta(angle - 2 * np.pi)
        return 0.0
