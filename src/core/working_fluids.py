"""Single-phase gas property adapters with consistent caloric derivatives.

Sources: fundamental relation du=cv*dT+(p-T*(dp/dT)_rho)/rho²*d(rho);
CoolProp low-level API, https://coolprop.org/coolprop/LowLevelAPI.html.
No reacting chemistry or two-phase control-volume closure is implemented.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class GasState:
    """Thermodynamic snapshot, all mass-specific quantities in SI.

    Attributes:
        pressure: p [Pa].
        temperature: T [K].
        density: rho [kg/m³].
        internal_energy: u [J/kg].
        enthalpy: h [J/kg].
        cv: (du/dT)_rho [J/(kg K)].
        cp: Isobaric heat capacity [J/(kg K)].
        dp_dT: (dp/dT)_rho [Pa/K].
        dp_drho: (dp/drho)_T [Pa m³/kg].
        du_drho: (du/drho)_T [J m³/kg²].
        entropy: s [J/(kg K)]; ideal-gas datum is arbitrary and consistent.
    """

    pressure: float
    temperature: float
    density: float
    internal_energy: float
    enthalpy: float
    cv: float
    cp: float
    dp_dT: float
    dp_drho: float
    du_drho: float
    entropy: float


class BaseWorkingFluid(ABC):
    """Property contract for a homogeneous single-phase gas."""

    @abstractmethod
    def evaluate(self, temperature: float, density: float) -> GasState:
        """Return properties at temperature [K] and density [kg/m³]."""

    @abstractmethod
    def from_pressure_temperature(
        self, pressure: float, temperature: float
    ) -> GasState:
        """Return properties at pressure [Pa] and temperature [K]."""


def _positive(*values: float) -> None:
    """Reject nonfinite or nonpositive SI state values."""
    if not np.all(np.isfinite(values)) or min(values) <= 0:
        raise ValueError("Gas state values must be finite and positive.")


@dataclass(frozen=True)
class IdealGas(BaseWorkingFluid):
    """Calorically perfect ideal gas: p=rho*R*T, u=cv*T, h=cp*T.

    Args:
        gas_constant: R [J/(kg K)], positive.
        gamma: cp/cv [-], strictly greater than one.
    """

    gas_constant: float = 287.0
    gamma: float = 1.35

    def __post_init__(self) -> None:
        """Validate thermodynamic constants [SI]."""
        _positive(self.gas_constant, self.gamma)
        if self.gamma <= 1:
            raise ValueError("Heat capacity ratio must exceed one.")

    def evaluate(self, temperature: float, density: float) -> GasState:
        """Return ideal-gas properties and analytic EOS derivatives [SI]."""
        _positive(temperature, density)
        r = self.gas_constant
        cv = r / (self.gamma - 1)
        cp = cv + r
        return GasState(
            density * r * temperature,
            temperature,
            density,
            cv * temperature,
            cp * temperature,
            cv,
            cp,
            density * r,
            r * temperature,
            0.0,
            float(cv * np.log(temperature) - r * np.log(density)),
        )

    def from_pressure_temperature(
        self, pressure: float, temperature: float
    ) -> GasState:
        """Return properties using rho=p/(RT) [kg/m³]."""
        _positive(pressure, temperature)
        return self.evaluate(
            temperature, pressure / (self.gas_constant * temperature)
        )


class CoolPropGas(BaseWorkingFluid):
    """CoolProp single-phase, nonreacting gas with variable caloric properties.

    Args:
        name: CoolProp fluid identifier, e.g. 'Air', 'Nitrogen' or 'CO2'.
        backend: CoolProp backend identifier, default 'HEOS'.

    Notes:
        Each simulator should own its adapter: AbstractState is mutable and
        must not be shared between concurrently executing integrations.
        Liquid, two-phase and out-of-range states are rejected, not clipped.
    """

    def __init__(self, name: str = "Air", backend: str = "HEOS") -> None:
        """Instantiate one reusable CoolProp state; require a valid gas EOS."""
        import CoolProp
        import CoolProp.CoolProp as cp

        self.name, self.backend = name, backend
        self._cp = cp
        try:
            self._state = CoolProp.AbstractState(backend, name)
        except (ValueError, RuntimeError) as error:
            raise ValueError(f"Cannot initialize CoolProp gas {name}: {error}")

    def _snapshot(self) -> GasState:
        """Return a validated SI snapshot after a CoolProp state update."""
        state, cp = self._state, self._cp
        if state.phase() not in (
            cp.iphase_gas,
            cp.iphase_supercritical,
            cp.iphase_supercritical_gas,
        ):
            raise ValueError("The 0D gas model requires a single gas phase.")
        temperature, density, pressure = state.T(), state.rhomass(), state.p()
        if not state.Tmin() <= temperature <= state.Tmax():
            raise ValueError(f"Temperature exceeds {self.name} EOS limits.")
        cv, capacity = state.cvmass(), state.cpmass()
        dp_dt = state.first_partial_deriv(cp.iP, cp.iT, cp.iDmass)
        dp_drho = state.first_partial_deriv(cp.iP, cp.iDmass, cp.iT)
        _positive(temperature, density, pressure, cv, capacity, dp_drho)
        return GasState(
            pressure,
            temperature,
            density,
            state.umass(),
            state.hmass(),
            cv,
            capacity,
            dp_dt,
            dp_drho,
            (pressure - temperature * dp_dt) / density**2,
            state.smass(),
        )

    def evaluate(self, temperature: float, density: float) -> GasState:
        """Return T-rho properties and homogeneous EOS derivatives [SI]."""
        _positive(temperature, density)
        try:
            self._state.update(self._cp.DmassT_INPUTS, density, temperature)
            return self._snapshot()
        except (ValueError, RuntimeError) as error:
            raise ValueError(f"Invalid {self.name} gas state: {error}")

    def from_pressure_temperature(
        self, pressure: float, temperature: float
    ) -> GasState:
        """Flash pressure [Pa] and temperature [K] into a gas snapshot."""
        _positive(pressure, temperature)
        try:
            self._state.update(self._cp.PT_INPUTS, pressure, temperature)
            return self._snapshot()
        except (ValueError, RuntimeError) as error:
            raise ValueError(f"Invalid {self.name} reservoir state: {error}")
