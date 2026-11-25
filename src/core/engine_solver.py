"""Composable homogeneous-gas thermodynamics, SI units and radians.

References:
    Heywood, Internal Combustion Engine Fundamentals (1988): slider-crank
    geometry, normalized Wiebe law, and open-system first-law balances.
    Woschni (1967), SAE 670931: empirical convection correlation; simplified
    here to piston-speed-driven velocity (no fitted combustion velocity).
    NASA Glenn, Mass Flow Choking:
    https://www.grc.nasa.gov/www/k-12/airplane/mflchk.html

The legacy OpposedPistonEngineCycle adapter retains its original assumptions:
synchronous identical pistons, constant ideal-gas R and gamma,
perfect mixing, infinite stationary intake/exhaust plenums, negligible kinetic
energy in cylinder, no leakage or species transport. Fuel is represented as
heat added to an already present charge; injection mass is neglected. AFR is
a prescribed charge property, not a computed scavenging/species balance.
A 720-degree run is a transient, not necessarily a periodic steady cycle.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from core.heat_models import (
    HeatReleaseModel,
    NoHeatRelease,
    OttoHeatRelease,
    WiebeCombustion,
    WoschniHeatTransfer,
)
from core.mechanics import (
    BaseMechanism,
    CrankMechanism,
    OpposedPistonMechanism,
    OttoCycleMechanism,
)
from core.valves_ports import (
    GasExchangeBoundary,
    PistonMaskedPort,
    PoppetValve,
    Reservoir,
    nozzle_mass_flow,
)
from core.working_fluids import (
    BaseWorkingFluid,
    CoolPropGas,
    GasState,
    IdealGas,
)

Array = NDArray[np.float64]


@dataclass(frozen=True)
class Port:
    """Rectangular piston-controlled port and stationary ideal-gas plenum.

    Args:
        top: Distance from piston TDC position to port top [m].
        height: Axial port height [m].
        width: Total circumferential width of all windows [m].
        pressure: Plenum stagnation pressure [Pa].
        temperature: Plenum stagnation temperature [K].
        discharge_coefficient: Effective nozzle coefficient [dimensionless].
    """

    top: float
    height: float
    width: float
    pressure: float
    temperature: float
    discharge_coefficient: float = 0.7

    def __post_init__(self) -> None:
        """Reject nonfinite geometry and invalid reservoir values [SI]."""
        values = (
            self.top,
            self.height,
            self.width,
            self.pressure,
            self.temperature,
            self.discharge_coefficient,
        )
        if not np.all(np.isfinite(values)):
            raise ValueError("Port parameters must be finite.")
        if self.top < 0 or min(values[1:]) <= 0:
            raise ValueError("Port dimensions and states must be positive.")
        if self.discharge_coefficient > 1:
            raise ValueError("Discharge coefficient must be in (0, 1].")


@dataclass(frozen=True)
class CycleResult:
    """Sampled trajectory: time [s], angle [rad], V [m³], p [Pa], T [K].

    mass is [kg]; pv is an (N, 2) array with columns [V, p]. The reported
    EOS residual compares integrated pressure with the injected fluid EOS.
    boundary_flows are inward-positive [kg/s]. internal_energy is total U [J];
    heat, work and enthalpy_exchange are cumulative Q/W/H [J]. Work is positive
    outward, heat/enthalpy inward. energy_relative_error audits Delta U=Q-W+H.
    """

    time: Array
    theta: Array
    volume: Array
    pressure: Array
    temperature: Array
    mass: Array
    pv: Array
    eos_relative_error: float
    boundary_flows: dict[str, Array] = field(default_factory=dict)
    internal_energy: Array | None = None
    heat: Array | None = None
    work: Array | None = None
    enthalpy_exchange: Array | None = None
    energy_relative_error: float = 0.0


@dataclass(frozen=True)
class EngineConfig:
    """Composition root for a prescribed-volume homogeneous gas simulator.

    Args:
        mechanism: Volume law [m³] and angular derivative [m³/rad].
        boundaries: Any combination of ports and valves; empty for closed gas.
        fluid: IdealGas or CoolPropGas (one adapter per concurrent simulation).
        combustion: Heat source dQ/dtheta [J/rad]; may include cooling.
        rpm: Constant process-coordinate speed [rev/min], positive.
        theta_start: Initial process angle [rad].
        theta_end: Final angle [rad], greater than theta_start.
        initial_state: Optional [p Pa, T K, mass kg] satisfying the chosen EOS.
        wall_heat_transfer: Optional callable(theta,p,T,V)->wall loss [W].
    """

    mechanism: BaseMechanism
    boundaries: tuple[GasExchangeBoundary, ...] = ()
    fluid: BaseWorkingFluid = field(default_factory=IdealGas)
    combustion: HeatReleaseModel = field(default_factory=NoHeatRelease)
    rpm: float = 2000.0
    theta_start: float = 0.0
    theta_end: float = 4 * np.pi
    initial_state: tuple[float, float, float] | None = None
    wall_heat_transfer: (
        Callable[[float, float, float, float], float] | None
    ) = None

    def __post_init__(self) -> None:
        """Validate integration bounds and unique exchange channel names."""
        if (
            not np.all(
                np.isfinite((self.rpm, self.theta_start, self.theta_end))
            )
            or self.rpm <= 0
            or self.theta_end <= self.theta_start
        ):
            raise ValueError("Invalid speed or angular integration interval.")
        if not isinstance(self.mechanism, BaseMechanism):
            raise ValueError("mechanism must implement BaseMechanism.")
        if not isinstance(self.fluid, BaseWorkingFluid):
            raise ValueError("fluid must implement BaseWorkingFluid.")
        if not callable(getattr(self.combustion, "heat_release", None)):
            raise ValueError("combustion must expose heat_release(theta).")
        if self.wall_heat_transfer is not None and not callable(
            self.wall_heat_transfer
        ):
            raise ValueError("wall_heat_transfer must be callable or None.")
        if any(
            not isinstance(b, GasExchangeBoundary) for b in self.boundaries
        ):
            raise ValueError("Invalid gas-exchange component.")
        if isinstance(self.combustion, OttoHeatRelease):
            if (
                not isinstance(self.mechanism, OttoCycleMechanism)
                or not isinstance(self.fluid, IdealGas)
                or self.wall_heat_transfer is not None
            ):
                raise ValueError(
                    "Otto benchmark requires ideal gas and dwells."
                )
            if (
                self.combustion.gamma != self.fluid.gamma
                or self.combustion.compression_ratio
                != self.mechanism.crank.compression_ratio
                or self.boundaries
            ):
                raise ValueError("Inconsistent air-standard Otto components.")
        object.__setattr__(self, "boundaries", tuple(self.boundaries))
        if self.initial_state is not None:
            object.__setattr__(
                self, "initial_state", tuple(self.initial_state)
            )
        names = [boundary.name for boundary in self.boundaries]
        if len(names) != len(set(names)):
            raise ValueError("Gas-exchange boundary names must be unique.")

    @classmethod
    def from_dict(cls, specification: Mapping[str, Any]) -> EngineConfig:
        """Build composition from SI dictionaries or component objects.

        Args:
            specification: EngineConfig field mapping. Nested dictionaries use
                type='crank'/'opposed'/'otto' for mechanism, 'ideal'/'coolprop'
                for fluid, 'none'/'wiebe'/'otto' for combustion, and
                'port'/'poppet' for boundaries. See examples/solve_modular.py.

        Returns:
            Validated component configuration without modifying the input map.
        """
        values = dict(specification)
        mechanism = values.get("mechanism")
        if isinstance(mechanism, Mapping):
            options = dict(mechanism)
            kind = options.pop("type", "crank")
            factories = {
                "crank": CrankMechanism,
                "opposed": OpposedPistonMechanism,
            }
            if kind == "otto":
                mechanism = OttoCycleMechanism(CrankMechanism(**options))
            elif kind in factories:
                mechanism = factories[kind](**options)
            else:
                raise ValueError(f"Unknown mechanism type: {kind}")
            values["mechanism"] = mechanism
        for key, factories, default in (
            ("fluid", {"ideal": IdealGas, "coolprop": CoolPropGas}, "ideal"),
            (
                "combustion",
                {
                    "none": NoHeatRelease,
                    "wiebe": WiebeCombustion,
                    "otto": OttoHeatRelease,
                },
                "none",
            ),
        ):
            item = values.get(key)
            if isinstance(item, Mapping):
                options = dict(item)
                kind = options.pop("type", default)
                if kind not in factories:
                    raise ValueError(f"Unknown {key} type: {kind}")
                values[key] = factories[kind](**options)
        boundaries = []
        for item in values.get("boundaries", ()):
            if not isinstance(item, Mapping):
                boundaries.append(item)
                continue
            options = dict(item)
            kind = options.pop("type")
            plenum = options.get("reservoir")
            if isinstance(plenum, Mapping):
                options["reservoir"] = Reservoir(**plenum)
            if kind == "port":
                if "piston_displacement" not in options:
                    if not hasattr(mechanism, "get_piston_displacement"):
                        raise ValueError(
                            "Masked ports need a displacement law."
                        )
                    options["piston_displacement"] = (
                        mechanism.get_piston_displacement
                    )
                boundaries.append(PistonMaskedPort(**options))
            elif kind == "poppet":
                boundaries.append(PoppetValve(**options))
            else:
                raise ValueError(f"Unknown exchange boundary type: {kind}")
        values["boundaries"] = tuple(boundaries)
        wall = values.get("wall_heat_transfer")
        if isinstance(wall, Mapping):
            values["wall_heat_transfer"] = WoschniHeatTransfer(**wall)
        return cls(**values)


class EngineSimulator:
    """Unified open/closed control-volume ODE engine, with injected components.

    For rho=m/V and u=u(T,rho), the first law gives
    Tdot=(Qdot-p*Vdot+Hdot-u*mdot-m*u_rho*rhodot)/(m*cv),
    pdot=p_T*Tdot+p_rho*rhodot. This remains valid for a nonideal gas.
    Signed inlet/outlet enthalpies always use the actual upstream state.
    """

    def __init__(self, config: EngineConfig | Mapping[str, Any]) -> None:
        """Resolve a typed/dictionary configuration and cache plenum states."""
        self.config = (
            EngineConfig.from_dict(config)
            if isinstance(config, Mapping)
            else config
        )
        self.mechanism = self.config.mechanism
        self.fluid = self.config.fluid
        self.omega = self.config.rpm * 2 * np.pi / 60
        self._reservoirs = {
            boundary.name: self.fluid.from_pressure_temperature(
                boundary.reservoir.pressure, boundary.reservoir.temperature
            )
            for boundary in self.config.boundaries
        }

    def flow_rates(
        self, theta: float, temperature: float, mass: float
    ) -> dict[str, float]:
        """Return channel rates [kg/s], all positive inward.

        Args:
            theta: Instantaneous process/crank angle [rad].
            temperature: Cylinder temperature [K].
            mass: Cylinder mass [kg].

        Returns:
            Unique channel names mapped to signed nozzle rates.
        """
        state = self.fluid.evaluate(
            temperature, mass / self.mechanism.get_volume(theta)
        )
        return self._flows(theta, state)

    def _flows(self, theta: float, state: GasState) -> dict[str, float]:
        """Evaluate effective areas once and calculate inward flows [kg/s]."""
        return {
            boundary.name: nozzle_mass_flow(
                boundary.get_effective_area(theta),
                state,
                self._reservoirs[boundary.name],
            )
            for boundary in self.config.boundaries
        }

    def _rates(
        self, t: float, y: Sequence[float], theta: float | None = None
    ) -> Array:
        """Return physical derivatives plus Q, W and enthalpy rates [SI]."""
        # Pressure is a redundant EOS diagnostic; implicit Jacobian probes may
        # perturb its independent coordinate without changing T or density.
        if not np.all(np.isfinite(y[:3])) or min(y[1:3]) <= 0:
            raise ValueError("Nonphysical cylinder state during integration.")
        angle = (
            self.config.theta_start + self.omega * t
            if theta is None
            else theta
        )
        _, temperature, mass = y[:3]
        volume = self.mechanism.get_volume(angle)
        if not np.isfinite(volume) or volume <= 0:
            raise ValueError("Mechanism volume must be finite and positive.")
        volume_rate = self.mechanism.get_dV_dtheta(angle) * self.omega
        state = self.fluid.evaluate(temperature, mass / volume)
        flows = self._flows(angle, state)
        mass_rate = sum(flows.values())
        enthalpy_rate = sum(
            flow
            * (self._reservoirs[name].enthalpy if flow > 0 else state.enthalpy)
            for name, flow in flows.items()
        )
        density_rate = (
            mass_rate / volume - state.density * volume_rate / volume
        )
        wall_loss = (
            0.0
            if self.config.wall_heat_transfer is None
            else self.config.wall_heat_transfer(
                angle, state.pressure, temperature, volume
            )
        )
        heat_rate = (
            self.omega * self.config.combustion.heat_release(angle) - wall_loss
        )
        work_rate = state.pressure * volume_rate
        temp_rate = (
            heat_rate
            - work_rate
            + enthalpy_rate
            - state.internal_energy * mass_rate
            - mass * state.du_drho * density_rate
        ) / (mass * state.cv)
        pressure_rate = state.dp_dT * temp_rate + state.dp_drho * density_rate
        return np.array(
            [
                pressure_rate,
                temp_rate,
                mass_rate,
                heat_rate,
                work_rate,
                enthalpy_rate,
            ]
        )

    def ode_system(
        self, t: float, y: Sequence[float], theta: float | None = None
    ) -> Array:
        """Return [pdot, Tdot, mdot] [Pa/s, K/s, kg/s] at time t [s].

        Args:
            t: Time since config.theta_start [s].
            y: Pressure [Pa], temperature [K], mass [kg].
            theta: Optional instantaneous angle [rad]; default advances with t.

        Returns:
            Derivatives from the open-system first law and chosen EOS.
        """
        if len(y) != 3 or not np.all(np.isfinite(y)) or min(y) <= 0:
            raise ValueError("Physical state requires positive p, T and mass.")
        return self._rates(t, y, theta)[:3]

    def solve_cycle(
        self,
        initial_state: Sequence[float] | None = None,
        samples: int = 1441,
        method: str = "Radau",
        rtol: float = 1e-7,
        atol: Sequence[float] = (1e-2, 1e-6, 1e-12),
        max_step_degrees: float = 1.0,
    ) -> CycleResult:
        """Integrate the configured interval and audit first-law energy.

        Args:
            initial_state: Optional [p Pa, T K, mass kg], overriding config.
                Default: 1 bar / 600 K; mass follows the chosen EOS.
            samples: Uniform output sample count, >=2.
            method: Implicit solve_ivp method: Radau or BDF.
            rtol: Relative tolerance [-].
            atol: Absolute tolerances [Pa, K, kg].
            max_step_degrees: Maximum process-angle step [degrees].

        Returns:
            SI trajectory, signed channel flows and cumulative Q/W/H [J].
            eos_relative_error compares p against the selected EOS.
            energy_relative_error audits Delta U=Q-W+H over the whole run.

        Raises:
            ValueError: Inconsistent EOS, invalid parameters or property state.
            RuntimeError: Failed numerical integration or nonphysical results.
        """
        if not isinstance(samples, int) or samples < 2:
            raise ValueError("samples must be an integer >= 2.")
        tolerances = np.asarray(atol, dtype=float)
        if (
            method not in {"Radau", "BDF"}
            or tolerances.shape != (3,)
            or not np.all(np.isfinite(tolerances))
            or np.any(tolerances <= 0)
            or not np.isfinite(rtol)
            or rtol <= 0
            or not np.isfinite(max_step_degrees)
            or max_step_degrees <= 0
        ):
            raise ValueError("Invalid method, tolerances or maximum step.")
        config = self.config
        volume0 = self.mechanism.get_volume(config.theta_start)
        given = (
            config.initial_state if initial_state is None else initial_state
        )
        if not np.isfinite(volume0) or volume0 <= 0:
            raise ValueError("Initial mechanism volume must be positive.")
        if given is None:
            default = self.fluid.from_pressure_temperature(1e5, 600.0)
            given = (1e5, 600.0, default.density * volume0)
        state = np.asarray(
            given,
            dtype=float,
        )
        if (
            state.shape != (3,)
            or not np.all(np.isfinite(state))
            or np.any(state <= 0)
        ):
            raise ValueError("Initial state requires three positive values.")
        eos = self.fluid.evaluate(state[1], state[2] / volume0)
        if not np.isclose(state[0], eos.pressure, rtol=1e-8, atol=0.0):
            raise ValueError(
                "Initial state must satisfy the EOS (ideal: pV=mRT)."
            )
        angles = np.linspace(config.theta_start, config.theta_end, samples)
        max_angle = np.deg2rad(max_step_degrees)
        if isinstance(config.combustion, WiebeCombustion):
            max_angle = min(max_angle, config.combustion.duration / 20)
        for boundary in config.boundaries:
            if isinstance(boundary, PoppetValve):
                duration = (
                    boundary.closing_theta - boundary.opening_theta
                ) % boundary.period
                max_angle = min(max_angle, duration / 20)
        solution = solve_ivp(
            self._rates,
            (0.0, (angles[-1] - angles[0]) / self.omega),
            np.concatenate((state, [0.0, 0.0, 0.0])),
            method=method,
            t_eval=(angles - angles[0]) / self.omega,
            rtol=rtol,
            atol=np.concatenate((tolerances, [1e-6, 1e-6, 1e-6])),
            max_step=max_angle / self.omega,
        )
        if not solution.success:
            raise RuntimeError(
                f"Engine integration failed: {solution.message}"
            )
        if not np.all(np.isfinite(solution.y)) or np.any(solution.y[:3] <= 0):
            raise RuntimeError("Integration produced a nonphysical state.")
        pressures, temperatures, masses = solution.y[:3]
        volumes = np.array([self.mechanism.get_volume(a) for a in angles])
        states = [
            self.fluid.evaluate(t, m / v)
            for t, m, v in zip(temperatures, masses, volumes)
        ]
        residual = max(
            abs(p / gas.pressure - 1) for p, gas in zip(pressures, states)
        )
        internal = masses * np.array([gas.internal_energy for gas in states])
        heat, work, enthalpy = solution.y[3:]
        change = internal[-1] - internal[0]
        energy_error = abs(change - heat[-1] + work[-1] - enthalpy[-1]) / max(
            abs(change), abs(heat[-1]), abs(work[-1]), abs(enthalpy[-1]), 1.0
        )
        flow_samples = [self._flows(a, gas) for a, gas in zip(angles, states)]
        return CycleResult(
            solution.t,
            angles,
            volumes,
            pressures,
            temperatures,
            masses,
            np.column_stack((volumes, pressures)),
            float(residual),
            {
                b.name: np.array([sample[b.name] for sample in flow_samples])
                for b in config.boundaries
            },
            internal,
            heat,
            work,
            enthalpy,
            float(energy_error),
        )


class OpposedPistonEngineCycle:
    """Constant-speed opposed-piston zero-dimensional research model.

    Args:
        bore: Cylinder bore [m].
        stroke: Stroke of EACH of the two identical pistons [m].
        con_rod: Connecting rod length [m], greater than stroke/2.
        rpm: Shaft speed [rev/min].
        compression_ratio: V_max/V_min for combined piston travel [-].
        fuel_mass: Fuel supplied per combustion event [kg].
        air_fuel_ratio: Prescribed air/fuel mass ratio [-]. Heat release
            is limited by min(1, AFR/stoichiometric_afr); no species tracked.
        theta_s: Start of combustion relative to TDC [rad], may be negative.
        delta_theta_b: Burn duration [rad], shorter than firing interval.
        m: Wiebe shape exponent [-], nonnegative.
        wiebe_a: Wiebe completeness coefficient [-].
        firing_interval: Combustion repeat interval [rad]; default 2*pi
            for a port-controlled two-stroke. 4*pi fires every other turn.
        lower_heating_value: Fuel LHV [J/kg].
        combustion_efficiency: Fraction of available fuel energy released [-].
        stoichiometric_afr: Stoichiometric air/fuel mass ratio [-].
        gas_constant: Working-fluid specific gas constant [J/(kg K)].
        gamma: Constant cp/cv [-], greater than one.
        wall_temperature: Fixed wall temperature [K].
        heat_transfer_multiplier: Woschni calibration factor [-]; zero
            disables convection for adiabatic verification.
        intake: Intake port; None creates a representative default port.
        exhaust: Exhaust port; None creates a representative default port.

    Notes:
        theta=0 is minimum volume, theta=pi maximum volume. Geometry and
        ports repeat every 2*pi regardless of firing_interval. Set port
        widths via a custom Port (or top >= stroke) to suppress exchange.
        The simplified Woschni coefficients require experimental calibration.
    """

    def __init__(
        self,
        bore: float = 0.08,
        stroke: float = 0.08,
        con_rod: float = 0.14,
        rpm: float = 2000.0,
        compression_ratio: float = 16.0,
        fuel_mass: float = 2e-5,
        air_fuel_ratio: float = 30.0,
        theta_s: float = np.deg2rad(-10.0),
        delta_theta_b: float = np.deg2rad(50.0),
        m: float = 2.0,
        wiebe_a: float = 6.908,
        firing_interval: float = 2 * np.pi,
        lower_heating_value: float = 42.5e6,
        combustion_efficiency: float = 0.98,
        stoichiometric_afr: float = 14.5,
        gas_constant: float = 287.0,
        gamma: float = 1.35,
        wall_temperature: float = 450.0,
        heat_transfer_multiplier: float = 1.0,
        intake: Port | None = None,
        exhaust: Port | None = None,
    ) -> None:
        """Validate configuration and derive SI geometric constants."""
        positive = (
            bore,
            stroke,
            con_rod,
            rpm,
            air_fuel_ratio,
            wiebe_a,
            firing_interval,
            lower_heating_value,
            stoichiometric_afr,
            gas_constant,
            wall_temperature,
            delta_theta_b,
        )
        other = (
            compression_ratio,
            fuel_mass,
            theta_s,
            m,
            gamma,
            combustion_efficiency,
            heat_transfer_multiplier,
        )
        if not np.all(np.isfinite(positive + other)):
            raise ValueError("All engine parameters must be finite.")
        if min(positive) <= 0 or compression_ratio <= 1 or gamma <= 1:
            raise ValueError("Invalid parameter, compression ratio or gamma.")
        if con_rod <= stroke / 2 or delta_theta_b >= firing_interval:
            raise ValueError("Invalid rod length or overlapping combustion.")
        if min(fuel_mass, m, heat_transfer_multiplier) < 0:
            raise ValueError("Fuel, m and heat factor must be nonnegative.")
        if not 0 <= combustion_efficiency <= 1:
            raise ValueError("Combustion efficiency must be in [0, 1].")
        self.bore, self.stroke, self.con_rod = bore, stroke, con_rod
        self.rpm, self.compression_ratio = rpm, compression_ratio
        self.fuel_mass, self.air_fuel_ratio = fuel_mass, air_fuel_ratio
        self.theta_s, self.delta_theta_b = theta_s, delta_theta_b
        self.m, self.wiebe_a = m, wiebe_a
        self.firing_interval = firing_interval
        self.gas_constant, self.gamma = gas_constant, gamma
        self.cv = gas_constant / (gamma - 1)
        self.cp = gamma * self.cv
        self.wall_temperature = wall_temperature
        self.heat_transfer_multiplier = heat_transfer_multiplier
        self.omega = rpm * 2 * np.pi / 60
        self.piston_area = np.pi * bore**2 / 4
        self.swept_volume = 2 * self.piston_area * stroke
        self.clearance_volume = self.swept_volume / (compression_ratio - 1)
        self.heat_per_event = (
            fuel_mass
            * lower_heating_value
            * combustion_efficiency
            * min(1.0, air_fuel_ratio / stoichiometric_afr)
        )
        self.intake = (
            intake
            if intake is not None
            else Port(0.82 * stroke, 0.12 * stroke, 0.08, 1.3e5, 320.0)
        )
        self.exhaust = (
            exhaust
            if exhaust is not None
            else Port(0.75 * stroke, 0.15 * stroke, 0.08, 1.0e5, 700.0)
        )

        self.mechanism = OpposedPistonMechanism(
            bore, stroke, con_rod, compression_ratio
        )
        self.fluid = IdealGas(gas_constant, gamma)
        self.combustion = WiebeCombustion(
            self.heat_per_event,
            theta_s,
            delta_theta_b,
            m,
            wiebe_a,
            firing_interval,
        )
        self.wall_model = WoschniHeatTransfer(
            bore, stroke, rpm, wall_temperature, heat_transfer_multiplier
        )
        boundaries = tuple(
            PistonMaskedPort(
                name=name,
                reservoir=Reservoir(port.pressure, port.temperature),
                discharge_coeff=port.discharge_coefficient,
                height=port.height,
                width_total=port.width,
                edge_position=port.top,
                piston_displacement=self.mechanism.get_piston_displacement,
            )
            for name, port in (
                ("intake", self.intake),
                ("exhaust", self.exhaust),
            )
        )
        self.simulator = EngineSimulator(
            EngineConfig(
                mechanism=self.mechanism,
                boundaries=boundaries,
                fluid=self.fluid,
                combustion=self.combustion,
                rpm=rpm,
                wall_heat_transfer=self.wall_model,
            )
        )

    def piston_motion(self, theta: float) -> tuple[float, float]:
        """Return single-piston displacement [m] and derivative [m/rad]."""
        return self.mechanism.piston_motion(theta)

    def volume(self, theta: float) -> tuple[float, float]:
        """Return volume [m³] and angular derivative [m³/rad]."""
        return (
            self.mechanism.get_volume(theta),
            self.mechanism.get_dV_dtheta(theta),
        )

    def combustion_heat_release(self, theta: float) -> float:
        """Return normalized periodic Wiebe heat release [J/rad]."""
        return self.combustion.heat_release(theta)

    def port_area(self, theta: float, port: Port) -> float:
        """Return uncovered rectangular geometric area [m²], excluding Cd."""
        return float(
            port.width
            * np.clip(
                self.mechanism.get_piston_displacement(theta) - port.top,
                0,
                port.height,
            )
        )

    def mass_flow(
        self, theta: float, p: float, temperature: float, port: Port
    ) -> float:
        """Return inward-positive signed nozzle flow [kg/s], including Cd."""
        return nozzle_mass_flow(
            port.discharge_coefficient * self.port_area(theta, port),
            self.fluid.from_pressure_temperature(p, temperature),
            self.fluid.from_pressure_temperature(
                port.pressure, port.temperature
            ),
        )

    def wall_heat_transfer(
        self, theta: float, p: float, temperature: float
    ) -> float:
        """Return simplified Woschni wall loss [W], positive outward."""
        return self.wall_model(theta, p, temperature, self.volume(theta)[0])

    def ode_system(
        self, t: float, y: Sequence[float], theta: float | None = None
    ) -> Array:
        """Return [pdot, Tdot, mdot] [Pa/s, K/s, kg/s] via EngineSimulator."""
        return self.simulator.ode_system(t, y, theta)

    def solve_cycle(
        self,
        initial_state: Sequence[float] | None = None,
        samples: int = 1441,
        method: str = "Radau",
        rtol: float = 1e-7,
        atol: Sequence[float] = (1e-2, 1e-6, 1e-12),
        max_step_degrees: float = 1.0,
    ) -> CycleResult:
        """Solve 0–720° with the shared configurable ODE integrator.

        Args:
            initial_state: [p Pa, T K, mass kg], satisfying pV=mRT.
            samples: Number of output angles, >=2.
            method: Radau or BDF implicit method.
            rtol: Relative integration tolerance [-].
            atol: Absolute tolerances [Pa, K, kg].
            max_step_degrees: Maximum crank-angle integration step [degrees].

        Returns:
            SI trajectory and conservation diagnostics. This adapter preserves
            the original public API while sharing the modular solver physics.
        """
        return self.simulator.solve_cycle(
            initial_state, samples, method, rtol, atol, max_step_degrees
        )
