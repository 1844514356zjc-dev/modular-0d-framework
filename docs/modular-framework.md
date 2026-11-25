# Component composition and extension guide

All core magnitudes are SI, angles are radians, and channel flow is positive
into the simulated control volume. This framework currently models one
homogeneous gas control volume, not an arbitrary hydraulic/thermal network.

## Typed configuration

`EngineConfig` supplies:

| Field | Contract |
| --- | --- |
| mechanism | `BaseMechanism.get_volume(theta)` [m³] and `get_dV_dtheta(theta)` [m³/rad] |
| boundaries | Tuple of `GasExchangeBoundary` objects with unique names; empty for closed gas |
| fluid | `BaseWorkingFluid.evaluate(T,rho)` and `from_pressure_temperature(p,T)` |
| combustion | Object exposing `heat_release(theta)` [J/rad], including negative heat removal |
| wall_heat_transfer | Optional callable `(theta,p,T,V) -> W`, positive heat loss |
| rpm | Constant process-coordinate rate [rev/min] |
| theta_start / theta_end | Integration bounds [rad], end strictly greater than start |
| initial_state | Optional `(pressure Pa, temperature K, mass kg)` satisfying the chosen EOS |

`EngineSimulator(config).ode_system(t,y)` exposes physical derivatives for
`y=[p,T,m]`. `solve_cycle()` also integrates Q, W and H for a first-law audit.
Pressure is a redundant diagnostic coordinate; thermodynamic properties are
computed from T and rho=m/V, and initial EOS consistency is checked explicitly.

`CycleResult.boundary_flows` preserves the actual configured channel names.
`internal_energy`, `heat`, `work`, and `enthalpy_exchange` are total/cumulative
energies in J, not mass-specific values. The energy residual checks
`Delta U - Q + W - H`. The EOS residual compares p against the injected EOS;
it does not assume pV=mRT when a CoolProp gas is selected.

## Dictionary configuration

Nested `type` selectors:

- mechanism: `crank`, `opposed`, `otto` (the latter wraps a single crank in an
  isochoric-dwell Otto volume schedule).
- fluid: `ideal` (gas_constant, gamma), `coolprop` (name, backend).
- combustion: `none`, `wiebe`, `otto`.
- boundary: `port`, `poppet`, with a nested reservoir `{pressure, temperature}`.

Already constructed objects may replace any nested dictionary. No input map is
mutated. Unknown component selectors are rejected. A wall-loss dictionary is
passed to `WoschniHeatTransfer`; a custom callable can replace it.

```python
from core.engine_solver import EngineSimulator

case = {
    "mechanism": {"type": "opposed", "bore": 0.08, "stroke": 0.08,
                  "con_rod": 0.14, "compression_ratio": 16.0},
    "fluid": {"type": "ideal", "gas_constant": 287.0, "gamma": 1.35},
    "combustion": {"type": "none"},
    "boundaries": [{
        "type": "port", "name": "scavenge",
        "height": 0.0096, "width_total": 0.08, "edge_position": 0.0656,
        "discharge_coeff": 0.65,
        "reservoir": {"pressure": 130000.0, "temperature": 320.0},
    }],
}
result = EngineSimulator(case).solve_cycle()
print(result.boundary_flows["scavenge"])
```

When creating a `PistonMaskedPort` directly, pass
`piston_displacement=mechanism.get_piston_displacement` (or a custom callable).
The dictionary factory automatically supplies this callback for the built-in
crank/opposed mechanisms. For a custom volume-only mechanism, the port requires
an explicit displacement callback: volume alone cannot determine edge masking.
Ports and valves can be mixed without changing `EngineSimulator`.

## Property and flow assumptions

`IdealGas` has constant R and gamma. `CoolPropGas` uses one mutable
`AbstractState` per adapter, evaluated at T,rho; do not share an adapter between
concurrently executing solvers. Independent simulations should construct their
own adapters. Named gases use the same EOS for cylinder and reservoirs.

The caloric first law includes `(du/drho)_T`; omitting this term while changing
to a real EOS is thermodynamically inconsistent. Liquid, two-phase, and
out-of-range temperature states are rejected. CoolProp does not turn the
prescribed Wiebe heat source into a reacting-mixture combustion model.

Nozzle flow is quasi-steady and locally idealized: upstream cp/cv and
`R_eff=p/(rho*T)` are used in the ideal-gas choking equation. It is not an exact
real-gas isentropic nozzle solver. Incoming enthalpy is the reservoir h;
outgoing enthalpy is cylinder h, using the same fluid reference convention.

## Preset interpretation

- **Four stroke:** theta=0 is firing TDC, expansion then exhaust/intake then
  compression. Geometry repeats every 360°; cam and heat laws repeat every
  720°. Default exhaust timing is 140–370°, intake 350–580°.
- **Opposed pistons:** two identical synchronous slider–cranks; stroke is per
  piston. Geometry and heat repeat every 360°. A novel single-stroke linkage
  requires its own motion-law implementation; none is inferred here.
- **Otto benchmark:** 0–180° compression, 180–360° isochoric heating,
  360–540° expansion, 540–720° isochoric cooling. Heat rejection matches the
  supplied ideal-gas gamma and geometry ratio. Inconsistent Otto components,
  open boundaries or wall losses are rejected rather than labeled ideal Otto.
- **Gas spring:** closed mass, no heat or wall loss; no defined thermal
  efficiency. Complete periods are reversible for the implemented idealization.

UI engine initial states use T=800 K and m=30 times fuel mass; pressure follows
the selected EOS. Core defaults are p=1 bar, T=600 K, with EOS-derived mass.
Neither default automatically represents a periodic steady engine charge.

UI fresh-charge tracing is a separate post-processing model with Y=0 at the
assessment start; it measures retention of the selected interval's fresh gas,
not a solved history of chemical residual species. Closed systems export N/A
metrics as JSON null. The full settings snapshot identifies active topology,
property backend and heat model. `heat_input_J` records source input energy;
`fuel_energy_J` is the legacy engine-only field and is null for closed cases.
Inactive preset defaults are not physical
components in the assembled simulation.

## Extension checklist

Implement the smallest appropriate contract, document formulas and SI units,
then verify independent limiting cases: volume derivatives, mass/energy
balance, entropy for a closed adiabatic gas, flow sign/choking, or an analytic
cycle efficiency. User-defined heat features must be resolved by the supplied
`max_step_degrees`; the built-in Wiebe model additionally limits steps to at
most one twentieth of its burn duration. Built-in cam intervals receive the
same one-twentieth resolution limit. Custom callbacks should be pure with
respect to solver state, because implicit solvers evaluate trial states and
numerical Jacobians.

The unchanged `OpposedPistonEngineCycle` constructor delegates to the same
mechanism, heat, boundary and ODE components. Existing `Port` callers remain
supported; new compositions should use the boundary classes directly.
