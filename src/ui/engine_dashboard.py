"""SI post-processing and Plotly figures for the opposed-piston laboratory.

Species post-processing uses the perfectly stirred control-volume balance
(Heywood, Internal Combustion Engine Fundamentals, gas exchange). It is a
passive fresh-air tracer, not a combustion chemistry or multi-zone model.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import json
from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from core.engine_solver import EngineConfig, EngineSimulator
from core.heat_models import (
    NoHeatRelease,
    OttoHeatRelease,
    WiebeCombustion,
    WoschniHeatTransfer,
)
from core.mechanics import (
    CrankMechanism,
    OpposedPistonMechanism,
    OttoCycleMechanism,
)
from core.scavenging_model import ScavengingEvaluator
from core.valves_ports import PistonMaskedPort, PoppetValve, Reservoir
from core.working_fluids import CoolPropGas, IdealGas

COLORS = {
    "Compression": "#798A99",
    "Combustion": "#ff7f0e",
    "Expansion": "#1f77b4",
    "Scavenging": "#2c3e50",
    "Gas exchange": "#2c3e50",
    "Heat addition": "#ff7f0e",
    "Heat rejection": "#6D7792",
}


@dataclass(frozen=True)
class EngineSettings:
    """UI configuration: lengths [mm], angles [deg], speed [rev/min].

    All topology presets feed EngineConfig and EngineSimulator. CoolProp
    identifiers are selected explicitly; Ideal gas uses R=287, gamma=1.35.
    """

    bore_mm: float = 80.0
    stroke_mm: float = 80.0
    rpm: float = 2000.0
    start_deg: float = -10.0
    duration_deg: float = 50.0
    wiebe_m: float = 2.0
    intake_height_mm: float = 9.6
    exhaust_height_mm: float = 12.0
    cd: float = 0.65
    fuel_mg: float = 20.0
    compression_ratio: float = 16.0
    topology: str = "opposed"
    fluid: str = "Ideal gas"
    valve_lift_mm: float = 8.0
    valve_diameter_mm: float = 30.0
    intake_open_deg: float = 350.0
    intake_close_deg: float = 580.0
    exhaust_open_deg: float = 140.0
    exhaust_close_deg: float = 370.0
    closed_model: str = "Air-standard Otto"
    heat_input_J: float = 500.0
    initial_pressure_bar: float = 1.0
    initial_temperature_K: float = 300.0


def build_simulator(settings: EngineSettings) -> EngineSimulator:
    """Assemble topology-specific components without UI or plotting state.

    Args:
        settings: Display-unit controls; conversion to SI occurs here.

    Returns:
        Shared EngineSimulator configured with real interchangeable components.

    Raises:
        ValueError: Unsupported topology, model, fluid or invalid geometry.
    """
    if settings.topology not in {"opposed", "four_stroke", "closed"}:
        raise ValueError("Unsupported system topology.")
    fluid = (
        IdealGas()
        if settings.fluid == "Ideal gas"
        else CoolPropGas(settings.fluid)
    )
    stroke, bore = settings.stroke_mm / 1000, settings.bore_mm / 1000
    mechanism_type = (
        OpposedPistonMechanism
        if settings.topology == "opposed"
        else CrankMechanism
    )
    crank = mechanism_type(
        bore, stroke, 1.75 * stroke, settings.compression_ratio
    )
    mechanism = crank
    boundaries = ()
    wall = None
    if settings.topology == "closed":
        if settings.closed_model == "Air-standard Otto":
            if not isinstance(fluid, IdealGas):
                raise ValueError(
                    "Air-standard Otto requires constant Ideal gas."
                )
            mechanism = OttoCycleMechanism(crank)
            heat = OttoHeatRelease(
                settings.heat_input_J, settings.compression_ratio, fluid.gamma
            )
        elif settings.closed_model == "Adiabatic gas spring":
            heat = NoHeatRelease()
        else:
            raise ValueError("Unknown closed-cycle heat model.")
        initial = fluid.from_pressure_temperature(
            settings.initial_pressure_bar * 1e5, settings.initial_temperature_K
        )
        mass = initial.density * mechanism.get_volume(0)
    else:
        mass = 30 * settings.fuel_mg * 1e-6
        if mass <= 0:
            raise ValueError("Engine fuel mass must be positive.")
        initial = fluid.evaluate(800.0, mass / mechanism.get_volume(0))
        heat = WiebeCombustion(
            settings.fuel_mg * 1e-6 * 42.5e6 * 0.98,
            np.deg2rad(settings.start_deg),
            np.deg2rad(settings.duration_deg),
            settings.wiebe_m,
            period=(
                2 * np.pi if settings.topology == "opposed" else 4 * np.pi
            ),
        )
        intake_reservoir = Reservoir(
            1.3e5 if settings.topology == "opposed" else 1.05e5, 320
        )
        exhaust_reservoir = Reservoir(1e5, 700)
        if settings.topology == "opposed":
            ports = []
            for name, edge, height, reservoir in (
                (
                    "intake",
                    0.82 * stroke,
                    settings.intake_height_mm / 1000,
                    intake_reservoir,
                ),
                (
                    "exhaust",
                    0.75 * stroke,
                    settings.exhaust_height_mm / 1000,
                    exhaust_reservoir,
                ),
            ):
                if edge + height > stroke:
                    raise ValueError(
                        "Port edge plus height exceeds piston stroke."
                    )
                ports.append(
                    PistonMaskedPort(
                        name=name,
                        reservoir=reservoir,
                        discharge_coeff=settings.cd,
                        height=height,
                        width_total=0.08,
                        edge_position=edge,
                        piston_displacement=crank.get_piston_displacement,
                    )
                )
            boundaries = tuple(ports)
        else:
            boundaries = tuple(
                PoppetValve(
                    name=name,
                    reservoir=reservoir,
                    discharge_coeff=settings.cd,
                    diameter=settings.valve_diameter_mm / 1000,
                    max_lift=settings.valve_lift_mm / 1000,
                    opening_theta=np.deg2rad(opening),
                    closing_theta=np.deg2rad(closing),
                )
                for name, opening, closing, reservoir in (
                    (
                        "intake",
                        settings.intake_open_deg,
                        settings.intake_close_deg,
                        intake_reservoir,
                    ),
                    (
                        "exhaust",
                        settings.exhaust_open_deg,
                        settings.exhaust_close_deg,
                        exhaust_reservoir,
                    ),
                )
            )
        wall = WoschniHeatTransfer(bore, stroke, settings.rpm)
    return EngineSimulator(
        EngineConfig(
            mechanism=mechanism,
            fluid=fluid,
            combustion=heat,
            boundaries=boundaries,
            wall_heat_transfer=wall,
            rpm=settings.rpm,
            initial_state=(initial.pressure, initial.temperature, mass),
        )
    )


def simulate(settings: EngineSettings) -> dict[str, Any]:
    """Solve the selected topology and post-process SI trajectory/metrics.

    Args:
        settings: Topology, model and display-unit input snapshot.

    Returns:
        Common report schema for engines and closed cycles. Closed systems
        return NaN gas-exchange efficiencies (N/A), never fictitious delivery.
        Work and heat budgets are integrated alongside the physical ODEs.
    """
    simulator = build_simulator(settings)
    config, mechanism = simulator.config, simulator.mechanism
    result = simulator.solve_cycle()
    angles = result.theta
    count = len(angles)
    flows = np.column_stack(
        [
            result.boundary_flows.get(name, np.zeros(count))
            for name in ("intake", "exhaust")
        ]
    )
    areas, lifts = {}, {}
    for boundary in config.boundaries:
        areas[boundary.name] = np.array(
            [boundary.get_geometric_area(a) for a in angles]
        )
        if isinstance(boundary, PoppetValve):
            lifts[boundary.name] = np.array(
                [boundary.get_lift(a) for a in angles]
            )
    if isinstance(mechanism, OttoCycleMechanism):
        displacement = (
            result.volume - mechanism.crank.clearance_volume
        ) / mechanism.crank.piston_area
    else:
        displacement = np.array(
            [mechanism.get_piston_displacement(a) for a in angles]
        )
    inlet_area = areas.get("intake", np.zeros(count))
    outlet_area = areas.get("exhaust", np.zeros(count))
    phases = []
    for index, angle in enumerate(angles):
        release = config.combustion.heat_release(angle)
        if settings.topology == "closed" and release != 0:
            phase = "Heat addition" if release > 0 else "Heat rejection"
        elif release > 0:
            phase = "Combustion"
        elif inlet_area[index] > 0 or outlet_area[index] > 0:
            phase = (
                "Scavenging"
                if settings.topology == "opposed"
                else "Gas exchange"
            )
        else:
            phase = (
                "Expansion"
                if mechanism.get_dV_dtheta(angle) >= 0
                else "Compression"
            )
        phases.append(phase)
    frame = pd.DataFrame(
        {
            "time_s": result.time,
            "theta_rad": angles,
            "theta_deg": np.rad2deg(angles),
            "volume_m3": result.volume,
            "pressure_Pa": result.pressure,
            "temperature_K": result.temperature,
            "mass_kg": result.mass,
            "piston_displacement_m": displacement,
            "intake_area_geo_m2": inlet_area,
            "exhaust_area_geo_m2": outlet_area,
            "intake_area_eff_m2": inlet_area * settings.cd,
            "exhaust_area_eff_m2": outlet_area * settings.cd,
            "intake_flow_kg_s": flows[:, 0],
            "exhaust_flow_kg_s": -flows[:, 1],
            "intake_lift_m": lifts.get("intake", np.zeros(count)),
            "exhaust_lift_m": lifts.get("exhaust", np.zeros(count)),
            "internal_energy_J": result.internal_energy,
            "net_heat_J": result.heat,
            "work_J": result.work,
            "enthalpy_exchange_J": result.enthalpy_exchange,
            "phase": phases,
        }
    )
    metrics = dict(
        DR=float("nan"),
        eta_tr=float("nan"),
        eta_sc=float("nan"),
        m_del=0.0,
        m_exhaust=0.0,
    )
    interval = None
    trapped, residual = None, None
    if settings.topology != "closed":
        if settings.topology == "opposed":
            opening = brentq(
                lambda a: (
                    mechanism.get_piston_displacement(a)
                    - min(b.edge_position for b in config.boundaries)
                ),
                0,
                np.pi,
            )
            start, end = 2 * np.pi + opening, 4 * np.pi - opening
        else:
            start = min(b.opening_theta for b in config.boundaries)
            end = max(b.closing_theta for b in config.boundaries)
        interval = list(np.rad2deg([start, end]))

        def tracer_rhs(time: float, fraction: np.ndarray) -> list[float]:
            """Return dY/dt [1/s] from the perfectly mixed species balance.

            Fresh intake carries Y=1, exhaust backflow Y=0; mixed outflow
            cancels in the fraction balance. Initial sweep inventory has Y=0.
            """
            angle = time * simulator.omega
            inflow = max(float(np.interp(angle, angles, flows[:, 0])), 0)
            backflow = max(float(np.interp(angle, angles, flows[:, 1])), 0)
            mass = float(np.interp(angle, angles, result.mass))
            return [
                (inflow * (1 - fraction[0]) - backflow * fraction[0]) / mass
            ]

        tracer = solve_ivp(
            tracer_rhs,
            (start / simulator.omega, end / simulator.omega),
            [0.0],
            rtol=1e-8,
            atol=1e-10,
            max_step=np.deg2rad(0.5) / simulator.omega,
        )
        if not tracer.success or not 0 <= tracer.y[0, -1] <= 1:
            raise RuntimeError("Perfect-mixing tracer integration failed.")
        mass_end = float(np.interp(end, angles, result.mass))
        trapped = mass_end * float(tracer.y[0, -1])
        residual = mass_end - trapped
        reference = simulator.fluid.from_pressure_temperature(
            config.boundaries[0].reservoir.pressure,
            config.boundaries[0].reservoir.temperature,
        )
        metrics = ScavengingEvaluator.evaluate_cycle_metrics(
            angles,
            flows[:, 0],
            -flows[:, 1],
            trapped,
            residual,
            mechanism.swept_volume,
            reference.density,
            rpm=settings.rpm,
            scavenging_interval=(start, end),
        )
    if isinstance(config.combustion, WiebeCombustion):
        heat_input = (
            4
            * np.pi
            / config.combustion.period
            * config.combustion.heat_per_event
            / 0.98
        )
    elif isinstance(config.combustion, OttoHeatRelease):
        heat_input = config.combustion.heat_input
    else:
        heat_input = 0.0
    work = float(result.work[-1])
    return {
        "schema_version": 2,
        "frame": frame,
        "settings": asdict(settings),
        "initial_state": dict(
            zip(
                ("pressure_Pa", "temperature_K", "mass_kg"),
                config.initial_state,
            )
        ),
        "metrics": metrics,
        "work_J": work,
        "thermal_efficiency": work / heat_input
        if heat_input
        else float("nan"),
        "heat_input_J": heat_input,
        "fuel_energy_J": heat_input if settings.topology != "closed" else None,
        "eos_error": result.eos_relative_error,
        "energy_error": result.energy_relative_error,
        "scavenging_interval_deg": interval,
        "m_trap_kg": trapped,
        "m_res_kg": residual,
        "model": f"{settings.topology}; {settings.fluid}; prescribed volume. "
        "Engine cases are start-up transients. CoolProp gas uses "
        "full caloric/EOS derivatives and a locally ideal nozzle. "
        "Exchange metrics use a perfectly mixed passive tracer; "
        "preexisting charge and exhaust backflow count as residual.",
    }


def style_figure(figure: go.Figure, x_title: str, y_title: str) -> go.Figure:
    """Style a scientific chart with labeled engineering-unit axes.

    Args:
        figure: Plotly figure to update in place.
        x_title: Horizontal quantity and unit.
        y_title: Vertical quantity and unit.

    Returns:
        Styled figure with a single aligned Cartesian grid.
    """
    figure.update_layout(
        template="plotly_white",
        height=410,
        margin=dict(l=60, r=60, t=48, b=55),
        font=dict(family="Arial, sans-serif", size=12, color="#2c3e50"),
        paper_bgcolor="#FFFFFF",
        plot_bgcolor="#FFFFFF",
        legend=dict(orientation="h", y=1.15, x=0, font_size=11),
        hovermode="x unified",
        hoverlabel=dict(font_size=12),
    )
    figure.update_xaxes(
        title=x_title,
        gridcolor="#E8EDF2",
        zeroline=False,
        showline=True,
        linecolor="#C9D3DC",
    )
    figure.update_yaxes(
        title=y_title,
        gridcolor="#E8EDF2",
        zeroline=False,
        showline=True,
        linecolor="#C9D3DC",
    )
    return figure


def thermodynamic_figure(
    frame: pd.DataFrame,
    temperature: bool = False,
    log_axes: bool = False,
) -> go.Figure:
    """Plot phase-colored P-V or T-angle paths with exact SI hover values.

    Args:
        frame: Simulation samples with explicitly named SI columns.
        temperature: Select T-angle rather than P-V.
        log_axes: Enable logarithmic volume/pressure axes for P-V only.

    Returns:
        Interactive figure; displayed scales are cm³/bar or degrees/K.
    """
    figure = go.Figure()
    x = frame.theta_deg if temperature else frame.volume_m3 * 1e6
    y = frame.temperature_K if temperature else frame.pressure_Pa / 1e5
    hover = frame[
        ["theta_rad", "volume_m3", "pressure_Pa", "temperature_K"]
    ].to_numpy()
    for phase, color in COLORS.items():
        mask = frame.phase.eq(phase).to_numpy()
        if not np.any(mask):
            continue
        # Include the sample before each transition without bridging gaps.
        mask[:-1] |= mask[1:]
        figure.add_trace(
            go.Scatter(
                x=x,
                y=np.where(mask, y, np.nan),
                customdata=hover,
                name=phase,
                mode="lines",
                line=dict(color=color, width=2.3),
                connectgaps=False,
                hovertemplate=(
                    "Angle: %{customdata[0]:.7f} rad<br>"
                    "Volume: %{customdata[1]:.8e} m³<br>"
                    "Pressure: %{customdata[2]:.5f} Pa<br>"
                    "Temperature: %{customdata[3]:.5f} K"
                    "<extra>%{fullData.name}</extra>"
                ),
            )
        )
    style_figure(
        figure,
        "Crank angle / deg" if temperature else "Cylinder volume / cm³",
        "Temperature / K" if temperature else "Pressure / bar",
    )
    if temperature:
        figure.update_xaxes(dtick=90, range=[0, 720])
    else:
        figure.update_layout(hovermode="closest")
        if log_axes:
            figure.update_xaxes(type="log")
            figure.update_yaxes(type="log")
    return figure


def area_figure(frame: pd.DataFrame, effective: bool = False) -> go.Figure:
    """Plot piston travel [mm] and port area [mm²] with SI hover values.

    Args:
        frame: Simulation samples in SI, plus theta_deg for display.
        effective: Include Cd in area if True; otherwise show geometric area.

    Returns:
        Dual-axis figure with synchronized ticks and one visible grid.
    """
    figure = make_subplots(specs=[[{"secondary_y": True}]])
    figure.add_trace(
        go.Scatter(
            x=frame.theta_deg,
            y=frame.piston_displacement_m * 1000,
            customdata=frame.piston_displacement_m,
            name="Piston displacement",
            line=dict(color="#8495A5", dash="dot"),
            hovertemplate=(
                "%{customdata:.8f} m<extra>Piston displacement</extra>"
            ),
        ),
        secondary_y=False,
    )
    kind = "eff" if effective else "geo"
    for name, column, color in (
        (f"Intake A_{kind}", f"intake_area_{kind}_m2", "#1f77b4"),
        (f"Exhaust A_{kind}", f"exhaust_area_{kind}_m2", "#ff7f0e"),
    ):
        figure.add_trace(
            go.Scatter(
                x=frame.theta_deg,
                y=frame[column] * 1e6,
                name=name,
                customdata=frame[column],
                line=dict(color=color, width=2.2),
                hovertemplate=(
                    "%{customdata:.8e} m²<extra>%{fullData.name}</extra>"
                ),
            ),
            secondary_y=True,
        )
    style_figure(figure, "Crank angle / deg", "Piston displacement / mm")
    figure.update_layout(
        yaxis2=dict(
            title="Effective area / mm²"
            if effective
            else "Geometric area / mm²",
            tickmode="sync",
            showgrid=False,
        )
    )
    figure.update_xaxes(dtick=90, range=[0, 720])
    return figure


def flow_figure(frame: pd.DataFrame) -> go.Figure:
    """Plot signed flows [g/s] with exact kg/s hover data.

    Args:
        frame: SI samples; intake is inward-positive, exhaust outward-positive.

    Returns:
        Signed mass-flow figure; negative flow denotes reversal.
    """
    figure = go.Figure()
    for label, column, color in (
        ("Intake / inward positive", "intake_flow_kg_s", "#1f77b4"),
        ("Exhaust / outward positive", "exhaust_flow_kg_s", "#ff7f0e"),
    ):
        figure.add_trace(
            go.Scatter(
                x=frame.theta_deg,
                y=frame[column] * 1000,
                name=label,
                customdata=frame[column],
                line=dict(color=color, width=2),
                hovertemplate="%{customdata:.9f} kg/s"
                "<extra>%{fullData.name}</extra>",
            )
        )
    style_figure(figure, "Crank angle / deg", "Mass flow / g s⁻¹")
    figure.add_hline(y=0, line_color="#9EADB9", line_width=1)
    figure.update_xaxes(dtick=90, range=[0, 720])
    return figure


def run_parameter_sweep(
    settings: EngineSettings,
    parameter: str,
    values: list[float],
) -> pd.DataFrame:
    """Solve independent transient cases, preserving failures as data gaps.

    Args:
        settings: Baseline configuration; all other fields stay constant.
        parameter: 'rpm' [rev/min] or 'compression_ratio' [dimensionless].
        values: Finite strictly increasing parameter values, 2 to 12 samples.

    Returns:
        Rows with parameter value, efficiency [-], work [J], peak pressure
        [Pa], trapping efficiency [-], and an English error message. Failed
        cases contain NaN outputs, never zero or interpolated replacements.

    Raises:
        ValueError: Unsupported parameter, invalid range or sample count.
    """
    if parameter not in {"rpm", "compression_ratio"}:
        raise ValueError("Select RPM or compression ratio for the sweep.")
    samples = np.asarray(values, dtype=float)
    if (
        samples.ndim != 1
        or not 2 <= samples.size <= 12
        or not np.all(np.isfinite(samples))
        or np.any(np.diff(samples) <= 0)
    ):
        raise ValueError("Supply 2–12 finite, increasing sweep values.")
    if parameter == "rpm" and (samples.min() < 1000 or samples.max() > 8000):
        raise ValueError("RPM sweep must remain within 1000–8000 rev/min.")
    if parameter == "compression_ratio" and (
        samples.min() < 8 or samples.max() > 24
    ):
        raise ValueError("Compression ratio sweep must remain within 8–24.")
    rows = []
    for value in samples:
        row = {
            "parameter_value": float(value),
            "thermal_efficiency": np.nan,
            "work_J": np.nan,
            "peak_pressure_Pa": np.nan,
            "trapping_efficiency": np.nan,
            "error": "",
        }
        try:
            case = simulate(replace(settings, **{parameter: float(value)}))
            row.update(
                thermal_efficiency=case["thermal_efficiency"],
                work_J=case["work_J"],
                peak_pressure_Pa=float(case["frame"].pressure_Pa.max()),
                trapping_efficiency=case["metrics"]["eta_tr"],
            )
        except (ValueError, RuntimeError, ArithmeticError) as error:
            row["error"] = str(error)
        rows.append(row)
    return pd.DataFrame(rows)


def sweep_figure(frame: pd.DataFrame, parameter: str) -> go.Figure:
    """Show efficiency [%] and work [J], leaving failed-case gaps visible.

    Args:
        frame: Sweep results from run_parameter_sweep.
        parameter: Swept field name, rpm or compression_ratio.

    Returns:
        Dual-axis chart with efficiency/work and raw SI hover values.
    """
    figure = make_subplots(specs=[[{"secondary_y": True}]])
    for name, column, scale, color, secondary in (
        (
            "Indicated thermal efficiency",
            "thermal_efficiency",
            100,
            "#1f77b4",
            False,
        ),
        ("Indicated work", "work_J", 1, "#ff7f0e", True),
    ):
        figure.add_trace(
            go.Scatter(
                x=frame.parameter_value,
                y=frame[column] * scale,
                customdata=frame[column],
                name=name,
                mode="lines+markers",
                connectgaps=False,
                line=dict(color=color, width=2.2),
                hovertemplate="%{customdata:.8f}"
                + (" J" if secondary else " (dimensionless)")
                + "<extra>%{fullData.name}</extra>",
            ),
            secondary_y=secondary,
        )
    style_figure(
        figure,
        "Speed / rev min⁻¹" if parameter == "rpm" else "Compression ratio / 1",
        "Thermal efficiency / %",
    )
    figure.update_layout(
        yaxis2=dict(
            title="Indicated work / J", tickmode="sync", showgrid=False
        )
    )
    return figure


def json_export(payload: dict[str, Any]) -> str:
    """Serialize SI reports, mapping every undefined metric to JSON null.

    Args:
        payload: Report returned by simulate, optionally containing a sweep.

    Returns:
        Strict JSON with topology, property model and full numerical budgets.
    """
    report = {key: value for key, value in payload.items() if key != "frame"}
    report["samples"] = payload["frame"].to_dict(orient="records")

    def clean(value: Any) -> Any:
        """Convert NumPy scalars and undefined numbers to strict JSON."""
        if isinstance(value, dict):
            return {key: clean(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [clean(item) for item in value]
        if isinstance(value, np.generic):
            return clean(value.item())
        if isinstance(value, float) and not np.isfinite(value):
            return None
        return value

    return json.dumps(
        clean(report), ensure_ascii=False, allow_nan=False, indent=2
    )


def energy_figure(frame: pd.DataFrame) -> go.Figure:
    """Plot integrated first-law terms [J] for the selected gas control volume.

    Args:
        frame: SI samples with cumulative heat/work/enthalpy and total U.

    Returns:
        Energy budget chart with all terms referenced to the initial state.
    """
    figure = go.Figure()
    for name, values, color in (
        (
            "Internal energy change",
            frame.internal_energy_J - frame.internal_energy_J.iloc[0],
            "#2c3e50",
        ),
        ("Net heat input", frame.net_heat_J, "#ff7f0e"),
        ("Work output", frame.work_J, "#1f77b4"),
        ("Enthalpy exchange", frame.enthalpy_exchange_J, "#8495A5"),
    ):
        figure.add_trace(
            go.Scatter(
                x=frame.theta_deg,
                y=values,
                name=name,
                line=dict(color=color, width=2),
                hovertemplate="%{y:.8f} J<extra>%{fullData.name}</extra>",
            )
        )
    return style_figure(figure, "Process angle / deg", "Energy / J")


def valve_lift_figure(frame: pd.DataFrame) -> go.Figure:
    """Plot the two cam-driven valve lifts [mm], with SI hover values.

    Args:
        frame: SI valve lift histories, sampled against crank angle.

    Returns:
        Intake/exhaust lift diagram over the 720-degree cam period.
    """
    figure = go.Figure()
    for name, column, color in (
        ("Intake valve", "intake_lift_m", "#1f77b4"),
        ("Exhaust valve", "exhaust_lift_m", "#ff7f0e"),
    ):
        figure.add_trace(
            go.Scatter(
                x=frame.theta_deg,
                y=frame[column] * 1000,
                customdata=frame[column],
                name=name,
                line=dict(color=color, width=2),
                hovertemplate=(
                    "%{customdata:.8f} m<extra>%{fullData.name}</extra>"
                ),
            )
        )
    return style_figure(figure, "Crank angle / deg", "Valve lift / mm")
