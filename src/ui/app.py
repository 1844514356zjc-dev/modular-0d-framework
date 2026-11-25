"""English research dashboard: streamlit run src/ui/app.py.

Engineering units are explicit in controls and axes; hover values and exported
State columns retain SI units. This module owns presentation and UI state.
"""

from __future__ import annotations

from dataclasses import asdict
from html import escape
import json
import logging
import os
from pathlib import Path
import sys
from typing import Any
from urllib.parse import urlparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
import streamlit as st  # noqa: E402

from ui.engine_dashboard import (  # noqa: E402
    EngineSettings,
    area_figure,
    energy_figure,
    flow_figure,
    json_export,
    run_parameter_sweep,
    simulate,
    sweep_figure,
    thermodynamic_figure,
    valve_lift_figure,
)

LOGGER = logging.getLogger(__name__)
CSS = """
<style>
.stApp {background:#F3F6F9; color:#2c3e50;}
.block-container {max-width:1600px; padding-top:4.5rem; padding-bottom:3rem;}
html,body,.stApp {font-family:"Segoe UI",Arial,sans-serif;}
h1,h2,h3 {color:#2c3e50; font-weight:600 !important;}
h2 {font-size:1.3rem !important;} h3 {font-size:1.05rem !important;}
[data-testid="stSidebar"] {background:#172B40; border-right:1px solid #243E55;}
[data-testid="stSidebar"] [data-testid="stSidebarContent"] {padding-top:1rem;}
[data-testid="stSidebar"] :is(h1,h2,h3,p,label,summary,span) {color:#E2EAF2;}
[data-testid="stSidebar"] [data-testid="stCaptionContainer"] p {color:#ADC0D1;}
[data-testid="stSidebar"] details {border-color:#375067; border-radius:6px;}
[data-testid="stSidebar"] summary {background:#223B52 !important;}
[data-testid="stSidebar"] summary:hover {background:#29475F !important;}
[data-testid="stSidebar"] svg {color:#CFDEEB;}
[data-testid="stSidebar"] :is(input,[data-baseweb="select"]>div) {
 background:#223B52; color:#F0F5FA; border-color:#496278;}
[data-testid="stSidebar"] [data-testid="stNumberInput"] button {
 background:#29445C; color:#EDF3F8;}
[data-testid="stSidebar"] [data-testid="stSliderTickBar"] {color:#C2CFDB;}
[data-testid="stSidebar"] button[kind="primary"] {
 background:#1f77b4; border:1px solid #459CD5; min-height:46px;}
[data-testid="stMetric"] {background:#FFF; padding:20px;
 border:1px solid #DEE6ED; border-radius:7px; border-top:3px solid #1f77b4;}
[data-testid="stMetricLabel"] {font-size:.8rem; color:#607487;}
[data-testid="stMetricValue"] {font-size:1.8rem; font-weight:600;
 font-variant-numeric:tabular-nums; letter-spacing:-.04em;}
[data-baseweb="tab-list"] {gap:24px; border-bottom:1px solid #D9E2EA;
 margin-bottom:18px;}
button[data-baseweb="tab"] {font-size:13px; padding-bottom:14px;}
[data-testid="stPlotlyChart"] {border:1px solid #E0E7EE; border-radius:7px;
 overflow:hidden; background:#FFF;}
.stButton button,.stDownloadButton button {border-radius:5px; font-weight:600;}
.lab-header {display:flex; align-items:flex-start; gap:24px;
 justify-content:space-between; padding-bottom:22px;
 border-bottom:1px solid #D9E2EA; margin-bottom:20px;}
.lab-header h1 {font-size:30px; line-height:1.2; letter-spacing:-1px;
 margin:0 0 12px; max-width:850px; padding:0;}
.lab-header p {font-size:14px; line-height:1.6; color:#64798A; margin:0;}
.repo-badge {display:inline-block; white-space:nowrap; color:#2c3e50;
 text-decoration:none; font-size:12px; background:#FFF;
 border:1px solid #CEDAE4;
 border-radius:5px; padding:8px 12px; margin-top:4px;}
.run-context {padding:12px 16px; background:#E8EFF5;
 border-left:3px solid #1f77b4;
 border-radius:3px; font-size:13px; color:#3C566C; margin:18px 0 22px;}
.empty-state {background:white; border:1px solid #DCE5ED; border-radius:8px;
 padding:44px; margin-top:24px; max-width:960px;}
.empty-state h2 {margin-bottom:16px;}
.empty-state p {color:#64798A; max-width:680px; line-height:1.7;}
@media(max-width:760px) {
 .block-container {padding:4.5rem 1rem 1rem;}
 .lab-header {flex-direction:column; gap:12px;}
 .lab-header h1 {font-size:25px;} .empty-state {padding:24px;}
 [data-testid="stMetric"] {padding:14px;}
}
</style>
"""


@st.cache_data(show_spinner=False, max_entries=8)
def cached_simulation(settings: EngineSettings) -> dict[str, Any]:
    """Cache results for a submitted engineering-unit configuration."""
    return simulate(settings)


def render_header() -> None:
    """Show the product identity and an optional verified repository link.

    Set GITHUB_REPOSITORY_URL to the project's public GitHub repository URL.
    Without it, display a static badge rather than inventing a repository.
    """
    repository = os.environ.get("GITHUB_REPOSITORY_URL", "")
    parsed = urlparse(repository)
    if (
        parsed.scheme == "https"
        and parsed.netloc == "github.com"
        and len(parsed.path.strip("/").split("/")) == 2
    ):
        badge = (
            f'<a class="repo-badge" href="{escape(repository)}" '
            'target="_blank" rel="noopener noreferrer">'
            "GitHub / Open source ↗</a>"
        )
    else:
        badge = '<span class="repo-badge">GitHub / Open source · MIT</span>'
    st.markdown(
        '<div class="lab-header"><div><h1>Modular 0D Thermo-Fluid '
        "&amp; Engine Framework</h1><p>Dynamic thermodynamics, "
        "modular mechanisms, gas exchange and closed-cycle analysis."
        "</p></div>" + badge + "</div>",
        unsafe_allow_html=True,
    )


def control_panel() -> tuple[EngineSettings, bool]:
    """Refresh model-specific panels from the selected system topology.

    Returns:
        Engineering-unit settings and the explicit simulation command.
    """
    with st.sidebar:
        st.header("Simulation workspace")
        st.caption("Control panel")
        with st.expander("System Topology", expanded=True):
            label = st.selectbox(
                "System topology",
                [
                    "Opposed-piston ports",
                    "Four-stroke poppet valves",
                    "Closed thermodynamic cycle",
                ],
                key="system_topology",
            )
            topology = {
                "Opposed-piston ports": "opposed",
                "Four-stroke poppet valves": "four_stroke",
                "Closed thermodynamic cycle": "closed",
            }[label]
            extra = dict(topology=topology)
            if topology == "closed":
                extra["closed_model"] = st.selectbox(
                    "Closed-cycle model",
                    [
                        "Air-standard Otto",
                        "Adiabatic gas spring",
                    ],
                )
            is_otto = (
                topology == "closed"
                and extra["closed_model"] == "Air-standard Otto"
            )
            fluids = (
                ["Ideal gas"]
                if is_otto
                else ["Ideal gas", "Air", "Nitrogen", "CO2"]
            )
            fluid = st.selectbox(
                "Working fluid / EOS",
                fluids,
                key=f"{topology}_fluid",
                help="Named fluids use CoolProp HEOS, single-phase gas only.",
            )
            extra["fluid"] = fluid
            st.caption(
                "Constant cp/cv"
                if fluid == "Ideal gas"
                else "CoolProp caloric properties and EOS derivatives"
            )
            if is_otto:
                st.caption(
                    "Air-standard Otto uses an isochoric-dwell "
                    "process coordinate, not physical cam timing."
                )
            elif topology == "four_stroke":
                st.caption("Single piston · 720° valve and firing period")
            elif topology == "opposed":
                st.caption("Synchronous pistons · 360° firing period")
        with st.expander("Engine / Cycle Parameters", expanded=True):
            bore = st.slider("Bore / mm", 50.0, 140.0, 80.0, 1.0)
            stroke = st.slider(
                "Stroke per piston / mm", 60.0, 160.0, 80.0, 1.0
            )
            rpm = st.slider("Engine speed / RPM", 1000, 8000, 2000, 100)
            ratio = st.slider(
                "Compression ratio / 1",
                8.0,
                24.0,
                16.0 if topology == "opposed" else 10.0,
                0.5,
                key=f"{topology}_ratio",
            )
            if topology == "closed":
                st.caption(
                    "RPM sets process duration; the ideal closed "
                    "benchmark has no rate-dependent losses."
                )
        model_label = (
            "Heat Source & Initial State"
            if topology == "closed"
            else "Combustion & Exchange Models"
        )
        with st.expander(model_label, expanded=topology == "closed"):
            if topology == "closed":
                extra["initial_pressure_bar"] = st.number_input(
                    "Initial pressure / bar", 0.5, 5.0, 1.0, 0.1
                )
                extra["initial_temperature_K"] = st.number_input(
                    "Initial temperature / K", 280.0, 800.0, 300.0, 10.0
                )
                if is_otto:
                    extra["heat_input_J"] = st.number_input(
                        "Isochoric heat input / J", 50.0, 2000.0, 500.0, 50.0
                    )
                else:
                    st.caption("No heat input, wall loss or gas exchange.")
            else:
                extra["start_deg"] = st.slider(
                    "Wiebe start / deg CA", -30.0, 10.0, -10.0, 1.0
                )
                extra["duration_deg"] = st.slider(
                    "Burn duration / deg CA", 20.0, 100.0, 50.0, 1.0
                )
                extra["wiebe_m"] = st.slider(
                    "Wiebe shape factor / 1", 0.0, 5.0, 2.0, 0.1
                )
                extra["fuel_mg"] = st.number_input(
                    "Fuel per firing / mg",
                    1.0,
                    100.0,
                    20.0 if fluid == "Ideal gas" else 10.0,
                    1.0,
                    key=f"{topology}_{fluid}_fuel",
                )
                extra["cd"] = st.slider(
                    "Discharge coefficient Cd / 1", 0.30, 1.00, 0.65, 0.01
                )
                if topology == "opposed":
                    extra["intake_height_mm"] = st.slider(
                        "Intake port height / mm", 1.0, 10.0, 9.6, 0.1
                    )
                    extra["exhaust_height_mm"] = st.slider(
                        "Exhaust port height / mm", 1.0, 14.0, 12.0, 0.1
                    )
                else:
                    extra["valve_diameter_mm"] = st.slider(
                        "Valve diameter / mm", 20.0, 45.0, 30.0, 1.0
                    )
                    extra["valve_lift_mm"] = st.slider(
                        "Maximum valve lift / mm", 2.0, 12.0, 8.0, 0.5
                    )
                    for key, title, low, high, default in (
                        ("intake_open_deg", "Intake opens", 300, 400, 350),
                        ("intake_close_deg", "Intake closes", 540, 650, 580),
                        ("exhaust_open_deg", "Exhaust opens", 110, 180, 140),
                        ("exhaust_close_deg", "Exhaust closes", 340, 420, 370),
                    ):
                        extra[key] = st.slider(
                            f"{title} / deg CA", low, high, default, 5
                        )
        run = st.button(
            "Run Simulation",
            type="primary",
            width="stretch",
            key="run_simulation",
        )
        st.caption("Radau integration · 0–720° · 1,441 samples")
        with st.expander("Boundary conditions & assumptions"):
            st.write(
                "Rod length = 1.75 × piston stroke. Homogeneous gas; "
                "no spatial resolution or reacting chemistry."
            )
            if topology != "closed":
                st.write("AFR = 30. Initial charge = AFR × fuel, at 800 K.")
                st.write(
                    "Intake: "
                    + ("1.30" if topology == "opposed" else "1.05")
                    + " bar / 320 K. Exhaust: 1.00 bar / 700 K."
                )
                st.write(
                    "Wall: 450 K. LHV: 42.5 MJ/kg. Burn efficiency: 0.98."
                )
                st.write("Fuel supplies heat; injection mass is neglected.")
                if fluid != "Ideal gas":
                    st.write(
                        "Nozzle flow uses a locally ideal approximation. "
                        "Out-of-range CoolProp states are rejected."
                    )
        return EngineSettings(
            bore_mm=bore,
            stroke_mm=stroke,
            rpm=rpm,
            compression_ratio=ratio,
            **extra,
        ), run


def sweep_panel(payload: dict[str, Any]) -> dict[str, Any] | None:
    """Run independent cases from the last successful simulation baseline.

    Args:
        payload: Successful baseline run with its immutable input snapshot.

    Returns:
        Matching sweep report, or None when no current report is available.
    """
    st.subheader("Parametric study")
    st.caption(
        "Each point solves a new 0–720° transient. Other parameters "
        "remain fixed at the last successful run. Maximum 12 cases."
    )
    with st.form("sweep_controls"):
        columns = st.columns([2, 2, 1])
        label = columns[0].selectbox(
            "Sweep parameter", ["Engine speed", "Compression ratio"]
        )
        parameter = "rpm" if label == "Engine speed" else "compression_ratio"
        # Both ranges remain visible because forms submit changes together.
        rpm_bounds = columns[1].slider(
            "Speed range / RPM", 1000, 8000, (1000, 5000), 100
        )
        ratio_bounds = columns[1].slider(
            "Compression ratio range / 1", 8.0, 24.0, (12.0, 20.0), 0.5
        )
        count = columns[2].number_input("Number of cases", 2, 12, 5)
        submitted = st.form_submit_button("Run Parameter Sweep")
    bounds = rpm_bounds if parameter == "rpm" else ratio_bounds
    signature = dict(
        baseline=payload["settings"],
        parameter=parameter,
        bounds=list(bounds),
        count=int(count),
    )
    if submitted:
        try:
            with st.spinner(f"Solving {count} independent operating points…"):
                table = run_parameter_sweep(
                    EngineSettings(**payload["settings"]),
                    parameter,
                    np.linspace(*bounds, int(count)).tolist(),
                )
            st.session_state["engine_sweep"] = dict(
                signature=signature,
                table=table,
            )
        except (ValueError, RuntimeError, ArithmeticError) as error:
            st.error(f"Sweep could not start: {error}")
        except Exception:
            LOGGER.exception("Unexpected parameter sweep failure")
            st.error(
                "An unexpected sweep error occurred. Previous results "
                "are preserved; check the server log for details."
            )
    sweep = st.session_state.get("engine_sweep")
    if sweep is None:
        st.info("Choose a parameter range and run a study to compare cases.")
        return None
    if sweep["signature"] != signature:
        st.warning(
            "Sweep inputs or baseline changed. Run the study again "
            "to update the comparison. Previous sweep is hidden."
        )
        return None
    table = sweep["table"]
    failed = table.error.ne("")
    if failed.any():
        st.warning(
            f"{int(failed.sum())} case(s) failed. Failed outputs are "
            "left empty and are not interpolated."
        )
    st.plotly_chart(
        sweep_figure(table, parameter),
        theme=None,
        width="stretch",
        key="sweep_plot",
    )
    st.dataframe(table, hide_index=True, width="stretch")
    metadata = dict(
        signature=signature, cases=json.loads(table.to_json(orient="records"))
    )
    left, right = st.columns(2)
    left.download_button(
        "Download sweep CSV",
        table.to_csv(index=False).encode("utf-8"),
        "engine_sweep.csv",
        "text/csv",
        width="stretch",
        on_click="ignore",
    )
    right.download_button(
        "Download sweep JSON",
        json.dumps(metadata, indent=2, allow_nan=False),
        "engine_sweep.json",
        "application/json",
        width="stretch",
        on_click="ignore",
    )
    return metadata


def render_results(payload: dict[str, Any]) -> None:
    """Render KPIs, three analysis tabs and snapshot-consistent downloads."""
    frame, metrics = payload["frame"], payload["metrics"]
    cards = st.columns(4)
    cards[0].metric("Indicated work", f"{payload['work_J']:.2f} J")
    efficiency = payload["thermal_efficiency"]
    cards[1].metric(
        "Thermal efficiency",
        f"{efficiency:.2%}" if np.isfinite(efficiency) else "N/A",
    )
    cards[2].metric(
        "Peak cylinder pressure", f"{frame.pressure_Pa.max() / 1e5:.2f} bar"
    )
    closed = payload["settings"]["topology"] == "closed"
    if closed:
        cards[3].metric(
            "Peak gas temperature", f"{frame.temperature_K.max():.0f} K"
        )
    else:
        cards[3].metric(
            "Fresh-charge trapping efficiency", f"{metrics['eta_tr']:.2%}"
        )
    saved = payload["settings"]
    st.markdown(
        f'<div class="run-context">Completed run &nbsp; | &nbsp; '
        f"{escape(saved['topology'])} / "
        f"{escape(saved['fluid'])} &nbsp; | &nbsp; "
        f"{saved['rpm']:g} RPM &nbsp; | &nbsp; "
        f"{saved['compression_ratio']:g}:1 compression &nbsp; | &nbsp; "
        f"0–720° process coordinate &nbsp; | &nbsp; Peak temperature "
        f"{frame.temperature_K.max():.0f} K</div>",
        unsafe_allow_html=True,
    )
    tabs = st.tabs(
        [
            "Thermodynamic & P-V Analysis",
            "Gas Exchange & Scavenging Dynamics",
            "Parametric Sweep & Data Export",
        ]
    )
    with tabs[0]:
        left, right = st.columns(2)
        with left:
            st.subheader("Pressure–volume diagram")
            logarithmic = st.checkbox("Logarithmic P–V axes", key="log_pv")
            st.plotly_chart(
                thermodynamic_figure(frame, log_axes=logarithmic),
                theme=None,
                key="pv",
                width="stretch",
            )
        with right:
            st.subheader("Temperature–crank angle")
            st.caption("Phase-resolved crank / process trajectory")
            st.plotly_chart(
                thermodynamic_figure(frame, True),
                theme=None,
                key="temperature",
                width="stretch",
            )
        st.caption(
            "Hover values use SI units. Phase priority: combustion, "
            "open-port scavenging, then expansion or compression."
        )
        with st.expander("Interpretation & initial state"):
            st.write(payload["model"])
            st.write(
                "Work and net heat are integrated alongside the ODE. "
                "Engine efficiency uses raw fuel LHV input; closed Otto "
                "efficiency uses the prescribed isochoric heat input. "
                "Adiabatic gas-spring efficiency is undefined (N/A)."
            )
            if closed:
                st.write(
                    "Otto uses compression, constant-volume heat addition, "
                    "expansion and constant-volume heat rejection. "
                    "The process coordinate includes volume dwells."
                )
            st.json(payload["initial_state"])
    with tabs[1]:
        if payload["settings"]["topology"] == "closed":
            st.subheader("Closed control-volume energy audit")
            st.info(
                "No valves or ports are attached. Mass is constant; "
                "delivery, trapping and scavenging metrics are not "
                "applicable. The first law is ΔU = Q − W + H."
            )
            st.plotly_chart(
                energy_figure(frame),
                theme=None,
                key="energy_audit",
                width="stretch",
            )
            st.caption(
                f"Relative first-law residual: "
                f"{payload['energy_error']:.3e}. Heat is positive inward; "
                "work is positive outward."
            )
        else:
            four_stroke = payload["settings"]["topology"] == "four_stroke"
            st.subheader(
                "Valve lift & flow area"
                if four_stroke
                else "Piston motion & port uncovering"
            )
            choice = st.radio(
                "Flow area convention",
                ["Geometric area A_geo", "Effective area Cd × A_geo"],
                horizontal=True,
            )
            st.plotly_chart(
                area_figure(frame, choice.startswith("Effective")),
                theme=None,
                key="area",
                width="stretch",
            )
            if four_stroke:
                st.plotly_chart(
                    valve_lift_figure(frame),
                    theme=None,
                    key="valve_lift",
                    width="stretch",
                )
                st.caption(
                    "Valve area is limited by curtain and throat areas. "
                    "The cam law repeats every 720°."
                )
            else:
                st.caption(
                    "Geometric area follows piston overlap. Effective "
                    "area applies the discharge coefficient once."
                )
            st.subheader("Transient mass exchange")
            st.plotly_chart(
                flow_figure(frame), theme=None, key="flow", width="stretch"
            )
            st.caption(
                "Intake is inward-positive; exhaust is outward-positive. "
                "Negative values indicate reverse flow."
            )
            st.subheader(
                "Fresh-charge retention"
                if four_stroke
                else "Single-sweep performance"
            )
            start, end = payload["scavenging_interval_deg"]
            st.caption(
                f"Assessment interval: {start:.2f}–{end:.2f}° CA. "
                "Reference density uses the selected EOS at intake "
                "conditions. Swept volume follows the mechanism."
            )
            cards = st.columns(3)
            cards[0].metric("Delivery ratio DR", f"{metrics['DR']:.3f}")
            cards[1].metric("Trapping efficiency", f"{metrics['eta_tr']:.2%}")
            cards[2].metric(
                "Fresh-charge purity"
                if four_stroke
                else "Scavenging efficiency",
                f"{metrics['eta_sc']:.2%}",
            )
            st.caption(
                f"Delivered fresh gas: {metrics['m_del'] * 1e6:.2f} mg. "
                f"Trapped fresh gas: {payload['m_trap_kg'] * 1e6:.2f} mg. "
                f"Residual charge: {payload['m_res_kg'] * 1e6:.2f} mg."
            )
            st.info(
                "Single-zone, perfect-mixing tracer estimate. Preexisting "
                "charge and exhaust backflow count as residual; only fresh "
                "intake from this exchange interval is tagged."
            )
    with tabs[2]:
        sweep = sweep_panel(payload)
        st.divider()
        st.subheader("Full trajectory & report")
        st.caption(
            "Exports contain the last successful run. Columns include "
            "SI units; theta_deg is provided for readability. JSON "
            "also includes inputs, assumptions and any matching sweep."
        )
        export = dict(payload)
        if sweep is not None:
            export["parametric_sweep"] = sweep
        left, right = st.columns(2)
        left.download_button(
            "Download trajectory CSV",
            frame.to_csv(index=False).encode("utf-8"),
            "engine_cycle.csv",
            "text/csv",
            width="stretch",
            on_click="ignore",
        )
        right.download_button(
            "Download simulation JSON",
            json_export(export),
            "engine_report.json",
            "application/json",
            width="stretch",
            on_click="ignore",
        )
        st.dataframe(frame, hide_index=True, width="stretch", height=480)
        with st.expander("Configuration & numerical diagnostics"):
            st.json(payload["settings"])
            st.json(payload["initial_state"])
            st.write(f"Maximum EOS residual: {payload['eos_error']:.3e}")
            st.write(
                f"Relative first-law residual: {payload['energy_error']:.3e}"
            )
            st.write("Radau / rtol = 1e-7 / output spacing = 0.5 deg CA")


def main() -> None:
    """Run on request, preserve successful results and surface safe errors."""
    st.set_page_config(
        page_title="0D Thermo-Fluid & Engine Framework", layout="wide"
    )
    st.markdown(CSS, unsafe_allow_html=True)
    settings, run = control_panel()
    render_header()
    if run:
        try:
            with st.spinner("Integrating cylinder states and gas exchange…"):
                candidate = cached_simulation(settings)
            st.session_state["engine_result"] = candidate
        except (ValueError, RuntimeError, ArithmeticError) as error:
            st.error(
                f"Simulation could not complete: {error} "
                "Check the parameters and run again."
            )
        except Exception:
            LOGGER.exception("Unexpected engine simulation failure")
            st.error(
                "An unexpected simulation error occurred. Previous "
                "Results are preserved; see the server log for details."
            )
    payload = st.session_state.get("engine_result")
    if payload is None:
        st.markdown(
            '<section class="empty-state"><h2>Your simulation workspace '
            "is ready</h2><p>Set the engine geometry, operating point and "
            "gas-exchange models in the control panel. Select "
            "<strong>Run Simulation</strong> to inspect the full "
            "pressure–volume path, mass exchange and scavenging performance."
            "</p><p>After a successful run, compare operating points with a "
            "parametric study or export the complete numerical trajectory."
            "</p></section>",
            unsafe_allow_html=True,
        )
        return
    if payload.get("schema_version") != 2:
        st.info("The model schema changed. Run Simulation to refresh results.")
        return
    if payload["settings"]["topology"] != settings.topology:
        st.info(
            "Topology changed. Configure this model and Run Simulation. "
            "Previous results remain preserved until a successful run."
        )
        return
    if payload["settings"] != asdict(settings):
        st.warning(
            "Parameters changed. Charts and exports still show the "
            "last successful run. Run Simulation to update them."
        )
    render_results(payload)


if __name__ == "__main__":
    main()
