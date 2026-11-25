"""End-to-end engine UI and physical post-processing checks."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.integrate import trapezoid
from streamlit.testing.v1 import AppTest

from ui.engine_dashboard import (
    EngineSettings,
    area_figure,
    json_export,
    run_parameter_sweep,
    simulate,
    sweep_figure,
    thermodynamic_figure,
)

APP = Path(__file__).resolve().parents[1] / "src/ui/app.py"


def test_metrics_and_exports() -> None:
    """Check raw-LHV denominator, species bounds and SI export round trip."""
    payload = simulate(EngineSettings())
    frame = payload["frame"]
    assert payload["fuel_energy_J"] == pytest.approx(2 * 20e-6 * 42.5e6)
    assert payload["work_J"] == pytest.approx(
        trapezoid(frame.pressure_Pa, frame.volume_m3), rel=5e-5
    )
    # Solver-integrated work avoids the coarse output-grid quadrature error.
    assert payload["work_J"] == frame.work_J.iloc[-1]
    assert payload["energy_error"] < 1e-6
    assert 0 <= payload["metrics"]["eta_tr"] <= 1
    assert 0 <= payload["metrics"]["eta_sc"] <= 1
    assert payload["scavenging_interval_deg"][0] > 360
    assert payload["scavenging_interval_deg"][1] < 720
    report = json.loads(json_export(payload))
    assert len(report["samples"]) == 1441
    assert report["settings"]["rpm"] == 2000
    assert report["samples"][0]["mass_kg"] == pytest.approx(30 * 20e-6)
    geometric = area_figure(frame)
    effective = area_figure(frame, True)
    np.testing.assert_allclose(
        effective.data[1].y, np.asarray(geometric.data[1].y) * 0.65
    )
    assert geometric.layout.yaxis2.tickmode == "sync"


def test_simulation_controls_and_stale_results() -> None:
    """Exercise run, three views, changed inputs and persistent exports."""
    app = AppTest.from_file(str(APP), default_timeout=30).run()
    assert not app.exception
    assert not app.metric
    app.button(key="run_simulation").click().run()
    assert not app.exception
    assert len(app.tabs) == 3
    assert len(app.metric) == 7
    assert len(app.dataframe[0].value) == 1441
    next(w for w in app.slider if w.label == "Engine speed / RPM").set_value(
        3000
    ).run()
    assert not app.exception
    assert any("Parameters changed" in item.value for item in app.warning)
    assert app.session_state["engine_result"]["settings"]["rpm"] == 2000
    app.button(key="run_simulation").click().run()
    assert not app.exception
    assert app.session_state["engine_result"]["settings"]["rpm"] == 3000


def test_failed_run_preserves_previous_result(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An integration error shows a friendly message without losing data."""
    import ui.engine_dashboard as dashboard

    app = AppTest.from_file(str(APP), default_timeout=30).run()
    app.button(key="run_simulation").click().run()

    def fail(settings: EngineSettings) -> dict:
        """Emulate a failed physical solve for a new parameter snapshot."""
        raise ValueError("invalid physical state")

    monkeypatch.setattr(dashboard, "simulate", fail)
    next(w for w in app.slider if w.label == "Engine speed / RPM").set_value(
        3100
    ).run()
    app.button(key="run_simulation").click().run()
    assert not app.exception
    assert any(
        "Simulation could not complete" in item.value for item in app.error
    )
    assert app.session_state["engine_result"]["settings"]["rpm"] == 2000


def test_sweep_retains_failed_case(monkeypatch: pytest.MonkeyPatch) -> None:
    """Failed points remain explicit gaps, with unchanged baseline inputs."""
    import ui.engine_dashboard as dashboard

    calls = []

    def solve_case(settings: EngineSettings) -> dict:
        """Record requested cases and emulate one rejected operating point."""
        calls.append(settings)
        if settings.rpm == 2000:
            raise ValueError("Rejected case")
        return {
            "thermal_efficiency": 0.3,
            "work_J": 400,
            "frame": pd.DataFrame({"pressure_Pa": [1e6]}),
            "metrics": {"eta_tr": 0.6},
        }

    monkeypatch.setattr(dashboard, "simulate", solve_case)
    baseline = EngineSettings(compression_ratio=18)
    table = run_parameter_sweep(baseline, "rpm", [1000, 2000, 3000])
    assert [item.rpm for item in calls] == [1000, 2000, 3000]
    assert all(item.compression_ratio == 18 for item in calls)
    assert baseline.rpm == 2000
    assert np.isnan(table.thermal_efficiency.iloc[1])
    assert table.error.iloc[1] == "Rejected case"
    assert not sweep_figure(table, "rpm").data[0].connectgaps
    with pytest.raises(ValueError):
        run_parameter_sweep(baseline, "unknown", [1000, 2000])


def test_english_figures_and_compression_ratio() -> None:
    """Verify new geometry input, English phases and raw SI hover metadata."""
    result = simulate(EngineSettings(compression_ratio=20))
    frame = result["frame"]
    assert frame.volume_m3.max() / frame.volume_m3.min() == pytest.approx(20)
    assert set(frame.phase) == {
        "Compression",
        "Combustion",
        "Expansion",
        "Scavenging",
    }
    figure = thermodynamic_figure(frame, log_axes=True)
    assert figure.layout.xaxis.type == "log"
    assert figure.layout.yaxis.type == "log"
    assert figure.data[0].customdata[0][2] == frame.pressure_Pa.iloc[0]
    assert " Pa" in figure.data[0].hovertemplate


def test_sweep_submission_and_baseline_invalidation() -> None:
    """Run a two-case study, then invalidate it when the baseline changes."""
    app = AppTest.from_file(str(APP), default_timeout=30).run()
    app.button(key="run_simulation").click().run()
    next(
        w for w in app.number_input if w.label == "Number of cases"
    ).set_value(2)
    next(
        w for w in app.button if w.label == "Run Parameter Sweep"
    ).click().run()
    assert not app.exception
    assert len(app.session_state["engine_sweep"]["table"]) == 2
    assert len(app.dataframe[0].value) == 2
    assert len(app.dataframe[1].value) == 1441
    next(w for w in app.slider if w.label == "Engine speed / RPM").set_value(
        3000
    )
    app.button(key="run_simulation").click().run()
    assert not app.exception
    assert any("baseline changed" in item.value for item in app.warning)
    assert len(app.dataframe[0].value) == 1441


def test_topology_switch_changes_active_components() -> None:
    """Switch real UI presets without reusing another topology's result."""
    app = AppTest.from_file(str(APP), default_timeout=30).run()
    app.button(key="run_simulation").click().run()
    app.selectbox(key="system_topology").select(
        "Four-stroke poppet valves"
    ).run()
    assert not app.exception
    labels = {widget.label for widget in app.slider}
    assert "Maximum valve lift / mm" in labels
    assert "Intake port height / mm" not in labels
    assert not app.metric
    app.button(key="run_simulation").click().run()
    assert not app.exception
    assert (
        app.session_state["engine_result"]["settings"]["topology"]
        == "four_stroke"
    )
    app.selectbox(key="system_topology").select(
        "Closed thermodynamic cycle"
    ).run()
    assert not app.exception
    labels = {widget.label for widget in app.slider}
    assert "Maximum valve lift / mm" not in labels
    assert "Wiebe start / deg CA" not in labels
    assert any(w.label == "Isochoric heat input / J" for w in app.number_input)
    app.button(key="run_simulation").click().run()
    assert not app.exception
    assert len(app.metric) == 4
    payload = app.session_state["engine_result"]
    expected = 1 - payload["settings"]["compression_ratio"] ** (-0.35)
    assert payload["thermal_efficiency"] == pytest.approx(expected, rel=2e-6)
    report = json.loads(json_export(payload))
    assert report["metrics"]["eta_tr"] is None
    assert report["settings"]["topology"] == "closed"
