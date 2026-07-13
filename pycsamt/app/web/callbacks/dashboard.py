# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/dashboard.py — Home-page survey dashboard.

Populates:
  - Per-line summary table
  - Data-quality indicators
  - Station analytics panel (tipper donut, freq bar, quality radar)
  - Profile page station dropdown
  - Map/Profile page "back to dashboard" navigation
"""

from __future__ import annotations

from dash import (
    Input,
    Output,
    State,
    clientside_callback,
    html,
    no_update,
)

from pycsamt.app.web.layout import IDs


def register_dashboard(app) -> None:

    # ── 1. Populate dashboard when store data changes ─────────────────────────
    @app.callback(
        Output("dash-empty-state", "style"),
        Output("dash-data-state", "style"),
        Output(IDs.DASH_QUALITY, "style"),
        Output("dash-launch-row", "style"),
        Output(IDs.DASH_LINE_TABLE, "children"),
        Output(IDs.DASH_QUALITY, "children"),
        Output("dash-survey-name", "children"),
        Output("dash-health-badge", "children"),
        Output("dash-health-badge", "className"),
        Input(IDs.STORE_DATA, "data"),
    )
    def update_dashboard(store_data):
        _hide = {"display": "none"}
        _show = {"display": "block"}
        _show_flex = {
            "display": "flex",
            "flexDirection": "column",
            "gap": "4px",
            "paddingTop": "6px",
        }

        if not store_data or not store_data.get("n_stations"):
            return (
                _show,
                _hide,
                _hide,
                _hide,
                [],
                [],
                "",
                "",
                "dash-health-badge",
            )

        n_sta = store_data.get("n_stations", 0)
        n_lines = store_data.get("n_lines", 1)
        data_dir = store_data.get("data_dir", "")
        line_counts: dict = store_data.get("line_counts", {})
        records = store_data.get("station_records", [])

        # ── Per-line table ────────────────────────────────────────────────
        if line_counts and n_lines > 1:
            header = html.Div(
                [
                    html.Span("Line", className="dlt-col dlt-head"),
                    html.Span(
                        "Stations", className="dlt-col dlt-head dlt-num"
                    ),
                    html.Span(
                        "Coverage", className="dlt-col dlt-head dlt-bar"
                    ),
                ],
                className="dlt-row dlt-header",
            )
            rows = [header]
            for line_name, cnt in line_counts.items():
                pct = int(cnt / n_sta * 100)
                rows.append(
                    html.Div(
                        [
                            html.Span(
                                str(line_name),
                                className="dlt-col dlt-name",
                                title=str(line_name),
                            ),
                            html.Span(str(cnt), className="dlt-col dlt-num"),
                            html.Div(
                                html.Div(
                                    className="dlt-fill",
                                    style={"width": f"{pct}%"},
                                ),
                                className="dlt-bar-track",
                            ),
                        ],
                        className="dlt-row",
                    )
                )
            line_table_content = html.Div(rows, className="dlt-table")
        else:
            line_table_content = html.Span(
                f"{n_sta} station(s) loaded from a single folder.",
                className="text-muted small",
            )

        # ── Data quality indicators ───────────────────────────────────────
        n_tipper = sum(1 for r in records if r.get("Tipper"))
        n_valid = sum(
            1 for r in records if r.get("Latitude") and r.get("Longitude")
        )
        n_nfreq = [r.get("N_freq", 0) for r in records if r.get("N_freq")]
        avg_freq = int(sum(n_nfreq) / len(n_nfreq)) if n_nfreq else 0

        def _qc_chip(icon: str, text: str, ok: bool) -> html.Div:
            col = "text-success" if ok else "text-warning"
            return html.Div(
                [
                    html.I(
                        className=f"bi bi-{'check-circle' if ok else 'exclamation-triangle'} me-1 {col}"
                    ),
                    html.Span(text, style={"fontSize": "11px"}),
                ],
                className="dash-qc-chip",
            )

        quality_content = [
            _qc_chip(
                "geo", f"{n_valid}/{n_sta} with coordinates", n_valid == n_sta
            ),
            _qc_chip(
                "freq", f"avg {avg_freq} frequencies/station", avg_freq >= 30
            ),
            _qc_chip("tip", f"{n_tipper} stations with tipper", n_tipper > 0),
        ]

        # ── Survey name badge ─────────────────────────────────────────────
        import os

        dir_name = (
            os.path.basename(data_dir.rstrip("/\\")) if data_dir else "survey"
        )
        survey_badge = dir_name or "survey"

        # ── Health badge ──────────────────────────────────────────────────
        all_ok = (n_valid == n_sta) and (avg_freq >= 30)
        badge_text = (
            "✓ All checks passed" if all_ok else "⚠ Review data quality"
        )
        badge_cls = (
            "dash-health-badge good" if all_ok else "dash-health-badge warn"
        )

        return (
            _hide,
            _show,
            _show,
            _show_flex,
            line_table_content,
            quality_content,
            survey_badge,
            badge_text,
            badge_cls,
        )

    # ── 2. Profile page station dropdown — populate options ───────────────────
    @app.callback(
        Output(IDs.PROFILE_PAGE_STATION, "options"),
        Output(IDs.PROFILE_PAGE_STATION, "value"),
        Input(IDs.STORE_DATA, "data"),
    )
    def populate_station_dropdown(store_data):
        if not store_data:
            return [], None
        records = store_data.get("station_records", [])
        opts = [{"label": r["ID"], "value": r["ID"]} for r in records]
        first = opts[0]["value"] if opts else None
        return opts, first

    # ── 3. Profile page station dropdown → write STORE_SELECTION ─────────────
    @app.callback(
        Output(IDs.STORE_SELECTION, "data", allow_duplicate=True),
        Input(IDs.PROFILE_PAGE_STATION, "value"),
        prevent_initial_call=True,
    )
    def profile_station_selected(station_id):
        if not station_id:
            return no_update
        return {"station_id": str(station_id)}

    # ── 4. STORE_SELECTION → sync profile page dropdown ──────────────────────
    @app.callback(
        Output(IDs.PROFILE_PAGE_STATION, "value", allow_duplicate=True),
        Input(IDs.STORE_SELECTION, "data"),
        State(IDs.PROFILE_PAGE_STATION, "options"),
        prevent_initial_call=True,
    )
    def sync_profile_dropdown(selection, opts):
        if not selection or not opts:
            return no_update
        sid = selection.get("station_id", "")
        valid = {o["value"] for o in opts}
        return sid if sid in valid else no_update

    # ── 4b. Profile page: station info header bar ────────────────────────────
    @app.callback(
        Output(IDs.PROF_INFO_BAR, "children"),
        Input(IDs.STORE_SELECTION, "data"),
        Input(IDs.STORE_DATA, "data"),
        prevent_initial_call=False,
    )
    def update_prof_info_bar(selection, store_data):
        from dash import html as _html

        def _placeholder():
            return [
                _html.Span(
                    [
                        _html.I(
                            className="bi bi-broadcast me-1",
                            style={"fontSize": "10px"},
                        ),
                        "No station selected — choose from the dropdown or click the table",
                    ],
                    className="prof-info-station",
                    style={"fontWeight": "400", "color": "var(--sub0)"},
                ),
            ]

        if not store_data:
            return _placeholder()

        sid = (selection or {}).get("station_id", "")
        records = store_data.get("station_records", [])
        match = next(
            (r for r in records if str(r.get("ID")) == str(sid)), None
        )

        if not match:
            if sid:
                return [
                    _html.Span(
                        [
                            _html.I(
                                className="bi bi-geo-alt-fill me-1",
                                style={"fontSize": "10px"},
                            ),
                            sid,
                        ],
                        className="prof-info-station",
                    ),
                ]
            return _placeholder()

        lat = match.get("Latitude", None)
        lon = match.get("Longitude", None)
        elev = match.get("Elevation", None)
        nf = match.get("N_freq", "—")
        tip = bool(match.get("Tipper"))
        line = match.get("Line", None)

        lat_str = f"{lat:.4f}°N" if isinstance(lat, (int, float)) else "—"
        lon_str = f"{lon:.4f}°E" if isinstance(lon, (int, float)) else "—"
        elev_str = f"{elev:.0f} m" if isinstance(elev, (int, float)) else ""

        children = [
            # Station name (bold)
            _html.Span(
                [
                    _html.I(
                        className="bi bi-geo-alt-fill me-1",
                        style={"fontSize": "10px", "color": "var(--blue)"},
                    ),
                    sid,
                ],
                className="prof-info-station",
            ),
            _html.Span("·", className="prof-info-sep"),
            # Coordinates
            _html.Span(
                f"{lat_str}, {lon_str}"
                + (f"  {elev_str}" if elev_str else ""),
                className="prof-info-coords",
            ),
            _html.Span("·", className="prof-info-sep"),
            # Frequency count badge
            _html.Span(
                [_html.I(className="bi bi-reception-4 me-1"), f"{nf} freq"],
                className="prof-info-badge prof-badge-freq",
            ),
            # Tipper badge
            _html.Span(
                (
                    [
                        _html.I(className="bi bi-check-circle-fill me-1"),
                        "Tipper",
                    ]
                    if tip
                    else [
                        _html.I(className="bi bi-dash-circle me-1"),
                        "No tipper",
                    ]
                ),
                className=f"prof-info-badge {'prof-badge-tipper' if tip else 'prof-badge-no-tipper'}",
            ),
        ]

        if line is not None:
            children.insert(2, _html.Span("·", className="prof-info-sep"))
            children.insert(
                3, _html.Span(f"Line {line}", className="prof-info-coords")
            )

        return children

    # ── 5. Map page: update station info card on map click ────────────────────
    @app.callback(
        Output(IDs.MAP_PAGE_INFO, "children"),
        Input(IDs.STORE_SELECTION, "data"),
        State(IDs.STORE_DATA, "data"),
    )
    def update_map_station_info(selection, store_data):
        if not selection:
            return html.Span(
                "Click a station on the map", className="text-muted small"
            )
        sid = selection.get("station_id", "")
        if not sid or not store_data:
            return html.Span(
                "Click a station on the map", className="text-muted small"
            )
        records = store_data.get("station_records", [])
        match = next(
            (r for r in records if str(r.get("ID")) == str(sid)), None
        )
        if not match:
            return html.Span(sid, className="fw-semibold small")
        lat = match.get("Latitude", "—")
        lon = match.get("Longitude", "—")
        elev = match.get("Elevation", "—")
        nf = match.get("N_freq", "—")
        tip = "Yes" if match.get("Tipper") else "No"
        line = match.get("Line", "—")
        return html.Div(
            [
                html.Span(sid, className="map-info-name"),
                html.Div(
                    [
                        html.Span("Line: ", className="map-info-key"),
                        html.Span(str(line)),
                    ],
                    className="map-info-row",
                ),
                html.Div(
                    [
                        html.Span("Lat: ", className="map-info-key"),
                        html.Span(
                            f"{lat:.4f}°"
                            if isinstance(lat, float)
                            else str(lat)
                        ),
                        html.Span("  Lon: ", className="map-info-key"),
                        html.Span(
                            f"{lon:.4f}°"
                            if isinstance(lon, float)
                            else str(lon)
                        ),
                    ],
                    className="map-info-row",
                ),
                html.Div(
                    [
                        html.Span("Elev: ", className="map-info-key"),
                        html.Span(
                            f"{elev:.1f} m"
                            if isinstance(elev, float)
                            else str(elev)
                        ),
                        html.Span("  N freq: ", className="map-info-key"),
                        html.Span(str(nf)),
                    ],
                    className="map-info-row",
                ),
                html.Div(
                    [
                        html.Span("Tipper: ", className="map-info-key"),
                        html.Span(tip),
                    ],
                    className="map-info-row",
                ),
            ]
        )

    # ── 6. Map page: populate survey-line filter checklist ───────────────────
    @app.callback(
        Output(IDs.MAP_PAGE_LINE_FILTER, "options"),
        Output(IDs.MAP_PAGE_LINE_FILTER, "value"),
        Output("map-line-filter-count", "children"),
        Input(IDs.STORE_DATA, "data"),
    )
    def populate_map_line_filter(store_data):
        if not store_data:
            return [], [], "(all)"
        line_counts: dict = store_data.get("line_counts", {})
        if not line_counts:
            return [], [], "(all)"
        opts = [
            {"label": f"{ln} ({cnt})", "value": ln}
            for ln, cnt in line_counts.items()
        ]
        values = [o["value"] for o in opts]
        count_label = f"({len(values)})" if values else "(all)"
        return opts, values, count_label

    # ── 7. Map page: stats strip ──────────────────────────────────────────────
    @app.callback(
        Output(IDs.MAP_PAGE_STATS, "children"),
        Input(IDs.STORE_DATA, "data"),
        Input(IDs.MAP_PAGE_LINE_FILTER, "value"),
    )
    def update_map_stats(store_data, line_filter):
        if not store_data:
            return ""
        records = store_data.get("station_records", [])
        if line_filter:
            records = [r for r in records if r.get("Line") in line_filter]
        n_sta = len(records)
        n_lines = len({r.get("Line") for r in records if r.get("Line")})
        n_tip = sum(1 for r in records if r.get("Tipper"))
        n_geo = sum(
            1 for r in records if r.get("Latitude") and r.get("Longitude")
        )

        def _stat(label, value):
            return html.Div(
                [
                    html.Span(str(value), className="map-stat-value"),
                    html.Span(label, className="map-stat-label"),
                ],
                className="map-stat-item",
            )

        return [
            _stat("stations", n_sta),
            _stat("lines", n_lines),
            _stat("with tipper", n_tip),
            _stat("with coordinates", n_geo),
        ]

    # ── 9. Analytics: 6 survey-level charts (slides 1 2 4 5 6 7) ────────────
    @app.callback(
        Output(IDs.DASH_CHART_DONUT, "figure"),
        Output(IDs.DASH_CHART_FREQ, "figure"),
        Output(IDs.DASH_CHART_GEO, "figure"),
        Output(IDs.DASH_CHART_RANK, "figure"),
        Output(IDs.DASH_CHART_HEATMAP, "figure"),
        Output(IDs.DASH_CHART_VIOLIN, "figure"),
        Input(IDs.STORE_DATA, "data"),
        Input(IDs.STORE_THEME, "data"),
    )
    def update_analytics_charts(store_data, theme):
        dark = (theme or "dark") == "dark"
        fc = "#cdd6f4" if dark else "#4c4f69"
        sc = "#a6adc8" if dark else "#7c7f93"  # subdued text
        gc = "rgba(205,214,244,0.09)" if dark else "rgba(76,79,105,0.11)"
        bg = "rgba(0,0,0,0)"
        sf = (
            "rgba(49,50,68,0.55)" if dark else "rgba(220,224,232,0.55)"
        )  # surface fill
        _e = _empty_analytics_fig(fc, bg)

        records = (store_data or {}).get("station_records", [])
        if not records:
            return _e, _e, _e, _e, _e, _e

        # ── Line palette (Catppuccin accent colours) ───────────────────────
        _LINE_PALETTE = [
            "#89b4fa",
            "#a6e3a1",
            "#fab387",
            "#f38ba8",
            "#cba6f7",
            "#f9e2af",
            "#94e2d5",
            "#89dceb",
        ]

        # Detect multi-line survey
        raw_lines = [str(r.get("Line", "") or "").strip() for r in records]
        line_names = list(dict.fromkeys(l for l in raw_lines if l))
        if not line_names:
            line_names = ["default"]
            raw_lines = ["default"] * len(records)
        is_multi = len(line_names) > 1
        line_color = {
            ln: _LINE_PALETTE[i % len(_LINE_PALETTE)]
            for i, ln in enumerate(line_names)
        }
        dict(
            zip(
                [r.get("ID", f"S{i}") for i, r in enumerate(records)],
                raw_lines,
            )
        )

        n_sta = len(records)
        n_tipper = sum(1 for r in records if r.get("Tipper"))
        n_no_tip = n_sta - n_tipper
        pct_tip = int(n_tipper / n_sta * 100) if n_sta else 0

        def _title(t):
            return {
                "text": t,
                "font": {"size": 10, "color": fc},
                "x": 0.5,
                "xanchor": "center",
                "pad": {"t": 2},
            }

        # ══ 1. Tipper coverage ════════════════════════════════════════════
        # Single-line: donut.  Multi-line: stacked bar per line.
        if is_multi:
            tip_yes = [
                sum(
                    1
                    for r in records
                    if str(r.get("Line", "") or "").strip() == ln
                    and r.get("Tipper")
                )
                for ln in line_names
            ]
            tip_no = [
                sum(
                    1
                    for r in records
                    if str(r.get("Line", "") or "").strip() == ln
                    and not r.get("Tipper")
                )
                for ln in line_names
            ]
            fig_donut = {
                "data": [
                    {
                        "type": "bar",
                        "name": "Tipper",
                        "x": line_names,
                        "y": tip_yes,
                        "marker": {
                            "color": "#89b4fa",
                            "line": {"color": bg, "width": 0},
                        },
                        "hovertemplate": "<b>%{x}</b><br>%{y} stations with tipper<extra></extra>",
                    },
                    {
                        "type": "bar",
                        "name": "No Tipper",
                        "x": line_names,
                        "y": tip_no,
                        "marker": {
                            "color": "#45475a",
                            "line": {"color": bg, "width": 0},
                        },
                        "hovertemplate": "<b>%{x}</b><br>%{y} stations without tipper<extra></extra>",
                    },
                ],
                "layout": {
                    "barmode": "stack",
                    "margin": {"l": 32, "r": 6, "t": 32, "b": 28},
                    "paper_bgcolor": bg,
                    "plot_bgcolor": bg,
                    "font": {"color": fc, "size": 8},
                    "showlegend": True,
                    "legend": {
                        "orientation": "h",
                        "y": -0.14,
                        "x": 0.5,
                        "xanchor": "center",
                        "font": {"size": 8},
                    },
                    "xaxis": {
                        "showgrid": False,
                        "zeroline": False,
                        "tickfont": {"size": 8},
                    },
                    "yaxis": {
                        "showgrid": True,
                        "gridcolor": gc,
                        "zeroline": False,
                        "tickfont": {"size": 7},
                        "title": {
                            "text": "Stations",
                            "font": {"size": 8},
                            "standoff": 4,
                        },
                    },
                    "title": _title(f"Tipper Coverage · {pct_tip}% overall"),
                    "bargap": 0.28,
                },
            }
        else:
            fig_donut = {
                "data": [
                    {
                        "type": "pie",
                        "values": [
                            max(n_tipper, 0.001),
                            max(n_no_tip, 0.001),
                        ],
                        "labels": ["With Tipper", "No Tipper"],
                        "hole": 0.68,
                        "marker": {
                            "colors": ["#89b4fa", "#45475a"],
                            "line": {"color": sf, "width": 2},
                        },
                        "textinfo": "none",
                        "hovertemplate": "<b>%{label}</b><br>%{value} stations · %{percent}<extra></extra>",
                        "direction": "clockwise",
                        "sort": False,
                        "rotation": 90,
                    }
                ],
                "layout": {
                    "margin": {"l": 10, "r": 10, "t": 32, "b": 32},
                    "paper_bgcolor": bg,
                    "plot_bgcolor": bg,
                    "font": {"color": fc, "size": 9},
                    "showlegend": True,
                    "legend": {
                        "orientation": "h",
                        "y": -0.08,
                        "x": 0.5,
                        "xanchor": "center",
                        "font": {"size": 9},
                        "itemwidth": 40,
                        "traceorder": "normal",
                    },
                    "annotations": [
                        {
                            "text": (
                                f"<b style='font-size:18px'>{pct_tip}%</b>"
                                f"<br><span style='font-size:9px;color:{sc}'>tipper</span>"
                            ),
                            "x": 0.5,
                            "y": 0.5,
                            "showarrow": False,
                            "font": {"size": 14, "color": "#89b4fa"},
                            "align": "center",
                            "xref": "paper",
                            "yref": "paper",
                        }
                    ],
                    "title": _title("Tipper Coverage"),
                },
            }

        # ══ 2. Freq / station ═════════════════════════════════════════════
        # Multi-line: one grouped-bar trace per line, colour-coded.
        # Single-line: sorted bar chart with gradient colours (original).
        nfreqs = [r.get("N_freq") or 0 for r in records]
        ids_raw = [r.get("ID", f"S{i}") for i, r in enumerate(records)]
        med_nf = sorted(nfreqs)[len(nfreqs) // 2] if nfreqs else 0

        if is_multi:
            freq_traces = []
            for ln in line_names:
                ln_recs = [
                    (r.get("ID", "?"), r.get("N_freq") or 0)
                    for r in records
                    if str(r.get("Line", "") or "").strip() == ln
                ]
                ln_recs.sort(key=lambda x: x[1])
                freq_traces.append(
                    {
                        "type": "bar",
                        "name": ln,
                        "x": [x[0] for x in ln_recs],
                        "y": [x[1] for x in ln_recs],
                        "marker": {
                            "color": line_color[ln],
                            "opacity": 0.88,
                            "line": {"color": bg, "width": 0},
                        },
                        "hovertemplate": f"<b>%{{x}}</b><br>%{{y}} frequencies · {ln}<extra></extra>",
                    }
                )
            # Overall median line
            all_ids_sorted = [
                r.get("ID", "?")
                for r in sorted(records, key=lambda r: r.get("N_freq") or 0)
            ]
            freq_traces.append(
                {
                    "type": "scatter",
                    "x": [all_ids_sorted[0], all_ids_sorted[-1]],
                    "y": [med_nf, med_nf],
                    "mode": "lines",
                    "line": {"color": "#f38ba8", "width": 1.2, "dash": "dot"},
                    "hoverinfo": "skip",
                    "showlegend": False,
                }
            )
            fig_freq = {
                "data": freq_traces,
                "layout": {
                    "barmode": "group",
                    "margin": {"l": 32, "r": 6, "t": 32, "b": 10},
                    "paper_bgcolor": bg,
                    "plot_bgcolor": bg,
                    "font": {"color": fc, "size": 8},
                    "xaxis": {"visible": False, "showgrid": False},
                    "yaxis": {
                        "showgrid": True,
                        "gridcolor": gc,
                        "zeroline": False,
                        "tickfont": {"size": 7},
                        "title": {
                            "text": "N freqs",
                            "font": {"size": 8},
                            "standoff": 4,
                        },
                    },
                    "legend": {
                        "orientation": "h",
                        "y": -0.12,
                        "x": 0.5,
                        "xanchor": "center",
                        "font": {"size": 8},
                        "traceorder": "normal",
                    },
                    "title": _title(
                        f"Frequencies / Station  ·  median {med_nf}"
                    ),
                    "bargap": 0.10,
                    "bargroupgap": 0.05,
                    "annotations": [
                        {
                            "x": 1,
                            "y": med_nf,
                            "xref": "paper",
                            "yref": "y",
                            "text": f" med={med_nf}",
                            "showarrow": False,
                            "font": {"size": 7, "color": "#f38ba8"},
                            "xanchor": "left",
                        }
                    ],
                },
            }
        else:
            pairs = sorted(zip(nfreqs, ids_raw, records), key=lambda x: x[0])
            nf_s = [p[0] for p in pairs]
            id_s = [p[1] for p in pairs]
            tip_s = [p[2].get("Tipper") for p in pairs]
            max_nf = max(nf_s) or 1

            def _bar_color(val, has_tip, max_v):
                t = val / max_v if max_v else 0
                if has_tip:
                    r2 = int(137 + (166 - 137) * (1 - t))
                    g2 = int(180 + (227 - 180) * (1 - t))
                    b2 = int(250 + (161 - 250) * (1 - t))
                    return f"rgb({r2},{g2},{b2})"
                return f"rgba(88,91,112,{0.45 + 0.5 * t:.2f})"

            bar_col = [_bar_color(v, t, max_nf) for v, t in zip(nf_s, tip_s)]
            fig_freq = {
                "data": [
                    {
                        "type": "bar",
                        "x": id_s,
                        "y": nf_s,
                        "marker": {
                            "color": bar_col,
                            "opacity": 0.92,
                            "line": {"color": bg, "width": 0},
                        },
                        "hovertemplate": "<b>%{x}</b><br>%{y} frequencies<extra></extra>",
                    },
                    {
                        "type": "scatter",
                        "x": [id_s[0], id_s[-1]] if id_s else [],
                        "y": [med_nf, med_nf],
                        "mode": "lines",
                        "line": {
                            "color": "#f38ba8",
                            "width": 1.2,
                            "dash": "dot",
                        },
                        "hoverinfo": "skip",
                        "showlegend": False,
                    },
                ],
                "layout": {
                    "margin": {"l": 32, "r": 6, "t": 32, "b": 14},
                    "paper_bgcolor": bg,
                    "plot_bgcolor": bg,
                    "font": {"color": fc, "size": 8},
                    "xaxis": {
                        "visible": False,
                        "showgrid": False,
                        "zeroline": False,
                    },
                    "yaxis": {
                        "showgrid": True,
                        "gridcolor": gc,
                        "zeroline": False,
                        "tickfont": {"size": 7},
                        "gridwidth": 0.5,
                        "title": {
                            "text": "N freqs",
                            "font": {"size": 8},
                            "standoff": 4,
                        },
                    },
                    "title": _title("Frequencies / Station"),
                    "bargap": 0.12,
                    "annotations": [
                        {
                            "x": 1,
                            "y": med_nf,
                            "xref": "paper",
                            "yref": "y",
                            "text": f" med={med_nf}",
                            "showarrow": False,
                            "font": {"size": 7, "color": "#f38ba8"},
                            "xanchor": "left",
                        }
                    ],
                },
            }

        # ══ 3. Survey Geometry — per-line colours ═════════════════════════
        geo_recs = [
            r for r in records if r.get("Latitude") and r.get("Longitude")
        ]
        if geo_recs:
            g_max = max((r.get("N_freq") or 1) for r in geo_recs) or 1
            geo_data = []
            for ln in line_names:
                ln_geo = [
                    r
                    for r in geo_recs
                    if str(r.get("Line", "") or "").strip() == ln
                ]
                if not ln_geo:
                    continue
                lons_l = [r["Longitude"] for r in ln_geo]
                lats_l = [r["Latitude"] for r in ln_geo]
                g_nf_l = [r.get("N_freq") or 1 for r in ln_geo]
                g_ids_l = [r.get("ID", "?") for r in ln_geo]
                g_sz_l = [max(5, min(18, v / g_max * 14 + 5)) for v in g_nf_l]
                col = line_color[ln]
                # Connecting line within this profile
                geo_data.append(
                    {
                        "type": "scatter",
                        "x": lons_l,
                        "y": lats_l,
                        "mode": "lines",
                        "line": {
                            "color": col.replace(")", ",0.25)").replace(
                                "rgb(", "rgba("
                            )
                            if col.startswith("rgb(")
                            else col + "40",
                            "width": 1.2,
                        },
                        "hoverinfo": "skip",
                        "showlegend": False,
                    }
                )
                # Station markers
                geo_data.append(
                    {
                        "type": "scatter",
                        "x": lons_l,
                        "y": lats_l,
                        "mode": "markers+text"
                        if len(ln_geo) <= 15
                        else "markers",
                        "name": ln,
                        "marker": {
                            "size": g_sz_l,
                            "color": col,
                            "opacity": 0.88,
                            "line": {"color": sf, "width": 1},
                        },
                        "text": [ln] * len(ln_geo)
                        if len(ln_geo) <= 15
                        else [],
                        "textposition": "top center",
                        "textfont": {"size": 7, "color": col},
                        "customdata": g_nf_l,
                        "hovertemplate": (
                            "<b>%{text}</b><br>"
                            "%{y:.4f}°N  %{x:.4f}°E<br>"
                            "%{customdata} freqs<extra>" + ln + "</extra>"
                        ),
                        "text": g_ids_l,
                    }
                )
            geo_layout = {
                "margin": {"l": 36, "r": 6, "t": 32, "b": 28},
                "paper_bgcolor": bg,
                "plot_bgcolor": bg,
                "font": {"color": fc, "size": 8},
                "showlegend": is_multi,
                "legend": {
                    "orientation": "h",
                    "y": -0.14,
                    "x": 0.5,
                    "xanchor": "center",
                    "font": {"size": 8},
                },
                "xaxis": {
                    "showgrid": True,
                    "gridcolor": gc,
                    "zeroline": False,
                    "tickfont": {"size": 7},
                    "title": {
                        "text": "Longitude (°E)",
                        "font": {"size": 8},
                        "standoff": 2,
                    },
                },
                "yaxis": {
                    "showgrid": True,
                    "gridcolor": gc,
                    "zeroline": False,
                    "tickfont": {"size": 7},
                    "title": {
                        "text": "Latitude (°N)",
                        "font": {"size": 8},
                        "standoff": 2,
                    },
                    "scaleanchor": "x",
                    "scaleratio": 1,
                },
                "title": _title(
                    f"Survey Geometry · {len(line_names)} line{'s' if len(line_names) != 1 else ''}"
                ),
            }
        else:
            geo_data = []
            geo_layout = {
                "paper_bgcolor": bg,
                "plot_bgcolor": bg,
                "xaxis": {"visible": False},
                "yaxis": {"visible": False},
                "margin": {"l": 0, "r": 0, "t": 32, "b": 0},
                "title": _title("Survey Geometry"),
                "annotations": [
                    {
                        "text": "No coordinate data available",
                        "x": 0.5,
                        "y": 0.5,
                        "showarrow": False,
                        "font": {"size": 10, "color": sc},
                        "xref": "paper",
                        "yref": "paper",
                    }
                ],
            }
        fig_geo = {"data": geo_data, "layout": geo_layout}

        # ══ 4. Station Ranking — colour-coded by line ═════════════════════
        id_to_line = {
            r.get("ID", f"S{i}"): str(r.get("Line", "") or "").strip()
            for i, r in enumerate(records)
        }
        rank_pairs = sorted(
            zip(nfreqs, ids_raw),
            key=lambda x: x[0],
            reverse=True,
        )
        rk_nf = [p[0] for p in rank_pairs]
        rk_ids = [p[1] for p in rank_pairs]
        n_rk = len(rk_nf)
        if is_multi:
            rk_col = [
                line_color.get(id_to_line.get(sid, ""), "#89b4fa")
                for sid in rk_ids
            ]
        else:
            rk_col = [
                f"rgba({int(137 + (250 - 137) * (1 - i / max(n_rk - 1, 1)))},"
                f"{int(180 + (179 - 180) * (1 - i / max(n_rk - 1, 1)))},"
                f"{int(250 + (135 - 250) * (1 - i / max(n_rk - 1, 1)))},0.88)"
                for i in range(n_rk)
            ]

        fig_rank = {
            "data": [
                {
                    "type": "bar",
                    "orientation": "h",
                    "y": rk_ids,
                    "x": rk_nf,
                    "marker": {
                        "color": rk_col,
                        "line": {"color": bg, "width": 0},
                    },
                    "customdata": [id_to_line.get(s, "") for s in rk_ids],
                    "hovertemplate": "<b>%{y}</b>  ·  %{customdata}<br>%{x} frequencies<extra></extra>",
                }
            ],
            "layout": {
                "margin": {"l": 54, "r": 16, "t": 32, "b": 20},
                "paper_bgcolor": bg,
                "plot_bgcolor": bg,
                "font": {"color": fc, "size": 8},
                "xaxis": {
                    "showgrid": True,
                    "gridcolor": gc,
                    "zeroline": False,
                    "tickfont": {"size": 7},
                    "title": {
                        "text": "N frequencies",
                        "font": {"size": 8},
                        "standoff": 2,
                    },
                },
                "yaxis": {
                    "showgrid": False,
                    "zeroline": False,
                    "tickfont": {"size": max(6, min(8, 120 // max(n_rk, 1)))},
                    "automargin": True,
                },
                "title": _title("Station Ranking by Coverage"),
                "bargap": 0.18,
            },
        }

        # ══ 5. Coverage Heatmap (new slide 6) ═════════════════════════════
        # Rows = survey lines, columns = station index, colour = N_freq.
        # Pads shorter lines with None so the grid stays rectangular.
        line_sta_nf = {}
        for ln in line_names:
            ln_recs = sorted(
                [
                    (r.get("ID", "?"), r.get("N_freq") or 0)
                    for r in records
                    if str(r.get("Line", "") or "").strip() == ln
                ],
                key=lambda x: x[1],
                reverse=True,
            )
            line_sta_nf[ln] = ln_recs

        max_sta = max((len(v) for v in line_sta_nf.values()), default=1)
        heat_z = []
        heat_text = []
        for ln in line_names:
            row = line_sta_nf[ln]
            vals = [x[1] for x in row]
            labels = [x[0] for x in row]
            # Pad to max_sta
            vals += [None] * (max_sta - len(vals))
            labels += [""] * (max_sta - len(labels))
            heat_z.append(vals)
            heat_text.append(
                [
                    f"{sid}<br>{nf} freqs" if sid else ""
                    for sid, nf in zip(labels, vals)
                ]
            )

        fig_heatmap = {
            "data": [
                {
                    "type": "heatmap",
                    "z": heat_z,
                    "x": list(range(max_sta)),
                    "y": line_names,
                    "text": heat_text,
                    "hovertemplate": "%{text}<extra></extra>",
                    "colorscale": [
                        [0, "rgba(69,71,90,0.6)"],
                        [0.3, "#313244"],
                        [0.6, "#89b4fa"],
                        [1, "#cba6f7"],
                    ],
                    "showscale": True,
                    "colorbar": {
                        "thickness": 10,
                        "len": 0.7,
                        "tickfont": {"size": 7, "color": sc},
                        "title": {
                            "text": "N freqs",
                            "font": {"size": 8, "color": fc},
                            "side": "right",
                        },
                    },
                    "zsmooth": "best",
                    "xgap": 1,
                    "ygap": 2,
                }
            ],
            "layout": {
                "margin": {"l": 38, "r": 60, "t": 32, "b": 22},
                "paper_bgcolor": bg,
                "plot_bgcolor": bg,
                "font": {"color": fc, "size": 8},
                "xaxis": {
                    "showgrid": False,
                    "zeroline": False,
                    "title": {
                        "text": "Station rank (high → low coverage)",
                        "font": {"size": 8},
                        "standoff": 2,
                    },
                    "tickfont": {"size": 6},
                },
                "yaxis": {
                    "showgrid": False,
                    "zeroline": False,
                    "tickfont": {"size": 9},
                    "automargin": True,
                },
                "title": _title(
                    "Coverage Heatmap  ·  frequency count per station"
                ),
            },
        }

        # ══ 6. Per-Line Frequency Distribution — box plot (new slide 7) ═══
        violin_traces = []
        for ln in line_names:
            ln_nf = [
                r.get("N_freq") or 0
                for r in records
                if str(r.get("Line", "") or "").strip() == ln
            ]
            if not ln_nf:
                continue
            col = line_color[ln]
            violin_traces.append(
                {
                    "type": "box",
                    "name": ln,
                    "y": ln_nf,
                    "marker": {
                        "color": col,
                        "size": 4,
                        "line": {"color": col, "width": 1},
                    },
                    "line": {"color": col, "width": 1.5},
                    "fillcolor": col + "28",
                    "boxpoints": "outliers",
                    "jitter": 0.3,
                    "whiskerwidth": 0.6,
                    "hovertemplate": (
                        f"<b>{ln}</b><br>"
                        "Q1: %{q1}<br>Median: %{median}<br>"
                        "Q3: %{q3}<extra></extra>"
                    ),
                }
            )

        fig_violin = {
            "data": violin_traces,
            "layout": {
                "margin": {"l": 36, "r": 8, "t": 32, "b": 28},
                "paper_bgcolor": bg,
                "plot_bgcolor": bg,
                "font": {"color": fc, "size": 8},
                "showlegend": False,
                "xaxis": {
                    "showgrid": False,
                    "zeroline": False,
                    "tickfont": {"size": 9},
                },
                "yaxis": {
                    "showgrid": True,
                    "gridcolor": gc,
                    "zeroline": False,
                    "tickfont": {"size": 7},
                    "title": {
                        "text": "N frequencies / station",
                        "font": {"size": 8},
                        "standoff": 4,
                    },
                },
                "title": _title("Per-Line Frequency Distribution"),
                "boxmode": "group",
            },
        }

        return fig_donut, fig_freq, fig_geo, fig_rank, fig_heatmap, fig_violin

    # ── 10. Analytics: quality radar or elevation profile (slide 3) ──────────
    @app.callback(
        Output(IDs.DASH_CHART_RADAR, "figure"),
        Input(IDs.STORE_SELECTION, "data"),
        Input(IDs.STORE_DATA, "data"),
        Input(IDs.STORE_THEME, "data"),
    )
    def update_station_radar(selection, store_data, theme):
        dark = (theme or "dark") == "dark"
        fc = "#cdd6f4" if dark else "#4c4f69"
        gc = "rgba(205,214,244,0.12)" if dark else "rgba(76,79,105,0.14)"
        bg = "rgba(0,0,0,0)"
        sc = "#a6adc8" if dark else "#7c7f93"

        records = (store_data or {}).get("station_records", [])
        if not records:
            return _empty_analytics_fig(fc, bg)

        station_id = (selection or {}).get("station_id")

        # ── No station selected ───────────────────────────────────────────
        if not station_id:
            _LINE_PAL = [
                "#89b4fa",
                "#a6e3a1",
                "#fab387",
                "#f38ba8",
                "#cba6f7",
                "#f9e2af",
                "#94e2d5",
                "#89dceb",
            ]
            raw_lines_r = [
                str(r.get("Line", "") or "").strip() for r in records
            ]
            line_names_r = list(dict.fromkeys(l for l in raw_lines_r if l))
            is_multi_r = len(line_names_r) > 1

            if is_multi_r:
                # Per-line average quality spider
                nfreqs_all = [r.get("N_freq") or 0 for r in records]
                max_nf_all = max(nfreqs_all) or 1
                med_nf_all = sorted(nfreqs_all)[len(nfreqs_all) // 2] or 1

                cats = [
                    "Freq. Range",
                    "Tipper",
                    "Coordinates",
                    "Data Fill",
                    "Elev. Rank",
                ]

                def _ring(r_val, color):
                    return {
                        "type": "scatterpolar",
                        "r": [r_val] * (len(cats) + 1),
                        "theta": cats + [cats[0]],
                        "mode": "lines",
                        "line": {"color": color, "width": 0.7, "dash": "dot"},
                        "hoverinfo": "skip",
                        "showlegend": False,
                    }

                spider_traces = [_ring(0.33, gc), _ring(0.67, gc)]
                for idx, ln in enumerate(line_names_r):
                    ln_recs = [
                        r
                        for r in records
                        if str(r.get("Line", "") or "").strip() == ln
                    ]
                    if not ln_recs:
                        continue
                    col = _LINE_PAL[idx % len(_LINE_PAL)]
                    col_fill = col + "22"

                    def _avg(vals):
                        v = [x for x in vals if x is not None]
                        return sum(v) / len(v) if v else 0

                    fs = _avg(
                        [
                            min((r.get("N_freq") or 0) / max_nf_all, 1)
                            for r in ln_recs
                        ]
                    )
                    ts = _avg(
                        [1.0 if r.get("Tipper") else 0.0 for r in ln_recs]
                    )
                    gs = _avg(
                        [
                            1.0
                            if (r.get("Latitude") and r.get("Longitude"))
                            else 0.0
                            for r in ln_recs
                        ]
                    )
                    ds = _avg(
                        [
                            min((r.get("N_freq") or 0) / med_nf_all, 1)
                            for r in ln_recs
                        ]
                    )

                    elevs_all_s = sorted(
                        r.get("Elevation") or 0 for r in records
                    )
                    n_e = len(elevs_all_s)
                    own_elevs = [r.get("Elevation") or 0 for r in ln_recs]
                    er_vals = []
                    for oe in own_elevs:
                        ci = min(
                            range(n_e), key=lambda i: abs(elevs_all_s[i] - oe)
                        )
                        er_vals.append((ci + 1) / n_e)
                    es = _avg(er_vals)

                    vals = [fs, ts, gs, ds, es]
                    spider_traces.append(
                        {
                            "type": "scatterpolar",
                            "r": vals + [vals[0]],
                            "theta": cats + [cats[0]],
                            "fill": "toself",
                            "fillcolor": col_fill,
                            "line": {"color": col, "width": 2},
                            "marker": {"size": 4, "color": col},
                            "name": ln,
                            "hovertemplate": f"<b>{ln}</b>  %{{theta}}: %{{r:.2f}}<extra></extra>",
                        }
                    )

                return {
                    "data": spider_traces,
                    "layout": {
                        "polar": {
                            "radialaxis": {
                                "visible": True,
                                "range": [0, 1],
                                "showticklabels": False,
                                "gridcolor": gc,
                                "linecolor": gc,
                                "gridwidth": 0.5,
                            },
                            "angularaxis": {
                                "tickfont": {"size": 8, "color": sc},
                                "gridcolor": gc,
                                "linecolor": gc,
                                "direction": "clockwise",
                            },
                            "bgcolor": bg,
                            "hole": 0.05,
                        },
                        "margin": {"l": 32, "r": 32, "t": 34, "b": 28},
                        "paper_bgcolor": bg,
                        "font": {"color": fc, "size": 9},
                        "showlegend": True,
                        "legend": {
                            "orientation": "h",
                            "y": -0.14,
                            "x": 0.5,
                            "xanchor": "center",
                            "font": {"size": 8},
                        },
                        "title": {
                            "text": "Per-Line Quality Comparison  —  click a station for detail",
                            "font": {"size": 9, "color": fc},
                            "x": 0.5,
                            "xanchor": "center",
                        },
                    },
                }

            # Single line or default: elevation profile
            ids = [r.get("ID", f"S{i}") for i, r in enumerate(records)]
            elevs = [r.get("Elevation") or 0 for r in records]
            return {
                "data": [
                    {
                        "type": "scatter",
                        "x": list(range(len(ids))),
                        "y": elevs,
                        "mode": "lines",
                        "fill": "tozeroy",
                        "fillcolor": "rgba(137,180,250,0.12)",
                        "line": {"color": "rgba(137,180,250,0)", "width": 0},
                        "hoverinfo": "skip",
                        "showlegend": False,
                    },
                    {
                        "type": "scatter",
                        "x": list(range(len(ids))),
                        "y": elevs,
                        "mode": "lines+markers",
                        "line": {
                            "color": "#89b4fa",
                            "width": 2,
                            "shape": "spline",
                            "smoothing": 0.8,
                        },
                        "marker": {
                            "size": 5,
                            "color": "#89b4fa",
                            "line": {"color": bg, "width": 0},
                        },
                        "hovertemplate": "<b>%{text}</b><br>%{y:.0f} m<extra></extra>",
                        "text": ids,
                        "showlegend": False,
                    },
                ],
                "layout": {
                    "margin": {"l": 38, "r": 8, "t": 32, "b": 14},
                    "paper_bgcolor": bg,
                    "plot_bgcolor": bg,
                    "font": {"color": fc, "size": 8},
                    "xaxis": {"visible": False, "showgrid": False},
                    "yaxis": {
                        "showgrid": True,
                        "gridcolor": gc,
                        "zeroline": False,
                        "tickfont": {"size": 7},
                        "gridwidth": 0.5,
                        "title": {
                            "text": "Elevation (m)",
                            "font": {"size": 8},
                            "standoff": 4,
                        },
                    },
                    "title": {
                        "text": "Elevation Profile  —  click a station for quality radar",
                        "font": {"size": 9, "color": fc},
                        "x": 0.5,
                        "xanchor": "center",
                    },
                },
            }

        # ── Station selected → 5-axis quality spider ──────────────────────
        match = next(
            (r for r in records if str(r.get("ID")) == str(station_id)), None
        )
        if not match:
            return _empty_analytics_fig(fc, bg)

        nfreqs = [r.get("N_freq") or 0 for r in records]
        max_nf = max(nfreqs) or 1
        med_nf = sorted(nfreqs)[len(nfreqs) // 2] or 1
        own_nf = match.get("N_freq") or 0

        freq_score = min(own_nf / max_nf, 1.0)
        tipper_score = 1.0 if match.get("Tipper") else 0.0
        geo_score = (
            1.0 if (match.get("Latitude") and match.get("Longitude")) else 0.0
        )
        fill_score = min(own_nf / med_nf, 1.0)

        elevs_s = sorted(r.get("Elevation") or 0 for r in records)
        own_e = match.get("Elevation") or 0
        c_idx = min(
            range(len(elevs_s)), key=lambda i: abs(elevs_s[i] - own_e)
        )
        elev_rank = (c_idx + 1) / len(elevs_s)

        cats = [
            "Freq. Range",
            "Tipper",
            "Coordinates",
            "Data Fill",
            "Elev. Rank",
        ]
        vals = [freq_score, tipper_score, geo_score, fill_score, elev_rank]

        # Background reference rings at 0.33 and 0.67
        def _ring(r_val, color):
            return {
                "type": "scatterpolar",
                "r": [r_val] * (len(cats) + 1),
                "theta": cats + [cats[0]],
                "mode": "lines",
                "line": {"color": color, "width": 0.8, "dash": "dot"},
                "hoverinfo": "skip",
                "showlegend": False,
            }

        return {
            "data": [
                _ring(0.33, gc),
                _ring(0.67, gc),
                {
                    "type": "scatterpolar",
                    "r": vals + [vals[0]],
                    "theta": cats + [cats[0]],
                    "fill": "toself",
                    "fillcolor": "rgba(137,180,250,0.20)",
                    "line": {"color": "#89b4fa", "width": 2.2},
                    "marker": {"size": 5, "color": "#89b4fa"},
                    "hovertemplate": "<b>%{theta}</b>: %{r:.2f}<extra></extra>",
                    "name": station_id,
                    "showlegend": False,
                },
            ],
            "layout": {
                "polar": {
                    "radialaxis": {
                        "visible": True,
                        "range": [0, 1],
                        "showticklabels": False,
                        "gridcolor": gc,
                        "linecolor": gc,
                        "gridwidth": 0.5,
                    },
                    "angularaxis": {
                        "tickfont": {"size": 8, "color": sc},
                        "gridcolor": gc,
                        "linecolor": gc,
                        "direction": "clockwise",
                    },
                    "bgcolor": bg,
                    "hole": 0.05,
                },
                "margin": {"l": 32, "r": 32, "t": 34, "b": 16},
                "paper_bgcolor": bg,
                "font": {"color": fc, "size": 9},
                "showlegend": False,
                "title": {
                    "text": f"Quality Radar  —  {station_id}",
                    "font": {"size": 9, "color": fc},
                    "x": 0.5,
                    "xanchor": "center",
                },
            },
        }

    # ── Navigation helpers ────────────────────────────────────────────────────
    clientside_callback(
        "function(n) { if (n) return 'home'; return window.dash_clientside.no_update; }",
        Output(IDs.NAV_SECTION, "data", allow_duplicate=True),
        Input("map-page-home-btn", "n_clicks"),
        prevent_initial_call=True,
    )
    clientside_callback(
        "function(n) { if (n) return 'home'; return window.dash_clientside.no_update; }",
        Output(IDs.NAV_SECTION, "data", allow_duplicate=True),
        Input("profile-page-home-btn", "n_clicks"),
        prevent_initial_call=True,
    )
    clientside_callback(
        "function(n) { if (n) return 'tab-rho-phi'; return window.dash_clientside.no_update; }",
        Output(IDs.PROFILE_TABS, "active_tab", allow_duplicate=True),
        Input("prof-pub-back-btn", "n_clicks"),
        prevent_initial_call=True,
    )

    # ── Survey Insights panel callbacks ───────────────────────────────────────

    from dash import ALL, callback

    # SI-1. Data health — auto-populated from STORE_DATA (no agent run needed)
    @callback(
        Output("si-health-dot", "className"),
        Output("si-health-text", "children"),
        Output("si-health-chips", "children"),
        Input(IDs.STORE_DATA, "data"),
        prevent_initial_call=False,
    )
    def _si_health(store_data):
        if not store_data or not store_data.get("n_stations"):
            return (
                "si-dot si-dot-empty",
                "No survey loaded",
                [],
            )

        n_stn = store_data.get("n_stations", 0)
        n_lines = store_data.get("n_lines", 0)
        lc = store_data.get("line_counts", {})
        data_dir = store_data.get("data_dir", "")

        dot_cls = "si-dot si-dot-ok"
        txt = f"{n_stn} station{'s' if n_stn != 1 else ''} loaded"

        chips: list = [
            html.Span(
                [
                    html.I(className="bi bi-broadcast me-1"),
                    f"{n_lines} profile{'s' if n_lines != 1 else ''}",
                ],
                className="si-chip si-chip-blue",
            ),
        ]

        # Per-profile breakdown (up to 4 shown)
        for _i, (line_id, cnt) in enumerate(list(lc.items())[:4]):
            chips.append(
                html.Span(
                    f"{line_id}: {cnt}", className="si-chip si-chip-sub"
                )
            )
        if len(lc) > 4:
            chips.append(
                html.Span(
                    f"+{len(lc) - 4} more", className="si-chip si-chip-sub"
                )
            )

        if data_dir:
            chips.append(
                html.Span(
                    [
                        html.I(className="bi bi-folder2 me-1"),
                        data_dir.split("/")[-1] or data_dir,
                    ],
                    className="si-chip si-chip-path",
                    title=data_dir,
                )
            )

        return dot_cls, txt, chips

    # SI-2. Last run — pulled from AGENTS_STORE (no re-execution)
    @callback(
        Output("si-last-run", "children"),
        Input(IDs.AGENTS_STORE, "data"),
        prevent_initial_call=False,
    )
    def _si_last_run(history):
        if not history:
            return html.Span("No runs yet", className="si-no-run")

        run = history[0]  # most recent first
        name = run.get("name", "Unknown")
        ts = run.get("timestamp", "")
        status = run.get("status", "success")
        summary = run.get("summary", "")

        icon_cls = (
            "bi bi-check-circle-fill si-run-ok"
            if status == "success"
            else "bi bi-x-circle-fill si-run-err"
        )

        summ_short = (
            (summary[:80] + "…") if summary and len(summary) > 80 else summary
        )

        return html.Div(
            [
                html.Div(
                    [
                        html.I(className=icon_cls),
                        html.Span(name, className="si-run-name"),
                        html.Span(ts, className="si-run-ts"),
                    ],
                    className="si-run-row",
                ),
                html.P(summ_short or "—", className="si-run-summary"),
                html.Button(
                    [
                        html.I(className="bi bi-arrow-right me-1"),
                        "View in Agents",
                    ],
                    id="si-view-agents-btn",
                    className="si-view-btn",
                    n_clicks=0,
                ),
            ]
        )

    # SI-3. Quick-workflow cards → navigate to Agents page
    clientside_callback(
        """function(clicks) {
            if (!clicks || clicks.every(c => !c)) return window.dash_clientside.no_update;
            return 'agents';
        }""",
        Output(IDs.NAV_SECTION, "data", allow_duplicate=True),
        Input({"type": "si-workflow", "index": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )

    # SI-4. "Open Agents" header button → navigate
    clientside_callback(
        "function(n) { if (n) return 'agents'; return window.dash_clientside.no_update; }",
        Output(IDs.NAV_SECTION, "data", allow_duplicate=True),
        Input("si-open-agents-btn", "n_clicks"),
        prevent_initial_call=True,
    )

    # SI-5. "View in Agents" button in Last Run → navigate
    # This button only appears after a run exists, so suppress missing-ID warnings.
    clientside_callback(
        "function(n) { if (n) return 'agents'; return window.dash_clientside.no_update; }",
        Output(IDs.NAV_SECTION, "data", allow_duplicate=True),
        Input("si-view-agents-btn", "n_clicks"),
        prevent_initial_call=True,
    )


# ── Private helpers ───────────────────────────────────────────────────────────


def _empty_analytics_fig(font_color: str, bg: str) -> dict:
    """Transparent empty figure shown before data is loaded."""
    return {
        "data": [],
        "layout": {
            "paper_bgcolor": bg,
            "plot_bgcolor": bg,
            "xaxis": {"visible": False},
            "yaxis": {"visible": False},
            "margin": {"l": 0, "r": 0, "t": 0, "b": 0},
            "annotations": [
                {
                    "text": "Load survey data to view",
                    "x": 0.5,
                    "y": 0.5,
                    "showarrow": False,
                    "font": {"size": 9, "color": "#585b70"},
                    "xref": "paper",
                    "yref": "paper",
                }
            ],
        },
    }
