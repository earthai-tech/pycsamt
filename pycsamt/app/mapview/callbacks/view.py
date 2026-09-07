# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""View-rail switching and the main canvas figure render."""

from __future__ import annotations

from dash import Input, Output, State, ctx

from pycsamt.app._borehole import (
    add_pcbh_to_figure,
    document_from_store,
    scene_borehole_traces,
    strip_log_figure,
)

from .._ids import IDs
from .._render import VIEW_TITLES, empty_figure, figure_for
from ..cache import get_view

_RAIL = {
    IDs.RAIL_MAP: "map",
    IDs.RAIL_3D: "map3d",
    IDs.RAIL_BH: "bh",
    IDs.RAIL_GEO: "geology",
}


def register_view(app) -> None:
    _register_rail(app)
    _register_render(app)


def _register_rail(app) -> None:
    @app.callback(
        Output(IDs.STORE_VIEW, "data"),
        Output(IDs.RAIL_MAP, "className"),
        Output(IDs.RAIL_3D, "className"),
        Output(IDs.RAIL_BH, "className"),
        Output(IDs.RAIL_GEO, "className"),
        Output(IDs.CANVAS_TITLE, "children"),
        Input(IDs.RAIL_MAP, "n_clicks"),
        Input(IDs.RAIL_3D, "n_clicks"),
        Input(IDs.RAIL_BH, "n_clicks"),
        Input(IDs.RAIL_GEO, "n_clicks"),
        prevent_initial_call=True,
    )
    def switch(*_clicks):
        active_id = ctx.triggered_id
        view = _RAIL.get(active_id, "map")

        def cls(bid):
            return (
                "mv-rail-btn mv-rail-active"
                if bid == active_id
                else "mv-rail-btn"
            )

        return (
            view,
            cls(IDs.RAIL_MAP),
            cls(IDs.RAIL_3D),
            cls(IDs.RAIL_BH),
            cls(IDs.RAIL_GEO),
            VIEW_TITLES.get(view, "Map"),
        )


def _borehole_view_figure(
    pcbh_store, mode, family, opacity, labels, as_tubes, radius, theme
):
    """Standalone Boreholes-view canvas: strip-log or free 3-D holes."""
    import plotly.graph_objects as go

    document = document_from_store(pcbh_store) if pcbh_store else None
    if document is None:
        fig = empty_figure(
            theme, "Open Borehole Studio to import or build boreholes"
        )
        return fig
    if (mode or "striplog") == "striplog":
        return strip_log_figure(document, family=family, theme=theme)

    # free-standing 3-D holes in their own PCBH coordinates
    from pycsamt.format.borehole import DisplayRadiusPolicy, build_render_model
    from pycsamt.map.borehole_align import (
        AlignedHole,
        AlignedSegment,
        SceneAlignment,
    )

    policy = (
        DisplayRadiusPolicy(mode="fixed", fixed_radius=float(radius))
        if radius
        else DisplayRadiusPolicy()
    )
    model = build_render_model(document, family=family, radius_policy=policy)
    holes = []
    for hole in model.boreholes:
        segs = tuple(
            AlignedSegment(
                borehole_id=hole.borehole_id,
                family=family,
                from_md=s.from_md,
                to_md=s.to_md,
                color=s.color,
                points=tuple((p.x, p.y, p.z) for p in s.points),
                metadata=dict(s.metadata),
            )
            for s in hole.interval_segments
        )
        holes.append(
            AlignedHole(
                borehole_id=hole.borehole_id,
                name=hole.collar.label,
                collar_scene=hole.collar.position,
                centerline=tuple(
                    (p.x, p.y, p.z) for p in hole.centerline.points
                ),
                segments=segs,
                display_radius=hole.display_radius,
                relation="native",
            )
        )
    alignment = SceneAlignment(holes=tuple(holes), family=family)
    fig = go.Figure()
    for trace in scene_borehole_traces(
        alignment, as_tubes=as_tubes, opacity=opacity, show_labels=labels
    ):
        fig.add_trace(trace)
    fig.update_layout(
        template="plotly_dark" if theme == "dark" else "plotly_white",
        scene={
            "xaxis_title": "Easting",
            "yaxis_title": "Northing",
            "zaxis_title": "Elevation (m)",
            "aspectmode": "data",
        },
        margin={"l": 0, "r": 0, "t": 0, "b": 0},
        paper_bgcolor="rgba(0,0,0,0)",
    )
    return fig


def _geology_view_figure(geo_store, struct_store, view_mode, theme):
    """Standalone Geology-view canvas: legend swatch bar or a structural
    profile-position / depth preview, chosen by ``GEO_VIEW_MODE``."""
    from pycsamt.app._geology import legend_from_store, legend_preview_figure
    from pycsamt.app._structure import (
        structure_from_store,
        structure_section_figure,
    )

    if (view_mode or "legend") == "structure":
        doc = structure_from_store(struct_store) if struct_store else None
        return structure_section_figure(
            doc.model if doc else None, theme=theme
        )
    legend = legend_from_store(geo_store) if geo_store else None
    return legend_preview_figure(legend)


def _register_render(app) -> None:
    @app.callback(
        Output(IDs.CANVAS_GRAPH, "figure"),
        Output(IDs.WELCOME, "style"),
        Input(IDs.STORE_DATA, "data"),
        Input(IDs.STORE_VIEW, "data"),
        Input(IDs.STORE_CONTROLS, "data"),
        Input(IDs.STORE_THEME, "data"),
        Input(IDs.STORE_LINES, "data"),
        Input(IDs.STORE_FIT, "data"),
        Input(IDs.STORE_MASKED, "data"),
        Input(IDs.PCBH_STORE, "data"),
        Input(IDs.PCBH_VISIBLE, "value"),
        Input(IDs.PCBH_LABELS, "value"),
        Input(IDs.PCBH_FAMILY, "value"),
        Input(IDs.PCBH_OPACITY, "value"),
        Input(IDs.PCBH_AS_TUBES, "value"),
        Input(IDs.PCBH_RADIUS, "value"),
        Input(IDs.PCBH_ON_MAP, "value"),
        Input(IDs.PCBH_IN_3D, "value"),
        Input(IDs.PCBH_3D_LEAN, "value"),
        Input(IDs.PCBH_3D_LEAN_DIR, "value"),
        Input(IDs.PCBH_3D_LABEL_ANGLE, "value"),
        Input(IDs.PCBH_3D_LABEL_SIZE, "value"),
        Input(IDs.PCBH_3D_COLLAR_SIZE, "value"),
        Input(IDs.PCBH_3D_DEPTH_TICKS, "value"),
        Input(IDs.PCBH_PATCH_GEOLOGY, "value"),
        Input(IDs.PCBH_PATCH_WIDTH, "value"),
        Input(IDs.PCPT_STORE, "data"),
        Input(IDs.PCPT_VISIBLE, "value"),
        Input(IDs.BH_VIEW_MODE, "value"),
        Input(IDs.GEO_STORE, "data"),
        Input(IDs.GEO_APPLY, "value"),
        Input(IDs.STRUCT_STORE, "data"),
        Input(IDs.STRUCT_APPLY, "value"),
        Input(IDs.GEO_VIEW_MODE, "value"),
        State(IDs.STORE_VIEWPORT, "data"),
        State(IDs.SESSION_ID, "data"),
        prevent_initial_call=False,
    )
    def render(
        store, view_name, controls, theme, lines, fit, masked, pcbh_store,
        pcbh_visible, pcbh_labels, pcbh_family, pcbh_opacity, pcbh_as_tubes,
        pcbh_radius, pcbh_on_map, pcbh_in_3d, pcbh_lean, pcbh_lean_dir,
        pcbh_label_angle, pcbh_label_size, pcbh_collar_size, pcbh_depth_ticks,
        pcbh_patch_geology, pcbh_patch_width,
        pcpt_store, pcpt_visible, bh_view_mode, geo_store, geo_apply,
        struct_store, struct_apply, geo_view_mode,
        viewport, session_id
    ):
        theme = theme or "light"
        geology_bands = None
        geology_patterns = None
        if geo_store and geo_apply:
            from pycsamt.app._geology import geology_bands_from_store

            try:
                geology_bands = geology_bands_from_store(geo_store)
            except Exception:  # noqa: BLE001
                geology_bands = None
            if (controls or {}).get("geology_fill") == "pattern":
                from pycsamt.app._geology import (
                    geology_pattern_stencils_from_store,
                )

                try:
                    geology_patterns = geology_pattern_stencils_from_store(
                        geo_store
                    )
                except Exception:  # noqa: BLE001
                    geology_patterns = None
        opacity = pcbh_opacity if pcbh_opacity is not None else 0.9
        family = pcbh_family or "lithology"
        radius = float(pcbh_radius) if pcbh_radius else None
        show_pcbh = bool(pcbh_visible)

        if (view_name or "map") == "bh":
            return _borehole_view_figure(
                pcbh_store if show_pcbh else None, bh_view_mode, family,
                opacity, bool(pcbh_labels), bool(pcbh_as_tubes), radius,
                theme,
            ), {"display": "none"}

        if view_name == "geology":
            return _geology_view_figure(
                geo_store, struct_store, geo_view_mode, theme,
            ), {"display": "none"}

        if not store or not store.get("n_stations"):
            fig = empty_figure(theme, "Upload PCBH or load survey lines")
            if pcbh_store and show_pcbh:
                fig = add_pcbh_to_figure(
                    fig,
                    pcbh_store,
                    visible=True,
                    show_labels=bool(pcbh_labels),
                    family=family,
                    opacity=opacity,
                )
            welcome = (
                {"display": "none"}
                if pcbh_store
                else {"display": "flex"}
            )
            return fig, welcome
        view = get_view(session_id)
        if view is None:
            return empty_figure(
                theme, "Session data unavailable — reload lines."
            ), {"display": "flex"}
        active = (lines or {}).get("active")
        borehole_opts = (
            {
                "store": pcbh_store,
                # 3-D scene insertion is opt-in (PCBH_IN_3D, default off)
                "visible": show_pcbh and bool(pcbh_in_3d),
                "labels": bool(pcbh_labels),
                "family": family,
                "opacity": opacity,
                "as_tubes": bool(pcbh_as_tubes),
                "radius": radius,
                "lean_deg": float(pcbh_lean) if pcbh_lean else 0.0,
                "lean_dir": pcbh_lean_dir or "N",
                "label_angle": (
                    float(pcbh_label_angle) if pcbh_label_angle else 0.0
                ),
                "label_size": (
                    float(pcbh_label_size) if pcbh_label_size else 11.0
                ),
                "collar_size": (
                    float(pcbh_collar_size) if pcbh_collar_size else 5.0
                ),
                "depth_ticks": (
                    float(pcbh_depth_ticks) if pcbh_depth_ticks else 0.0
                ),
                "patch_geology": bool(pcbh_patch_geology),
                "patch_width": (
                    float(pcbh_patch_width) if pcbh_patch_width else 20.0
                ),
                "show_on_map": show_pcbh and bool(pcbh_on_map),
                "points_store": pcpt_store,
                "points_visible": bool(pcpt_visible),
            }
            if (pcbh_store or pcpt_store)
            else None
        )
        structure_opts = (
            {"store": struct_store}
            if struct_store and struct_apply
            else None
        )
        fig = figure_for(
                view_name or "map",
                view,
                controls,
                theme=theme,
                active_lines=active,
                masked=masked,
                fit=int(fit or 0),
                boreholes=borehole_opts,
                geology=geology_bands,
                geology_patterns=geology_patterns,
                structure=structure_opts,
                viewport=viewport,
            )
        return fig, {"display": "none"}
