# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Clientside navigation — page switching, active-button highlight, sidebar toggle."""
from __future__ import annotations

from dash import ALL, ClientsideFunction, Input, Output, State, clientside_callback

from dash import Input, Output

from pycsamt.app.web.layout import IDs, _NAV_ENTRIES

_PAGE_IDS = [pid for pid, _, _ in _NAV_ENTRIES]


def register_navigation(app) -> None:

    # ── Page navigation ──────────────────────────────────────────────────────

    # Each nav-btn click → write its page_id into NAV_SECTION store
    for page_id, _, _ in _NAV_ENTRIES:
        clientside_callback(
            f"""
            function(n) {{
                if (n) return "{page_id}";
                return window.dash_clientside.no_update;
            }}
            """,
            Output(IDs.NAV_SECTION, "data", allow_duplicate=True),
            Input(f"nav-btn-{page_id}", "n_clicks"),
            prevent_initial_call=True,
        )

    # NAV_SECTION → toggle page visibility (display flex / none)
    _pages_json = str(_PAGE_IDS)
    clientside_callback(
        f"""
        function(section) {{
            var pages = {_pages_json};
            var styles = [];
            for (var i = 0; i < pages.length; i++) {{
                if (pages[i] === section) {{
                    styles.push({{
                        display: "flex",
                        flexDirection: "column",
                        flex: "1",
                        overflow: "hidden",
                        minHeight: "0",
                        animationName: "page-in"
                    }});
                }} else {{
                    styles.push({{display: "none"}});
                }}
            }}
            return styles;
        }}
        """,
        [Output(f"page-{pid}", "style") for pid in _PAGE_IDS],
        Input(IDs.NAV_SECTION, "data"),
    )

    # NAV_SECTION → active class on nav buttons
    clientside_callback(
        f"""
        function(section) {{
            var pages = {_pages_json};
            var classes = [];
            for (var i = 0; i < pages.length; i++) {{
                classes.push(pages[i] === section ? "nav-btn active" : "nav-btn");
            }}
            return classes;
        }}
        """,
        [Output(f"nav-btn-{pid}", "className") for pid in _PAGE_IDS],
        Input(IDs.NAV_SECTION, "data"),
    )

    # NAV_SECTION → update navbar breadcrumb chip
    _labels_json = str(
        {pid: lbl for pid, _, lbl in _NAV_ENTRIES}
    ).replace("'", '"')
    clientside_callback(
        f"""
        function(section) {{
            var labels = {_labels_json};
            return labels[section] || section;
        }}
        """,
        Output("nav-page-chip", "children"),
        Input(IDs.NAV_SECTION, "data"),
    )

    # ── Logo click → home ───────────────────────────────────────────────────
    clientside_callback(
        "function(n) { if (n) return 'home'; return window.dash_clientside.no_update; }",
        Output(IDs.NAV_SECTION, "data", allow_duplicate=True),
        Input("nav-brand", "n_clicks"),
        prevent_initial_call=True,
    )

    # ── Sidebar collapse toggle ──────────────────────────────────────────────

    # Toggle button click → flip "expanded" / "collapsed" in the store
    clientside_callback(
        """
        function(n, current) {
            if (!n) return window.dash_clientside.no_update;
            return current === "expanded" ? "collapsed" : "expanded";
        }
        """,
        Output(IDs.STORE_SIDEBAR, "data"),
        Input(IDs.BTN_SIDEBAR_TOGGLE, "n_clicks"),
        State(IDs.STORE_SIDEBAR, "data"),
        prevent_initial_call=True,
    )

    # Sidebar store → sidebar CSS class (adds / removes "collapsed")
    clientside_callback(
        """
        function(state) {
            return state === "collapsed" ? "sidebar collapsed" : "sidebar";
        }
        """,
        Output(IDs.SIDEBAR_DIV, "className"),
        Input(IDs.STORE_SIDEBAR, "data"),
    )

    # Sidebar store → flip the chevron icon direction
    clientside_callback(
        """
        function(state) {
            return state === "collapsed"
                ? "bi bi-chevron-right"
                : "bi bi-chevron-left";
        }
        """,
        Output("sidebar-chevron", "className"),
        Input(IDs.STORE_SIDEBAR, "data"),
    )

    # Sidebar store → update title attr on toggle button
    clientside_callback(
        """
        function(state) {
            return state === "collapsed" ? "Expand sidebar" : "Collapse sidebar";
        }
        """,
        Output(IDs.BTN_SIDEBAR_TOGGLE, "title"),
        Input(IDs.STORE_SIDEBAR, "data"),
    )

    # ── Welcome overlay ──────────────────────────────────────────────────────

    # STORE_DATA populated → animate overlay out, then hide it permanently.
    #
    # CRITICAL: must return the actual className string (not no_update).
    # Returning no_update leaves Dash's VirtualDOM thinking className is
    # still "wlc-overlay".  When the many server callbacks triggered by
    # STORE_DATA all settle and React reconciles the full component tree, it
    # resets the element to the VirtualDOM value — stripping wlc-hiding and
    # restoring pointer-events, which makes the welcome cards clickable again
    # and re-opens the modal on accidental clicks.
    #
    # By returning the real class string, Dash records it in the prop store
    # so React reconciliation keeps it.  The final display:none is applied
    # via inline style (outside React's control) after the animation ends.
    clientside_callback(
        """
        function(store_data) {
            if (!store_data || !store_data.n_stations)
                return window.dash_clientside.no_update;
            // Force full removal via inline style after the 0.9 s animation.
            // Inline style is not managed by React (no style prop on this div),
            // so it survives subsequent reconciliation passes.
            var ov = document.getElementById('welcome-overlay');
            if (ov) {
                setTimeout(function() { ov.style.display = 'none'; }, 1000);
            }
            return 'wlc-overlay wlc-hiding';
        }
        """,
        Output(IDs.WELCOME_OVERLAY, "className"),
        Input(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )

    # ── Theme toggle ─────────────────────────────────────────────────────────

    # BTN_THEME click → flip "dark" / "light" in the store
    clientside_callback(
        """
        function(n, current) {
            if (!n) return window.dash_clientside.no_update;
            return current === "dark" ? "light" : "dark";
        }
        """,
        Output(IDs.STORE_THEME, "data"),
        Input(IDs.BTN_THEME, "n_clicks"),
        State(IDs.STORE_THEME, "data"),
        prevent_initial_call=True,
    )

    # Theme store → set data-theme attribute on <html> for CSS variable swap
    clientside_callback(
        """
        function(theme) {
            document.documentElement.setAttribute("data-theme", theme || "dark");
            return window.dash_clientside.no_update;
        }
        """,
        Output(IDs.STORE_THEME, "data", allow_duplicate=True),
        Input(IDs.STORE_THEME, "data"),
        prevent_initial_call=True,
    )
