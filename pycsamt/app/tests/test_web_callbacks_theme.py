# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.theme (dark / light theme toggle).

Trivial, pure-function callback -- no data, no mocks needed.
"""

from __future__ import annotations

from pycsamt.app.web.layout import IDs

# ── Dash callback-lookup helpers (shared pattern across web-callback tests) ──


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _cb(web_app, output_id_prop):
    return _unwrap(web_app.callback_map[output_id_prop])


def _cb_multi(web_app, *substrings):
    key = next(k for k in web_app.callback_map if all(s in k for s in substrings))
    return _unwrap(web_app.callback_map[key])


def _cb_by_input(web_app, output_substr, input_id):
    for k, v in web_app.callback_map.items():
        if output_substr not in k:
            continue
        if any(input_id in str(i.get("id")) for i in v.get("inputs", [])):
            return _unwrap(v)
    raise AssertionError(
        f"no callback found for output~={output_substr!r} input={input_id!r}"
    )


# ── toggle_theme ─────────────────────────────────────────────────────────────


class TestToggleTheme:
    def _fn(self, web_app):
        return _cb_multi(
            web_app, f"{IDs.STORE_THEME}.data", f"{IDs.BTN_THEME}.children"
        )

    def test_dark_to_light(self, web_app):
        new_theme, label = self._fn(web_app)(1, "dark")
        assert new_theme == "light"
        assert label == "☀ Theme"

    def test_light_to_dark(self, web_app):
        new_theme, label = self._fn(web_app)(1, "light")
        assert new_theme == "dark"
        assert label == "☾ Theme"

    def test_none_current_theme_defaults_to_dark_then_toggles_light(self, web_app):
        # current_theme=None -> treated as "dark" -> toggles to "light"
        new_theme, label = self._fn(web_app)(1, None)
        assert new_theme == "light"
        assert label == "☀ Theme"

    def test_unknown_current_theme_treated_as_light_toggles_dark(self, web_app):
        # Any value other than "dark" is treated as "light" by the ternary,
        # so an unrecognised string still toggles deterministically to dark.
        new_theme, label = self._fn(web_app)(1, "sepia")
        assert new_theme == "dark"
        assert label == "☾ Theme"
