"""Accessible collapsible literal includes for long documentation examples.

``code-dropdown`` deliberately delegates file slicing, ``pyobject`` lookup,
syntax highlighting, line numbering, dependency tracking, and diagnostics to
Sphinx's own :class:`~sphinx.directives.code.LiteralInclude` implementation.
Only the presentation wrapper is custom.
"""

from __future__ import annotations

from docutils import nodes
from docutils.parsers.rst import directives
from sphinx.directives.code import LiteralInclude


class code_dropdown(nodes.General, nodes.Element):
    """Semantic disclosure wrapper around a Sphinx literal include."""


class code_dropdown_summary(nodes.General, nodes.Element):
    """Summary row for :class:`code_dropdown`."""


class CodeDropdownDirective(LiteralInclude):
    """Include source code inside an accessible, collapsed disclosure panel."""

    option_spec = dict(LiteralInclude.option_spec)
    option_spec.update(
        {
            "title": directives.unchanged,
            "open": directives.flag,
        }
    )

    def run(self):
        included = super().run()
        if not included:
            return included

        title = self.options.get("title", "View source code")
        wrapper = code_dropdown()
        wrapper["open"] = "open" in self.options
        wrapper["title"] = title
        wrapper["classes"].append("pyc-code-dropdown")

        summary = code_dropdown_summary()
        summary += nodes.Text(title)
        wrapper += summary
        wrapper.extend(included)
        return [wrapper]


def visit_dropdown_html(translator, node):
    opened = " open" if node.get("open") else ""
    translator.body.append(
        f'<details class="pyc-code-dropdown"{opened}>'
    )


def depart_dropdown_html(translator, node):
    translator.body.append("</details>")


def visit_summary_html(translator, node):
    title = translator.encode(node.astext())
    translator.body.append(
        '<summary class="pyc-code-dropdown__summary">'
        '<span class="pyc-code-dropdown__icon" aria-hidden="true">&lt;/&gt;</span>'
        '<span class="pyc-code-dropdown__labels">'
        f'<span class="pyc-code-dropdown__title">{title}</span>'
        '<span class="pyc-code-dropdown__hint" '
        'data-open="Click to inspect and copy the complete code" '
        'data-close="Click to collapse the source code">'
        'Click to inspect and copy the complete code'
        '</span></span>'
        '<span class="pyc-code-dropdown__chevron" aria-hidden="true"></span>'
        '</summary>'
    )
    raise nodes.SkipNode


def depart_summary_html(translator, node):  # pragma: no cover - SkipNode
    pass


def visit_dropdown_fallback(translator, node):
    """Keep code visible for PDF, man-page, text, and other builders."""


def depart_dropdown_fallback(translator, node):
    pass


def visit_summary_fallback(translator, node):
    raise nodes.SkipNode


def depart_summary_fallback(translator, node):  # pragma: no cover - SkipNode
    pass


def setup(app):
    fallback = (visit_dropdown_fallback, depart_dropdown_fallback)
    summary_fallback = (visit_summary_fallback, depart_summary_fallback)
    app.add_node(
        code_dropdown,
        html=(visit_dropdown_html, depart_dropdown_html),
        latex=fallback,
        text=fallback,
        man=fallback,
        texinfo=fallback,
    )
    app.add_node(
        code_dropdown_summary,
        html=(visit_summary_html, depart_summary_html),
        latex=summary_fallback,
        text=summary_fallback,
        man=summary_fallback,
        texinfo=summary_fallback,
    )
    app.add_directive("code-dropdown", CodeDropdownDirective)
    return {
        "version": "1.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
