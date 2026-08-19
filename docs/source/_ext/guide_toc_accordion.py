"""Collapsible sub-branches for the numbered ``pycsamt-guide-toc`` outline.

The ``.pycsamt-guide-toc`` toctrees (see ``user_guide/index.rst`` and every
section's own ``index.rst``) render every in-page heading down to
``:maxdepth: 3`` -- intentional, it's what gives the outline its
scikit-learn-style completeness (see the ``:titlesonly:`` note in
``user_guide/index.rst``). Rendered flat, that same completeness makes a page
pulling in a whole section's worth of sub-headings unreadably long.

This module runs after Sphinx has resolved a page's toctree HTML and rewrites
any list item that carries its own nested list -- other than the outline's
top level (the chapters, which stay always visible) -- into a native
``<details>``/``<summary>`` disclosure, closed by default. Siblings under the
same parent share a ``name`` attribute, so opening one natively closes the
others: the browser handles the mutual exclusion, no JS is shipped, and
browsers that predate the ``name`` attribute (Safari <17.4, Firefox <131)
simply fall back to independent, non-exclusive disclosures.
"""

from __future__ import annotations

from bs4 import BeautifulSoup


def _collapse_children(soup, ul, group_seq):
    """Wrap every ``<li>`` in *ul* that has its own nested ``<ul>`` in a
    closed ``<details>`` disclosure. Wrapped siblings in *ul* share one
    ``name`` so the browser treats them as a single accordion group."""
    group_seq[0] += 1
    group_name = f"pyc-toc-{group_seq[0]}"

    for li in ul.find_all("li", recursive=False):
        nested = li.find("ul", recursive=False)
        if nested is None:
            continue

        # Collapse grandchildren first so each nesting level gets its own
        # independent accordion group rather than sharing its parent's.
        _collapse_children(soup, nested, group_seq)

        link = li.find("a", recursive=False)
        title = link.get_text(strip=True) if link else "subsections"

        toggle = soup.new_tag(
            "span", attrs={"class": "pyc-toc-toggle", "aria-hidden": "true"}
        )
        toggle.append(soup.new_tag("i", attrs={"class": "fa-solid fa-chevron-down"}))

        summary = soup.new_tag(
            "summary", attrs={"aria-label": f"Show subsections of {title}"}
        )
        summary.append(toggle)

        details = soup.new_tag("details", attrs={"name": group_name})
        details.append(summary)
        nested.extract()
        details.append(nested)

        li["class"] = [*li.get("class", []), "has-children"]
        li.append(details)


def _collapse_guide_toctrees(app, pagename, templatename, context, doctree):
    body = context.get("body")
    if not body or "pycsamt-guide-toc" not in body:
        return

    soup = BeautifulSoup(body, "html.parser")
    wrappers = soup.find_all("div", class_="pycsamt-guide-toc")
    if not wrappers:
        return

    group_seq = [0]
    for wrapper in wrappers:
        top_ul = wrapper.find("ul", recursive=False)
        if top_ul is None:
            continue
        # The outline's first level (the chapters) always stays visible;
        # only recurse into what each chapter nests underneath it.
        for chapter in top_ul.find_all("li", recursive=False):
            nested = chapter.find("ul", recursive=False)
            if nested is not None:
                _collapse_children(soup, nested, group_seq)

    context["body"] = str(soup)


def setup(app):
    app.connect("html-page-context", _collapse_guide_toctrees)
    return {
        "version": "1.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
