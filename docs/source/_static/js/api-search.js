/* Flat, searchable and sortable public API index. */
(function () {
  "use strict";

  function initialise() {
    const controls = document.getElementById("api-search");
    if (!controls) return;

    const input = document.getElementById("api-search-input");
    const limitSelect = document.getElementById("api-limit-filter");
    const status = document.getElementById("api-search-status");
    const pagination = document.getElementById("api-pagination");
    const content = controls.closest("main") || document;
    const table = content.querySelector("table.autosummary");
    if (!table) return;

    table.classList.add("api-object-table");
    const body = table.tBodies[0];
    const entries = Array.from(body.rows).map(function (row, originalIndex) {
      const link = row.cells[0].querySelector("a");
      const qualified = link ? link.textContent.trim() : row.cells[0].textContent.trim();
      const parts = qualified.split(".");
      const objectName = parts.pop();
      const owner = parts.join(".");
      if (link) link.textContent = objectName;

      const ownerLine = document.createElement("div");
      ownerLine.className = "api-object-owner";
      ownerLine.textContent = owner || "pycsamt";
      row.cells[1].appendChild(ownerLine);

      return {
        row: row,
        name: objectName.toLowerCase(),
        description: row.cells[1].textContent.toLowerCase(),
        search: (qualified + " " + row.textContent).toLowerCase(),
        originalIndex: originalIndex
      };
    });

    const head = table.createTHead();
    const headerRow = head.insertRow();
    ["Object", "Description"].forEach(function (label, column) {
      const th = document.createElement("th");
      th.scope = "col";
      const button = document.createElement("button");
      button.type = "button";
      button.className = "api-sort";
      button.dataset.column = String(column);
      button.innerHTML = "<span>" + label + "</span><span class=\"api-sort__icon\" aria-hidden=\"true\">◆</span>";
      th.appendChild(button);
      headerRow.appendChild(th);
    });

    let page = 1;
    let sortColumn = -1;
    let sortDirection = 1;

    function pageButton(label, target, disabled, current) {
      const button = document.createElement("button");
      button.type = "button";
      button.textContent = label;
      button.disabled = disabled;
      if (current) button.setAttribute("aria-current", "page");
      button.addEventListener("click", function () {
        page = target;
        refresh();
      });
      return button;
    }

    function refresh() {
      const terms = input.value.trim().toLowerCase().split(/\s+/).filter(Boolean);
      const filtered = entries.filter(function (entry) {
        return terms.every(function (term) { return entry.search.includes(term); });
      });

      filtered.sort(function (a, b) {
        if (sortColumn < 0) return a.originalIndex - b.originalIndex;
        const left = sortColumn === 0 ? a.name : a.description;
        const right = sortColumn === 0 ? b.name : b.description;
        return left.localeCompare(right) * sortDirection;
      });
      filtered.forEach(function (entry) { body.appendChild(entry.row); });

      const limit = limitSelect.value === "all" ? Math.max(filtered.length, 1) : Number(limitSelect.value);
      const pages = Math.max(1, Math.ceil(filtered.length / limit));
      page = Math.min(page, pages);
      const start = (page - 1) * limit;
      const end = Math.min(start + limit, filtered.length);
      const visible = new Set(filtered.slice(start, end));
      entries.forEach(function (entry) { entry.row.hidden = !visible.has(entry); });

      status.textContent = filtered.length
        ? "Showing " + (start + 1) + " to " + end + " of " + filtered.length + " entries"
        : "No matching API objects";

      pagination.replaceChildren();
      if (pages > 1) {
        pagination.appendChild(pageButton("Previous", page - 1, page === 1, false));
        const first = Math.max(1, page - 2);
        const last = Math.min(pages, first + 4);
        for (let number = first; number <= last; number += 1) {
          pagination.appendChild(pageButton(String(number), number, false, number === page));
        }
        pagination.appendChild(pageButton("Next", page + 1, page === pages, false));
      }
    }

    table.querySelectorAll(".api-sort").forEach(function (button) {
      button.addEventListener("click", function () {
        const column = Number(button.dataset.column);
        sortDirection = sortColumn === column ? -sortDirection : 1;
        sortColumn = column;
        table.querySelectorAll(".api-sort").forEach(function (item) {
          item.removeAttribute("data-direction");
        });
        button.dataset.direction = sortDirection === 1 ? "ascending" : "descending";
        page = 1;
        refresh();
      });
    });
    input.addEventListener("input", function () { page = 1; refresh(); });
    limitSelect.addEventListener("change", function () { page = 1; refresh(); });
    refresh();
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initialise);
  } else {
    initialise();
  }
}());
