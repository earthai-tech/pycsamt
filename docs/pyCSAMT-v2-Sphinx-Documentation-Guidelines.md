# pyCSAMT-v2 Sphinx Documentation Guidelines

These guidelines define how the **pyCSAMT-v2** Sphinx documentation should be written, structured, illustrated, and maintained. The objective is to produce documentation that is scientifically rigorous, reproducible, easy to follow, and useful both to new users and experienced geophysicists.

---

## 1. General Writing Style

Documentation pages should read naturally and progressively. Avoid fragmented explanations or sections that appear mechanically generated.

Each page should guide the reader from one concept to the next through clear transitions. Introduce the motivation for a method before explaining how it works, and explain the practical meaning of results after presenting them.

Avoid creating unnecessary sections merely to separate small ideas. In particular, do not systematically create sections such as:

- "Mathematical Formulation"
- "Analysis"
- "Interpretation"
- "Output"
- "Discussion"

unless a separate section is genuinely necessary for the scientific or pedagogical structure of the page.

Instead, mathematical explanations, output interpretation, and scientific discussion should normally be integrated naturally into the surrounding text.

The documentation must explain **why** a method is used, **what** it computes, **how** it works, and **how users should interpret the result**.

---

## 2. Terminology and Glossary References

Important technical terms should be explicitly defined.

When a Sphinx glossary term is used through:

```rst
:term:`term`
```

the corresponding term must exist in `glossary.rst` and contain a sufficiently clear definition.

Do not introduce glossary references merely for formatting. Glossary entries should be reserved for terminology that users may reasonably need to consult independently, such as:

- electromagnetic quantities;
- inversion terminology;
- processing concepts;
- statistical quantities;
- machine-learning concepts;
- pyCSAMT-specific concepts.

When appropriate, connect the glossary definition to the relevant documentation page so that users can move from a short definition to a complete scientific explanation.

---

## 3. Scientific and Mathematical Reproducibility

Scientific methods must be documented with enough detail for a knowledgeable user to understand and reproduce the calculation.

Mathematical formulations should therefore be added wherever they improve clarity or reproducibility.

Do **not** create a dedicated mathematical section automatically. Instead, introduce equations naturally when the scientific explanation requires them.

For example, a page may transition as follows:

> The apparent resistivity alone does not describe the complete electromagnetic response because the phase carries complementary information about the relative timing of the electric and magnetic fields. For a measured impedance \(Z\), the apparent resistivity is therefore obtained from ...

followed directly by the relevant equation.

### Mathematical notation

Use consistent mathematical notation throughout a page and, whenever possible, throughout related documentation pages.

A symbol should not silently change meaning between equations.

Before introducing a complex equation:

1. explain the physical or computational quantity being represented;
2. present the equation;
3. define every important symbol;
4. explain its practical interpretation;
5. connect it to the corresponding pyCSAMT implementation when relevant.

Avoid equations that are included only for appearance. Every equation should help the reader understand the implementation or scientific method.

---

## 4. Equation Labels and Cross-References

Theory and scientific-background pages should use labeled equations whenever an equation may be referenced later.

For example:

```rst
.. math::
   :label: eq-apparent-resistivity

   \rho_a =
   \frac{1}{\mu_0 \omega}
   \left| Z \right|^2
```

The equation can then be referenced in the surrounding text using the appropriate Sphinx equation reference.

For example:

```rst
As shown in :eq:`eq-apparent-resistivity`, the apparent resistivity
depends on both the impedance magnitude and angular frequency.
```

Equation labels should be:

- descriptive;
- unique;
- stable;
- written consistently.

Prefer labels such as:

```text
eq-apparent-resistivity
eq-phase-tensor
eq-static-shift-factor
eq-occam-objective
eq-ai-inversion-loss
```

rather than generic labels such as:

```text
eq1
formula2
math3
```

Equation numbering will be handled globally through the Sphinx configuration.

---

## 5. Code Examples

Interactive or executable examples should preferably use the `pycon` syntax rather than a standard Python code block.

Prefer:

```rst
.. code-block:: pycon

   >>> from pycsamt import CSAMT
   >>> data = CSAMT(...)
   >>> data.fit(...)
```

instead of:

```rst
.. code-block:: python

   from pycsamt import CSAMT
   data = CSAMT(...)
   data.fit(...)
```

The `pycon` format makes the documentation closer to an interactive Python session and clearly distinguishes user input from displayed results.

Use ordinary Python source formatting only when presenting complete scripts or code that is not intended to resemble an interactive session.

---

## 6. Code Examples Must Show Their Results

Whenever practical, code examples should be accompanied by their resulting output.

A code example should not leave the user guessing what they are expected to obtain.

The general workflow is:

1. run the example externally or in a temporary environment;
2. capture the actual output;
3. add the relevant output immediately after the example;
4. explain what the result means in the surrounding documentation.

For example:

```rst
.. code-block:: pycon

   >>> result = processor.fit(data)
   >>> print(result.status)
   completed
```

For structured outputs, reproduce the important portion exactly enough for users to compare with their own execution.

Do not fabricate output. Display output obtained from the actual documented example whenever execution is possible.

If an output is already present and remains correct, it does not need to be regenerated unnecessarily.

---

## 7. Figures Produced by Code Examples

If an example generates a figure, run the example and preserve the generated figure in the documentation image hierarchy.

The preferred locations are:

```text
docs/source/images/user_guide/
```

or:

```text
docs/source/images/theory/
```

depending on the page type.

Use additional subdirectories when useful for organization.

For example:

```text
docs/source/images/user_guide/processing/
docs/source/images/user_guide/inversion/
docs/source/images/theory/csamt/
docs/source/images/theory/ai_inversion/
```

The generated figure should then be referenced directly below or near the corresponding example.

For example:

```rst
.. figure:: /images/user_guide/processing/static_shift_example.png
   :align: center
   :width: 85%

   Static-shift correction obtained from the example above.
```

The surrounding text should explain what is visible in the figure and how it should be interpreted.

Do not create a separate "Interpretation" section simply for this purpose. Integrate the explanation into the natural flow of the page.

---

## 8. Scientific Interpretation of Figures

Figures should not appear without explanation.

After displaying a figure, explain the scientifically important features that the user should observe.

The explanation should answer questions such as:

- What pattern is visible?
- Why does that pattern occur?
- What physical quantity is represented?
- What does a change in color, amplitude, phase, resistivity, loss, uncertainty, or depth imply?
- What result would indicate a potential problem?
- How does the figure relate to the preceding calculation?

For example, rather than writing:

> The corrected data are shown below.

prefer an explanation such as:

> The corrected curves retain the frequency-dependent shape of the original
> response while reducing the approximately frequency-independent offset
> between neighbouring stations. This behaviour is consistent with the
> expected effect of galvanic static shift: the apparent-resistivity
> amplitude changes, whereas the phase response remains comparatively
> unaffected.

The objective is to help the user understand the result rather than merely display it.

---

## 9. Creating Plots When No Plotting Function Exists

Some API examples may produce numerical results without providing a dedicated `plot_*` function.

If visualization would significantly improve understanding, the documentation may include additional plotting code.

For short examples, the plotting instructions may be shown directly on the page.

For longer plotting implementations, place the plotting code in `docs/scripts/` and expose it through a code dropdown.

The generated figure should still be saved under the appropriate documentation image directory and referenced from the page.

Documentation-specific plotting code should remain clearly separated from the pyCSAMT public API unless the plotting functionality itself belongs in the library.

---

## 10. Multiple Figures

Avoid presenting many related figures vertically when they can be compared more effectively together.

When an example naturally produces several panels, prefer a single figure created using Matplotlib axes, for example:

```python
fig, axes = plt.subplots(...)
```

This is particularly appropriate because many pyCSAMT plotting functions already support an `ax` parameter.

If several related plots need to be compared, arrange them as:

- horizontal panels;
- vertical panels where scientifically appropriate;
- a compact grid.

The goal is to make visual comparison straightforward rather than forcing users to scroll through many nearly independent figures.

If separate image files already exist, they may also be arranged through suitable Sphinx layout mechanisms, but generating one coherent multi-panel figure is generally preferable for closely related outputs.

---

## 11. Long Code Examples

Long examples should not dominate a documentation page.

If an example becomes too large for comfortable inline reading, move the complete executable code into:

```text
docs/scripts/
```

and expose it with the project's `code-dropdown` directive.

For example:

```rst
.. code-dropdown:: /../../scripts/generate_ai_inversion_figures.py
   :language: python
   :pyobject: make_losses_regularization_tradeoff
   :linenos:
   :title: View regularization-trade-off source code
```

The documentation page should contain the scientific explanation, important usage commands, outputs, and resulting figures, while the complete implementation remains available through the expandable source-code block.

---

## 12. Do Not Import Documentation Scripts into Examples

Scripts stored in `docs/scripts/` are documentation resources and should not be imported into user-facing code examples.

Avoid examples such as:

```python
from docs.scripts.generate_ai_inversion_figures import (
    make_losses_regularization_tradeoff,
)
```

This may work inside the repository but fail when users copy the code into a normal Python environment or Jupyter notebook.

Instead:

- place long supporting code in `docs/scripts/`;
- expose it through `.. code-dropdown::`;
- keep the actual user-facing example based only on public pyCSAMT APIs and normal scientific Python dependencies.

A documentation example should remain meaningful when copied outside the pyCSAMT repository.

---

## 13. Prefer Public APIs in Examples

Examples should demonstrate interfaces that users are expected to call.

Avoid relying unnecessarily on:

- private functions;
- undocumented internal classes;
- repository-relative imports;
- testing utilities;
- documentation-only modules;
- temporary development interfaces.

If an internal function must be discussed for architectural or developer documentation, clearly identify it as internal.

User-guide examples should normally rely on stable public APIs.

---

## 14. Reproducible Example Inputs

Examples should make the origin of their input data clear.

When using packaged sample data, explicitly identify the dataset and explain what it contains.

When using synthetic data, describe how the synthetic data are constructed and set random seeds whenever stochastic generation is involved.

For example:

```rst
.. code-block:: pycon

   >>> import numpy as np
   >>> rng = np.random.default_rng(42)
```

If an algorithm itself is stochastic, document any random-state or seed parameter required to reproduce the result.

Users should be able to understand why their result might differ from the documentation.

---

## 15. Numerical Results

Important numerical outputs should include sufficient context.

For example, do not simply write:

```text
0.947
```

Explain the quantity:

```text
Confidence ratio: 0.947
```

and then interpret the value in the text.

Where appropriate, specify:

- units;
- normalization;
- expected range;
- threshold interpretation;
- numerical tolerance;
- dependence on the input dataset.

Avoid reporting excessive decimal precision when it has no scientific meaning.

---

## 16. Theory Pages

Pages under the theory or scientific-background documentation require particular attention to scientific rigor.

A theory page should normally progress through:

1. physical or methodological motivation;
2. relevant quantities and assumptions;
3. mathematical relationships;
4. interpretation of those relationships;
5. connection to the implementation in pyCSAMT;
6. practical implications or limitations.

These elements do not necessarily need separate headings.

The objective is a coherent narrative rather than a rigid template.

Assumptions must be stated explicitly when they influence the validity of an equation or algorithm.

For example:

- 1-D earth assumptions;
- isotropy;
- stationarity;
- plane-wave approximation;
- smoothness assumptions;
- boundary conditions;
- regularization assumptions;
- statistical independence;
- Gaussian noise assumptions.

---

## 17. Units and Symbols

Physical quantities should use consistent units.

Where appropriate, use SI units and state them explicitly:

- frequency \(f\): Hz;
- angular frequency \(\omega\): rad s\(^{-1}\);
- resistivity \(\rho\): \(\Omega\,\mathrm{m}\);
- conductivity \(\sigma\): S m\(^{-1}\);
- electric field: V m\(^{-1}\);
- magnetic field: A m\(^{-1}\);
- distance or depth: m or km.

If the API accepts different units, document the required input unit and the unit of the returned value.

Mathematical notation and API variable names do not always need to be identical, but their relationship should be clear.

---

## 18. Parameters and Defaults

Important parameters should not merely be listed.

Explain their effect on the calculation.

For example, for an inversion regularization parameter, describe:

- what term it controls;
- what happens when its value increases;
- what happens when it decreases;
- the default used by pyCSAMT;
- whether the default is data-dependent;
- when users should consider changing it.

Whenever possible, connect parameter explanations to the underlying mathematical expression.

---

## 19. Warnings, Limitations, and Scientific Preconditions

Important scientific limitations should appear close to the relevant method.

Use Sphinx admonitions where appropriate:

```rst
.. note::

   ...
```

```rst
.. warning::

   ...
```

```rst
.. important::

   ...
```

Warnings should explain genuine scientific, computational, or reproducibility concerns rather than obvious statements.

Examples include:

- invalid dimensional assumptions;
- insufficient frequency coverage;
- poorly conditioned inversion;
- unreliable extrapolation;
- missing optional dependencies;
- interpretation limits of AI-derived models;
- sensitivity to initialization or regularization;
- preprocessing requirements.

---

## 20. Relationship Between Explanation, Code, and Result

A strong documentation example should follow a natural sequence:

**scientific motivation → method → executable example → observed result → interpretation**

For example:

```text
Why static-shift correction is needed
        ↓
Physical explanation
        ↓
Relevant mathematical relationship
        ↓
pyCSAMT implementation
        ↓
Executable example
        ↓
Captured numerical output or figure
        ↓
Interpretation of the result
        ↓
Practical recommendation or limitation
```

This sequence should emerge naturally in the prose. It should not be reproduced as a rigid set of headings on every page.

---

## 21. Existing Outputs and Figures

Before generating new outputs, check whether the documentation already contains:

- the expected console output;
- the required figure;
- an equivalent validated figure;
- an existing script producing the same result.

If the existing material is still accurate, reuse it.

Do not regenerate figures or duplicate content without a reason.

Likewise, avoid creating documentation galleries simply to collect generated figures. Figures should normally remain close to the examples and scientific explanations they support.

---

## 22. Documentation Scripts

The `docs/scripts/` directory should contain scripts that support documentation generation, especially when:

- an example is too long for the page;
- several figures must be reproduced;
- a complex scientific demonstration requires many steps;
- deterministic generation of documentation assets is desirable.

Scripts should be written so they can be executed independently when possible.

A documentation script should preferably:

- contain a clear module docstring;
- use public pyCSAMT interfaces;
- expose logically named functions;
- use deterministic seeds where required;
- save figures explicitly;
- avoid machine-specific absolute paths;
- create required output directories safely.

For example:

```text
docs/scripts/
├── generate_static_shift_figures.py
├── generate_edi_processing_example.py
├── generate_ai_inversion_figures.py
└── generate_resistivity_section.py
```

---

## 23. Image Naming

Use descriptive and stable filenames.

Prefer:

```text
static_shift_before_after.png
occam_1d_convergence.png
edi_phase_resistivity_example.png
ai_inversion_regularization_tradeoff.png
```

Avoid:

```text
fig1.png
test2.png
new_plot.png
output_final2.png
```

The filename should convey the scientific content even when viewed outside the documentation page.

---

## 24. Cross-References

Take advantage of Sphinx cross-referencing rather than repeating the same explanation throughout the documentation.

Use references to:

- glossary terms;
- equations;
- figures;
- sections;
- API objects;
- classes;
- methods;
- related tutorials.

For example:

```rst
See :term:`static shift` for the general definition.
```

or:

```rst
The correction factor follows the relationship introduced in
:eq:`eq-static-shift-factor`.
```

or:

```rst
See :ref:`user-guide-edi-processing` for the complete processing workflow.
```

Cross-references should help users navigate the documentation without forcing duplicated explanations across multiple pages.

---

## 25. Final Quality Check for Each Page

Before considering a documentation page complete, verify that:

- technical terminology is defined or linked appropriately;
- important scientific claims are explained;
- mathematical expressions are reproducible;
- symbols are defined consistently;
- important equations on theory pages have labels;
- code uses public and realistic APIs;
- interactive examples use `pycon` where appropriate;
- executable examples show their actual outputs;
- generated plots are included when they improve understanding;
- figures are scientifically interpreted;
- related figures are combined efficiently when appropriate;
- long source code has been moved to `docs/scripts/`;
- long scripts are exposed using `.. code-dropdown::`;
- documentation scripts are not imported in user-facing examples;
- units and assumptions are explicit;
- parameters are explained rather than merely listed;
- relevant limitations or warnings are stated;
- existing valid outputs are not unnecessarily duplicated;
- the prose flows naturally from concept to implementation and interpretation.

The final documentation should allow a user not only to reproduce a pyCSAMT workflow, but also to understand the scientific reasoning behind the workflow and correctly interpret the results.