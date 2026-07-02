Configuration Commands
======================

``pycsamt config`` is the command group for making the CLI behave the
way you want every time you open a terminal. It stores user preferences
in ``~/.pycsamt.toml`` and loads them when pyCSAMT starts.

Use it for settings that are not part of the survey itself:

* terminal behaviour, such as verbosity, colors, default output format,
  and worker count;
* figure export defaults, such as DPI, file format, save directory, and
  transparent backgrounds;
* plotting and interpretation style presets;
* pipeline runtime defaults, such as error policy and progress display;
* AI-agent provider preferences and budget limits.

Configuration is intentionally separate from survey data. Changing
``~/.pycsamt.toml`` does not edit EDI files, AVG files, inversion
directories, or generated results.

When to use config
------------------

Use ``pycsamt config`` when a preference should follow you across many
commands and projects. For example, if every plot for a report should be
saved as PDF at 300 DPI, configure that once:

.. code-block:: console

   pycsamt config set plot.fmt pdf
   pycsamt config set plot.dpi 300

Use command-line options when a value is specific to one run:

.. code-block:: console

   pycsamt pipe run --config workflow.yaml --dpi 150 --plot-fmt png

Use environment variables when a value is secret, machine-specific, or
better controlled outside pyCSAMT. API keys belong in environment
variables, not in ``~/.pycsamt.toml``.

How values are resolved
-----------------------

The configuration system combines several sources:

1. Package defaults provide safe baseline behaviour.
2. ``~/.pycsamt.toml`` stores your persistent user overrides.
3. Environment variables override or complement file settings where
   supported.
4. Explicit command-line options override behaviour for the current
   command invocation.

Two commands are useful because they answer different questions:

``pycsamt config list``
    Shows what is persisted in ``~/.pycsamt.toml``. If a value is coming
    from an environment variable or from package defaults, it may not
    appear here.

``pycsamt config show``
    Shows the current effective configuration after pyCSAMT has loaded
    defaults, the TOML file, and environment variables.

Start with ``show`` when you want to understand how the CLI will behave
right now:

.. code-block:: console

   pycsamt config show
   pycsamt config show plot
   pycsamt config show agent --format json

The config file
---------------

User configuration is stored at:

.. code-block:: text

   ~/.pycsamt.toml

You normally do not need to edit this file by hand. The CLI validates and
writes it for you:

.. code-block:: console

   pycsamt config set plot.dpi 300

The resulting TOML is simple. A file may look like this:

.. code-block:: toml

   [plot]
   fmt = "pdf"
   dpi = 300
   savedir = "figures"

   [style]
   preset = "publication"

   [pipe]
   on_step_error = "warn"

   [agent]
   provider = "openai"
   budget_usd = 5.0

Key syntax
----------

The CLI uses dot notation:

.. code-block:: console

   pycsamt config set plot.dpi 300
   pycsamt config get plot.dpi

The first part is the section name. Everything after it is the setting.
Nested settings are also written with dots:

.. code-block:: console

   pycsamt config set control.rho.view log10
   pycsamt config set control.phase.wrap true

Internally, nested TOML keys are flattened with double underscores. The
user-facing command remains the dotted form.

Value types
-----------

``config set`` receives values as text and converts simple scalar types:

.. list-table::
   :header-rows: 1
   :widths: 28 36 36

   * - Value typed
     - Stored as
     - Example
   * - ``true``, ``yes``, ``on``
     - Boolean ``True``
     - ``control.phase.wrap true``
   * - ``false``, ``no``, ``off``
     - Boolean ``False``
     - ``cli.no_color false``
   * - ``300``
     - Integer
     - ``plot.dpi 300``
   * - ``5.0``
     - Float
     - ``agent.budget_usd 5.0``
   * - ``pdf``
     - String
     - ``plot.fmt pdf``

Use ``--dry-run`` when you want pyCSAMT to validate a value without
writing it:

.. code-block:: console

   pycsamt config set plot.dpi 300 --dry-run

Sections
--------

Configuration is organized by section. The section name is always the
first part of a key, as in ``plot.dpi`` or ``pipe.on_step_error``.

``cli``
    Terminal behaviour: verbosity, color, output format, default output
    directory, and parallel worker count. These are general defaults used
    across command groups.

``control``
    Plot-axis and display controls used by geophysical plots, such as
    apparent-resistivity scaling, phase handling, and x-axis views.

``style``
    Visual style for MT components, line/marker choices, and common
    plotting presets.

``section_view``
    Section-plot layout choices, including figure geometry, colorbar
    behaviour, and axis direction.

``station``
    Station tick and marker rendering: density, rotation, label style,
    and station display controls.

``interp``
    Geological interpretation display defaults, including interpretation
    figure styles and color palettes.

``plot``
    Figure export settings: file format, DPI, save directory, bounding
    box, and transparency.

``agent``
    AI-agent preferences such as provider, model, and budget. API keys
    are not stored here; use environment variables for secrets.

``pipe``
    Processing pipeline runtime settings, including error behaviour,
    progress display, plotting defaults, and output policy.

``view``
    Tabular view backend, for example pyCSAMT's native view or pandas.

Inspecting configuration
------------------------

List only values that are persisted in the TOML file:

.. code-block:: console

   pycsamt config list
   pycsamt config list plot
   pycsamt config list agent --format json

Read one effective value:

.. code-block:: console

   pycsamt config get plot.dpi
   pycsamt config get agent.provider
   pycsamt config get control.rho.view

Show the full effective state:

.. code-block:: console

   pycsamt config show
   pycsamt config show pipe
   pycsamt config show plot --format json

If ``config get`` says a key is not found, check the section name and the
actual attribute path. ``config show SECTION`` is the easiest way to see
what that section exposes.

Editing configuration
---------------------

Set a value:

.. code-block:: console

   pycsamt config set plot.dpi 300
   pycsamt config set plot.fmt pdf
   pycsamt config set plot.savedir figures

Remove one user override:

.. code-block:: console

   pycsamt config unset plot.dpi

Reset one section:

.. code-block:: console

   pycsamt config reset plot --yes

Reset all user configuration:

.. code-block:: console

   pycsamt config reset --yes

Resetting removes values from ``~/.pycsamt.toml``. It does not unset
environment variables. If an environment variable still exists, it can
continue to affect the effective configuration.

Style presets
-------------

Style presets are shortcuts for common plotting looks. They are easier
and safer than setting many visual keys one by one.

.. code-block:: console

   pycsamt config style pycsamt
   pycsamt config style publication
   pycsamt config style dark
   pycsamt config style modem

Use ``publication`` when preparing figures for reports or papers. Use
``dark`` for dark-background notebooks and presentations. Use ``modem``
when you want plotting conventions that sit well beside ModEM-style
outputs.

To test a preset without saving it:

.. code-block:: console

   pycsamt config style publication --no-persist

Interpretation presets
----------------------

Interpretation presets configure geological interpretation displays:

.. code-block:: console

   pycsamt config interp default
   pycsamt config interp publication
   pycsamt config interp dark
   pycsamt config interp accessible

Use ``accessible`` when figures must remain readable for color-vision
deficiency. Use ``publication`` when exporting interpretation products to
reports.

As with style presets, ``--no-persist`` applies the preset only for the
current process:

.. code-block:: console

   pycsamt config interp accessible --no-persist

Environment variables
---------------------

Environment variables are best for machine-level defaults and secrets.
Display all recognized variables with:

.. code-block:: console

   pycsamt config env
   pycsamt config env --section plot
   pycsamt config env --section agent

Common CLI and plotting variables:

.. list-table::
   :header-rows: 1
   :widths: 38 62

   * - Variable
     - Meaning
   * - ``PYCSAMT_VERBOSE``
     - Verbosity level: ``0`` quiet, ``1`` info, ``2`` debug.
   * - ``PYCSAMT_NO_COLOR``
     - Disable ANSI color when set.
   * - ``PYCSAMT_OUTPUT``
     - Default output format: ``text``, ``json``, or ``csv``.
   * - ``PYCSAMT_OUTPUT_DIR``
     - Default output directory.
   * - ``PYCSAMT_JOBS``
     - Default parallel worker count.
   * - ``PYCSAMT_FMT``
     - Figure export format, such as ``png``, ``pdf``, or ``svg``.
   * - ``PYCSAMT_DPI``
     - Raster figure DPI.
   * - ``PYCSAMT_SAVEDIR``
     - Default directory for saved figures.
   * - ``PYCSAMT_TRANSPARENT``
     - Transparent figure background when set.
   * - ``PYCSAMT_BBOX``
     - Figure bounding-box mode, for example ``tight``.

Agent key variables:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Variable
     - Provider
   * - ``OPENAI_API_KEY``
     - OpenAI
   * - ``PYCSAMT_OPENAI_API_KEY``
     - OpenAI alias
   * - ``ANTHROPIC_API_KEY``
     - Claude / Anthropic
   * - ``PYCSAMT_CLAUDE_API_KEY``
     - Claude alias
   * - ``GOOGLE_API_KEY``
     - Gemini
   * - ``GOOGLE_GENERATIVEAI_API_KEY``
     - Gemini alias
   * - ``PYCSAMT_GEMINI_API_KEY``
     - Gemini alias

Agent configuration
-------------------

The ``agent`` configuration section stores provider preferences, not API
keys. This is deliberate: keys should stay in environment variables so
they do not end up in dotfiles, commits, or logs.

Check agent status:

.. code-block:: console

   pycsamt config agent status
   pycsamt config agent status --format json

Configure a provider preference:

.. code-block:: console

   pycsamt config set agent.provider openai
   pycsamt config set agent.budget_usd 5.0

Get key setup instructions:

.. code-block:: console

   pycsamt config agent set-key openai
   pycsamt config agent set-key claude
   pycsamt config agent set-key gemini

The ``set-key`` command prints the environment variable to use. It does
not write the secret into ``~/.pycsamt.toml``.

Common recipes
--------------

Publication plotting defaults:

.. code-block:: console

   pycsamt config style publication
   pycsamt config set plot.fmt pdf
   pycsamt config set plot.dpi 300
   pycsamt config set plot.savedir figures

Quiet scripting defaults:

.. code-block:: console

   pycsamt config set cli.no_color true
   pycsamt config set cli.output json
   pycsamt config set cli.verbose 0

Pipeline debugging defaults:

.. code-block:: console

   pycsamt config set pipe.on_step_error warn
   pycsamt config set pipe.show_progress true
   pycsamt config show pipe

Agent setup:

.. code-block:: console

   pycsamt config agent set-key openai
   pycsamt config set agent.provider openai
   pycsamt config set agent.budget_usd 5.0
   pycsamt config agent status

Troubleshooting
---------------

``config list`` is empty but pyCSAMT still has settings
    The values may be package defaults or environment variables. Use
    ``pycsamt config show`` to inspect the effective configuration.

Changing a value did not affect one command
    Check whether that command was called with an explicit option, such
    as ``--dpi`` or ``--format``. Command-line options are intended to
    override user defaults for that run.

Resetting did not remove a value
    Look for an environment variable. ``pycsamt config reset`` edits only
    ``~/.pycsamt.toml``.

An API key is not detected
    Run ``pycsamt config env --section agent`` and confirm that the
    correct variable is set in the shell where you run pyCSAMT.

You are unsure which key name to use
    Run ``pycsamt config show SECTION`` first, then set one value at a
    time with ``--dry-run`` before persisting it.

