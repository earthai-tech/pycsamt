.. _applications-hosted:

Hosted Applications
====================

Three of the four pyCSAMT application surfaces also run as hosted
instances — the same Dash apps described in :doc:`overview`, already
running on a server, reachable straight from a browser. No install, no
Python environment, no local server to keep open.

Each card below previews its app with the same walkthrough animation used
in its guide. The secured, per-user hosted rollout is still in progress
(see the note below the table) — click a card for its current status.

.. grid:: 1 1 3 3
   :gutter: 3
   :class-container: pyc-hosted-grid

   .. grid-item-card:: Agent Master
      :img-top: /_static/applications/agent_master/agent-master-walkthrough.gif
      :class-card: sd-card-hover
      :link: https://agent-master.pycsamt.org
      :link-type: url

      Conversational workflow surface — load EDI data, chat through a
      pipeline, and collect reports and figures.

      .. container:: pyc-hosted-install

         :doc:`agent_master/installation`

      +++
      Open Agent Master →

   .. grid-item-card:: Web App
      :img-top: /_static/applications/web/web-walkthrough.gif
      :class-card: sd-card-hover
      :link: https://web.pycsamt.org
      :link-type: url

      Browser dashboard for loading, processing, and visualizing EDI
      surveys end to end.

      .. container:: pyc-hosted-install

         :doc:`web/installation`

      +++
      Open Web App →

   .. grid-item-card:: MapView
      :img-top: /_static/applications/mapview/mapview-walkthrough.gif
      :class-card: sd-card-hover
      :link: https://mapview.pycsamt.org
      :link-type: url

      Dedicated map workbench for viewing survey lines and stations in
      space.

      .. container:: pyc-hosted-install

         :doc:`mapview/installation`

      +++
      Open MapView →

The desktop GUI has no hosted equivalent — it is a native application
installed on your own machine; see :doc:`desktop/index`.

Choosing Hosted Or Local
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * -
     - Hosted instance
     - Your own install / server
   * - Setup
     - None — open the link
     - ``pip install "pycsamt[app]"``, then launch
   * - Where data lives
     - The shared host running the instance
     - Your machine, or a server you control
   * - Who can see loaded data
     - Anyone using the same running instance
     - Only you, or your trusted network / team
   * - Good for
     - A first look, a demo, trying a workflow
     - Real survey data, repeated work, private or team use

If you decide local or self-hosted is the better fit, each guide's
*Installation And Launch* page covers it: :doc:`agent_master/installation`,
:doc:`web/installation`, :doc:`mapview/installation`. For running your own
instance behind a reverse proxy and authentication, see
:doc:`web/deployment`.

.. caution::

   Hosted instances are shared, unauthenticated demo servers, not a
   private workspace.

   * **Do not load sensitive or confidential survey data** through them —
     treat anything you load as visible to anyone else using the same
     running instance at the same time.
   * There is no login and no per-user isolation: loaded surveys,
     generated figures, and reports are held in server-side state shared
     by whoever is currently using that instance, the same way a single
     local install would be shared on a trusted LAN (see
     :doc:`web/deployment`).
   * AI API keys, when a hosted app asks for one, are kept in your own
     browser (``localStorage``) and are never sent to or stored on the
     server — see *Where State Lives* in :doc:`web/deployment` for the
     full split between browser-held and server-held state.
   * Availability is best-effort. For anything you depend on, run your
     own instance instead.

.. note::

   A secured deployment — authentication and per-user data isolation — is
   in progress for these hosted instances. Until it ships, treat them as
   shared demo servers per the caution above, and check back on this page
   for the update.

Next Steps
----------

* :doc:`overview` -- what each application surface is for.
* :doc:`web/deployment` -- host, port, and security considerations for
  running your own instance.
