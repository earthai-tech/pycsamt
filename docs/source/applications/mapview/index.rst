.. _applications-mapview:

MapView
=======

MapView is the dedicated map workbench of pyCSAMT: a Dash/Plotly browser
application focused entirely on **seeing a survey in space** — basemap station
maps, profile overlays, Surfer-style contours, and interactive 3-D resistivity
scenes with fence sections, blocks, depth slices, iso-surfaces, and draped
topography. It is the point-and-click layer over the :doc:`Python mapping API
</user_guide/map/index>`: everything the app draws comes from the same
:class:`~pycsamt.map.MapView` façade you can drive from code, so a scene you
build here can be reproduced exactly in a script.

.. raw:: html

   <a class="pycsamt-video-card"
      href="https://youtu.be/DcDA3THoFtw" target="_blank" rel="noopener">
     <span class="pycsamt-video-thumb">
       <img src="/_static/applications/mapview/mapview-demo-cover.jpg"
            alt="MapView video demo thumbnail" loading="lazy">
       <span class="pycsamt-video-play" aria-hidden="true">
         <svg viewBox="0 0 68 48" xmlns="http://www.w3.org/2000/svg">
           <rect x="0" y="0" width="68" height="48" rx="12"
                 fill="rgba(20,20,20,0.75)"></rect>
           <path d="M27 15l18 9-18 9V15z" fill="#fff"></path>
         </svg>
       </span>
     </span>
     <span class="pycsamt-video-caption">
       <strong>Watch the demo</strong> &mdash; pyCSAMT v2.6 as an
       integrated geoscience interpretation platform: EM processing and
       inversion, 3-D visualization, boreholes, geological information,
       and the USGS AI Pack in one workflow, illustrated on a real
       Cu&ndash;Mo (copper&ndash;molybdenum) mineralization case study.
       Opens on YouTube in a new tab.
     </span>
   </a>

The case study shown is the audio-magnetotelluric survey and 2-D
inversion of [Kouabena2025]_ (bundled as ``data/AMT/WILLY_DATA/``; see
:doc:`/user_guide/models/occam2d` and the dataset's own ``README.md``).
Prefer a shorter, silent tour instead? See the animated walkthrough on
the :doc:`overview` page.

A good first pass, in order, is **Overview → Installation → Loading &
Sessions → Views & Controls → Exports → Troubleshooting**.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   overview
   installation
   loading_and_sessions
   views
   exports
   troubleshooting
