.. _api:

API reference
=============

This is the public class and function reference for :mod:`pycsamt`.
Use the search field to find an object by name, package, or description.  For
conceptual guidance and worked examples, see the :doc:`../user_guide/index`.
Runtime settings and the ``configure_*`` / ``reset_*`` pattern are documented
separately in :doc:`../api_guide/index`.  Private objects and test modules are
not included.  Each entry links to the object's canonical generated reference;
public facade paths are used when an object is implemented in a private module.
The owning public module is shown beneath its description and is searchable;
only supported public import paths are displayed.  The catalogue is rebuilt
from each subpackage's current public exports, including the IoT field stack.

.. raw:: html

   <div id="api-search" class="api-search" role="search" aria-label="Search the public API">
     <div class="api-search__limit">
       <label for="api-limit-filter" class="visually-hidden">Entries per page</label>
       <select id="api-limit-filter">
         <option value="all">All</option>
         <option value="25">25</option>
         <option value="50">50</option>
         <option value="100">100</option>
       </select>
       <span>entries per page</span>
     </div>
     <div class="api-search__field">
       <label for="api-search-input">Search:</label>
       <input id="api-search-input" type="search" autocomplete="off" spellcheck="false">
     </div>
     <p id="api-search-status" class="api-search__status" aria-live="polite"></p>
     <nav id="api-pagination" class="api-pagination" aria-label="API table pages"></nav>
   </div>

.. public-api-index::

   pycsamt.api
   pycsamt.core
   pycsamt.io
   pycsamt.metadata
   pycsamt.site
   pycsamt.seg
   pycsamt.z
   pycsamt.jones
   pycsamt.zonge
   pycsamt.stratagem
   pycsamt.emtools
   pycsamt.transformers
   pycsamt.interp
   pycsamt.tdem
   pycsamt.utils
   pycsamt.forward
   pycsamt.inversion
   pycsamt.models
   pycsamt.backends
   pycsamt.ai
   pycsamt.agents
   pycsamt.pipeline
   pycsamt.session
   pycsamt.topo
   pycsamt.map
   pycsamt.gis
   pycsamt.iot
   pycsamt.app

.. toctree::
   :maxdepth: 1
   :hidden:

   pycsamt
   api
   core
   io
   metadata
   site
   seg
   z
   jones
   zonge
   stratagem
   emtools
   transformers
   interp
   tdem
   utils
   forward
   inversion
   models
   backends
   ai
   agents
   pipeline
   session
   topo
   map
   gis
   iot
   app
