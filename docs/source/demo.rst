Demo 
===== 

Available filters
-----------------

Electromagnetic Array Profiling (EMAP) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Trimming moving average (TMA) mostly used by `Zonge International Engineering <http://zonge.com/>`_. 
* Fixed-length-dipole moving average (FLMA) also implemented by `Zonge international Engineering <https://zonge.com.au/>`_. 
* Adaptative moving-average (AMA) based on idea of `Torres-Verdin <https://sci-hub.se/http://dx.doi.org/10.1190/1.1443273>`_.

Magnetotellurics
^^^^^^^^^^^^^^^^^

* MT Removal distorsion (dist)
* Static shift removal (ss)

For example, to remove the static shift effect into a corrupted `SEG <https://seg.org/Default.aspx?TabId=176&language=en-US>`_ EDI data, 
(e.g.,FLMA) user can apply one of above filters FLMA) to correct the EDI data located in data/edi directory as the below command line (CLI) indicates.

.. code-block:: bash 

   $ correctedi data/edi -ft flma --ndipole 5 --dipole-length=50  
   
   
Plot inversion and geostratigraphy(G) misfits
----------------------------------------------
For this quickstart, we assume the user has already the forward modeling files (``*.resp``, ``*.dat``, ``*.mesh``, ``*.iter`` and ``*.logfile``(optional)),
 and the supplement data (i.e. boreholes/wells and geological data) collected in the survey area.
 
1.Plotting some fitting curves of resistivity and phase inversion after applying on observed data the static shift correction.
  For instance, we visualize the fitting curves of four survey lines with their
  corresponding RMS. The stations (e.g. ``S00``,``S04``, ``s08``, and ``S12``) are randomly chosen of each survey line: 
 
 .. code-block:: 
    
	$ fitforward  data/inversionFiles --stations S00 S04 s08 S12 --rms 1.013 1.451 1. 1.069 
	
 Click `here<https://github.com/WEgeophysics/pyCSAMT/blob/develop/examples/examplefitcurves.png>`_ to see the reference output.

2. Plotting the ``misfit`` of the model response from the forward modeling (FM) of line ``K1``. Set kind argument to phase for the error in phase.

 .. code-block:: bash 
 
    $ misfit2d data/inversionFiles/K1.dat data/inversionFiles/K1.resp --kind=rho 
	
 To see the output, click on the   `reference output<https://github.com/WEgeophysics/pyCSAMT/blob/develop/examples/misfit.png>`_.
3. Construction of the new resistivity model (NM) also called the geostratigraphy model. The latter is built from
 the additional data (e.g. borehole/well data and/or geological data) collected in the exploration area combined with 
 the FM. The best approach to fastly build the NM model is to gather all the data collected in the exploration area into a 
 single configuration file in ``*.json`` or ``*.yml`` format. For illustrating, we will use ``mysurveyarea_data.yml`` file like:

 .. code-block:: bash 
	
	# mysurveyarea_data.yml
	input_layers:
		- river water
		- sedimentary rocks
		- fracture zone
		- gravel
		- granite
		- igneous rocks
		- basement rocks
	input_resistivities:
		- 10
		- 66
		- 70
		- 180
		- 1000
		- 3000
		- 7000
	data_fn: data/occam2D\K1.dat
	iter_fn: data/occam2D\K1.iter
	mesh_fn: data/occam2D\Occam2DMesh
	model_fn: data/occam2D\Occam2DModel
	ptol: 0.1
	beta: 5
	n_epochs: 100
	build: true

 where ``(data_fn, iter_fn, mesh_fn ,model_fn)`` and ``(ptol, beta, n_epochs, build)`` are OCCAM2D inversion files resulting from FM (CRM) 
 and the constructor parameters respectively. From the CLI below, NM is created.
 
 .. code-block:: bash 
 
    $ nm -c mysurveyarea_data.json --show --misfit
	
 Indeed, the ``Misfit G`` helps to determine whether the different layers with their corresponding resistivity values 
 are misclassified or not. An example of misfit G map can be visualized here.

 **Note**: For CSAMT data processing and an example of a deep implementation, please refer to our 
 `wiki page<https://github.com/WEgeophysics/pyCSAMT/wiki/How-pyCSAMT-works-%3F>`_.


Plot the pseudostratigraphic log
--------------------------------

Once the NM is built, user just needs to query the software parameters memory (SPM) to fech at each station the
 corresponding pseudostratigraphic log. An example of the command to do this task
 (e.g. extracting the log of station ``S00``) is:

.. code-block:: bash 

   $ pseudostratigraphic -s S00 
   
The output below with layer thicknesses estimation are displayed.

.. code-block:: bash 

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ PseudoStratigraphic Details: Station = S00 ]~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	------------------------------------------------------------------------------------------------------
	|      Rank |            Stratum             |         Thick-range(m)         |     Thickness(m)     |
	------------------------------------------------------------------------------------------------------
	|        1. |         fracture zone          |         0.0 ----- 6.0          |         6.0          |
	|        2. |             gravel             |         6.0 ----- 13.0         |         7.0          |
	|        3. |            granite             |        13.0 ----- 29.0         |         16.0         |
	|        4. |         igneous rocks          |        29.0 ----- 49.0         |         20.0         |
	|        5. |         basement rocks         |        49.0 ----- 249.0        |        200.0         |
	|        6. |         igneous rocks          |       249.0 ----- 289.0        |         40.0         |
	|        7. |            granite             |       289.0 ----- 529.0        |        240.0         |
	|        8. |         igneous rocks          |       529.0 ----- 699.0        |        170.0         |
	|        9. |         basement rocks         |       699.0 ----- 999.0        |        300.0         |
	------------------------------------------------------------------------------------------------------
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Survey Line: Occam2D files properties ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	|model = Occam2DModel     |iter  = ITER17.iter      |mesh  = Occam2DMesh      |data  = OccamDataFile.dat|
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
Click `here <https://github.com/WEgeophysics/pyCSAMT/blob/develop/examples/pseudostratigraphic_log.PNG>`_ to see the predicted log.
Obviously, it does not make sense to expect to drill until to reach ``1km`` depth. Therefore, another feature is implemented to help 
the user to only fetch from the SPM the most interesting part of the predicted log for a specific purpose. To do such task,
 one needs to fiddle with the ``zoom`` parameter. For instance, the CLI below with ``zoom=25%`` only displays the first ``250m`` assuming 
 that the investigation depth is ``1000m`` maximum.
 
.. code-block:: bash 
   
   $ pseudostratigraphic --station=S00 --zoom 25%
   
Check the following `ouput<https://github.com/WEgeophysics/pyCSAMT/blob/develop/examples/zoom25.PNG>`_ to see the new log. Futhermore,
 it's also possible to provide the top (e.g. ``10m``) and the bottom(e.g. ``120m``) of the log for visualization as:
 
.. code-block:: bash 

   $ pseudostratigraphic --station S00 --zoom 10 120 --fontsize 12
   
   
