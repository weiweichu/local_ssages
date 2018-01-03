.. _umbrella-sampling:

Umbrella Sampling
-----------------

Introduction
^^^^^^^^^^^^

Calculations of thermodynamic data and other properties rely on proper sampling of the configurational space. 
However, the presence of energy barriers can prevent certain configurations from being sampled properly or even sampled at 
all. Umbrella sampling is a simulation technique that helps to overcome those barriers and improve sampling 
by applying a bias along a specified collective variable. The bias takes the form of a harmonic potential and is typically constant throughout a simulation. 
Usually, a series of umbrella-sampled simulations are performed and analyzed together using the weighted histogram analysis method 
(WHAM).

The functional form of the artificial bias is

.. math::

    U_{umbrella} = \frac{1}{2} k \left(s - s_0\right)^2

where :math:`k` is the spring constant, :math:`s` is the current value of the collective variable and :math:`s_0` is the desired value of the collective variable.
The success of the umbrella sampling method depends on the correct choice of :math:`k` and :math:`s_0` for the different simulations. Suitable values of :math:`k` and :math:`s_0` depend on the distances between adjacent umbrellas, and the gradient of the free energy surface at :math:`s_0`.

For more information about the umbrella sampling method, the interested reader is referred  to Ref. [1]_.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

The following parameters need to be set under `"method"` in the JSON input file:

.. code-block:: javascript

    "type" : "Umbrella"
	
The following options are available for Umbrella Sampling:

centers (required)
	Array of target values for each CV. This can either be an array of numbers, or 
	an array of an array of numbers, for multiple walkers. If only a single array is 
	defined, it is used across all walkers.

ksprings (required)
	Array of spring constants for each CV. This can either be an array of numbers, or 
	an array of an array of numbers, for multiple walkers. If only a single array is 
	defined, it is used across all walkers.

output_file (required)
	Output file name for umbrella sampling data. This can be a string or 
	an array of strings for multiple walkers.

output_frequency (optional)
	Frequency of writing to output file (default: 1)

append (optional)
	Boolean value which will cause umbrella sampling to either append to 
	an existing file or override it entirely. 

The umbrella sampling method can also be run with time dependent centers. 
This is equivalent to the "steered MD" method. To do so, omit `"centers"` 
and use the following options instead: 

centers0 (required)
	Array of initial target values for each CV. This can either be an array of numbers, or 
	an array of an array of numbers, for multiple walkers. If only a single array is 
	defined, it is used across all walkers. 

centers1 (required)
	Array of final target values for each CV. This can either be an array of numbers, or 
	an array of an array of numbers, for multiple walkers. If only a single array is 
	defined, it is used across all walkers. 

timesteps (required)
	The number of timesteps over which to scale the umbrella centers 

By setting the above options, the umbrella method will linearly interpolate 
between `centers0` and `centers1` in `timesteps` iterations. 

.. _Umbrella_tutorial:

Tutorial
^^^^^^^^

This tutorial will go through running Umbrella Sampling on an atomistic model of butane using LAMMPS as the 
MD engine. 
Umbrella sampling will be performed on the torsional CV of the butane C atoms. 

The butane implementation in LAMMPS requires several modules to be added before being linked to SSAGES.
To do this, return to your LAMMPS ``src`` directory and issue the following commands

.. code-block:: bash

    make yes-molecule
    make yes-kspace
    make

which add the ``MOLECULE`` and ``KSPACE`` LAMMPS packages. Additional information about building
LAMMPS with optional packages can be found in the
`LAMMPS documentation <http://lammps.sandia.gov/doc/Section_start.html#start-3>`_.


The files for running this example can 
be found in ``Examples/User/Umbrella`` and consist of the following files:

``Butane_SSAGES.in``
	LAMMPS input file

``Butane.data``
	LAMMPS data file describing butane molecule.

``umbrella_input.json``
	Template JSON input containing information for one Umbrella Sampling simulation. 

``umbrella_multiwalker_gen.py``
    Python script for creating `multiwalker_umbrella.json` input file. The total number of
    simulations and the 'centers' values are controlled in this file.
	
Once in the directory, the appropriate `.json` file needs to be generated. A `.json` file
is already in the directory, ``umbrella_input.json``, which contains the CV information
and specifies the LAMMPS input files to be used. A single-walker umbrella simulation can
be run directly using

.. code-block:: bash

    ssages umbrella_input.json

The simulation will create an output file named `umbrella.dat1` containing the value of
the CV and the target value (the center) every 100 timesteps. From this histogram, the
local free energy can be calculated.

While it is possible to run Umbrella sampling using a simle walker, typically multiple
walker (multiple umbrellas) are simulated. For multiwalker Umbrella sampling of butane,
you can generate an input file using the `umbrella_multiwalker_gen.py` script via

.. code-block:: bash

	python umbrella_multiwalker_gen.py

This will generate an input file called ``multiwalker_umbrella.json`` containing the
information from ``umbrella_input.json`` duplicated 12 times with varying values of
``centers``. These values correspond to the target values of the torsional CV. 

To run multiwalker SSAGES issue the command:

.. code-block:: bash 

	mpiexec -np 12 /path/to/SSAGES/build/ssages multiwaler_umbrella.json
	
This will run 12 different umbrella sampling simulations simultaneously.
Ideally, this example will be run in computing environment where each process can run
on a different processor. The example will still work if run on a users local desktop
or laptop machine, but the runtime of the code will be very large.

During the simulation 12 different output files will be generated, each containing the
iteration, target value of the corresponding 'center' CV,  and the value of the CV at
the iteration number. 

These output files can then be used to construct a complete free energy surface using
the WHAM algorithm [2]_. Though SSAGES does not currently contain its own implementation
of WHAM, there are many implementations available, such as that provided by the
Grossfield Lab [3]_.

References
^^^^^^^^^^

.. [1] Kästner, J. (2011). *Umbrella sampling*. Wiley Interdiscip Rev Comput Mol Sci, 1(6), 932–942. 

.. [2] Kumar, S., Rosenberg, J., & Bouzida, D. (1992). The weighted histogram analysis method for free‐energy calculations on biomolecules. I. The method. Journal of Computational Chemistry, 13(8), 1011–1021. 

.. [3] Grossfield, A. WHAM: the weighted histogram analysis method. `http://membrane.urmc.rochester.edu/content/wham <http://membrane.urmc.rochester.edu/content/wham>`_

