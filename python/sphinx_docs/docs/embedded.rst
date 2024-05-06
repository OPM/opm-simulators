OPM Embedded Python Documentation
=================================

The PYACTION keyword is a flow specific keyword which allows for executing embedded Python
code in the SCHEDULE section. The embedded Python code will then be executed at the end of each successful timestep.

The PYACTION keyword is inspired
by the ACTIONX keyword, but instead of a .DATA formatted condition you
are allowed to implement the condition with a general Python script. The
ACTIONX keywords are very clearly separated in a condition part and an
action part in the form of a list of keywords which are effectively injected in
the SCHEDULE section when the condition evaluates to true.
This is not so for PYACTION where there is one Python script in which both
conditions can be evaluated and changes applied.

In order to enable the PYACTION keyword:

1. OPM flow must be compiled with the cmake switches -DOPM ENABLE EMBEDDED PYTHON=ON and -DOPM ENABLE PYTHON=ON, the default is to build with these switches set to OFF.

2. The keyword PYACTION must be added to the SCHEDULE section:

.. code-block:: python

	<PYACTION\_NAME>  <SINGLE/UNLIMITED> /
	<pythonscript> /

3. You need to provide the Python script.


To interact with the simulator in the embedded Python code, you can access four variables from the simulator:

.. code-block:: python

	# Python module opm_embedded
	import opm_embedded
	# The current EclipseState
	ecl_state = opm_embedded.current_ecl_state
	# The current Schedule
	schedule = opm_embedded.current_schedule
	# The current SummaryState
	summary_state = opm_embedded.current_summary_state
	# The current report step
	report_step = opm_embedded.current_report_step 


- current_ecl_state: An instance of the EclipseState class - this is a representation of all static properties in the model, ranging from porosity to relperm tables. The content of the ecl state is immutable - you are not allowed to change the static properties at runtime.

- current_schedule: An instance of the Schedule class - this is a representation of all the content from the SCHEDULE section, notably all well and group information and the timestepping.

- current_report_step: This is an integer for the report step we are currently working on. Observe that the PYACTION is called for every simulator timestep, i.e. it will typically be called multiple times with the same value for the report step argument.

- current_summary_state: An instance of the SummaryState class, this is where the current summary results of the simulator are stored. The SummaryState class has methods to get hold of well, group and general variables.
