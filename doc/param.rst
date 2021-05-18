**************
Parameter File
**************
  
``trains`` reads two parameter files to configure its running.
A typical primary parameter file looks like::
  
  % parameter file
  % lines beginning with '%' are regarded as comments and are neglected
  %========================================================================

  FileDir         .                           % directory for where all the output files are stored
  FlagRec         1                           % 0, generate mock data; 1, reconstruct model from data
  ParamPT         param/param_pt              % parameter file for model

A typical secondary parameter file looks like::

  % parameter file
  % lines beginning with '%' are regarded as comments and are neglected
  % 
  % this is a full list. 
  % in a run, only some parameters are required to specify
  %========================================================================
  % generic
  %
  FlagMethod                  0                           % 0, search GW source by FreePhase method
                                                          % 1, search GW source by MxPhase method
                                                          % 2, search GW source by AvPhase method
  FlagEvolve                  0                           % Is binary orbit evolve or not?
  %========================================================================
  % data file
  %
  PulsarCatalog               data/pulsar_catalog.txt     % file for pulsar catalog
  PTRFile                     data/sim_ptr.txt            % file for pulsar timing residuals
  PTRConstructFileOut         data/pptr.txt               % output filename for ptr reconstruction
  %=========================================================================
  % pulsar timing array generic configuration
  %

  SourceNumber                1                           % Number of GW source in simulation and free phase search
  TimingNumber                130                         % Number of points in timing residuals of a pulsar in simulation
  %========================================================================
  % set fixed GW source parameters and their fixed values
  % do not put space in the strings
  % 1: fixed; 0: not fixed; values are separated by ":"
  %

  SourceParFix                00000000
  SourceParFixVal             0.0：0.0

.. note::
  
  In the subdirectory ``example/``, some examples of parameter file are provided. 
  Users can choose appropriate parameter files with their purposes.