****************
Data
****************

Data Input
==========

The input data in ``trains`` depend on the option ``FlagRec`` specified in the parameter file (see :ref:`Parameter File`). 

* ``FlagRec = 0``

  create simulation data for given pulsar timing array.
  ``PulsarCatalog`` need to be specified.


* ``FlagRec= 1``
  
  fit the signal model to the given pulsar timing residuals.
  ``PulsarCatalog`` and `PTRFile` need to be specified.


Data Format
===========

* pulsar catalog should contain five columns, 
  which represent right ascension, declination, distance, distance uncertainty and timing uncertainty respectively.
  Units of right ascension and declination should be radian.
  Units of distance and distance uncertainty should be lightyear.
  Units of timing uncertainty should be second.

* data of pulsar timing residuals contain two columns,
  which represent observation time and timing residual respectively.
  Units of observation time and timing residual should be year and second respectively.
  
  .. note::

    In current version, the number of data points in timing residuals :math:`N_{\rm t}` of each pulsar should be the same.
    The number of rows in the data file should be the product of :math:`N_{\rm t}` and the total number of pulsars :math:`N_{\rm p}`.
    Thus, row :math:`1` to row :math:`N_{\rm t}` are timing residuals of the first pulsar, 
    row :math:`N_{\rm t} + 1` to row :math:`2N_{\rm t}` are timing residuals of the second pulsar, and so on.
    The order of pulsars should be the same as that of the pulsar catalog.