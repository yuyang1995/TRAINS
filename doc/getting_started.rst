.. _getting_started:

***************
Getting Started
***************

Third-party Package Dependence
===============================

  * **MPICH** --- an MPI implementation library, available at http://www-unix.mcs.anl.gov/mpi/mpich

  * **GSL** --- the GNU Scientific Library, downloaded at http://www.gnu.org/software/gsl

  * **CDNest** --- Diffusive nested sampling, downloaded at https://github.com/LiyrAstroph/CDNest

.. note::
  In Linux systems, there are package managers that can install the above libraries conveniently (except for CDNest).
  For example, DNF package managers in Fedora distribution.  
  If so, use them. In such a case, the libraries usually are installed in standard environment path. Otherwise, 
  if any of the above libraries is not installed in standard locations on your system, the :ref:`Makefile` provided 
  with the code may need slight adjustments. 

Compiling
=============================

Edit the configurations in :ref:`Makefile` to be consistent with your system if necessary. Then compile the package with the command

.. code:: bash

   make

This creates an executable file ``trains``.

Running
=============================

To run the package in a parallel computer/cluster, use the following command: 

.. code:: bash

   mpiexec -n np ./trains param/param

where ``np`` is the number of cores and ``param`` is the parameter file (see :ref:`Parameter File`), which specifies configurations for ``trains``.

An exemplary pulsar timing residuals dataset is provided in the subdirectory ``data/``, including:

.. code-block:: bash

   pulsar_catalog.txt            # location and timing uncertainties of pulsars
   sim_ptr.txt                   # timing residuals of pulsars
  

One can try to run the above command to test ``trains`` with the provided dataset.
See :ref:`Data Format` for the information of data format.

Command-line Options
======================

``trains`` also adimits several command-line options:

.. code-block:: bash

    -h or --help
        print help information.
    -p or --post_proc
        only do posterior processing.
    -r or --restart
        restart from the backup. 
    -t or --temperature or --temp
        specify tempering temperature in posterior processing, e.g., ``-t 2``, ``--temperature 2``, or ``--temp 2``.
    -s or --seed 
        set a seed for the random number generator, e.g., ``-s 100`` or ``--seed 100``
    -c or --recalc_info
        only do posterior processing, but recalculate the posterior sample information.
    -e or --exam_prior
        examine the priors.
    -v or --version
        print version.
    -n or --para_name
        print parameter name.

For example, if one wants to redo posterior processing with a different temperature, say 10 (the default is 1), one may use the command

.. code:: bash

   ./trains src/param -pt10



MCMC Sampling
=============

The output Markov chain is stored in ``data/posterior_sample_pt.txt``.

The parameter names and prior ranges are stored in ``data/para_names_pt.txt``. 
The last column of those files indicates the prior type of the parameter with ``1`` means Gaussian and ``2`` means uniform.

One need to tune the corresponding option files ``OPTIONSPT``, which specify configurations for nested sampling.
