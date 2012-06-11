Linux (Debian/Ubuntu)
=====================

This section details how to install the Python version and the command line
version of Psignifit from a source download.

.. note::

    Currently (Jun 2012) we provide experimental support for packages via
    `Neurodebian <http://neuro.debian.net/>`_ for `all current releases of
    Debian *and* Ubuntu
    <http://neuro.debian.net/pkgs/python-pypsignifit.html>`_. If you have
    configured the NeuroDebian repositories correctly, the following command
    will install Psignifit3 for Python::

        $ sudo apt-get install python-pypsignifit

    This option is the easiest route to getting python-psignifit installed.

.. note::

    If you wish (or need) to install from source, first look at:

    * `Obtaining the Sources`_
    * `Installing Dependencies`_

    And then continue to install the desired version:

    * `Installing Psignifit for Python (python-psignifit)`_

    Or:

    * `Installing the Command Line Interface (cli)`_

.. note::

    If you need to build the documentation on your local machine, you may want to
    look at:

    * `Building the documentation`_

Obtaining the Sources
---------------------

You will want to download the most recent version of psignifit from:
`<http://sourceforge.net/projects/psignifit/files/>`_. You will want the
``zip`` file, for example ``psignifit_3.0_beta.20120207.1.zip``.

Extract the file, and enter the directory by typing::

    unzip psignifit_3.0_beta.<date of the snapshot>.1.zip
    cd psignifit_3.0_beta.<date of the snapshot>.1

replacing ``<date of the snapshot>`` by the date string in the file
name.

Installing Dependencies
-----------------------

If you are using `Debian <http://www.debian.org/>`_ or `Ubuntu
<http://www.ubuntu.com/>`_ the following packages need to be installed.

* ``make``
* ``gcc``

If  you additionally wish to install the Python version, you must install:

* ``python``
* ``python-dev``
* ``python-numpy (provides python-numpy-dev)``
* ``python-scipy``
* ``python-matplotlib``
* ``python-nose``
* ``swig``

In order to check whether or not you have the packages already installed, type::

    aptitude search make gcc python python-dev\
        python-numpy python-scipy python-matplotlib python-nose swig

Packages that are installed on your machine are listed with a leading ``<i>``

In order to install missing packages, type::

    sudo aptitude install make gcc python python-dev\
        python-numpy python-scipy python-matplotlib python-nose

If you are using a different Linux distribution, please refer to the
documentation of that system for more information on how to install the required
software.

Installing Psignifit for Python (python-psignifit)
--------------------------------------------------

System-wide installation
~~~~~~~~~~~~~~~~~~~~~~~~

Installing system wide is simple, type::

    make install

as ``root`` and everything will be installed to the right place.

Install into a Custom Location
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you do not have root privileges, you can use the ``DESTDIR`` variable to
install into an alternative location, for example the directory
``psignifit-install`` in your home-directory::

    mkdir $HOME/psignifit-install
    make DESTDIR=$HOME/psignifit-install install

The first command first creates the directory ``$HOME/psignifit-install``
where ``$HOME`` is automatically replaced by the name of your own
home-directory.  The second command will install Psignifit into that directory.

To use psignifit from, you will also have to set the ``$PYTHONPATH`` and
``$LD_LIBRARY_PATH`` environment variables. Either invoke the Python
interpreter from the command line by typing::

   PYTHONPATH=psignifit-install/lib/python/ LD_LIBRARY_PATH=psignifit-install/usr/lib/ ipython

Or, you set the variables in your ``.bashrc`` (or equivalent) file
by adding the line::

   export PYTHONPATH=psignifit-install/lib/python/
   export LD_LIBRARY_PATH=psignifit-install/usr/lib/

Testing your Installation
~~~~~~~~~~~~~~~~~~~~~~~~~

Open a Python interpreter and type::

    import pypsignifit
    pypsignifit.version

If you can see a version string, such as ``snap-2011-10-28-59-gda5adc6`` you
installation was probably successful.

Installing the Command Line Interface (cli)
-------------------------------------------

System-wide installation
~~~~~~~~~~~~~~~~~~~~~~~~
Installing system wide is simple, type::

    make cli-install

as ``root`` and everything will be installed to the right place.

Install into a Custom Location
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you do not have root privileges, you can use the ``DESTDIR`` variable to
install into an alternative location, for example the directory
``psignifit-install`` in your home-directory::

    mkdir $HOME/psignifit-install
    make DESTDIR=$HOME/psignifit-install cli-install

The first command first creates the directory ``$HOME/psignifit-install``
where ``$HOME`` is automatically replaced by the name of your own
home-directory.  The second command will install Psignifit into that directory.

Additionally you must add the installation directory to you ``$PATH``, by
typing::

    export PATH=$PATH:$HOME/psignifit-install/usr/bin

Or, for a more permanent solution and the above line to ``.bashrc`` (or
equivalent) file.

Testing your Installation
~~~~~~~~~~~~~~~~~~~~~~~~~

Now, you should be able to call::

    psignifit-mcmc -h
    psignifit-diagnostics -h
    psignifit-bootstrap -h
    psignifit-mapestimate -h

And see some usage messages after each call.

Building the documentation
--------------------------

We recommend you to use the documentation we supply on the website. But in case
you would like to build it locally, you must first install the following
packages:

* ``python-sphinx``
* ``doxygen``

And then do::

    make doc

The documentation can be found in ``doc-html/``. Point your browser at
``doc-html/index.html`` to begin.
