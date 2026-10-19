Installation
============

CRYSTALClear is on `PyPI <https://pypi.org/project/CRYSTALClear/>`_. It needs
**Python 3**.

It should be installed in an environment of its own rather than in the system
Python. Either a conda environment or a virtual environment may be used; conda
also provides the Python version, whereas a virtual environment uses the
``python`` already installed.

.. tab-set::

   .. tab-item:: conda

      .. code-block:: sh

         conda create --name crystalclear python=3.9
         conda activate crystalclear
         pip install CRYSTALClear

   .. tab-item:: venv

      .. code-block:: sh

         python -m venv ~/.venvs/crystalclear
         source ~/.venvs/crystalclear/bin/activate    # Windows: ~\.venvs\crystalclear\Scripts\activate
         pip install CRYSTALClear

Then import it:

.. code-block:: python

   import CRYSTALClear
   print(CRYSTALClear.__version__)

:doc:`usage` shows how a calculation is read and plotted, and the
:doc:`example notebooks <notebooks/index>` go through each kind of property.

Dependencies
------------

These are installed automatically if not already present:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Package
     - Why
   * - pymatgen
     - structures, lattices and symmetry; it brings numpy, scipy and matplotlib
       with it, for the arithmetic and the plots
   * - mendeleev
     - element symbols from atomic numbers
   * - ase
     - ASE ``Atoms`` objects, built from TOPOND data

Installation from source
------------------------

.. code-block:: sh

   git clone https://github.com/crystaldevs/CRYSTALClear.git
   cd CRYSTALClear
   pip install -e .

``-e`` installs it in editable mode: changes to the source take effect without
installing again.
