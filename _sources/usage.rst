Usage
=====

.. rst-class:: cc-lead

Every use of CRYSTALClear follows the same path: read a file, extract what you
need from it, and plot it.

.. raw:: html

   <div class="cc-flow">
     <div class="cc-step"><span class="cc-step-label">1 · Read</span>
       a file into an object of <code>crystal_io</code></div>
     <div class="cc-arrow" aria-hidden="true"></div>
     <div class="cc-step"><span class="cc-step-label">2 · Extract</span>
       a property with one of its methods</div>
     <div class="cc-arrow" aria-hidden="true"></div>
     <div class="cc-step"><span class="cc-step-label">3 · Plot</span>
       it with a function of <code>plot</code></div>
     <div class="cc-arrow" aria-hidden="true"></div>
     <div class="cc-step"><span class="cc-step-label">Then</span>
       restyle and save it as any matplotlib figure</div>
   </div>

Files to try it on
------------------

If you have no CRYSTAL calculation to hand, take these:
:download:`crystalclear-examples.zip <examples/crystalclear-examples.zip>`
(76 kB), real calculations from the example notebooks:

- ``spectra/``: a frequency calculation of thiourea, with IR intensities;
- ``electronic/``: the band structure and density of states of MgO, written
  by PROPERTIES;
- ``elastic/``: an elastic constants calculation of coesite.

Unpack it and start Python in ``crystalclear-examples/``: every example on this
page then runs as it is written, and the last one ends with the figure of MgO's
bands beside its density of states. For coesite, ``plot_cry_ela`` draws the
directional Young's modulus:

.. code-block:: python

    from CRYSTALClear.crystal_io import Crystal_output
    import CRYSTALClear.plot as CCplt

    co = Crystal_output("elastic/coesite_ela.out")
    co.get_elatensor()
    CCplt.plot_cry_ela(co, "young", ndeg=50)

Read a file
-----------

The ``crystal_io`` module has one class for each kind of file CRYSTAL reads or
writes:

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: cc-classes

   * - Class
     - Reads
   * - ``Crystal_input``
     - a CRYSTAL input deck (``.d12``), to read, edit and write
   * - ``Crystal_output``
     - a CRYSTAL output (``.out``)
   * - ``Properties_input``
     - a PROPERTIES input deck (``.d3``)
   * - ``Properties_output``
     - the data files a PROPERTIES run writes (``.DAT``, ``.f25``)
   * - ``Crystal_gui``
     - a geometry file (``.gui``, ``fort.34``), to read or to write
   * - ``Crystal_density``
     - a density matrix (``.f98``)
   * - ``External_unit``
     - other files CRYSTAL writes: ``IRSPEC.DAT``, ``RAMSPEC.DAT``, phonon bands
       and densities of states

An object is made by passing it the path of the file:

.. code-block:: python

    from CRYSTALClear.crystal_io import Crystal_output

    co = Crystal_output("spectra/thiourea_pbed3_Ahlrichs-pVTZ_freqcalc.out")

Extract what you need
---------------------

Its methods read one property each, and keep it on the object:

.. code-block:: python

    co.get_IR()          # the IR intensities of a frequency calculation
    print(co.IR_HO_0K)   # ...now an attribute of the object

Plot it
-------

The functions of the ``plot`` module draw it, and return the matplotlib
``fig`` and ``ax``, to be adjusted with the usual matplotlib commands:

.. code-block:: python

    from CRYSTALClear.plot import plot_cry_spec

    fig, ax = plot_cry_spec(co.IR_HO_0K)

    ax.set_title("Infrared spectrum")
    ax.set_xlim(1200, 2500)                 # cm⁻¹
    fig.savefig("IR_1200-2500.png", dpi=300)

A complete example
------------------

The band structure of MgO beside its density of states, from the files
PROPERTIES writes to a saved figure:

.. code-block:: python

    from CRYSTALClear.crystal_io import Properties_output
    from CRYSTALClear.plot import plot_electron_banddos

    # 1 · Read and 2 · Extract: one method each, for the two files
    bands = Properties_output().read_electron_band("electronic/mgo_band.BAND")
    doss = Properties_output().read_electron_dos("electronic/mgo_doss_totalao.DOSS")

    # 3 · Plot
    fig, ax = plot_electron_banddos(bands, doss,
                                    k_labels=["Gamma", "X", "W", "L", "Gamma"],
                                    energy_range=[-10, 15])

    # Then save
    fig.savefig("mgo_bands_dos.png", dpi=300)
