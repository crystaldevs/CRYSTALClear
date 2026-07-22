Usage
=====

Once you have installed **CRYSTALClear**, you can easily parse and visualize 
results from **CRYSTAL** calculations.

The general workflow is:

1. **Create a parser object** using the classes in the ``crystal_io`` module.
2. **Extract data** from the CRYSTAL output files by calling the appropriate methods.
3. **Plot and customize** physical properties using the functions in the ``plot`` module.

----

1. Creating a Parser Object
---------------------------

The ``crystal_io`` module provides seven main classes:

- ``Crystal_input``
- ``Crystal_output``
- ``Properties_input``
- ``Properties_output``
- ``Crystal_gui``
- ``Crystal_density``
- ``External_unit``

Each class is designed to handle specific CRYSTAL input or output files.

**Example:**

.. code-block:: python

    from CRYSTALClear.crystal_io import Crystal_output

    # Create an object by passing the path to a CRYSTAL output file
    co = Crystal_output("path/to/CRYSTAL_output_file.out")


----

2. Extracting Data
------------------

Once the object is created, you can call its **data-extraction methods** to 
populate attributes with relevant information from the CRYSTAL output.

**Example:**

.. code-block:: python

    # Extract IR spectrum
    co.get_IR()

    # Now the object contains attributes with the parsed data
    print(co.IR_HO_0K)

----

3. Plotting Results
-------------------

The ``plot`` module provides advanced, publication-quality plotting utilities 
for physical properties. These functions return standard **Matplotlib** objects 
(``fig``, ``ax``), so you can further customize them with familiar Matplotlib 
commands.

**Example:**

.. code-block:: python

    from CRYSTALClear.plot import plot_cry_spec
    import matplotlib.pyplot as plt

    # Generate the plot
    fig, ax = plot_cry_spec(co.IR_HO_0K)

    # Customize the plot using Matplotlib
    ax.set_title("Infrared Spectrum")
    ax.set_ylim(1200, 2500)

    # Save the figure
    fig.savefig("IR_1200-2500.png", dpi=300)

----

4. Summary Workflow
-------------------

1. **Load** CRYSTAL output with the appropriate ``crystal_io`` class.
2. **Extract** desired properties via dedicated methods.
3. **Visualize** results using ``plot`` module functions.
4. **Customize** plots with Matplotlib as needed.

----

Example Session
---------------

.. code-block:: python

    from CRYSTALClear.crystal_io import External_unit
    from CRYSTALClear import plot as CCplt

    # Step 1: Read output file
    co = External_unit()

    # Step 2: Extract band structure
    co.read_phonon_band("my_bands.f25")

    # Step 3: Plot
    fig, ax = plot_phonon_band(co.band_structure)

    # Step 4: Customize
    ax.set_title("Example Phonon Band Structure")
    fig.tight_layout()
    fig.savefig("phonon_band_structure.png")

With these simple steps, **CRYSTALClear** lets you go from raw CRYSTAL output files 
to beautiful, customizable plots of your results.

