.. CRYSTALClear documentation master file

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: User Guide
   :glob:

   introduction
   installation
   usage
   notebooks/index
   cite

.. toctree::
   :hidden:
   :maxdepth: 3
   :caption: API Documentation
   :glob:

   convert
   crystal_io
   execute
   flatmaps
   plot
   units

CRYSTALClear
============

.. raw:: html

   <div class="cc-hero">
     <img class="cc-logo cc-logo-light" src="_static/CRYSTAL_logo.png" alt="CRYSTALClear">
     <img class="cc-logo cc-logo-dark" src="_static/CRYSTAL_logo-dark.png" alt="CRYSTALClear">
     <p class="cc-tagline">A Python framework for the <strong>CRYSTAL</strong>
       quantum-chemistry code: it reads what a calculation writes and plots
       the properties it returns.</p>
     <pre class="cc-install"><code>pip install CRYSTALClear</code></pre>
     <div class="cc-buttons">
       <a class="cc-btn cc-btn-primary" href="installation.html">Get started</a>
       <a class="cc-btn" href="notebooks/index.html">Example notebooks</a>
       <a class="cc-btn" href="crystal_io.html">API reference</a>
     </div>
   </div>

A quick look
------------

Read an output, ask it for what you need, and plot it. The IR spectrum of a
frequency calculation, as the vibrational spectra notebook draws it:

.. code-block:: python

   from CRYSTALClear.crystal_io import Crystal_output
   import CRYSTALClear.plot as CCplt

   CO = Crystal_output("spectra/thiourea_pbed3_Ahlrichs-pVTZ_freqcalc.out")
   CO.get_IR()
   CCplt.plot_cry_spec(CO.IR_HO_0K)

Every plotting function returns matplotlib objects, to be restyled, combined or
saved like any other figure.

Examples
--------

Each notebook works through one kind of property on a real calculation.

.. raw:: html

   <div class="cc-gallery">
     <a class="cc-card" href="notebooks/Electronic_Structure/Electronic_Structure_CRYSTALClear.html">
       <img src="_static/gallery/electronic-structure.png" alt="Band structure and density of states of MgO">
       <span class="cc-card-title">Electronic structure</span>
       <span class="cc-card-text">Band structures and densities of states, side by side or compared across functionals.</span>
     </a>
     <a class="cc-card" href="notebooks/Equation_Of_State/EOS.html">
       <img src="_static/gallery/equation-of-state.png" alt="Equations of state of two phases">
       <span class="cc-card-title">Equation of state</span>
       <span class="cc-card-text">Energy against volume for two phases, and the pressure at which one gives way to the other.</span>
     </a>
     <a class="cc-card" href="notebooks/elastic_properties/elastic_properties_tutorial.html">
       <img src="_static/gallery/elastic-properties.png" alt="Directional elastic property of coesite">
       <span class="cc-card-title">Elastic properties</span>
       <span class="cc-card-text">Young and shear moduli, linear compressibility and Poisson ratio of coesite, in every direction.</span>
     </a>
     <a class="cc-card" href="notebooks/electron_transport_properties/tutorial_electron_transport_properties.html">
       <img src="_static/gallery/transport.png" alt="Power factor at several temperatures">
       <span class="cc-card-title">Transport properties</span>
       <span class="cc-card-text">Seebeck coefficient, conductivity, power factor and ZT across temperatures.</span>
     </a>
     <a class="cc-card" href="notebooks/vibrational_spectra/vibrational_spectra_tutorial.html">
       <img src="_static/gallery/vibrational-spectra.png" alt="Spectra at three temperatures">
       <span class="cc-card-title">Vibrational spectra</span>
       <span class="cc-card-text">IR and Raman spectra of a frequency calculation, with Gaussian or Lorentzian broadening.</span>
     </a>
   </div>

More about the project, where it comes from and what it covers, is in the
:doc:`introduction <introduction>`. If CRYSTALClear helps your work, here is
:doc:`how to cite it <cite>`.
