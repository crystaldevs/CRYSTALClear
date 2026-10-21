Introduction
============

.. raw:: html

   <p class="cc-lead">CRYSTALClear turns what a <a href="https://www.crystal.unito.it">CRYSTAL</a>
   calculation writes into numbers you can use and figures you can publish.
   Point it at an output and ask: the geometry, the energy, the band gap, the
   elastic tensor, the spectrum. Hand the answer to a plotting function and
   you have a matplotlib figure.</p>

How it fits together
--------------------

.. raw:: html

   <div class="cc-flow">
     <div class="cc-step">
       <span class="cc-step-label">CRYSTAL writes</span>
       <code>.out</code> outputs, and the <code>.DAT</code> and <code>.f25</code>
       files of a PROPERTIES run
     </div>
     <div class="cc-arrow" aria-hidden="true"></div>
     <div class="cc-step">
       <span class="cc-step-label">crystal_io reads</span>
       <code>Crystal_output</code>, <code>Properties_output</code>,
       <code>External_unit</code>
     </div>
     <div class="cc-arrow" aria-hidden="true"></div>
     <div class="cc-step">
       <span class="cc-step-label">plot draws</span>
       one function per kind of plot: <code>plot_electron_band</code>,
       <code>plot_cry_ela</code>, <code>plot_cry_spec</code>…
     </div>
     <div class="cc-arrow" aria-hidden="true"></div>
     <div class="cc-step">
       <span class="cc-step-label">you get</span>
       matplotlib figures and axes, to restyle, combine or save
     </div>
   </div>

:doc:`usage` goes through the three reading classes, and the
:doc:`example notebooks <notebooks/index>` through each kind of property.

What it covers
--------------

.. raw:: html

   <div class="cc-tiles">
     <div class="cc-tile">
       <h3>Electronic structure</h3>
       <p>Band structures and densities of states, alone or side by side,
       projected or total; charge density maps and their differences; X-ray
       diffraction patterns.</p>
       <p class="cc-tile-api"><a href="plot.html#CRYSTALClear.plot.plot_electron_banddos"><code>plot_electron_banddos</code></a> <a href="plot.html#CRYSTALClear.plot.plot_cry_contour"><code>plot_cry_contour</code></a> <a href="plot.html#CRYSTALClear.plot.plot_cry_xrd"><code>plot_cry_xrd</code></a></p>
     </div>
     <div class="cc-tile">
       <h3>Lattice dynamics</h3>
       <p>Phonon bands and densities of states; IR and Raman spectra, harmonic
       and anharmonic; VCI levels, cuts through the potential energy surface,
       anharmonic scans.</p>
       <p class="cc-tile-api"><a href="plot.html#CRYSTALClear.plot.plot_phonon_banddos"><code>plot_phonon_banddos</code></a> <a href="plot.html#CRYSTALClear.plot.plot_cry_spec"><code>plot_cry_spec</code></a> <a href="plot.html#CRYSTALClear.plot.plot_cry_vci"><code>plot_cry_vci</code></a> <a href="plot.html#CRYSTALClear.plot.plot_cry_pes_1D"><code>plot_cry_pes_1D</code></a></p>
     </div>
     <div class="cc-tile">
       <h3>Elastic properties</h3>
       <p>Young and shear moduli, linear compressibility and Poisson ratio in
       every direction, as surfaces or as sections through a plane.</p>
       <p class="cc-tile-api"><a href="plot.html#CRYSTALClear.plot.plot_cry_young"><code>plot_cry_young</code></a> <a href="plot.html#CRYSTALClear.plot.plot_cry_comp"><code>plot_cry_comp</code></a> <a href="plot.html#CRYSTALClear.plot.plot_cry_ela_2D"><code>plot_cry_ela_2D</code></a></p>
     </div>
     <div class="cc-tile">
       <h3>Equations of state</h3>
       <p>Energy against volume, Birch–Murnaghan fits, and the enthalpy
       difference that marks a transition between phases.</p>
       <p class="cc-tile-api"><a href="crystal_io.html#CRYSTALClear.crystal_io.Crystal_output.get_EOS"><code>get_EOS</code></a> <a href="plot.html#CRYSTALClear.plot.plot_cry_EOS"><code>plot_cry_EOS</code></a></p>
     </div>
     <div class="cc-tile">
       <h3>Topology of the density</h3>
       <p>TOPOND results: maps of the density and its Laplacian, profiles
       along a line, magnetisation and spin-current vector fields.</p>
       <p class="cc-tile-api"><a href="crystal_io.html#CRYSTALClear.crystal_io.Properties_output.read_topond_tlap"><code>read_topond_tlap</code></a> <a href="plot.html#CRYSTALClear.plot.plot_cry_lapl_profile"><code>plot_cry_lapl_profile</code></a> <a href="plot.html#CRYSTALClear.plot.plot_vecfield2D_m"><code>plot_vecfield2D_m</code></a></p>
     </div>
     <div class="cc-tile">
       <h3>Transport</h3>
       <p>Seebeck coefficient, electrical conductivity, power factor and ZT,
       against chemical potential or carrier concentration.</p>
       <p class="cc-tile-api"><a href="plot.html#CRYSTALClear.plot.plot_cry_seebeck_potential"><code>plot_cry_seebeck_potential</code></a> <a href="plot.html#CRYSTALClear.plot.plot_cry_powerfactor_carrier"><code>plot_cry_powerfactor_carrier</code></a> <a href="plot.html#CRYSTALClear.plot.plot_cry_zt"><code>plot_cry_zt</code></a></p>
     </div>
     <div class="cc-tile cc-tile-wide">
       <h3>Everything else in an output</h3>
       <p>Geometry and symmetry, lattice and reciprocal lattice, energies and
       convergence, band gap and Fermi energy, Mulliken charges, forces, the
       dielectric tensor: each a single method of <code>Crystal_output</code> away.</p>
       <p class="cc-tile-api"><a href="crystal_io.html#CRYSTALClear.crystal_io.Crystal_output.get_geometry"><code>get_geometry</code></a> <a href="crystal_io.html#CRYSTALClear.crystal_io.Crystal_output.get_band_gap"><code>get_band_gap</code></a> <a href="crystal_io.html#CRYSTALClear.crystal_io.Crystal_output.get_mulliken_charges"><code>get_mulliken_charges</code></a> <a href="crystal_io.html#CRYSTALClear.crystal_io.Crystal_output.get_dielectric_tensor"><code>get_dielectric_tensor</code></a></p>
     </div>
   </div>

Where it comes from
-------------------

CRYSTALClear started as a fork of `CRYSTALpytools <cpt_>`_, and builds on the
work done there: its functionality has been expanded, its interface improved
and its documentation made more detailed. It is also what
`CRYSTALLine <https://crystaldevs.github.io/CRYSTALLine/>`_, the desktop
program for CRYSTAL structures, reads files and draws plots with.

It grows with the people who use it. Bug reports, ideas and pull requests are
welcome on `GitHub <https://github.com/crystaldevs/CRYSTALClear>`_.

.. _cpt: https://github.com/crystal-code-tools/CRYSTALpytools
