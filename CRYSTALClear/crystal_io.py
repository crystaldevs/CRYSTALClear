#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Objects of input / output files of CRYSTAL. Methods to edit or substract data
from corresponding files are provided.
"""
from CRYSTALClear import units
from CRYSTALClear.base.crysd12 import Crystal_inputBASE


class Crystal_input(Crystal_inputBASE):
    """
    Crystal input object inherited from the :ref:`Crystal_inputBASE <ref-base-crysd12>`
    object. For the basic set-ups of keywords, please refer to manuals there.
    """

    def __init__(self):
        super(Crystal_input, self).__init__()

    def geom_from_cif(self, file, zconv=None, keyword='EXTERNAL',
                      pbc=[True, True, True], gui_name='fort.34', **kwargs):
        """
        Read geometry from cif file and put infomation to geom block, either as
        "EXTERNAL" or "CRYSTAL". CIF files with a single geometry only.

        CIF with Symmetry is required.

        .. note::

            CIF files should have the ".cif" extension.

        Args:
            file (str): CIF file.
            keyword (str): "EXTERNAL" or "CRYSTAL".
            zconv (list[list[int, int]]): 1st element: The **index** of atom;
                2nd element: The new conventional atomic number. Atoms of the
                irreducible unit is required.
            pbc (list[bool]): Valid only if ``keyword = EXTERNAL``. Force to remove
                periodic boundary conditions along x, y, z axis.
            gui_name (str): Valid only if ``keyword = EXTERNAL``. Gui file's name.
            **kwargs: Passed to Pymatgen `SpacegroupAnalyzer <https://pymatgen.org/pymatgen.symmetry.html#pymatgen.symmetry.analyzer.SpacegroupAnalyzer>`_ object.
        """
        import re

        from pymatgen.core.lattice import Lattice
        from pymatgen.core.structure import IStructure, Structure

        struc_3d = IStructure.from_file(file)
        if keyword == 'CRYSTAL':
            struc = struc_3d
        else:
            latt_mx = struc_3d.lattice.matrix
            latt = Lattice(latt_mx, pbc=pbc)
            struc = Structure(latt,
                              species=list(struc_3d.atomic_numbers),
                              coords=struc_3d.cart_coords,
                              coords_are_cartesian=True)
        self.geom_from_pmg(struc, zconv, keyword, gui_name, **kwargs)

        return self

    def geom_from_pmg(self, struc, zconv=None, keyword='EXTERNAL',
                      gui_name='fort.34', **kwargs):
        """
        Read geometry defined by PyMatGen structure object and put infomation
        into geom block, either as "EXTERNAL" or "CRYSTAL".

        See Crystal_input.geom_from_cif for definition of arguments.

        .. note::

            Coordinates of corresponding atoms may not consistent with the
            original CIF file if "CRYSTAL" is used, in which case
            coordinates of another symmetry equivalent atom is used.
        """
        import re

        from CRYSTALClear.convert import cry_pmg2gui
        from CRYSTALClear.geometry import refine_geometry

        if re.match(r'^EXTERNAL$', keyword, re.IGNORECASE):
            super(Crystal_input, self).geom.external()
            gui = cry_pmg2gui(struc, gui_file=gui_name, symmetry=True,
                              zconv=zconv, **kwargs)
        elif re.match(r'^CRYSTAL$', keyword, re.IGNORECASE):
            sg, _, latt, natom, atom = refine_geometry(struc, **kwargs)
            if zconv != None:
                z_atom_index = [i[0] for i in zconv]
                for i in range(natom):
                    try:
                        atom_to_sub = z_atom_index.index(i)
                        z_input = zconv[atom_to_sub][1]
                    except ValueError:
                        z_input = atom[i][0]
                    atom[i][0] = z_input
            super(Crystal_input, self).geom.crystal(
                IGR=sg, latt=latt, atom=atom)
        else:
            raise ValueError("Input keyword format error: {}".format(keyword))

        return self

    def bs_from_bse(self, name, element, zconv=None):
        """
        Download basis set definitions from BSE. A wrapper of BasisSetBASE.from_bse.

        Args:
            name (str): Basis set's name.
            element (list[str] | list[int]): List of elements.
            zconv (list[int]): If not none, use the conventional atomic number.
                Its length must be the same as element. Its sequence must be
                consistent with basis set's.
        """
        from CRYSTALClear.base.basisset import BasisSetBASE

        return BasisSetBASE.from_bse(name=name, element=element, zconv=zconv)

    def bs_from_string(bs_str, fmt='crystal'):
        """
        Define basis set from a string. A wrapper of BasisSetBASE.from_string.

        Args:
            bs_str (str)
            fmt (str): Format string. Consistent with BSE python API.
        """
        from CRYSTALClear.base.basisset import BasisSetBASE

        return BasisSetBASE.from_string(bs_str=bs_str, fmt=fmt)

    def bs_from_file(file, fmt='crystal'):
        """
        Define a basis set from a file. A wrapper of BasisSetBASE.from_file.
        """
        from CRYSTALClear.base.basisset import BasisSetBASE

        return BasisSetBASE.from_file(file=file, fmt=fmt)


class Crystal_output:
    """This class reads a CRYSTAL output and generates an object."""

    def __init__(self, output=None):
        """
        Initialize the Crystal_output.

        Args:
            output (str): Filename
        """
        if output != None:
            self._read_output(output)

    def _read_output(self, output_name):
        """
        Reads a CRYSTAL output file.

        Args:
            output_name (str): Name of the output file.
        Returns:
            self (Crystal_output)
        """
        import re
        import sys

        self.name = output_name

        # Check if the file exists
        try:
            if output_name[-3:] != 'out' and output_name[-4:] != 'outp':
                output_name = output_name+'.out'
            file = open(output_name, 'r', errors='ignore')
            self.data = file.readlines()
            file.close()
        except:
            raise FileNotFoundError(
                'EXITING: a .out file needs to be specified')

        # Check the calculation terminated correctly
        self.terminated = False

        for i, line in enumerate(self.data[::-1]):
            if re.match(r'^ EEEEEEEEEE TERMINATION', line):
                self.terminated = True
                # This is the end of output
                self.eoo = len(self.data)-1-i
                break

        if self.terminated == False:
            self.eoo = len(self.data)

        return self

    def read_cry_output(self, output_name):
        """
        Deprecated
        """
        import warnings

        warnings.warn('Deprecated. Define output during initialization.')
        return self._read_output(output_name)

    def get_dielectric_tensor(self):
        """Extracts the dielectric tensor from the output.

        Returns:
            self.dielectric_tensor (list[float]): Dielectric tensor values.
        """
        import re

        for i, line in enumerate(self.data):
            if re.match(r'^ TENSOR IN PRINCIPAL AXES SYSTEM', line):
                # This is the end of output
                self.dielectric_tensor = [
                    float(x) for x in self.data[i+1].split()[1::2]]
                return self.dielectric_tensor
        return None

    def get_eigenvectors(self):
        """Extracts eigenvectors from the output.

        Returns:
            self.num_ao (int): Number of atomic orbitals.
            self.num_k (int): Number of k points.
        """
        import re

        for i, line in enumerate(self.data):
            if re.match(r'\s NUMBER OF AO', line) != None:
                self.num_ao = int(line.split()[3])

            if re.match(r'\s SHRINK. FACT.(MONKH.)', line) != None:
                self.num_k = int(line.split()[13])

            if re.match(r'\s SHRINK. FACT.(MONKH.)', line) != None:
                self.num_k = int(line.split()[13])

    #### Geometry ####

    def get_dimensionality(self):
        """
        Gets the dimensionality of the system.

        Returns:
            self.dimensionality (int): Dimensionality of the system.
        """
        import re

        self.dimensionality = None
        for line in self.data:
            if re.match(r'\sGEOMETRY FOR WAVE FUNCTION - DIMENSIONALITY OF THE SYSTEM', line):
                self.dimensionality = int(line.split()[9])
                break

        if self.dimensionality == None:
            raise Exception('Invalid file. Dimension information not found.')

        return self.dimensionality

    def get_symmops(self):
        """
        Return the symmetry operators

        Returns:
            self.symmops (numpy.ndarray): Symmetry operators.
        """
        import re

        import numpy as np

        self.n_symmops = 0
        self.symmops = np.array([])

        symmops = []
        for i, line in enumerate(self.data):
            if re.match(r'^ \*\*\*\*   \d+ SYMMOPS - TRANSLATORS IN FRACTIONAL UNITS', line):
                self.n_symmops = int(line.split()[1])
                for j in range(0, self.n_symmops):
                    symmops.append(self.data[i+3+j].split()[2:])
                break

        if self.n_symmops > 0:
            self.symmops = np.reshape(np.array(symmops, dtype=float), [
                                      self.n_symmops, 4, 3])

        return self.symmops

    def get_geometry(self, initial=True, write_gui=False,
                     gui_name=None, symmetry='pymatgen', **kwargs):
        """
        Get the geometry.

        Args:
            initial (bool): Read the initial or last gemetry. Useful in
                case of geometry optimization.
            write_gui (bool): If True, write .gui file
            gui_name (str): Valid only if ``write_gui = True``. Gui file
                is named as ``gui_name``. If None, use ``basename.gui``. The
                basename is the same as output file.
            symmetry (str): Valid only if ``write_gui = True``. "pymatgen"
                to use symmetry info from a pymatgen SpacegroupAnalyzer;
                "initial" to use symmetry information on output file. If
                "None", no symmetry. Otherwise it is taken from the existing
                gui file.
            **kwargs: Valid only if ``write_gui = True`` and
                ``symmetry = "pymatgen"``.  Passed to Pymatgen
                SpacegroupAnalyzer object.

        Returns:
            self.geometry (Structure | Molecule): A pymatgen Structure or
                molecule object.
        """
        import os
        import re
        import warnings

        import numpy as np

        from CRYSTALClear.base.crysout import GeomBASE
        from CRYSTALClear.convert import cry_pmg2gui
        from CRYSTALClear.crystal_io import Crystal_gui
        from CRYSTALClear.geometry import rotate_lattice

        # Get geometry
        bg_line = -1
        if initial == True:
            for nline, line in enumerate(self.data[:self.eoo]):
                # Use atom coords to read molecule geometries. Go 4 lines up for periodic systems
                if re.match(r'^\s*ATOMS IN THE ASYMMETRIC UNIT', line):
                    bg_line = nline - 4
                    break
        else:
            for nline, line in enumerate(self.data[self.eoo::-1]):
                # Use atom coords to read molecule geometries. Go 4 lines up for periodic systems
                if re.match(r'^\s*ATOMS IN THE ASYMMETRIC UNIT', line):
                    bg_line = len(self.data[:self.eoo]) - nline - 4
                    break

        if bg_line < 0:
            raise Exception('Geometry information not found.')

        output = GeomBASE.read_geom(self.data, bg_line)
        struc = output[1]

        # Get the last lattice matrix: structure obtained by GeomBASE might be rotated.
        ndimen = self.get_dimensionality()
        lattice_line = -1
        if ndimen != 0:
            if initial == True:
                for nline, line in enumerate(self.data[:self.eoo]):
                    if re.match(r'^ DIRECT LATTICE VECTORS CARTESIAN', line):
                        lattice_line = nline
                        break
            else:
                for nline, line in enumerate(self.data[self.eoo::-1]):
                    if re.match(r'^ DIRECT LATTICE VECTORS CARTESIAN', line):
                        lattice_line = len(self.data[:self.eoo]) - nline
                        break

        if lattice_line != -1:
            latt_crys = [
                self.data[lattice_line + 2].strip().split(),
                self.data[lattice_line + 3].strip().split(),
                self.data[lattice_line + 4].strip().split()
            ]
            latt_crys = np.array(latt_crys, dtype=float)
            latt_pmg = struc.lattice.matrix
            # latt_crys = latt_pmg @ rot, rotation matrix obtained from the last config
            rot = np.linalg.inv(latt_pmg) @ latt_crys
            struc = rotate_lattice(struc, rot)

        self.geometry = struc

        # Write gui files
        if write_gui == True:
            # Conventional atomic numbers
            zconv = [[i, self.atom_numbers[i]] for i in range(self.n_atoms)]
            if gui_name == None:
                gui_name = os.path.splitext(self.name)[0]
                gui_name = '{}.gui'.format(gui_name)

            if symmetry == 'pymatgen':
                # zconv needs to be None, given SpacegroupAnalyzer usage
                gui = cry_pmg2gui(struc, gui_file=gui_name,
                                  symmetry=True, zconv=None, **kwargs)
            elif symmetry == None:
                gui = cry_pmg2gui(struc, gui_file=gui_name,
                                  symmetry=False, zconv=zconv)
            elif symmetry == 'initial':
                self.get_symmops()
                gui = cry_pmg2gui(struc, gui_file=None,
                                  symmetry=False, zconv=zconv)
                gui.symmops = self.symmops
                # Translation vector of the symmop to cartesian
                for n, symmops in enumerate(gui.symmops):
                    symmops[3][0] = symmops[3][0]*struc.lattice.a
                    symmops[3][1] = symmops[3][1]*struc.lattice.b
                    symmops[3][2] = symmops[3][2]*struc.lattice.c
                    gui.symmops[n][3] = symmops[3]
                gui.n_symmops = self.n_symmops
                gui.space_group = self.sg_number
                gui.write_gui(gui_name, symm=True)
            else:
                warnings.warn('Symmetry adapted from reference geometry. Make sure that is desired.',
                              stacklevel=2)
                gui_ref = Crystal_gui().read_gui(symmetry)
                gui = cry_pmg2gui(struc, gui_file=None,
                                  symmetry=False, zconv=zconv)
                # Replace the symmops with the reference file
                gui.symmops = gui_ref.symmops
                gui.n_symmops = gui_ref.n_symmops
                gui.space_group = gui_ref.space_group
                gui.write_gui(gui_list[idx_s], symm=True)

        return self.geometry

    def get_last_geom(self, write_gui_file=False, symm_info='initial'):
        """ Return the last optimised geometry.

        Returns:
            self.last_geom (list): To be documented.

        """
        struc = self.get_geometry(
            initial=False, write_gui=write_gui_file, symmetry=symm_info)
        if 'Molecule' in str(type(struc)):
            self.last_geom = [[[500., 0., 0.], [0., 500., 0.], [0., 0., 500.]],
                              self.atom_numbers,
                              self.atom_positions_cart.tolist()]
        else:
            self.last_geom = [struc.lattice.abc,
                              struc.lattice.angles,
                              struc.lattice.volume,
                              struc.lattice.matrix.tolist(),
                              self.atom_numbers,
                              self.atom_positions_cart.tolist()]

        return self.last_geom

    def get_lattice(self, initial=True):
        """
        Returns the lattice of the system. Unit: Angstrom.

        Args:
            initial (bool): Read the initial or last lattice. Useful in
                case of geometry optimization.
        Returns:
            self.lattice (np.ndarray): Lattice of the system.
        """
        import re
        import warnings

        import numpy as np

        ndimen = self.get_dimensionality()
        self.lattice = None
        if ndimen == 0:
            warnings.warn('0D system. No lattice.')
            return self.lattice

        self.get_geometry(initial=initial, write_gui=False)
        self.lattice = self.geometry.lattice.matrix

        return self.lattice

    def get_reciprocal_lattice(self, initial=True):
        """
        Returns the reciprocal lattice of the system. Unit: Angstrom^-1.
        Reciprocal lattice is consistent with CRYSTAL: 2pi is added.

        Args:
            initial (bool): Read the initial or last lattice. Useful in
                case of geometry optimization.
        Returns:
            self.reciprocal_lattice (np.ndarray): Lattice of the system.
        """
        import warnings

        ndimen = self.get_dimensionality()
        self.reciprocal_lattice = None
        if ndimen == 0:
            warnings.warn('0D system. No lattice.')
            return self.reciprocal_lattice

        self.get_lattice(initial=initial)
        self.reciprocal_lattice = self.geometry.lattice.reciprocal_lattice.matrix

        return self.reciprocal_lattice

    @property
    def sg_number(self):
        """
        CRYSTAL 0~3D space group number. Before geometry editing.
        """
        import re

        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        ndimen = self.get_dimensionality()
        sg = 'unknown'
        if ndimen == 0:
            rexp = r'^\s*POINT  GROUP  N\.'
        elif ndimen == 1:
            rexp = r'^POLYMER GROUP N\.'
        elif ndimen == 2:
            rexp = r'^\s*TWO\-SIDED PLANE GROUP N\.'

        if ndimen < 3:
            for nline, line in enumerate(self.data):
                if re.match(rexp, line):
                    if ndimen == 0:
                        sg = int(line.strip().split()[3])
                    elif ndimen == 1:
                        sg = int(line.strip()[16:19].strip())
                    elif ndimen == 2:
                        sg = int(line.strip().split()[3])
                    break
        else:
            struc = self.get_geometry(initial=True, write_gui=False)
            sg = SpacegroupAnalyzer(struc).get_space_group_number()
        return sg

    @property
    def sg_symbol(self):
        """
        CRYSTAL 0~3D 0~3D space group symbol. Before geometry editing.
        """
        import re
        import warnings

        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        ndimen = self.get_dimensionality()
        sg = 'unknown'
        if ndimen == 0:
            rexp = r'^\s*POINT  GROUP  N\.'
        elif ndimen == 1:
            rexp = r'^POLYMER GROUP N\.'
        elif ndimen == 2:
            rexp = r'^\s*TWO\-SIDED PLANE GROUP N\.'
        elif ndimen == 3:
            rexp = r'^\s*SPACE GROUP \(CENTROSYMMETRIC\)'

        for nline, line in enumerate(self.data):
            if re.match(rexp, line):
                if ndimen == 0:
                    sg = line.strip().split(':')[1].split('OR')[1].strip()
                elif ndimen == 1:
                    sg = line.strip().split(':')[1].strip()
                elif ndimen == 2:
                    sg = line.strip().split(':')[1].strip()
                else:
                    sg = line.strip().split(':')[1].strip()
                break

        if sg == 'unknown':
            warnings.warn('Symmstry information lost. Trying to get from pymatgen...',
                          stacklevel=2)
            struc = self.get_geometry(initial=True, write_gui=False)
            sg = SpacegroupAnalyzer(struc).get_space_group_symbol()
        return sg

    @property
    def n_atoms(self):
        """
        Number of atoms. After geometry editing.
        """
        struc = self.get_geometry(initial=True, write_gui=False)
        return struc.num_sites

    @property
    def atom_symbols(self):
        """
        Atom symbols. After geometry editing.
        """
        from mendeleev import element

        symbol = []
        for i in self.atom_numbers:
            symbol.append(element(int(i % 100)).symbol)

        return symbol

    @property
    def atom_numbers(self):
        """
        Conventional atom numbers. After geometry editing.
        """
        from CRYSTALClear.base.crysout import GeomBASE

        _, conv_z = GeomBASE.read_conv_z(self.data, 0)
        return conv_z

    @property
    def atom_positions(self):
        """
        Composite fractional / Cartesian atomic coordinates. Consistent
        with CRYSTAL definitions. 3D: Fractional; 2D: Frac, Frac, Cart; 1D
        Frac, Cart, Cart; 0D: Cart, Cart, Cart. After geometry editing
        (before optimization).
        """
        import numpy as np

        ndimen = self.get_dimensionality()
        struc = self.get_geometry(initial=True, write_gui=False)

        if ndimen == 0:
            composite_coord = struc.cart_coords.tolist()
        elif ndimen == 3:
            composite_coord = struc.frac_coords.tolist()
        else:
            cart_coord = struc.cart_coords.tolist()
            frac_coord = struc.frac_coords.tolist()
            composite_coord = []
            for i in range(struc.num_sites):
                composite_coord.append(
                    frac_coord[i][:ndimen] + cart_coord[i][ndimen:])

        return np.array(composite_coord, dtype=float)

    @property
    def atom_positions_cart(self):
        """
        Cartesian atomic coordinates. After geometry editing (before optimization).
        """
        struc = self.get_geometry(initial=True, write_gui=False)
        return struc.cart_coords

    @property
    def atom_positions_frac(self):
        """
        Fractional atomic coordinates. After geometry editing (before optimization).
        """
        import warnings

        ndimen = self.get_dimensionality()
        if ndimen != 3:
            warnings.warn("Low dimension systems. Property 'atom_positions' is called instead.",
                          stacklevel=2)
            return self.atom_positions
        else:
            struc = self.get_geometry(initial=True, write_gui=False)
            return struc.frac_coords

    def get_trans_matrix(self):
        """
        Get cell transformation matrix

        Returns:
            self.trans_matrix (np.ndarray): 3\*3 array of supercell
                expansion matrix
        """
        import re

        import numpy as np

        ndimen = self.get_dimensionality()
        self.trans_matrix = np.eye(3, dtype=float)
        for i, line in enumerate(self.data[:self.eoo]):
            if re.match(r'^\s+EXPANSION MATRIX OF PRIMITIVE CELL', line):
                mx = np.array([
                    self.data[i + 1].strip().split()[1:],
                    self.data[i + 2].strip().split()[1:],
                    self.data[i + 3].strip().split()[1:],
                ], dtype=float)
                self.trans_matrix[:ndimen, :ndimen] = mx[:ndimen, :ndimen]
                break

        return self.trans_matrix

    def get_primitive_geometry(self, initial=True, write_gui=False, gui_name=None,
                               symmetry='pymatgen', **kwargs):
        """
        Get the primitive geometry, reduced by cell transformation matrix
        inversed.

        .. note::
            This is not the standard primitive cell. This method returns
            geometry before supercell expansion keywords such as "SUPERCEL"
            or "SCELPHONO".

            Conventional atomic numbers are not available.

        Args:
            initial (bool): Read the initial or last geometry. Useful in
                case of geometry optimization.
            write_gui (bool): If True, write .gui file
            gui_name (str): Valid only if ``write_gui = True``. Gui file
                is named as "gui_name". If None, use "basename.gui". The
                basename is the same as output file.
            symmetry (str): Valid only if ``write_gui = True``. "pymatgen"
                to use symmetry info from a pymatgen SpacegroupAnalyzer;
                "initial" to use symmstry information on output file. If
                None, no symmetry. Otherwise it is taken from the existing
                gui file.
            **kwargs: Valid only if ``write_gui = True`` and
                ``symmetry = 'pymatgen'``.  Passed to Pymatgen
                SpacegroupAnalyzer object.

        Returns:
            self.primitive_geometry (Structure | Molecule): A pymatgen
                Structure or molecule object.
        """
        import os
        import re
        import warnings

        import numpy as np

        from CRYSTALClear.convert import cry_pmg2gui
        from CRYSTALClear.crystal_io import Crystal_gui
        from CRYSTALClear.geometry import get_pcel

        ndimen = self.get_dimensionality()
        self.get_geometry(initial=initial, write_gui=False)
        self.get_trans_matrix()

        if ndimen == 0:
            warnings.warn('0D system. Nothing to reduce.', stacklevel=2)
            self.primitive_geometry = self.geometry
            return self.primitive_geometry

        shrink_mx = np.linalg.inv(self.trans_matrix)
        pstruc = get_pcel(self.geometry, self.trans_matrix)

        self.primitive_geometry = pstruc
        # Write gui files
        if write_gui == True:
            if gui_name == None:
                gui_name = os.path.splitext(self.name)[0]
                gui_name = '{}.gui'.format(gui_name)

            if symmetry == 'pymatgen':
                gui = cry_pmg2gui(pstruc, gui_file=gui_name,
                                  symmetry=True, **kwargs)
            elif symmetry == None:
                gui = cry_pmg2gui(pstruc, gui_file=gui_name, symmetry=False)
            elif symmetry == 'initial':
                self.get_symmops()
                gui = cry_pmg2gui(pstruc, gui_file=None, symmetry=False)
                gui.symmops = self.symmops
                gui.n_symmops = self.n_symmops
                gui.space_group = self.sg_number
                gui.write_gui(gui_name, symm=True)
            else:
                warnings.warn('Symmetry adapted from reference geometry. Make sure that is desired.',
                              stacklevel=2)
                gui_ref = Crystal_gui().read_gui(symmetry)
                gui = cry_pmg2gui(pstruc, gui_file=None, symmetry=False)
                # Replace the symmops with the reference file
                gui.symmops = gui_ref.symmops
                gui.n_symmops = gui_ref.n_symmops
                gui.space_group = gui_ref.space_group
                gui.write_gui(gui_name, symm=True)

        return self.primitive_geometry

    def get_primitive_lattice(self, initial=True):
        """
        Returns the primitive lattice of the system reduced by cell
        transformation matrix inverse. Unit: Angstrom.

        Args:
            initial (bool): Read the initial or last lattice. Useful in
                case of geometry optimization.
        Returns:
            self.primitive_lattice (np.ndarray): Lattice of the system.
        """
        import warnings

        ndimen = self.get_dimensionality()
        self.primitive_lattice = None
        if ndimen == 0:
            warnings.warn('0D system. No lattice.')
            return self.primitive_lattice

        self.get_primitive_geometry(initial=initial, write_gui=False)
        self.primitive_lattice = self.primitive_geometry.lattice.matrix

        return self.primitive_lattice

    def get_primitive_reciprocal_lattice(self, initial=False):
        """
        Returns the primitive reciprocal lattice of the system before
        expansion by cell transformation matrix inverse. Unit: Angstrom^-1.

        Args:
            initial (bool): Read the initial or last lattice. Useful in
                case of geometry optimization.
        Returns:
            self.primitive_reciprocal_lattice (np.ndarray): Lattice of the system.
        """
        import warnings

        ndimen = self.get_dimensionality()
        self.primitive_reciprocal_lattice = None
        if ndimen == 0:
            warnings.warn('0D system. No lattice.')
            return self.primitive_reciprocal_lattice

        self.get_primitive_lattice(initial=initial)
        self.primitive_reciprocal_lattice = self.primitive_geometry.lattice.reciprocal_lattice.matrix

        return self.primitive_reciprocal_lattice

    def get_config_analysis(self, return_multiplicity=False):
        """
        Return the configuration analysis for solid solutions (CONFCON keyword in input)

        Args:
            return_multiplicity (bool, optional): Return multiplicity
            information (default is False).
        Returns:
            list or str: Configuration analysis if available, or a warning message
        """
        import re

        import numpy as np

        # Check this is a configuration analysis calculation
        try:
            begin = self.data.index(
                '                             CONFIGURATION ANALYSIS\n')
        except:
            return "WARNING: this is not a CONFCNT analysis."

        for line in self.data[::-1]:
            if '----' in line:
                dash_line = line.rstrip().lstrip()
                break

        for i, line in enumerate(self.data[begin:]):
            if re.match(r'^ COMPOSITION', line):
                self.n_classes = line.split()[9]
                original_atom = str(line.split()[2])
                begin = begin+i
        config_list = []

        # Read all the configurations
        for line in self.data[begin:]:
            if not re.match(r'^   WARNING', line):
                config_list.extend(line.split())

        '''multiplicity = []

        for i in range(len(config_list)):
            if config_list[i] == 'MULTIPLICITY' and i < len(config_list) - 1:
                number = re.findall(r'\d+', config_list[i+1])
                if number:
                    config_list.append(int(number[0]))'''

        config_list = np.array(config_list)
        warning = np.where(config_list == 'WARNING')
        config_list = np.delete(config_list, warning)
        atom1_begin = np.where(config_list == original_atom)[0]
        atom1_end = np.where(
            config_list == dash_line)[0]
        atom2_begin = np.where(config_list == 'XX')[0]
        atom2_end = np.where(
            config_list == '<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')[0]
        end = np.where(
            config_list == '===============================================================================')[0][-1]
        atom2_end = np.append(atom2_end, end)
        atom_type1 = []
        atom_type2 = []
        if return_multiplicity == True:
            input_string = ' '.join(config_list.tolist())
            matches = re.findall(r'MULTIPLICITY\s*:\s*(\d+)', input_string)

            # multiplicity_tmp = config_list[np.where(config_list == 'MULTIPLICITY')[0]+1]
            multiplicity = [int(x) for x in matches]
            self.multiplicity = multiplicity
        config_list = config_list.tolist()
        for i in range(len(atom1_end)):
            atom_type1.append(
                [int(x) for x in config_list[atom1_begin[i+1]+1:atom1_end[i]]])
            atom_type2.append(
                [int(x) for x in config_list[atom2_begin[i]+1:atom2_end[i]]])

        self.atom_type1 = atom_type1
        self.atom_type2 = atom_type2

        if return_multiplicity == True:
            return [self.atom_type1, self.atom_type2, multiplicity]
        else:
            return [self.atom_type1, self.atom_type2]

    #### SCF / OPT Convergence ####

    def get_convergence(self, history=False):
        """
        The upper level of ``get_scf_convergence`` and ``get_opt_convergence``.
        For analysing the geometry and energy convergence.

        .. note::

            It might not work well with SCF / OPT cycles of multiple systems
            such as PREOPTGEOM + EOS calculations

        Args:
            history (bool): If true, the convergence history of optimisation
                (energy,gradient, displacement) / SCF is returned.

        Returns:
            self (Crystal_output): New attributes listed below
            self.final_energy (float) The converged energy of SCF / Opt. Unit: eV

        For other attributes, see ``get_scf_convergence`` and
        ``get_opt_convergence`` methods on the same page.
        """
        import re
        import warnings

        self.final_energy = None
        for line in self.data[self.eoo::-1]:
            if re.match(r'\s\W OPT END - CONVERGED', line) != None:
                data = line.strip().split()
                self.final_energy = units.H_to_eV(float(data[7]))
                self.opt_cycles = int(data[-2])
                if re.search('CONVERGED', line):
                    self.opt_status = 'converged'
                elif re.search('FAILED', line):
                    self.opt_status = 'failed'
                else:
                    warnings.warn('Unknown termination.', stacklevel=2)
                    self.opt_status = 'unknown'
                is_scf = False
                break

            elif re.match(r'^ == SCF ENDED', line) != None:
                data = line.strip().split()
                self.final_energy = units.H_to_eV(float(data[8]))
                self.scf_cycles = int(data[-1])
                if re.search('CONVERGENCE', line):
                    self.scf_status = 'converged'
                elif re.search('TOO MANY CYCLES', line):
                    self.scf_status = 'too many cycles'
                else:
                    warnings.warn('Unknown termination.', stacklevel=2)
                    self.scf_status = 'unknown'
                is_scf = True
                break

        if self.final_energy == None:
            warnings.warn('No final energy found in the output file. self.final_energy = None',
                          stacklevel=2)

        if history == True:
            if is_scf == True:
                self.get_scf_convergence(all_cycles=False)
            else:
                self.get_opt_convergence()

        return self

    #### SCF ####

    def get_scf_convergence(self, all_cycles=False):
        """
        Returns the scf convergence energy and energy difference. A wrapper of
        ``CRYSTALClear.base.SCFBASE.read_convergence``.

        Args:
            all_cycles (bool, optional): Return all SCF steps for a geometry
                opt. The "ONELOG'" CRYSTAL keyword is needed.

        Returns:
            self (Crystal_output): New attributes listed below
            self.scf_cycles (int | array): Number of cycles. Array if
                ``all_cycles=True``.
            self.scf_status (str | list): "terminated", "converged",
                "too many cycles" and "unknown". List if ``all_cycles=True``.
            self.scf_energy (array): SCF energy convergence. Unit: eV
            self.scf_deltae (array): Energy difference. Unit: eV
        """
        import numpy as np

        from CRYSTALClear.base.crysout import SCFBASE

        if all_cycles == True:
            self.scf_cycles = []
            self.scf_status = []
            self.scf_energy = []
            self.scf_deltae = []

        countline = 0
        while countline < self.eoo:
            output = SCFBASE.read_convergence(self.data[:self.eoo], countline)
            if all_cycles == False:
                self.scf_cycles = output[1]
                self.scf_status = output[2]
                self.scf_energy = output[3]
                self.scf_deltae = output[4]
                break
            else:
                countline = output[0]
                self.scf_cycles.append(output[1])
                self.scf_status.append(output[2])
                self.scf_energy.append(output[3])
                self.scf_deltae.append(output[4])
                countline += 1

        if all_cycles == True:
            self.scf_cycles = np.array(self.scf_cycles, dtype=int)
            self.scf_energy = np.array(self.scf_energy)
            self.scf_deltae = np.array(self.scf_deltae)

        return self

    def get_fermi_energy(self, history=False):
        """
        Returns the system Fermi energy.

        Args:
            history (bool): Whether to read the convergence history of Fermi energy.

        Returns:
            self.fermi_energy (float | array): Fermi energy of the system. For
                spin-polarized insulating systems, ``self.fermi_energy`` would
                be either a 2\*1 array (``history=False``) or a nCYC\*2 array
                (``history=True``).
        """
        from CRYSTALClear.base.crysout import SCFBASE

        output = SCFBASE.read_fermi_energy(
            self.data[:self.eoo], self.eoo - 1, history=history)
        self.spin_pol = output[1]
        self.fermi_energy = output[2]

        return self.fermi_energy

    def get_band_gap(self, history=False):
        """
        Returns the system band gap.

        Args:
            history (bool): Whether to read the convergence history of band gap.

        Returns:
            self.band_gap (float | array): Band gap of the system. For
                spin-polarized systems, ``self.band_gap`` would be either a
                2\*1 array (``history=False``) or a nCYC\*2 array
                (``history=True``).
        """
        from CRYSTALClear.base.crysout import SCFBASE

        output = SCFBASE.read_band_gap(
            self.data[:self.eoo], self.eoo - 1, history=history)
        self.spin_pol = output[1]
        self.band_gap = output[2]

        return self.band_gap

    def get_mulliken_charges(self):
        """
        Return the atomic Mulliken charges (PPAN keyword in input).

        Returns:
            self.mulliken_charges (array): natom\*1 for non spin-polarised systems.
                natom\*3 for spin-polarised systems. [total, :math:`alpha`, :math:`beta`].
        """
        import re
        import warnings

        import numpy as np

        mulliken = []  # empty, 1*1 or 2*1 list
        countline = 0
        countm = 0
        while countline < self.eoo:
            line = self.data[countline]
            if re.match(r'\s*MULLIKEN POPULATION ANALYSIS', line):
                mulliken_charge = []  # natom*1
                countline += 4
                line2 = self.data[countline]
                while len(line2.strip()) != 0:
                    if re.match(r'^\s+[0-9]+\s+[A-Z, a-z]+\s+[0-9+]', line2):
                        data = line2.strip().split()
                        mulliken_charge.append(data[3])
                    countline += 1
                    line2 = self.data[countline]

                mulliken.append(mulliken_charge)
                continue
            else:
                countline += 1

        if len(mulliken) == 0:
            warnings.warn('Mulliken analysis not found.', stacklevel=2)
            self.mulliken_charges = np.array([], dtype=float)
        elif len(mulliken) == 1:
            self.mulliken_charges = np.array(mulliken[0], dtype=float)
            self.spin_pol = False
        else:
            apb = np.array(mulliken[0], dtype=float)
            amb = np.array(mulliken[1], dtype=float)
            self.mulliken_charges = np.array(
                [apb, (apb + amb) / 2, (apb - amb) / 2])
            self.spin_pol = True

        return self.mulliken_charges

    def get_final_energy(self):
        """
        Get the final energy of the system. A wrapper of ``self.get_convergence``.

        Returns:
            self.final_energy (float): The final energy of the system.
        """
        self.get_convergence(history=False)

        return self.final_energy

    def get_num_cycles(self):
        """Deprecated

        Returns:
            self.scf_cycles (int): Number of SCF cycles.
        """
        import warnings

        warnings.warn('Deprecated. Use get_scf_convergence.', stacklevel=2)
        self.get_scf_convergence(all_cycles=False)

        return self.scf_cycles

    #### Optimization ####

    def get_opt_convergence_energy(self):
        """Deprecated. Returns the energy for each opt step.

        Returns:
            self.opt_energy (array): Energy for each optimization step.
        """
        import warnings

        warnings.warn('Deprecated. Use get_opt_convergence.', stacklevel=2)
        self.get_opt_convergence()

        return self.opt_energy

    def get_opt_convergence(self, primitive=False, scf_history=False,
                            write_gui=False, gui_name=None,
                            symmetry='pymatgen', **kwargs):
        """
        Returns optimisation convergence. A wrapper of
        ``CRYSTALClear.base.OptBASE.read_convergence``.

        Args:
            primitive (bool): Restore the primitive cell (multiply by the
                cell transform matrix inversed)
            scf_history (bool): Read SCF history of each optimisation step.
                Keyword "ONELOG" is needed. Please refer to
                ``self.get_scf_convergence(all_cycles=True)`` method.
            write_gui (bool): If True, write .gui file of each step
            gui_name (str): Valid only if ``write_gui = True``. Gui file
                is named as "gui_name-optxxx.gui". If "None", use
                "basename-optxxx.gui". The basename is the same as output.
            symmetry (str): Valid only if ``write_gui = True``. "pymatgen"
                to use symmetry info from a pymatgen SpacegroupAnalyzer;
                "initial" to use symmetry information on output file. If
                "None", no symmetry. Otherwise it is taken from the existing
                gui file.
            **kwargs: Valid only if ``write_gui = True`` and
                ``symmetry = 'pymatgen'``.  Passed to Pymatgen
                SpacegroupAnalyzer object.

        Returns:
            self (Crystal_output): New attributes listed below
            self.opt_cycles (int): Number of cycles.
            self.opt_status (str): "terminated", "converged", "failed" and
            "unknown"
            self.opt_energy (array): Total energy convergence. Unit: eV
            self.opt_deltae (array): Total energy difference. Unit: eV
            self.opt_geometry (list): Pymatgen structure at each step.
            self.opt_maxgrad (array): Maximum gradient convergence. Unit: Hartree/Bohr
            self.opt_rmsgrad (array): RMS gradient convergence. Unit: Hartree/Bohr
            self.opt_maxdisp (array): Maximum displacement convergence. Unit: Bohr
            self.opt_rmsdisp (array): RMS displacement convergence. Unit: Bohr
        """
        import os
        import re

        import numpy as np
        from pymatgen.core.lattice import Lattice

        from CRYSTALClear.base.crysout import OptBASE
        from CRYSTALClear.convert import cry_pmg2gui
        from CRYSTALClear.crystal_io import Crystal_gui
        from CRYSTALClear.geometry import get_pcel, rotate_lattice

        ndimen = self.get_dimensionality()
        countline = 0

        # Initial geometry
        struc0 = self.get_geometry(initial=True, write_gui=False)

        # Initial SCF energy. No SCF when optimisation is restarted
        e0 = None
        self.opt_cycles = 0
        lattice_line = -1
        while countline < self.eoo:
            line = self.data[countline]
            # Initial step SCF
            if re.match(r'^\s*== SCF ENDED', line):
                line_data = line.strip().split()
                e0 = float(line_data[8])
                countline += 1
            # Initial step SCF empirical corrections
            elif re.match(r'^\s*TOTAL ENERGY \+', line):
                line_data = line.strip().split()
                e0 = float(line_data[-1])
                countline += 1
            # Initial Gradient
            elif re.match(r'^\s+MAX GRADIENT', line):
                line_data = line.strip().split()
                maxg0 = float(line_data[2])
                countline += 1
            elif re.match(r'^\s+RMS GRADIENT', line):
                line_data = line.strip().split()
                rmsg0 = float(line_data[2])
                countline += 1
            # Enter the Opt block
            elif re.match(r'^\s*OPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPT', line):
                output = OptBASE.read_optblock(self.data[:self.eoo], countline)
                self.opt_cycles = output[1]
                self.opt_status = output[2]
                self.opt_energy = output[3]
                self.opt_deltae = output[4]
                self.opt_geometry = output[5]
                self.opt_maxgrad = output[6]
                self.opt_rmsgrad = output[7]
                self.opt_maxdisp = output[8]
                self.opt_rmsdisp = output[9]
                break
            # Lattice matrix rotation issue
            elif re.match(r'^ DIRECT LATTICE VECTORS CARTESIAN', line):
                lattice_line = countline
                countline += 1
            else:
                countline += 1

        # Converged at initial step, or wrong file
        if self.opt_cycles == 0:
            if e0 == None:
                raise Exception('Valid data not found.')
            else:
                self.opt_cycles += 1
                self.opt_status = 'converged at initial step'
                self.opt_energy = np.array([e0,])
                self.opt_deltae = np.array([e0,])
                self.opt_geometry = [struc0,]
                self.opt_maxgrad = np.array([maxg0,])
                self.opt_rmsgrad = np.array([rmsg0,])
                self.opt_maxdisp = np.array([])
                self.opt_rmsdisp = np.array([])
        else:
            # restarted from previous opt
            if e0 == None:
                pass
            else:
                self.opt_cycles += 1
                self.opt_energy = np.concatenate([[e0,], self.opt_energy])
                self.opt_deltae = np.concatenate([[e0,], self.opt_deltae])
                self.opt_geometry = [struc0,] + self.opt_geometry
                self.opt_maxgrad = np.concatenate([[maxg0,], self.opt_maxgrad])
                self.opt_rmsgrad = np.concatenate([[rmsg0,], self.opt_rmsgrad])
                self.opt_maxdisp = np.concatenate([[0.,], self.opt_maxdisp])
                self.opt_rmsdisp = np.concatenate([[0.,], self.opt_rmsdisp])

        # Lattice matrix rotation issue
        if ndimen != 0 and lattice_line != -1:
            mx_crys = [
                self.data[lattice_line + 2].strip().split(),
                self.data[lattice_line + 3].strip().split(),
                self.data[lattice_line + 4].strip().split()
            ]
            mx_crys = np.array(mx_crys, dtype=float)
            # Reference lattice
            latt_ref = Lattice(mx_crys)
            # Pmg lattice
            latt_pmg = Lattice.from_parameters(
                a=latt_ref.a, b=latt_ref.b, c=latt_ref.c,
                alpha=latt_ref.alpha, beta=latt_ref.beta, gamma=latt_ref.gamma
            )
            mx_pmg = latt_pmg.matrix
            # mx_crys = mx_pmg @ rot
            rot = np.linalg.inv(mx_pmg) @ mx_crys
            for idx_s, s in enumerate(self.opt_geometry):
                self.opt_geometry[idx_s] = rotate_lattice(s, rot)

        # Get primitive cell
        if primitive == True:
            ndimen = self.get_dimensionality()
            if ndimen == 0:
                warnings.warn('0D system. Nothing to reduce.', stacklevel=2)
            else:
                self.get_trans_matrix()
                shrink_mx = np.linalg.inv(self.trans_matrix)
                for idx_s, s in enumerate(self.opt_geometry):
                    self.opt_geometry[idx_s] = get_pcel(s, self.trans_matrix)

        # SCF history
        if scf_history == True:
            self.get_scf_convergence(all_cycles=True)

        # Write gui files
        if write_gui == True:
            if gui_name == None:
                gui_name = os.path.splitext(self.name)[0]
            gui_list = [
                '{}-opt{:0=3d}.gui'.format(gui_name, i+1) for i in range(self.opt_cycles)]

            if symmetry == 'pymatgen':
                for idx_s, s in enumerate(self.opt_geometry):
                    gui = cry_pmg2gui(s, gui_file=gui_list[idx_s],
                                      symmetry=True, **kwargs)
            elif symmetry == None:
                for idx_s, s in enumerate(self.opt_geometry):
                    gui = cry_pmg2gui(
                        s, gui_file=gui_list[idx_s], symmetry=False)
            elif symmetry == 'initial':
                self.get_symmops()
                for idx_s, s in enumerate(self.opt_geometry):
                    gui = cry_pmg2gui(s, gui_file=None, symmetry=False)
                    gui.symmops = self.symmops
                    gui.n_symmops = self.n_symmops
                    gui.space_group = self.sg_number
                    gui.write_gui(gui_list[idx_s], symm=True)
            else:
                warnings.warn('Symmetry adapted from reference geometry. Make sure that is desired.',
                              stacklevel=2)
                gui_ref = Crystal_gui().read_gui(symmetry)
                for idx_s, s in enumerate(self.opt_geometry):
                    gui = cry_pmg2gui(s, gui_file=None, symmetry=False)
                    # Replace the symmops with the reference file
                    gui.symmops = gui_ref.symmops
                    gui.n_symmops = gui_ref.n_symmops
                    gui.space_group = gui_ref.space_group
                    gui.write_gui(gui_list[idx_s], symm=True)

        return self

    def get_forces(self, initial=True, grad=False):
        """
        Read forces.

        Args:
            initial (bool): Return forces from the initial calculation. If
                ``initial=False``, return to the last forces, which is valid
                only for geometry optimizations and the keyword "ONELOG" is
                needed in .d12 file.
            grad (bool): Return gradient convergence history. For optimizations
                only.

        Returns:
            self (Crystal_output): New attributes listed below
            self.forces_atoms (array): natom\*3 array. Atomic forces. Unit: Hartree/Bohr
            self.forces_cell (array): 3\*3 array. Cell forces, 3D only. Unit: Hartree/Bohr
            self.opt_maxgrad (array): Maximum gradient convergence. Unit: Hartree/Bohr
            self.opt_rmsgrad (array): RMS gradient convergence. Unit: Hartree/Bohr
        """
        import re
        import warnings

        import numpy as np

        if initial == False or grad == True:
            if ' OPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPT\n' not in self.data:
                warnings.warn(
                    'Not a geometry optimisation: Set initial = True and grad = False', stacklevel=2)
                initial = True
                grad = False

        self.forces_atoms = []
        self.forces_cell = []

        if grad == True:
            self.get_opt_convergence()

        if initial == True:
            for i, line in enumerate(self.data):
                if re.match(r'^ CARTESIAN FORCES IN HARTREE/BOHR \(ANALYTICAL\)', line):
                    for j in range(i+2, i+2+self.n_atoms):
                        self.forces_atoms.append(
                            self.data[j].strip().split()[2:])
                    self.forces_atoms = np.array(
                        self.forces_atoms, dtype=float)

                if re.match(r'^ GRADIENT WITH RESPECT TO THE CELL PARAMETER IN HARTREE/BOHR', line):
                    for j in range(i+4, i+7):
                        self.forces_cell.append(self.data[j].strip().split())
                    self.forces_cell = np.array(self.forces_cell, dtype=float)

        else:
            for i, line in enumerate(self.data[::-1]):
                if re.match(r'^ GRADIENT WITH RESPECT TO THE CELL PARAMETER IN HARTREE/BOHR', line):
                    for j in range(len(self.data)-i+3, len(self.data)-i+6):
                        self.forces_cell.append(self.data[j].strip().split())
                    self.forces_cell = np.array(self.forces_cell, dtype=float)

                if re.match(r'^ CARTESIAN FORCES IN HARTREE/BOHR \(ANALYTICAL\)', line):
                    for j in range(len(self.data)-i+1, len(self.data)-i+1+self.n_atoms):
                        self.forces_atoms.append(
                            self.data[j].strip().split()[2:])
                    self.forces_atoms = np.array(
                        self.forces_atoms, dtype=float)

        return self

    #### Lattice Dynamics ####

    def get_phonon(self, read_eigvt=False, rm_imaginary=True, rm_overlap=True,
                   imaginary_tol=-1e-4, q_overlap_tol=1e-4, eigvt_amplitude=1.):
        """
        Read phonon-related properties from output file.

        Args:
            read_eigvt (bool): Whether to read phonon eigenvectors and
                normalize it to 1.
            rm_imaginary (bool): Remove the modes with negative frequencies and
                set all the related properties to NaN.
            rm_overlap (bool): *For dispersion calculations* Remove repeated q
                points and recalculate their weights.
            imaginary_tol (float): Valid only if ``rm_imaginary = True``.  The threshold
                of negative frequencies.
            q_overlap_tol (float): Valid only if ``rm_overlap = True``.  The threshold of
                overlapping points, defined as the 2nd norm of the difference
                of fractional q vectors
            eigvt_amplitude (float | str): Valid if ``read_eigvt = True``.
                Amplitude of normalization, or "classical", "classical-rev",
                classical amplitude and revmove classical amplitude.


        .. note::

            In QHA calculations, the "q point" dimension refer to harmonic
            phonons computed. In other cases it refers to points in reciprocal
            space.

        Returns:
            self (Crystal_output): New attributes listed below
            self.edft (array[float]): :math:`E_{0}` Energy with empirical
                correction. Unit: kJ/mol.
            self.nqpoint (int): Number of q points
            self.qpoint (list[list[array[float], float]]): A nqpoint\*1 list of
                2\*1 list whose first element is a 3\*1 array of fractional
                coordinates and the second is its weight.
            self.nmode (array[int]): Number of modes at q point. nqpoint\*1
                array.
            self.frequency (array[float]): nqpoint\*nmode array ofvibrational
                frequency. Unit: THz
            self.intens (array[float]): nqpoint\*nmode array of harmonic
                intensiy. Unit: km/mol
            self.IR (array[bool]): nqpoint\*nmode array of boolean values
                specifying whether the mode is IR active
            self.Raman (array[bool]): nqpoint\*nmode array of boolean values
                specifying whether the mode is Raman active
            self.eigenvector (array[complex]): Valid if ``read_eigvt = True``.
                nqpoint\*nmode\*natom\*3 array of eigenvectors. Normalized to 1.
        """

        import re
        import numpy as np

        from CRYSTALClear.base.crysout import PhononBASE
        from CRYSTALClear.units import H_to_kjmol

        is_freq = False
        found_anti = True
        self.edft = []
        self.nqpoint = 0
        self.qpoint = []
        self.nmode = []
        self.frequency = []
        self.intens = []
        self.IR = []
        self.Raman = []
        self.eigenvector = []

        countline = 0
        while countline < self.eoo:
            line = self.data[countline]
            # Whether is a frequency file
            if re.match(r'^\s*\+\+\+\sSYMMETRY\sADAPTION\sOF\sVIBRATIONAL\sMODES\s\+\+\+', line):
                is_freq = True
                countline += 1
                continue
            # E_0 with empirical corrections
            # fix
            elif (re.match(r'^\s+CENTRAL POINT', line) and
                  re.match(r'^\sATOM', self.data[countline-1])):
                # fix
                self.edft.append(float(line.strip().split()[2]))
                countline += 1
                continue
            # Q point info + frequency
            # Dispersion
            elif re.match(r'^.+EXPRESSED IN UNITS\s+OF DENOMINATOR', line):
                shrink = int(line.strip().split()[-1])
                countline += 1
                continue
            elif re.match(r'\s+DISPERSION K POINT NUMBER', line):
                coord = np.array(line.strip().split()[7:10], dtype=float)
                weight = float(line.strip().split()[-1])
                self.qpoint.append([coord / shrink, weight])
                self.nqpoint += 1
                countline += 2
                # Read phonons
                phonon = PhononBASE.readmode_basic(
                    self.data[:self.eoo], countline)
                countline = phonon[0]
                self.frequency.append(phonon[1])
                self.intens.append(phonon[2])
                self.IR.append(phonon[3])
                self.Raman.append(phonon[4])
            # Gamma point
            elif re.match(r'^\s+MODES\s+EIGV\s+FREQUENCIES\s+IRREP', line) and self.nqpoint == 0:
                countline += 2
                # Read phonons
                phonon = PhononBASE.readmode_basic(
                    self.data[:self.eoo], countline)
                countline = phonon[0]
                self.frequency.append(phonon[1])
                self.intens.append(phonon[2])
                self.IR.append(phonon[3])
                self.Raman.append(phonon[4])
            # Phonon eigenvector
            # Gamma point: real numbers. Imaginary = 0
            elif re.match(r'^\s+NORMAL MODES NORMALIZED', line):
                if read_eigvt == False:
                    countline += 1
                    continue
                countline += 2
                eigvt = PhononBASE.readmode_eigenvector(
                    self.data[:self.eoo], countline)
                countline = eigvt[0]
                self.eigenvector.append(eigvt[1] + 0.j)
            # Dispersion: complex numbers
            elif re.match(r'^\s+MODES IN PHASE', line):
                if read_eigvt == False:
                    countline += 1
                    continue
                if found_anti == False:  # Real k point
                    self.eigenvector.append(tmp_eigvt)
                countline += 2
                found_anti = False
                eigvt = PhononBASE.readmode_eigenvector(
                    self.data[:self.eoo], countline)
                countline = eigvt[0]
                tmp_eigvt = eigvt[1] + 0.j
            elif re.match(r'^\s+MODES IN ANTI\-PHASE', line):
                if read_eigvt == False:
                    countline += 1
                    continue
                countline += 2
                found_anti = True
                eigvt_anti = PhononBASE.readmode_eigenvector(
                    self.data[:self.eoo], countline)
                countline = eigvt_anti[0]
                self.eigenvector.append(tmp_eigvt + eigvt_anti[1] * 1.j)
            # Other data
            else:
                countline += 1
                continue

        if is_freq == False:
            raise Exception('Not a frequency calculation.')
        if found_anti == False and read_eigvt == True:  # The last real k point
            self.eigenvector.append(tmp_eigvt)

        # Format data
        # HA/QHA Gamma point calculation
        if self.nqpoint == 0:
            self.nqpoint = len(self.edft)
            self.qpoint = [[np.zeros([3,]), 1.] for i in range(self.nqpoint)]
            if len(self.edft) == 1:
                self.edft = [self.edft[0] for i in range(self.nqpoint)]
        # Dispersion
        else:
            self.qpoint = [[i[0], i[1] / self.nqpoint] for i in self.qpoint]
            self.edft = [self.edft[0] for i in range(self.nqpoint)]

        self.edft = H_to_kjmol(np.array(self.edft))

        self.frequency = np.array(self.frequency)
        self.nmode = np.array([len(i) for i in self.frequency], dtype=int)
        if self.intens[0] == []:
            self.intens = []
            self.IR = []
            self.Raman = []
        else:
            self.intens = np.array(self.intens)

        if self.eigenvector != []:
            self.eigenvector = np.array(self.eigenvector)
            # already normalised to classical amplitude
            if str(eigvt_amplitude).lower() == 'classical':
                pass
            # remove classical amplitude
            elif str(eigvt_amplitude).lower() == 'classical-rev':
                struc = self.get_geometry(initial=False, write_gui=False)

            # To a specific value
            else:
                for idx_q in range(self.nqpoint):
                    self.eigenvector[idx_q] = PhononBASE.normalize_eigenvector(
                        self.eigenvector[idx_q],
                        amplitude=float(eigvt_amplitude),
                    )

        if rm_imaginary == True:
            self = PhononBASE.clean_imaginary(self, threshold=imaginary_tol)

        if rm_overlap == True and self.nqpoint > 1:
            self = PhononBASE.clean_q_overlap(self, threshold=q_overlap_tol)

        return self

    def get_q_info(self):
        """
        Deprecated.
        """
        import warnings

        warnings.warn(
            'This method is deprecated. Use `get_phonon`.', stacklevel=2)
        return self

    def get_mode(self):
        """
            Deprecated.
            """
        return self.get_q_info()

    def get_phonon_eigenvector(self):
        """
        Deprecated.
        """
        return self.get_q_info()

    def get_IR(self):
        """
        Extract harmonic IR spectrum from CRYSTAL output.
        
        Returns:
            self.IR_HO_0K (np.array): 2D array containing harmonic IR frequency and intensities computed at 0 K.
        """

        import numpy as np
        from CRYSTALClear.units import thz_to_cm

        self.get_phonon()

        self.IR_HO_0K = np.stack([thz_to_cm(self.frequency), self.intens], axis=2).reshape(-1, 2)

        return self

    def get_Raman(self):
        """
        Extract harmonic Raman spectrum from CRYSTAL output.

        Returns:
            self (Crystal_output): New attributes listed below.
            self.Ram_HO_0K_tot (np.array): 2D array containing harmonic Raman frequency and intensities (total) computed at 0 K.
            self.Ram_HO_0K_par (np.array): 2D array containing harmonic Raman frequency and intensities (parallel) computed at 0 K.
            self.Ram_HO_0K_per (np.array): 2D array containing harmonic Raman frequency and intensities (perpendicular) computed at 0 K.
            self.Ram_HO_0K_comp_xx (np.array): 2D array containing harmonic Raman frequency and intensities (xx component) computed at 0 K.

        Note:
            Please, note that for the sake of brevity, only the xx component attributes have been listed here, but the yy, zz, xy, xz, yz components are available as well.
        """

        import re
        import numpy as np

        savePC = False
        saveSC = False

        self.Ram_HO_0K_tot = []
        self.Ram_HO_0K_par = []
        self.Ram_HO_0K_per = []
        self.Ram_HO_0K_comp_xx = []
        self.Ram_HO_0K_comp_xy = []
        self.Ram_HO_0K_comp_xz = []
        self.Ram_HO_0K_comp_yy = []
        self.Ram_HO_0K_comp_yz = []
        self.Ram_HO_0K_comp_zz = []

        for i, line in enumerate(self.data):

            if re.match(r'\s*POLYCRYSTALLINE ISOTROPIC INTENSITIES', line):
                savePC = True

            if re.match(r'\s*SINGLE CRYSTAL DIRECTIONAL INTENSITIES', line):
                savePC = False
                saveSC = True

            if re.match(r'\s*NORMAL MODES NORMALIZED TO CLASSICAL AMPLITUDES', line):
                break

            if re.match(r'^\s*\d', line):

                if(savePC):
                    ndeg = int(line.split()[1].replace("-", "")) \
                            - int(line.split()[0].replace("-", "")) \
                            + 1
                    freq = float(line.split()[2])
                    intens_tot = float(line.split()[-3]) / ndeg 
                    intens_par = float(line.split()[-2]) / ndeg 
                    intens_per = float(line.split()[-1]) / ndeg 
                    for k in range(ndeg):
                        self.Ram_HO_0K_tot.append([freq, intens_tot]) 
                        self.Ram_HO_0K_par.append([freq, intens_par]) 
                        self.Ram_HO_0K_per.append([freq, intens_per]) 

                if(saveSC):
                    ndeg = int(line.split()[1].replace("-", "")) \
                            - int(line.split()[0].replace("-", "")) \
                            + 1
                    freq = float(line.split()[2])
                    intens_xx = float(line.split()[-6]) / ndeg 
                    intens_xy = float(line.split()[-5]) / ndeg 
                    intens_xz = float(line.split()[-4]) / ndeg 
                    intens_yy = float(line.split()[-3]) / ndeg 
                    intens_yz = float(line.split()[-2]) / ndeg 
                    intens_zz = float(line.split()[-1]) / ndeg 
                    for k in range(ndeg):
                        self.Ram_HO_0K_comp_xx.append([freq, intens_xx])
                        self.Ram_HO_0K_comp_xy.append([freq, intens_xy])
                        self.Ram_HO_0K_comp_xz.append([freq, intens_xz])
                        self.Ram_HO_0K_comp_yy.append([freq, intens_yy])
                        self.Ram_HO_0K_comp_yz.append([freq, intens_yz])
                        self.Ram_HO_0K_comp_zz.append([freq, intens_zz])

        self.Ram_HO_0K_tot = np.array(self.Ram_HO_0K_tot)
        self.Ram_HO_0K_par = np.array(self.Ram_HO_0K_par)
        self.Ram_HO_0K_per = np.array(self.Ram_HO_0K_per)
        self.Ram_HO_0K_comp_xx = np.array(self.Ram_HO_0K_comp_xx)
        self.Ram_HO_0K_comp_xy = np.array(self.Ram_HO_0K_comp_xy)
        self.Ram_HO_0K_comp_xz = np.array(self.Ram_HO_0K_comp_xz)
        self.Ram_HO_0K_comp_yy = np.array(self.Ram_HO_0K_comp_yy)
        self.Ram_HO_0K_comp_yz = np.array(self.Ram_HO_0K_comp_yz)
        self.Ram_HO_0K_comp_zz = np.array(self.Ram_HO_0K_comp_zz)

        return self


    def get_anh_const(self):
        """
        Extract anharmonic terms of the PES (ANHAPES).

        Returns:
            self (Crystal_output): New attributes listed below.
            self.PES_single(np.array): 2D array containing single-mode anharmonic force constants (col 1: mode; cols 2-3: cubic and quartic single-mode terms of the PES, respectively).
            self.PES_couple(np.array): 2D array contaning anharmonic force constants coupling two phonon modes (col 1: mode I; col 2: mode J; cols 3-7: two-mode anharmonic force constants of type IIJ, IJJ, IIIJ, IJJJ and IIJJ, respectively.)
            self.PES_triplet(np.array): 2D array contaning anharmonic force constant coupling three phonon modes (col 1: mode I; col 2: mode J; col 3: mode K; cols 4-7: three-mode anharmonic force constants of type IJK, IIJK, IJJK and IJKK, respectively)
        """

        import re
        import numpy as np

        save = False
        PES_single = []
        PES_couple = []
        PES_triplet = []

        for i, line in enumerate(self.data):

            if re.match(r'\s*CALCULATION OF CUBIC AND QUARTIC TERMS', line):
                save = True

            if re.match(r'\s*TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SCANPES', line):
                save = False

            if (re.match(r'\s*MODE:', line) and save):
                a = int(line.split()[1])
                for k, line2 in enumerate(self.data[i:]):
                    if re.match(r'\s*ETA', line2):
                        b = float((self.data[i+k]).split()[6])
                        c = float((self.data[i+k+1]).split()[7])
                        break
                PES_single.append([a, b, c])

            if (re.match(r'\s*COUPLE', line) and save):
                a = int(line.split()[3])
                b = int(line.split()[4])
                for k, line2 in enumerate(self.data[i:]):
                    if re.match(r'\s*ETA', line2):
                        c = float((self.data[i+k]).split()[6])
                        d = float((self.data[i+k+1]).split()[6])
                        e = float((self.data[i+k+2]).split()[7])
                        f = float((self.data[i+k+3]).split()[7])
                        g = float((self.data[i+k+4]).split()[7])
                        break
                PES_couple.append([a, b, c, d, e, f, g])

            if (re.match(r'\s*TRIPLET', line) and save):
                a = int(line.split()[3])
                b = int(line.split()[4])
                c = int(line.split()[5])
                for k, line2 in enumerate(self.data[i:]):
                    if re.match(r'\s*ETA', line2):
                        d = float((self.data[i+k]).split()[6])
                        e = float((self.data[i+k+1]).split()[7])
                        f = float((self.data[i+k+2]).split()[7])
                        g = float((self.data[i+k+3]).split()[7])
                        break
                PES_triplet.append([a, b, c, d, e, f, g])

        self.PES_single = np.array(PES_single)
        self.PES_couple = np.array(PES_couple)
        self.PES_triplet = np.array(PES_triplet)

        return self


    def get_anh_spectra(self):
        """
        Extract anharmonic (VSCF and VCI) IR and Raman spectra (in development).

        Returns:
            self (Crystal_output): New attributes listed below.
            self.IR_HO_0K (np.array): 2D array containing harmonic IR frequency and intensities computed at 0 K.
            self.IR_HO_T (np.array): 2D array containing harmonic IR frequency and intensities computed at temperature T.
            self.IR_VSCF_0K (np.array): 2D array containing VSCF IR frequency and intensities computed at 0 K.
            self.IR_VSCF_T (np.array): 2D array containing VSCF IR frequency and intensities computed at temperature T.
            self.IR_VCI_0K (np.array): 2D array containing VCI IR frequency and intensities computed at 0 K.
            self.IR_VCI_T (np.array): 2D array containing VCI IR frequency and intensities computed at temperature T.

            self.Ram_HO_0K_tot (np.array): 2D array containing harmonic Raman frequency and intensities (total) computed at 0 K.
            self.Ram_HO_0K_per (np.array): 2D array containing harmonic Raman frequency and intensities (perpendicular component ) computed at temperature 0 K.
            self.Ram_HO_0K_par (np.array): 2D array containing harmonic Raman frequency and intensities (parallel component ) computed at temperature 0 K.
            self.Ram_HO_T_tot (np.array): 2D array containing harmonic Raman frequency and intensities (total) computed at temperature T.
            self.Ram_HO_T_per (np.array): 2D array containing harmonic Raman frequency and intensities (perpendicular component ) computed at temperature T.
            self.Ram_HO_T_par (np.array): 2D array containing harmonic Raman frequency and intensities (parallel component ) computed at temperature T.

            self.Ram_VSCF_0K_tot (np.array): 2D array containing VSCF Raman frequency and intensities (total) computed at 0 K. self.Ram_VSCF_0K_per (array[float]): 2D numpy array containing VSCF Raman frequency and intensities (perpendicular component) computed at 0 K.
            self.Ram_VSCF_0K_par (np.array): 2D array containing VSCF Raman frequency and intensities (parallel component) computed at 0 K.
            self.Ram_VSCF_T_tot (np.array): 2D array containing VSCF Raman frequency and intensities (total) computed at temperature T. self.Ram_VSCF_T_per (array[float]): 2D numpy array containing VSCF Raman frequency and intensities (perpendicular component) computed at temperature T.
            self.Ram_VSCF_T_par (np.array): 2D array containing VSCF Raman frequency and intensities (parallel component) computed at temperature T.

            self.Ram_VCI_0K_tot (np.array): 2D array containing VCI Raman frequency and intensities (total) computed at 0 K.
            self.Ram_VCI_0K_per (np.array): 2D array containing VCI Raman frequency and intensities (perpendicular component) computed at 0 K.
            self.Ram_VCI_0K_par (np.array): 2D array containing VCI Raman frequency and intensities (parallel component) computed at 0 K.
            self.Ram_VCI_T_tot (np.array): 2D array containing VCI Raman frequency and intensities (total) computed at temperature T.
            self.Ram_VCI_T_per (np.array): 2D array containing VCI Raman frequency and intensities (perpendicular component) computed at temperature T.
            self.Ram_VCI_T_par (np.array): 2D array containing VCI Raman frequency and intensities (parallel component) computed at temperature T.

            self.Ram_HO_0K_comp_xx (np.array): 2D array containing harmonic Raman frequency and intensities (xx component) computed at 0 K.
            self.Ram_HO_T_comp_xx (np.array): 2D array containing harmonic Raman frequency and intensities (xx component) computed at temperature T.
            self.Ram_VSCF_0K_comp_xx (np.array): 2D array containing VSCF Raman frequency and intensities (xx component) computed at 0 K.
            self.Ram_VSCF_T_comp_xx (np.array): 2D array containing VSCF Raman frequency and intensities (xx component) computed at temperature T.
            self.Ram_VCI_0K_comp_xx (np.array): 2D array containing VCI Raman frequency and intensities (xx component) computed at 0 K.
            self.Ram_VCI_T_comp_xx (np.array): 2D array containing VCI Raman frequency and intensities (xx component) computed at temperature T.

        Note:
            Please, note that for the sake of brevity, only the xx Raman component attributes have been listed here, but the yy, zz, xy, xz, yz components are available as well.
        """

        import re
        import numpy as np

        # Initialize some logical variables
        save = False
        HO_IR = False
        HO_Ram_0K_tot = False
        HO_Ram_T_tot = False
        HO_Ram_0K_comp = False
        HO_Ram_T_comp = False
        VSCF_IR = False
        VSCF_Ram_0K_tot = False
        VSCF_Ram_0K_comp = False
        VSCF_Ram_T_tot = False
        VSCF_Ram_T_comp = False
        VCI_IR = False
        VCI_Ram_0K_tot = False
        VCI_Ram_0K_comp = False
        VCI_Ram_T_tot = False
        VCI_Ram_T_comp = False
        # anscan
        imode_IR = []
        imode_Ram = []
        anscan_IR = False
        anscan_Ram_0K_tot = False
        anscan_Ram_0K_comp = False
        anscan_Ram_T_tot = False
        anscan_Ram_T_comp = False
        # anscan

        # Initialize some member variables
        IR_HO = []
        IR_VSCF = []
        IR_VCI = []
        Ram_HO_0K_tot = []
        Ram_HO_T_tot = []
        Ram_HO_0K_comp = []
        Ram_HO_T_comp = []
        Ram_VSCF_0K_tot = []
        Ram_VSCF_T_tot = []
        Ram_VSCF_0K_comp = []
        Ram_VSCF_T_comp = []
        Ram_VCI_0K_tot = []
        Ram_VCI_T_tot = []
        Ram_VCI_0K_comp = []
        Ram_VCI_T_comp = []
        # anscan
        IR_anscan = []
        # anscan

        # Initialize some buffers
        bufferHO_IR = []
        bufferHO_Ram_0K_tot = []
        bufferHO_Ram_T_tot = []
        bufferHO_Ram_0K_comp = []
        bufferHO_Ram_T_comp = []
        bufferVSCF_IR = []
        bufferVSCF_Ram_0K_tot = []
        bufferVSCF_Ram_T_tot = []
        bufferVSCF_Ram_0K_comp = []
        bufferVSCF_Ram_T_comp = []
        bufferVCI_IR = []
        bufferVCI_Ram_0K_tot = []
        bufferVCI_Ram_T_tot = []
        bufferVCI_Ram_0K_comp = []
        bufferVCI_Ram_T_comp = []
        # anscan
        bufferanscan_IR = []
        bufferanscan_Ram_0K_tot = []
        bufferanscan_Ram_T_tot = []
        bufferanscan_Ram_0K_comp = []
        bufferanscan_Ram_T_comp = []
        # anscan

        # Big loop over lines of CRYSTAL output file -->
        for i, line in enumerate(self.data):

            if re.match(r'\s*HARMONIC IR SPECTRUM', line):
                HO_IR = True
                save = True

            if re.match(r'\s*ANHARMONIC IR SPECTRUM', line):
                if re.match(r'\s*VSCF', self.data[i-1]):
                    VSCF_IR = True
                elif re.match(r'\s*VCI*', self.data[i-1]):
                    VCI_IR = True
                # anscan
                elif re.match(r'\s*MODE*', self.data[i-1]):
                    anscan_IR = True
                    imode_IR.append(int(self.data[i-1].split()[1]))
                # anscan
                save = True

            if re.match(r'\s*HARMONIC RAMAN SPECTRUM', line):
                if re.match(r'\s*\[ 0 K \]', self.data[i+1]):
                    if re.match(r'\s*I_TOT', self.data[i+4]):
                        HO_Ram_0K_tot = True
                    else:
                        HO_Ram_0K_comp = True
                else:
                    if re.match(r'\s*I_TOT', self.data[i+4]):
                        HO_Ram_T_tot = True
                    else:
                        HO_Ram_T_comp = True
                save = True

            if re.match(r'\s*ANHARMONIC RAMAN SPECTRUM', line):
                if re.match(r'\s*VSCF', self.data[i-1]):
                    if re.match(r'\s*\[ 0 K \]', self.data[i+1]):
                        if re.match(r'\s*I_TOT', self.data[i+4]):
                            VSCF_Ram_0K_tot = True
                        else:
                            VSCF_Ram_0K_comp = True
                    else:
                        if re.match(r'\s*I_TOT', self.data[i+4]):
                            VSCF_Ram_T_tot = True
                        else:
                            VSCF_Ram_T_comp = True
                save = True
                if re.match(r'\s*VCI*', self.data[i-1]):
                    if re.match(r'\s*\[ 0 K \]', self.data[i+1]):
                        if re.match(r'\s*I_TOT', self.data[i+4]):
                            VCI_Ram_0K_tot = True
                        else:
                            VCI_Ram_0K_comp = True
                    else:
                        if re.match(r'\s*I_TOT', self.data[i+4]):
                            VCI_Ram_T_tot = True
                        else:
                            VCI_Ram_T_comp = True
                # anscan
                if re.match(r'\s*MODE*', self.data[i-1]):
                    if re.match(r'\s*\[ 0 K \]', self.data[i+1]):
                        if re.match(r'\s*I_TOT', self.data[i+4]):
                            imode_Ram.append(int(self.data[i-1].split()[1]))
                            anscan_Ram_0K_tot = True
                        else:
                            anscan_Ram_0K_comp = True
                    else:
                        if re.match(r'\s*I_TOT', self.data[i+4]):
                            anscan_Ram_T_tot = True
                        else:
                            anscan_Ram_T_comp = True
                # anscan

                save = True

            if re.match(r'\s*HHHHHHHHHHHHH', line):
                save = False
                HO_IR = False
                HO_Ram_0K_tot = False
                HO_Ram_T_tot = False
                HO_Ram_0K_comp = False
                HO_Ram_T_comp = False
                VSCF_IR = False
                VSCF_Ram_0K_tot = False
                VSCF_Ram_0K_comp = False
                VSCF_Ram_T_tot = False
                VSCF_Ram_T_comp = False
                VCI_IR = False
                VCI_Ram_0K_tot = False
                VCI_Ram_0K_comp = False
                VCI_Ram_T_tot = False
                VCI_Ram_T_comp = False
                # anscan
                anscan_IR = False
                anscan_Ram_0K_tot = False
                anscan_Ram_0K_comp = False
                anscan_Ram_T_tot = False
                anscan_Ram_T_comp = False
                # anscan

            if save:
                if HO_IR:
                    bufferHO_IR.append(line)
                if HO_Ram_0K_tot:
                    bufferHO_Ram_0K_tot.append(line)
                if HO_Ram_T_tot:
                    bufferHO_Ram_T_tot.append(line)
                if HO_Ram_0K_comp:
                    bufferHO_Ram_0K_comp.append(line)
                if HO_Ram_T_comp:
                    bufferHO_Ram_T_comp.append(line)
                if VSCF_IR:
                    bufferVSCF_IR.append(line)
                if VSCF_Ram_0K_tot:
                    bufferVSCF_Ram_0K_tot.append(line)
                if VSCF_Ram_T_tot:
                    bufferVSCF_Ram_T_tot.append(line)
                if VSCF_Ram_0K_comp:
                    bufferVSCF_Ram_0K_comp.append(line)
                if VSCF_Ram_T_comp:
                    bufferVSCF_Ram_T_comp.append(line)
                if VCI_IR:
                    bufferVCI_IR.append(line)
                if VCI_Ram_0K_tot:
                    bufferVCI_Ram_0K_tot.append(line)
                if VCI_Ram_T_tot:
                    bufferVCI_Ram_T_tot.append(line)
                if VCI_Ram_0K_comp:
                    bufferVCI_Ram_0K_comp.append(line)
                if VCI_Ram_T_comp:
                    bufferVCI_Ram_T_comp.append(line)
                # anscan
                if anscan_IR:
                    bufferanscan_IR.append([])
                    idx = len(imode_IR) - 1
                    bufferanscan_IR[idx].append(line)
                if anscan_Ram_0K_tot:
                    bufferanscan_Ram_0K_tot.append([])
                    idx = len(imode_Ram) - 1
                    bufferanscan_Ram_0K_tot[idx].append(line)
                if anscan_Ram_T_tot:
                    bufferanscan_Ram_T_tot.append([])
                    idx = len(imode_Ram) - 1
                    bufferanscan_Ram_T_tot[idx].append(line)
                if anscan_Ram_0K_comp:
                    bufferanscan_Ram_0K_comp.append([])
                    idx = len(imode_Ram) - 1
                    bufferanscan_Ram_0K_comp[idx].append(line)
                if anscan_Ram_T_comp:
                    bufferanscan_Ram_T_comp.append([])
                    idx = len(imode_Ram) - 1
                    bufferanscan_Ram_T_comp[idx].append(line)
                # anscan

        # <-- Big loop over lines of CRYSTAL output file

        # debug
        # for i in range(len(imode_IR)):
        #    print(bufferanscan_IR[i])
        # for i in range(len(imode_Ram)):
        #    print(bufferanscan_Ram_0K_tot[i])
        # debug

        # Save and parse VSCF data for IR spectrum
        n_VSCF_ir = len(bufferVSCF_IR)
        if n_VSCF_ir > 0:
            for i, line in enumerate(bufferVSCF_IR[5:n_VSCF_ir-1]):
                IR_VSCF.append(line.split()[3:6])
                for j in range(3):
                    IR_VSCF[i][j] = float(IR_VSCF[i][j])
            IR_VSCF = np.array(IR_VSCF)
            self.IR_VSCF_T = IR_VSCF[:, 0:3:2]
            self.IR_VSCF_0K = IR_VSCF[:, 0:2]

        # Save and parse VCI data for IR spectrum
        n_VCI_ir = len(bufferVCI_IR)
        if n_VCI_ir > 0:
            for i, line in enumerate(bufferVCI_IR[5:n_VCI_ir-1]):
                IR_VCI.append(line.split()[3:6])
                for j in range(3):
                    IR_VCI[i][j] = float(IR_VCI[i][j])
            IR_VCI = np.array(IR_VCI)
            self.IR_VCI_T = IR_VCI[:, 0:3:2]
            self.IR_VCI_0K = IR_VCI[:, 0:2]

        # Save and parse HO data for IR spectrum
        n_HO_ir = len(bufferHO_IR)
        if n_HO_ir > 0:
            for i, line in enumerate(bufferHO_IR[5:n_HO_ir-1]):
                IR_HO.append(line.split()[3:6])
                for j in range(3):
                    IR_HO[i][j] = float(IR_HO[i][j])
            IR_HO = np.array(IR_HO)
            self.IR_HO_T = IR_HO[:, 0:3:2]
            self.IR_HO_0K = IR_HO[:, 0:2]

        # anscan
        self.IR_anscan_T = []
        self.IR_anscan_0K = []
        # Save and parse anscan data for IR spectrum
        for k in range(len(imode_IR)):
            n_anscan_IR = len(bufferanscan_IR[k])
            if n_anscan_IR > 0:
                IR_anscan = []
                for i, line in enumerate(bufferanscan_IR[k][5:n_anscan_IR-1]):
                    IR_anscan.append(line.split()[3:6])
                    for j in range(3):
                        IR_anscan[i][j] = float(IR_anscan[i][j])
                IR_anscan = np.array(IR_anscan)
                self.IR_anscan_T.append(IR_anscan[:, 0:3:2])
                self.IR_anscan_0K.append(IR_anscan[:, 0:2])
            # anscan

       # Save and parse HO data for Raman spectrum (0K, tot)
        n_HO_Ram_0K_tot = len(bufferHO_Ram_0K_tot)
        if n_HO_Ram_0K_tot > 0:
            for i, line in enumerate(bufferHO_Ram_0K_tot[6:n_HO_Ram_0K_tot-1]):
                Ram_HO_0K_tot.append(line.split()[3:7])
                for j in range(4):
                    Ram_HO_0K_tot[i][j] = float(Ram_HO_0K_tot[i][j])
            Ram_HO_0K_tot = np.array(Ram_HO_0K_tot)
            self.Ram_HO_0K_tot = Ram_HO_0K_tot[:, 0:2]
            self.Ram_HO_0K_per = Ram_HO_0K_tot[:, 0:3:2]
            self.Ram_HO_0K_par = Ram_HO_0K_tot[:, 0:4:3]

       # Save and parse HO data for Raman spectrum (T, tot)
        n_HO_Ram_T_tot = len(bufferHO_Ram_T_tot)
        if n_HO_Ram_T_tot > 0:
            for i, line in enumerate(bufferHO_Ram_T_tot[6:n_HO_Ram_T_tot-1]):
                Ram_HO_T_tot.append(line.split()[3:7])
                for j in range(4):
                    Ram_HO_T_tot[i][j] = float(Ram_HO_T_tot[i][j])
            Ram_HO_T_tot = np.array(Ram_HO_T_tot)
            self.Ram_HO_T_tot = Ram_HO_T_tot[:, 0:2]
            self.Ram_HO_T_per = Ram_HO_T_tot[:, 0:3:2]
            self.Ram_HO_T_par = Ram_HO_T_tot[:, 0:4:3]

       # Save and parse HO data for Raman spectrum (0K, comp)
        n_HO_Ram_0K_comp = len(bufferHO_Ram_0K_comp)
        if n_HO_Ram_0K_comp > 0:
            for i, line in enumerate(bufferHO_Ram_0K_comp[6:n_HO_Ram_0K_comp-1]):
                Ram_HO_0K_comp.append(line.split()[3:10])
                for j in range(7):
                    Ram_HO_0K_comp[i][j] = float(Ram_HO_0K_comp[i][j])
            Ram_HO_0K_comp = np.array(Ram_HO_0K_comp)
            self.Ram_HO_0K_comp_xx = Ram_HO_0K_comp[:, 0:2]
            self.Ram_HO_0K_comp_xy = Ram_HO_0K_comp[:, 0:3:2]
            self.Ram_HO_0K_comp_xz = Ram_HO_0K_comp[:, 0:4:3]
            self.Ram_HO_0K_comp_yy = Ram_HO_0K_comp[:, 0:5:4]
            self.Ram_HO_0K_comp_yz = Ram_HO_0K_comp[:, 0:6:5]
            self.Ram_HO_0K_comp_zz = Ram_HO_0K_comp[:, 0:7:6]

       # Save and parse HO data for Raman spectrum (T, comp)
        n_HO_Ram_T_comp = len(bufferHO_Ram_T_comp)
        if n_HO_Ram_T_comp > 0:
            for i, line in enumerate(bufferHO_Ram_T_comp[6:n_HO_Ram_T_comp-1]):
                Ram_HO_T_comp.append(line.split()[3:10])
                for j in range(7):
                    Ram_HO_T_comp[i][j] = float(Ram_HO_T_comp[i][j])
            Ram_HO_T_comp = np.array(Ram_HO_T_comp)
            self.Ram_HO_T_comp_xx = Ram_HO_T_comp[:, 0:2]
            self.Ram_HO_T_comp_xy = Ram_HO_T_comp[:, 0:3:2]
            self.Ram_HO_T_comp_xz = Ram_HO_T_comp[:, 0:4:3]
            self.Ram_HO_T_comp_yy = Ram_HO_T_comp[:, 0:5:4]
            self.Ram_HO_T_comp_yz = Ram_HO_T_comp[:, 0:6:5]
            self.Ram_HO_T_comp_zz = Ram_HO_T_comp[:, 0:7:6]

       # Save and parse VSCF data for Raman spectrum (0K, tot)
        n_VSCF_Ram_0K_tot = len(bufferVSCF_Ram_0K_tot)
        if n_VSCF_Ram_0K_tot > 0:
            for i, line in enumerate(bufferVSCF_Ram_0K_tot[6:n_VSCF_Ram_0K_tot-1]):
                Ram_VSCF_0K_tot.append(line.split()[3:7])
                for j in range(4):
                    Ram_VSCF_0K_tot[i][j] = float(Ram_VSCF_0K_tot[i][j])
            Ram_VSCF_0K_tot = np.array(Ram_VSCF_0K_tot)
            self.Ram_VSCF_0K_tot = Ram_VSCF_0K_tot[:, 0:2]
            self.Ram_VSCF_0K_per = Ram_VSCF_0K_tot[:, 0:3:2]
            self.Ram_VSCF_0K_par = Ram_VSCF_0K_tot[:, 0:4:3]

       # Save and parse VSCF data for Raman spectrum (T, tot)
        n_VSCF_Ram_T_tot = len(bufferVSCF_Ram_T_tot)
        if n_VSCF_Ram_T_tot > 0:
            for i, line in enumerate(bufferVSCF_Ram_T_tot[6:n_VSCF_Ram_T_tot-1]):
                Ram_VSCF_T_tot.append(line.split()[3:7])
                for j in range(4):
                    Ram_VSCF_T_tot[i][j] = float(Ram_VSCF_T_tot[i][j])
            Ram_VSCF_T_tot = np.array(Ram_VSCF_T_tot)
            self.Ram_VSCF_T_tot = Ram_VSCF_T_tot[:, 0:2]
            self.Ram_VSCF_T_per = Ram_VSCF_T_tot[:, 0:3:2]
            self.Ram_VSCF_T_par = Ram_VSCF_T_tot[:, 0:4:3]

       # Save and parse VSCF data for Raman spectrum (0K, comp)
        n_VSCF_Ram_0K_comp = len(bufferVSCF_Ram_0K_comp)
        if n_VSCF_Ram_0K_comp > 0:
            for i, line in enumerate(bufferVSCF_Ram_0K_comp[6:n_VSCF_Ram_0K_comp-1]):
                Ram_VSCF_0K_comp.append(line.split()[3:10])
                for j in range(7):
                    Ram_VSCF_0K_comp[i][j] = float(Ram_VSCF_0K_comp[i][j])
            Ram_VSCF_0K_comp = np.array(Ram_VSCF_0K_comp)
            self.Ram_VSCF_0K_comp_xx = Ram_VSCF_0K_comp[:, 0:2]
            self.Ram_VSCF_0K_comp_xy = Ram_VSCF_0K_comp[:, 0:3:2]
            self.Ram_VSCF_0K_comp_xz = Ram_VSCF_0K_comp[:, 0:4:3]
            self.Ram_VSCF_0K_comp_yy = Ram_VSCF_0K_comp[:, 0:5:4]
            self.Ram_VSCF_0K_comp_yz = Ram_VSCF_0K_comp[:, 0:6:5]
            self.Ram_VSCF_0K_comp_zz = Ram_VSCF_0K_comp[:, 0:7:6]

       # Save and parse VSCF data for Raman spectrum (T, comp)
        n_VSCF_Ram_T_comp = len(bufferVSCF_Ram_T_comp)
        if n_VSCF_Ram_T_comp > 0:
            for i, line in enumerate(bufferVSCF_Ram_T_comp[6:n_VSCF_Ram_T_comp-1]):
                Ram_VSCF_T_comp.append(line.split()[3:10])
                for j in range(7):
                    Ram_VSCF_T_comp[i][j] = float(Ram_VSCF_T_comp[i][j])
            Ram_VSCF_T_comp = np.array(Ram_VSCF_T_comp)
            self.Ram_VSCF_T_comp_xx = Ram_VSCF_T_comp[:, 0:2]
            self.Ram_VSCF_T_comp_xy = Ram_VSCF_T_comp[:, 0:3:2]
            self.Ram_VSCF_T_comp_xz = Ram_VSCF_T_comp[:, 0:4:3]
            self.Ram_VSCF_T_comp_yy = Ram_VSCF_T_comp[:, 0:5:4]
            self.Ram_VSCF_T_comp_yz = Ram_VSCF_T_comp[:, 0:6:5]
            self.Ram_VSCF_T_comp_zz = Ram_VSCF_T_comp[:, 0:7:6]

       # Save and parse VCI data for Raman spectrum (0K, tot)
        n_VCI_Ram_0K_tot = len(bufferVCI_Ram_0K_tot)
        if n_VCI_Ram_0K_tot > 0:
            for i, line in enumerate(bufferVCI_Ram_0K_tot[6:n_VCI_Ram_0K_tot-1]):
                Ram_VCI_0K_tot.append(line.split()[3:7])
                for j in range(4):
                    Ram_VCI_0K_tot[i][j] = float(Ram_VCI_0K_tot[i][j])
            Ram_VCI_0K_tot = np.array(Ram_VCI_0K_tot)
            self.Ram_VCI_0K_tot = Ram_VCI_0K_tot[:, 0:2]
            self.Ram_VCI_0K_per = Ram_VCI_0K_tot[:, 0:3:2]
            self.Ram_VCI_0K_par = Ram_VCI_0K_tot[:, 0:4:3]

       # Save and parse VCI data for Raman spectrum (T, tot)
        n_VCI_Ram_T_tot = len(bufferVCI_Ram_T_tot)
        if n_VCI_Ram_T_tot > 0:
            for i, line in enumerate(bufferVCI_Ram_T_tot[6:n_VCI_Ram_T_tot-1]):
                Ram_VCI_T_tot.append(line.split()[3:7])
                for j in range(4):
                    Ram_VCI_T_tot[i][j] = float(Ram_VCI_T_tot[i][j])
            Ram_VCI_T_tot = np.array(Ram_VCI_T_tot)
            self.Ram_VCI_T_tot = Ram_VCI_T_tot[:, 0:2]
            self.Ram_VCI_T_per = Ram_VCI_T_tot[:, 0:3:2]
            self.Ram_VCI_T_par = Ram_VCI_T_tot[:, 0:4:3]

       # Save and parse VCI data for Raman spectrum (0K, comp)
        n_VCI_Ram_0K_comp = len(bufferVCI_Ram_0K_comp)
        if n_VCI_Ram_0K_comp > 0:
            for i, line in enumerate(bufferVCI_Ram_0K_comp[6:n_VCI_Ram_0K_comp-1]):
                Ram_VCI_0K_comp.append(line.split()[3:10])
                for j in range(7):
                    Ram_VCI_0K_comp[i][j] = float(Ram_VCI_0K_comp[i][j])
            Ram_VCI_0K_comp = np.array(Ram_VCI_0K_comp)
            self.Ram_VCI_0K_comp_xx = Ram_VCI_0K_comp[:, 0:2]
            self.Ram_VCI_0K_comp_xy = Ram_VCI_0K_comp[:, 0:3:2]
            self.Ram_VCI_0K_comp_xz = Ram_VCI_0K_comp[:, 0:4:3]
            self.Ram_VCI_0K_comp_yy = Ram_VCI_0K_comp[:, 0:5:4]
            self.Ram_VCI_0K_comp_yz = Ram_VCI_0K_comp[:, 0:6:5]
            self.Ram_VCI_0K_comp_zz = Ram_VCI_0K_comp[:, 0:7:6]

       # Save and parse VCI data for Raman spectrum (T, comp)
        n_VCI_Ram_T_comp = len(bufferVCI_Ram_T_comp)
        if n_VCI_Ram_T_comp > 0:
            for i, line in enumerate(bufferVCI_Ram_T_comp[6:n_VCI_Ram_T_comp-1]):
                Ram_VCI_T_comp.append(line.split()[3:10])
                for j in range(7):
                    Ram_VCI_T_comp[i][j] = float(Ram_VCI_T_comp[i][j])
            Ram_VCI_T_comp = np.array(Ram_VCI_T_comp)
            self.Ram_VCI_T_comp_xx = Ram_VCI_T_comp[:, 0:2]
            self.Ram_VCI_T_comp_xy = Ram_VCI_T_comp[:, 0:3:2]
            self.Ram_VCI_T_comp_xz = Ram_VCI_T_comp[:, 0:4:3]
            self.Ram_VCI_T_comp_yy = Ram_VCI_T_comp[:, 0:5:4]
            self.Ram_VCI_T_comp_yz = Ram_VCI_T_comp[:, 0:6:5]
            self.Ram_VCI_T_comp_zz = Ram_VCI_T_comp[:, 0:7:6]

        # anscan
        # Save and parse anscan data for Raman spectrum (0K, tot)
        self.Ram_anscan_0K_tot = []
        self.Ram_anscan_0K_per = []
        self.Ram_anscan_0K_par = []
        for k in range(len(imode_Ram)):
            n_anscan_Ram_0K_tot = len(bufferanscan_Ram_0K_tot[k])
            if n_anscan_Ram_0K_tot > 0:
                Ram_anscan = []
                for i, line in enumerate(bufferanscan_Ram_0K_tot[k][6:n_anscan_Ram_0K_tot-1]):
                    Ram_anscan.append(line.split()[3:7])
                    for j in range(4):
                        Ram_anscan[i][j] = float(Ram_anscan[i][j])
                Ram_anscan = np.array(Ram_anscan)
                self.Ram_anscan_0K_tot.append(Ram_anscan[:, 0:2])
                self.Ram_anscan_0K_per.append(Ram_anscan[:, 0:3:2])
                self.Ram_anscan_0K_par.append(Ram_anscan[:, 0:4:3])

        # Save and parse anscan data for Raman spectrum (T, tot)
        self.Ram_anscan_T_tot = []
        self.Ram_anscan_T_per = []
        self.Ram_anscan_T_par = []
        for k in range(len(imode_Ram)):
            n_anscan_Ram_T_tot = len(bufferanscan_Ram_T_tot[k])
            if n_anscan_Ram_T_tot > 0:
                Ram_anscan = []
                for i, line in enumerate(bufferanscan_Ram_T_tot[k][6:n_anscan_Ram_T_tot-1]):
                    Ram_anscan.append(line.split()[3:7])
                    for j in range(4):
                        Ram_anscan[i][j] = float(Ram_anscan[i][j])
                Ram_anscan = np.array(Ram_anscan)
                self.Ram_anscan_T_tot.append(Ram_anscan[:, 0:2])
                self.Ram_anscan_T_per.append(Ram_anscan[:, 0:3:2])
                self.Ram_anscan_T_par.append(Ram_anscan[:, 0:4:3])

        # Save and parse anscan data for Raman spectrum (0K, comp)
        self.Ram_anscan_0K_xx = []
        self.Ram_anscan_0K_xy = []
        self.Ram_anscan_0K_xz = []
        self.Ram_anscan_0K_yy = []
        self.Ram_anscan_0K_yz = []
        self.Ram_anscan_0K_zz = []
        for k in range(len(imode_Ram)):
            n_anscan_Ram_0K_comp = len(bufferanscan_Ram_0K_comp[k])
            if n_anscan_Ram_0K_comp > 0:
                Ram_anscan = []
                for i, line in enumerate(bufferanscan_Ram_0K_comp[k][6:n_anscan_Ram_0K_comp-1]):
                    Ram_anscan.append(line.split()[3:10])
                    for j in range(7):
                        Ram_anscan[i][j] = float(Ram_anscan[i][j])
                Ram_anscan = np.array(Ram_anscan)
                self.Ram_anscan_0K_xx.append(Ram_anscan[:, 0:2])
                self.Ram_anscan_0K_xy.append(Ram_anscan[:, 0:3:2])
                self.Ram_anscan_0K_xz.append(Ram_anscan[:, 0:4:3])
                self.Ram_anscan_0K_yy.append(Ram_anscan[:, 0:5:4])
                self.Ram_anscan_0K_yz.append(Ram_anscan[:, 0:6:5])
                self.Ram_anscan_0K_zz.append(Ram_anscan[:, 0:7:6])

        # Save and parse anscan data for Raman spectrum (T, comp)
        self.Ram_anscan_T_xx = []
        self.Ram_anscan_T_xy = []
        self.Ram_anscan_T_xz = []
        self.Ram_anscan_T_yy = []
        self.Ram_anscan_T_yz = []
        self.Ram_anscan_T_zz = []
        for k in range(len(imode_Ram)):
            n_anscan_Ram_T_comp = len(bufferanscan_Ram_T_comp[k])
            if n_anscan_Ram_T_comp > 0:
                Ram_anscan = []
                for i, line in enumerate(bufferanscan_Ram_T_comp[k][6:n_anscan_Ram_T_comp-1]):
                    Ram_anscan.append(line.split()[3:10])
                    for j in range(7):
                        Ram_anscan[i][j] = float(Ram_anscan[i][j])
                Ram_anscan = np.array(Ram_anscan)
                self.Ram_anscan_T_xx.append(Ram_anscan[:, 0:2])
                self.Ram_anscan_T_xy.append(Ram_anscan[:, 0:3:2])
                self.Ram_anscan_T_xz.append(Ram_anscan[:, 0:4:3])
                self.Ram_anscan_T_yy.append(Ram_anscan[:, 0:5:4])
                self.Ram_anscan_T_yz.append(Ram_anscan[:, 0:6:5])
                self.Ram_anscan_T_zz.append(Ram_anscan[:, 0:7:6])
        # anscan

        return self

    # ANSCAN+DWELL
    def get_anscan(self, anscanwf):
        """
        Work in progress for ANSCAN stuff (in development).
        """

        import re

        import numpy as np

        ndispl = 0
        storeE = False
        storeFC = False
        self.energy = []
        self.force_const = []

        for i, line in enumerate(self.data):
            if re.match(r'\s*SCAN ALONG NORMAL MODES', line):
                ndispl = (int(self.data[i+1].split()[5]) -
                          int(self.data[i+1].split()[2]))
                if (ndispl < 1):
                    raise Exception('Number of displacements too low.')
                # Save range
                first = float(self.data[i+1].split()[2])
                last = float(self.data[i+1].split()[5])
                step = float(self.data[i+1].split()[7])
                self.rangescan = [first*step, last*step]
            if re.match(r'\s*\[DISPLAC\]\s*\[\s*SCAN POTENTIAL\s*\]', line):
                V = np.zeros([ndispl+1, 2])
                for k in range(ndispl+1):
                    if (re.match(r'^\s+0.0000*', self.data[i+2+k])):
                        continue
                    V[k, 0] = float(self.data[i+2+k].split()[2])
                    V[k, 1] = float(self.data[i+2+k].split()[4])
                self.anhpot = V[:, 0]
                self.harmpot = V[:, 1]
                # Save index of anscan mode
                strtmp = self.data[i-1].split()[1]
                strtmp = strtmp[:strtmp.find('(')]
                anhmode = int(strtmp)
            if re.match(r'\s*ANHARMONIC VIBRATIONAL STATES', self.data[i-3]):
                storeE = True
            if re.match(r'\s*POTENTIAL ENERGY DERIVATIVES', line):
                storeE = False
            if (storeE and (len(self.data[i].split()) != 0)):
                self.energy.append(float(self.data[i].split()[2]))
            if re.match(r'\s*POTENTIAL ENERGY DERIVATIVES', self.data[i-3]):
                storeFC = True
            if (storeFC and (re.match(r'^(?!\s*\d)', line))):
                storeFC = False
            if (storeFC):
                self.force_const.append(float(line.split()[2]))

        # Call get_phonon() method and save harmonic freq
        self.get_phonon(rm_imaginary=False)
        self.harm_freq = self.frequency[0, anhmode-1]*33.333333

        # Read ANSCANWF.DAT
        try:
            file = open(anscanwf, 'r', errors='ignore')
            data = file.readlines()
            file.close()
        except:
            raise FileNotFoundError(
                'EXITING: a .anscanwf file needs to be specified')

        self.wf = np.zeros([len(self.energy), 10])
        self.alpha = None

        for i, line in enumerate(data):
            if re.match(r'^\s*$', line):
                break
            if re.match(r'.ALP.*', line):
                self.alpha = float(line.split()[2])
                break
            if re.match(r'.ANH.*', line):
                continue
            for j in range(10):
                self.wf[i-1, j] = float(line.split()[j])

    # ANSCAN+DWELL

    def get_elatensor(self):
        """
        Extracts the elastic tensor from CRYSTAL output.

        Returns:
            self.elatensor(np.array): Symmetrized 6x6 elastic tensor.
        """

        import numpy as np

        startstring = " SYMMETRIZED ELASTIC"
        stopstring = " ELASTIC MODULI"
        self.elatensor = []
        buffer = []
        strtensor = []
        copy = False

        # Search for elastic tensor and save it into buffer
        for line in self.data:
            if line.startswith(startstring):
                copy = True
            elif line.startswith(stopstring):
                copy = False
            elif copy:
                buffer.append(line)

        # Build tensor
        for i in range(6):
            # Clean buffer and copy it in strtensor
            strtensor.append(
                buffer[i + 1].replace(" |", " ").replace("\n", ""))
            # Split strtensor strings and copy them in tensor
            self.elatensor.append(strtensor[i].split())
            # Conversion str -> float
            for j in range(6 - i):
                self.elatensor[i][j] = float(self.elatensor[i][j])
            # Add zeros
            for k in range(i):
                self.elatensor[i].insert(0, 0)
        buffer.clear()

        # Symmetrize tensor
        for i in range(6):
            for j in range(6):
                self.elatensor[j][i] = self.elatensor[i][j]

        self.elatensor = np.array(self.elatensor)

        return self.elatensor

    def get_EOS(self):
        """
        Extracts the Equation of state output data and corresponding fittings

        Returns:
            self (Crystal_output): New attributes listed below.
            self.VvsE (np.array): numpy array containing the computed volumes and
                                  the corresponding energies in the first and second column respectively
            self.murnaghan (np.array): numpy array containing the fitted thermodynamics
                                       functions with the Murnaghan Equation
            self.bmurnaghan (np.array): numpy array containing the fitted thermodynamics
                                        functions with the Birch-Murnaghan Equation
            self.pt (np.array): numpy array containing the fitted thermodynamics
                                functions with the Poirier-Tarantola Equation
            self.vinet (np.array): numpy array containing the fitted thermodynamics
                                   functions with the Vinet Equation
        """

        import re
        import numpy as np

        # Definition of the delimiters regex -->
        volume_start = r'\s+SORTING VOLUMES/ENERGIES'
        volume_end = r'\s+\++\sFITTING\sUSING\sALL\sPOINTS\s\++\s'
        murnaghan_start = r'\s+THERMODYNAMIC\s\w+\s\w+\s\w+\s\d\sK\s\w+\sEOS:\sMURNAGHAN\s\d+'
        bmurnaghan_start = r'\s+THERMODYNAMIC\s\w+\s\w+\s\w+\s\d\sK\s\w+\sEOS:\sBIRCH-MURNAGHAN\s\d+'
        poiriertarantola_start = r'\s+THERMODYNAMIC\s\w+\s\w+\s\w+\s\d\sK\s\w+\sEOS:\sPOIRIER-TARANTOLA\s\d+'
        vinet_start = r'\s+THERMODYNAMIC\s\w+\s\w+\s\w+\s\d\sK\s\w+\sEOS:\sVINET\s\d+'
        eos_end = r' EOSEOSEOSEOSEOSEOSEOSEOSEOSEOSEOSEOSEOSEOSEOSEOSEOSEOSEOSEOSEOSEOSEOSEOSEOSEOS'
        eos_end_logic = False
        # <--

        # Search of delimiters index -->
        for index, line in enumerate(self.data):
            if re.match(volume_start, line):
                volume_start_idx = index
                eos_end_logic = True

            if re.match(volume_end, line):
                volume_end_idx = index

            if re.match(murnaghan_start, line):
                murnaghan_start_idx = index

            if re.match(bmurnaghan_start, line):
                bmurnaghan_start_idx = index

            if re.match(poiriertarantola_start, line):
                poiriertarantola_start_idx = index

            if re.match(vinet_start, line):
                vinet_start_idx = index

            if (re.match(eos_end, line) and eos_end_logic):
                eos_end_idx = index
        # <--

        # Creation of the numpy array for volume and energy -->
        no_volumes = len(self.data[volume_start_idx+4:volume_end_idx-1])
        self.VvsE = np.ones((no_volumes, 2))
        for i, line in enumerate(self.data[(volume_start_idx + 4): (volume_end_idx - 1)]):
            lsplit = line.split()
            self.VvsE[i, 0] = float(lsplit[0])
            self.VvsE[i, 1] = float(lsplit[1])
        # <--

        # Creation of the Murnaghan numpy array -->
        no_points = len(
            self.data[murnaghan_start_idx + 7: bmurnaghan_start_idx - 3])
        self.murnaghan = np.ones((no_points, 5))
        for i, line in enumerate(self.data[murnaghan_start_idx + 7: bmurnaghan_start_idx - 3]):
            lsplit = line.split()
            self.murnaghan[i, 0] = float(lsplit[0])
            self.murnaghan[i, 1] = float(lsplit[1])
            self.murnaghan[i, 2] = float(lsplit[2])
            self.murnaghan[i, 3] = float(lsplit[3])
            self.murnaghan[i, 4] = float(lsplit[4])
        # <--

        # Creation of Birch-Murnaghan numpy array -->
        no_points = len(
            self.data[bmurnaghan_start_idx + 7: poiriertarantola_start_idx - 3])
        self.bmurnaghan = np.ones((no_points, 5))
        for i, line in enumerate(self.data[bmurnaghan_start_idx + 7: poiriertarantola_start_idx - 3]):
            lsplit = line.split()
            self.bmurnaghan[i, 0] = float(lsplit[0])
            self.bmurnaghan[i, 1] = float(lsplit[1])
            self.bmurnaghan[i, 2] = float(lsplit[2])
            self.bmurnaghan[i, 3] = float(lsplit[3])
            self.bmurnaghan[i, 4] = float(lsplit[4])
        # <--

        # Creation of Poirier-Tarantola numpy array -->
        no_points = len(
            self.data[poiriertarantola_start_idx + 7: vinet_start_idx - 3])
        self.pt = np.ones((no_points, 5))
        for i, line in enumerate(self.data[poiriertarantola_start_idx + 7: vinet_start_idx - 3]):
            lsplit = line.split()
            self.pt[i, 0] = float(lsplit[0])
            self.pt[i, 1] = float(lsplit[1])
            self.pt[i, 2] = float(lsplit[2])
            self.pt[i, 3] = float(lsplit[3])
            self.pt[i, 4] = float(lsplit[4])
        # <--

        # Creation of Vinet numpy array -->
        no_points = len(
            self.data[vinet_start_idx + 7: eos_end_idx - 1])
        self.vinet = np.ones((no_points, 5))
        for i, line in enumerate(self.data[vinet_start_idx + 7: eos_end_idx - 1]):
            lsplit = line.split()
            self.vinet[i, 0] = float(lsplit[0])
            self.vinet[i, 1] = float(lsplit[1])
            self.vinet[i, 2] = float(lsplit[2])
            self.vinet[i, 3] = float(lsplit[3])
            self.vinet[i, 4] = float(lsplit[4])
        # <--

        return self


class Properties_input:
    """Create a properties_input object"""

    def __init__(self, input_name=None):
        """Initialise the object"""

        self.is_newk = False

    def from_file(self, input_name):
        """
        Read the properties input from a file.

        Args:
            input_name (str): Name of the input file.
        Returns:
            self (Properties_input): Properties_input object.
        """
        import sys

        self.name = input_name
        if input_name is not None:
            try:
                if input_name[-3:] != 'd12':
                    input_name = input_name+'.d12'
                file = open(input_name, 'r')
                self.data = file.readlines()
                file.close()
            except:
                print('EXITING: a .d3 file needs to be specified')
                sys.exit(1)

            # Check if NEWK is in the input
            if 'NEWK\n' in self.data:
                self.is_newk = True
                self.newk_block = self.data[0:2]
                self.property_block = self.data[2:]
            else:
                self.is_newk = False
                self.property_block = self.data

        return self

    def make_newk_block(self, shrink1, shrink2, Fermi=1, print_option=0):
        """
        Returns the newk block.

        Args:
            shrink1 (int): First newk shrinking factor.
            shrink2 (int): Second newk shrinking factor.
            Fermi (int): Fermi recalculation option (default is 1).
            print_option (int): Properties printing option (default is 0).
        """

        self.is_newk = True

        self.newk_block = ['NEWK\n', '%s %s\n' % (shrink1, shrink2),
                           '%s %s\n' % (Fermi, print_option)]

    def make_bands_block(self, k_path, n_kpoints, first_band, last_band, print_eig=0, print_option=1,
                         precision=5, title='BAND STRUCTURE CALCULATION'):
        """
        Returns the bands block for a bands calculation.

        Args:
            k_path (list or HighSymmKpath): k-path for bands calculation.
            n_kpoints (int): Number of k-points along the path.
            first_band (int): Index of the first band.
            last_band (int): Index of the last band.
            print_eig (int): Printing options for eigenvalues (default is 0).
            print_option (int): Properties printing options (default is 1).
            precision (int): Number of zeros in the calculation of the gcd
            title (str): Title of the calculation (default is "BAND STRUCTURE CALCULATION").
        """

        import sys
        import numpy as np

        bands_block = []

        # path from a pymatgen k_path object
        if 'HighSymmKpath' in str(type(k_path)):
            k_path_flat = [item for sublist in k_path.kpath['path']
                           for item in sublist]
            k_path_pmg = []
            for i in k_path_flat:
                # This is a pmg HighSymmKpath object
                k_path_pmg.append(k_path.kpath['kpoints'][i].tolist())
            k_path = np.array(k_path_pmg)

        elif type(k_path[0]) == list:
            # This is a list of lists
            k_path = np.array(k_path)

        else:
            print('EXITING: k_path type must be a list of list (k coordinates) or\
                a pymatgen HighSymmKpath object. %s selected' % type(k_path))
            sys.exit(1)

        k_unique = np.unique(k_path)

        # Find the shrinking factor
        k_unique = np.array(np.around(k_unique, precision)
                            * 10**precision, dtype=int)
        if len(k_unique) > 2:
            gcd = np.gcd.reduce(k_unique)
        else:
            gcd = np.gcd(k_unique[0], k_unique[1])
        k_path = np.array((k_path/gcd)*10**precision, dtype=int)
        shrink = int(10**precision/gcd)

        bands_block.append('BAND\n')
        bands_block.append(title+'\n')

        bands_block.append(str(len(k_path)-1)+' '+str(shrink)+' '+str(n_kpoints) +
                           ' '+str(first_band)+' '+str(last_band)+' ' +
                           str(print_option)+' '+str(print_eig)+'\n')

        # Add the symmetry line
        for i in range(len(k_path[:-1])):
            bands_block.append(' '.join([str(x) for x in k_path[i]])+'  ' +
                               ' '.join([str(x) for x in k_path[i+1]])+'\n')

        bands_block.append('END\n')

        self.property_block = bands_block

        return self

    def make_doss_block(self, n_points=200, band_range=None, e_range=None, plotting_option=2,
                        poly=12, print_option=1):
        """
        Returns the doss block for a doss calculation.

        Args:
            n_points (int): Number of points in the DOS plot (default is 200).
            band_range (list or tuple): Range of bands to include in the DOS calculation (default is None).
            e_range (list or tuple): Energy range for the DOS calculation (default is None).
            plotting_option (int): DOS plotting options (default is 2).
            poly (int): Degree of the polynomial for smearing (default is 12).
            print_option (int): Properties printing options (default is 1).
        """

        import sys

        # either band_range or e_range needs to be specified
        doss_block = []
        if band_range == None and e_range == None:
            print('EXITING: please specify either band_range or e_range. None selected')
            sys.exit(1)
        elif band_range != None and e_range != None:
            print('EXITING: please specify either band_range or e_range. Both selected')
            sys.exit(1)
        elif type(band_range) == list and len(band_range) == 2:
            doss_range = band_range
        elif type(e_range) == list and len(e_range) == 2:
            doss_range = [-1, -1]

        else:
            print('EXITING: either the band_range argument or the e_range argument\
                do not match the required format (2 item list)')
            sys.exit(1)

        doss_block.append('DOSS\n')
        doss_block.append(str(0)+' '+str(n_points)+' '+str(doss_range[0])+' ' +
                          str(doss_range[1])+' '+str(plotting_option)+' '+str(poly)+' ' +
                          str(print_option)+'\n')

        if doss_range == [-1, -1]:
            doss_block.append(
                str(units.eV_to_H(e_range[0]))+' '+str(units.eV_to_H(e_range[1]))+'\n')

        doss_block.append('END\n')

        self.property_block = doss_block

        return self

    def make_pdoss_block(self, projections, proj_type='atom', output_file=None, n_points=200, band_range=None,
                         e_range=None, plotting_option=2, poly=12, print_option=1):
        """
        Returns the pdoss block for a pdoss calculation.

        Args:
            projections (dict): Dictionary specifying the projections for the pdoss calculation.
            proj_type (str): Type of projection ("atom" or "site", default is "atom").
            output_file (str): Output file name (default is None).
            n_points (int): Number of points in the DOS plot (default is 200).
            band_range (list or tuple): Range of bands to include in the DOS calculation (default is None).
            e_range (list or tuple): Energy range for the DOS calculation (default is None).
            plotting_option (int): DOS plotting options (default is 2).
            poly (int): Degree of the polynomial for smearing (default is 12).
            print_option (int): Properties printing options (default is 1).
        """

        import sys

        pdoss_block = []
        if band_range == None and e_range == None:
            print('EXITING: please specify either band_range or e_range. None selected')
            sys.exit(1)
        elif band_range != None and e_range != None:
            print('EXITING: please specify either band_range or e_range. Both selected')
            sys.exit(1)
        elif type(band_range) == list and len(band_range) == 2:
            pdoss_range = band_range
            range_is_bands = True
        elif type(e_range) == list and len(e_range) == 2:
            pdoss_range = [-1, -1]
            range_is_bands = False

        else:
            print('EXITING: either the band_range argument or the e_range argument\
                do not match the required format (2 item list)')
            sys.exit(1)

        pdoss_block.append('DOSS\n')
        pdoss_block.append(str(len(projections))+' '+str(n_points)+' '+str(pdoss_range[0])+' ' +
                           str(pdoss_range[1])+' '+str(plotting_option)+' '+str(poly)+' ' +
                           str(print_option)+'\n')

        if range_is_bands == False:
            pdoss_block.append(
                str(round(units.eV_to_H(e_range[0]), 6))+' '+str(round(units.eV_to_H(e_range[1]), 6))+'\n')

        flat_proj = [x for sublist in projections for x in sublist]

        if all(isinstance(x, int) for x in flat_proj):
            if proj_type == 'atom':
                for proj in projections:
                    pdoss_block.append(str(-len(proj))+' ' +
                                       ' '.join([str(x) for x in proj])+'\n')
            if proj_type == 'ao':
                for proj in projections:
                    pdoss_block.append(str(len(proj))+' ' +
                                       ' '.join([str(x) for x in proj])+'\n')
            elif proj_type != 'atom' and proj_type != 'ao':
                print(
                    'EXITING: please specify either atom or ao projection. %s selected' % proj_type)
                sys.exit(1)
        elif all(isinstance(x, str) for x in flat_proj):
            if output_file == None:
                print(
                    'EXITING: please specify an output file to use the atoms projection.')
                sys.exit(1)
            else:
                output = Crystal_output(output_file)
                output.get_last_geom()
                atoms_symbols = output.atom_symbols
                atoms_symbols.insert(0, 0)

                for proj in projections:
                    atom_positions_list = []
                    for element in proj:
                        index = [i for i, ele in enumerate(
                            atoms_symbols) if ele == element.upper()]
                        atom_positions_list.append([str(x) for x in index])
                    pdoss_block.append(
                        str(-len(index))+' '+' '.join([str(x) for x in index])+'\n')

        pdoss_block.append('END\n')

        self.property_block = pdoss_block

        return self

    def write_properties_input(self, input_name):
        """
        Writes the properties input to a file.

        Args:
            input_name (str): Name of the output file.
        """

        import itertools
        import sys

        if self.is_newk == False:
            property_input_list = self.property_block
        if self.is_newk == True:
            property_input_list = list(itertools.chain(
                self.newk_block, self.property_block))

        with open(input_name, 'w') as file:
            for line in property_input_list:
                file.writelines(line)


class Properties_output:
    """Creates a Properties_output object."""

    def __init__(self):
        """Initialize the Properties_output object."""

        pass

    def read_file(self, properties_output):
        """Parse the properties output file.

        Args:
            properties_output (str): Properties output file.
        Returns:
            Properties_output: Updated Properties_output object.
        """
        import os

        self.file_name = properties_output

        try:
            file = open(self.file_name, 'r')
            self.data = file.readlines()
            file.close()

            # directory
            dir_name = os.path.split(properties_output)[0]
            self.abspath = os.path.join(dir_name)

            # title (named "title" only to distinguish from "file_name" which means another thing)
            self.title = os.path.split(properties_output)[1]

        except:
            raise FileNotFoundError(
                'EXITING: a CRYSTAL properties file needs to be specified')

    def read_electron_band(self, band_file, output=None):
        """
        Generate bands object from CRYSTAL BAND.DAT or fort.25 file.
        Energy unit: eV. E Fermi is aligned to 0.

        Args:
            band_file (str): Name of BAND.DAT or fort.25 file
            output (str): Properties output file (.outp or .out). For 3D, k
                coordinates and geometry information.

        Returns:
            self.bands (BandsBASE): A Bands base object
        """
        from CRYSTALClear.base.propout import BandsBASE, OutBASE

        self.read_file(band_file)
        if '-%-' in self.data[0]:  # fort.25 file format
            self.bands = BandsBASE.f25_parser(self.data)
        else:  # BAND.DAT file format
            self.bands = BandsBASE.BAND_parser(self.data)

        if output != None:
            self.bands.geometry = OutBASE.get_geometry(output)
            self.bands.tick_pos3d, self.bands.k_point_pos3d = OutBASE.get_3dkcoord(
                output)
            self.bands.reciprocal_latt = self.bands.geometry.lattice.reciprocal_lattice.matrix

        return self.bands

    def read_electron_dos(self, properties_output):
        """
        Generate doss object from CRYSTAL DOSS.DAT or fort.25 file.
        Energy unit: eV. E Fermi is aligned to 0.

        Args:
            properties_output (str): File name

        Returns:
            self.doss (DOSBASE): A DOS base object
        """
        from CRYSTALClear.base.propout import DOSBASE

        self.read_file(properties_output)
        if '-%-' in self.data[0]:  # fort.25 file format
            self.doss = DOSBASE.f25_parser(self.data)
        else:  # DOSS.DAT file format
            self.doss = DOSBASE.DOSS_parser(self.data)

        return self.doss

    def read_cry_bands(self, properties_output):
        """
        Deprecated.
        """
        import warnings

        warnings.warn('Deprecated. Use read_electron_band instead.')
        return self.read_electron_band(properties_output)

    def read_cry_doss(self, properties_output):
        """
        Deprecated.
        """
        import warnings

        warnings.warn('Deprecated. Use read_electron_dos instead.')
        return self.read_electron_dos(properties_output)

    def read_cry_contour(self, properties_output):
        """Read the CRYSTAL contour files (SURFRHOO, SURFLAPP, SURFLAPM, SURFGRHO,
           SURFELFB, SURFVIRI, SURFGKIN, SURFKKIN) to create the contour objects.

        Args:
            properties_output (str): Properties output file.
        Returns:
            Properties_output: Updated Properties_output object.
        """
        import re
        import sys

        import numpy as np
        import pandas as pd

        self.read_file(properties_output)

        filename = str(properties_output)

        tipo = ''

        if (filename.endswith('.SURFRHOO')):
            self.tipo = 'SURFRHOO'
            self.path = filename
        elif (filename.endswith('.SURFLAPP')):
            self.tipo = 'SURFLAPP'
            self.path = filename
        elif (filename.endswith('.SURFLAPM')):
            self.tipo = 'SURFLAPM'
            self.path = filename
        elif (filename.endswith('.SURFGRHO')):
            self.tipo = 'SURFGRHO'
            self.path = filename
        elif (filename.endswith('.SURFELFB')):
            self.tipo = 'SURFELFB'
            self.path = filename
        elif (filename.endswith('.SURFVIRI')):
            self.tipo = 'SURFVIRI'
            self.path = filename
        elif (filename.endswith('.SURFGKIN')):
            self.tipo = 'SURFGKIN'
            self.path = filename
        elif (filename.endswith('.SURFKKIN')):
            self.tipo = 'SURFKKIN'
            self.path = filename
        else:
            sys.exit('Please choose a valid file')

        l_dens = self.data

        n_punti_x = int(l_dens[1].strip().split()[0])
        n_punti_y = int(l_dens[1].strip().split()[1])

        self.npx = n_punti_x

        x_min = units.au_to_angstrom(float(l_dens[2].strip().split()[0]))
        x_max = units.au_to_angstrom(float(l_dens[2].strip().split()[1]))
        x_step = units.au_to_angstrom(float(l_dens[2].strip().split()[2]))

        y_min = units.au_to_angstrom(float(l_dens[3].strip().split()[0]))
        y_max = units.au_to_angstrom(float(l_dens[3].strip().split()[1]))
        y_step = units.au_to_angstrom(float(l_dens[3].strip().split()[2]))

        l_dens = l_dens[5:]

        m_dens = []
        for i in l_dens:
            m_dens.append(re.sub("\s\s+", " ", i))

        n_dens = []
        for i in m_dens:
            n_dens.append(i.replace('\n', '').split())

        self.df = pd.DataFrame(n_dens)

        self.x_points = np.linspace(x_min, x_max, n_punti_x)
        self.y_points = np.linspace(y_min, y_max, n_punti_y)

        a = x_max - x_min
        b = y_max - y_min
        r = a/b

        self.x_graph_param = 10
        self.y_graph_param = 10 / r

        ctr1 = np.array([0.002, 0.004, 0.008, 0.02, 0.04,
                         0.08, 0.2, 0.4, 0.8, 2, 4, 8, 20])
        colors1 = ['r', 'r', 'r', 'r', 'r', 'r',
                   'r', 'r', 'r', 'r', 'r', 'r', 'r']
        ls1 = ['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']

        ctr2 = np.array([-8, -4, -2, -0.8, -0.4, -0.2, -0.08, -0.04, -0.02, -0.008, -0.004, -0.002, 0.002, 0.004, 0.008, 0.02, 0.04, 0.08,
                         0.2, 0.4, 0.8, 2, 4, 8])
        colors2 = ['b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b',
                   'b', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
        ls2 = ['--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--',
               '--', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']

        ctr3 = np.array([0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90,
                         0.95, 1])
        colors3 = ['k', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b',
                   'b', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
        ls3 = ['dotted', '--', '--', '--', '--', '--', '--', '--', '--',
               '--', '--', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']

        if (self.tipo == 'SURFRHOO') or (self.tipo == 'SURFGRHO') or (self.tipo == 'SURFGKIN'):
            self.levels = ctr1
            self.colors = colors1
            self.linestyles = ls1
            self.fmt = '%1.3f'
        elif (self.tipo == 'SURFLAPP') or (self.tipo == 'SURFLAPM') or (self.tipo == 'SURFVIRI') or (self.tipo == 'SURFKKIN'):
            self.levels = ctr2
            self.colors = colors2
            self.linestyles = ls2
            self.fmt = '%1.3f'
        elif (self.tipo == 'SURFELFB'):
            self.levels = ctr3
            self.colors = colors3
            self.linestyles = ls3
            self.fmt = '%1.2f'

        return self

    def read_cry_xrd_spec(self, properties_output):
        """
        Read XRD spectrum data from a file.

        Args:
            properties_output (str): Path to the properties output file.
        Returns:
            self: The modified object with extracted XRD spectrum data.
        """
        import re
        import sys

        import pandas as pd

        self.read_file(properties_output)

        data = self.data
        filename = self.abspath
        title = self.title

        if filename.endswith('.outp'):
            pass
        else:
            sys.exit('please, choose a valid file or rename it properly')

        spectrum = re.compile(
            '2THETA    INTENS  INTENS-LP INTENS-LP-DW', re.DOTALL)

        match = []

        a = 0

        for line in data:
            if spectrum.search(line):
                match.append('WRITE LINE:' + line)
                a = 1
            else:
                match.append('WRONG LINE:' + line)

        if (a == 0):
            sys.exit('please, choose a valid file or rename it properly')

        df = pd.DataFrame(match)

        num_riga = (df[df[0].str.contains(u'WRITE')].index.values)

        num_riga = num_riga[0]

        match = match[num_riga:]

        pattern = re.compile(
            '\s+ \d+\.\d+ \s+ \d+\.\d+ \s+ \d+\.\d+ \s+ \d+\.\d+\n', re.DOTALL)

        match_2 = []

        for line in match:
            if pattern.search(line):
                # pulisco dalle scritte di prima
                line = line.replace('WRONG LINE:', '')
                match_2.append(line)

        df = pd.DataFrame([i.strip().split() for i in match_2])

        for i in range(0, 4):
            df[i] = df[i].astype(float)

        df = df.rename(columns={0: '2THETA', 1: 'INTENS',
                                2: 'INTENS-LP', 3: 'INTENS-LP-DW'})

        self.x = df['2THETA']
        self.y = df['INTENS-LP']

        self.title = title[:-1]

        return self

    def read_cry_rholine(self, properties_output):
        """
        Read density line data from a file.

        Args:
            properties_output (str): Path to the properties output file.
        Returns:
            self: Modified object with extracted density line data.
        """
        import re
        import sys

        import pandas as pd

        self.read_file(properties_output)

        l_dens = self.data
        filename = self.abspath
        title = self.title

        if filename.endswith('.RHOLINE'):
            pass
        else:
            sys.exit('please, choose a valid file or rename it properly')

        m_dens = []
        for i in l_dens:
            m_dens.append(re.sub("\s\s+", " ", i))

        n_dens = []
        for i in m_dens:
            n_dens.append(i.replace('\n', '').split())

        df_dens = pd.DataFrame(n_dens)
        df_dens = df_dens.dropna()

        for i in range(0, len(df_dens.columns)):
            df_dens[i] = pd.to_numeric(df_dens[i])

        self.x = units.au_to_angstrom((df_dens[0]-5.55))
        self.y = df_dens[1]/0.148184743

        self.title = title[:-4]

        return self

    def read_cry_seebeck(self, properties_output):
        """
        Read Seebeck coefficient data from a file.

        Args:
            properties_output (str): Path to the properties output file.
        Returns:
            self: Modified object with extracted Seebeck coefficient data.
        """
        import re
        import sys

        import pandas as pd

        self.read_file(properties_output)

        data = self.data
        filename = self.abspath
        title = self.title

        spectrum = re.compile('Npoints', re.DOTALL)

        match = []

        for line in data:
            if spectrum.search(line):
                match.append('RIGHT LINE:' + line)
            else:
                match.append('WRONG LINE:' + line)

        df = pd.DataFrame(match)
        indx = list(df[df[0].str.contains("RIGHT")].index)

        lin = []
        for i in indx:
            lin.append(i+1)

        diffs = [abs(x - y) for x, y in zip(lin, lin[1:])]

        length = diffs[0] - 1

        lif = []
        for i in lin:
            lif.append(i+length)

        c = []
        for i in range(len(lin)):
            c.append(lin[i])
            c.append(lif[i])

        d = [c[i:i + 2] for i in range(0, len(c), 2)]

        l = []
        for i in range(0, len(d)):
            pd.DataFrame(l.append(df[d[i][0]:d[i][1]]))

        right = df[df[0].str.contains("RIGHT")]
        right = right.reset_index().drop('index', axis=1)

        self.temp = []

        for i in range(0, len(right)):
            self.temp.append(float(str(right[0][i])[20:24]))

        ll = []
        for k in range(0, len(l)):
            ll.append(l[k].reset_index().drop('index', axis=1))

        self.all_data = []
        for k in range(0, len(ll)):
            for i in ll[k]:
                self.all_data.append(ll[k][i].apply(
                    lambda x: x.replace('WRONG LINE:', '')))

        self.volume = (float(str(match[2:3])[-13:-4]))

        self.title = title

        return self

    def read_cry_sigma(self, properties_output):
        """
        Read electrical conductivity data from a file.

        Args:
            properties_output (str): Path to the properties output file.
        Returns:
            self: Modified object with extracted electrical conductivity data.
        """
        import re
        import sys

        import pandas as pd

        self.read_file(properties_output)

        data = self.data
        filename = self.abspath
        title = self.title

        spectrum = re.compile('Npoints', re.DOTALL)

        match = []

        for line in data:
            if spectrum.search(line):
                match.append('RIGHT LINE:' + line)
            else:
                match.append('WRONG LINE:' + line)

        df = pd.DataFrame(match)
        indx = list(df[df[0].str.contains("RIGHT")].index)

        lin = []
        for i in indx:
            lin.append(i+1)

        diffs = [abs(x - y) for x, y in zip(lin, lin[1:])]

        length = diffs[0] - 1

        lif = []
        for i in lin:
            lif.append(i+length)

        c = []
        for i in range(len(lin)):
            c.append(lin[i])
            c.append(lif[i])

        d = [c[i:i + 2] for i in range(0, len(c), 2)]

        l = []
        for i in range(0, len(d)):
            pd.DataFrame(l.append(df[d[i][0]:d[i][1]]))

        right = df[df[0].str.contains("RIGHT")]
        right = right.reset_index().drop('index', axis=1)

        self.temp = []

        for i in range(0, len(right)):
            self.temp.append(float(str(right[0][i])[20:24]))

        ll = []
        for k in range(0, len(l)):
            ll.append(l[k].reset_index().drop('index', axis=1))

        self.all_data = []
        for k in range(0, len(ll)):
            for i in ll[k]:
                self.all_data.append(ll[k][i].apply(
                    lambda x: x.replace('WRONG LINE:', '')))

        self.volume = (float(str(match[2:3])[-13:-4]))

        self.title = title

        return self

    def read_cry_lapl_profile(self, properties_output):
        """
        Read Laplacian profile data from a file.

        Args:
            properties_output (str): Path to the properties output file.
        Returns:
            self: Modified object with extracted Laplacian profile data.
        """
        import re

        import numpy as np
        import pandas as pd

        data = self.data
        filename = self.abspath
        title = self.title

        self.read_file(properties_output)

        spectrum = re.compile('PROFILE ALONG THE POINTS', re.DOTALL)

        match = []

        for line in data:
            if spectrum.search(line):
                match.append('RIGHT LINE: ' + line)
            else:
                match.append('WRONG LINE: ' + line)

        df = pd.DataFrame(match)
        num_riga = (df[df[0].str.contains(u'RIGHT')].index.values)
        num_in = num_riga + 8

        spectrum_fin = re.compile('EEEEEEEEEE TERMINATION  DATE', re.DOTALL)

        match_fin = []

        for line in data:
            if spectrum_fin.search(line):
                match_fin.append('RIGHT LINE: ' + line)
            else:
                match_fin.append('WRONG LINE: ' + line)

        df_fin = pd.DataFrame(match_fin)

        num_fin = (df_fin[df_fin[0].str.contains(u'RIGHT')].index.values)
        match = match[num_in[0]:num_fin[0]-2]

        match_2 = []
        for line in match:
            line = line.replace('WRONG LINE:', '')
            match_2.append(line)

        df = pd.DataFrame([i.strip().split() for i in match_2])
        df = df.rename({0: 'x', 1: 'y', 2: 'z', 3: 'dist',
                        4: 'rho', 5: 'lapl', 6: 'L3', 7: 'ellip'}, axis=1)
        df = df.drop(df[df.lapl == '********'].index)
        df = df.reset_index().drop('index', axis=1)
        df['lapl'] = df['lapl'].apply(pd.to_numeric)
        df['dist'] = df['dist'].apply(pd.to_numeric)
        self.datax = df.dist
        self.datay = df.lapl

        return self

    def read_cry_density_profile(self, properties_output):
        """
        Read density profile data from a file.

        Args:
            properties_output (str): Path to the properties output file.
        Returns:
            self: Modified object with extracted density profile data.
        """
        import re

        import numpy as np
        import pandas as pd

        self.read_file(properties_output)

        data = self.data
        filename = self.abspath
        title = self.title

        spectrum = re.compile('PROFILE ALONG THE POINTS', re.DOTALL)

        match = []

        for line in data:
            if spectrum.search(line):
                match.append('RIGHT LINE: ' + line)
            else:
                match.append('WRONG LINE: ' + line)

        df = pd.DataFrame(match)
        num_riga = (df[df[0].str.contains(u'RIGHT')].index.values)
        num_in = num_riga + 8

        spectrum_fin = re.compile('EEEEEEEEEE TERMINATION  DATE', re.DOTALL)

        match_fin = []

        for line in data:
            if spectrum_fin.search(line):
                match_fin.append('RIGHT LINE: ' + line)
            else:
                match_fin.append('WRONG LINE: ' + line)

        df_fin = pd.DataFrame(match_fin)

        num_fin = (df_fin[df_fin[0].str.contains(u'RIGHT')].index.values)
        match = match[num_in[0]:num_fin[0]-2]

        match_2 = []
        for line in match:
            line = line.replace('WRONG LINE:', '')
            match_2.append(line)

        df = pd.DataFrame([i.strip().split() for i in match_2])
        df = df.rename({0: 'x', 1: 'y', 2: 'z', 3: 'dist',
                        4: 'rho', 5: 'lapl', 6: 'L3', 7: 'ellip'}, axis=1)
        df = df.drop(df[df.lapl == '********'].index)
        df = df.reset_index().drop('index', axis=1)
        df['rho'] = df['rho'].apply(pd.to_numeric)
        df['dist'] = df['dist'].apply(pd.to_numeric)
        self.datax = df.dist
        self.datay = df.rho

        return self

    def read_vecfield(self, properties_output, which_prop):
        """Reads the fort.25 file to return data arrays containing one or more vectiorial density properties.

        Args:
            properties_output (str): Properties output file.
            which_prop (str): Density property selected by the user: "m" (magnetization), "j" (spin current), "J" (spin current density)
        Returns:
            Properties_output (str): fort.25 output file.
        """

        import numpy as np

        self.read_file(properties_output)

        data = self.data

        # Reads the header information
        self.nrow = int(data[0].split()[1])
        self.ncol = int(data[0].split()[2])
        self.stepx = float(data[0].split()[3])
        self.stepy = float(data[0].split()[4])
        self.cosxy = float(data[0].split()[5])

        self.a = np.array([float(data[1].split()[0]), float(data[1].split()[1]),
                           float(data[1].split()[2])])
        self.b = np.array([float(data[1].split()[3]), float(data[1].split()[4]),
                           float(data[1].split()[5])])
        self.c = np.array([float(data[2].split()[0]), float(data[2].split()[1]),
                           float(data[2].split()[2])])
        self.naf = int(data[2].split()[3])
        self.ldim = int(data[2].split()[4])

        # Elaborates the header data

        gridsize = (self.nrow*self.ncol)
        lines = int(gridsize/6)
        if (gridsize % 6) != 0:
            lines += 1

        # Reads the types of density properties requested by the user and initializes the data arrays
        check = np.zeros(3, dtype=int)
        if 'm' in which_prop:
            check[0] = 1
            dens_m_array = np.zeros(gridsize*3, dtype=float)
            self.dens_m = np.zeros((self.nrow, self.ncol, 3), dtype=float)
        if 'j' in which_prop:
            check[1] = 1
            dens_j_array = np.zeros(gridsize*3, dtype=float)
            self.dens_j = np.zeros((self.nrow, self.ncol, 3), dtype=float)
        if 'J' in which_prop:
            check[2] = 1
            dens_JX_array = np.zeros(gridsize*3, dtype=float)
            self.dens_JX = np.zeros((self.nrow, self.ncol, 3), dtype=float)
            dens_JY_array = np.zeros(gridsize*3, dtype=float)
            self.dens_JY = np.zeros((self.nrow, self.ncol, 3), dtype=float)
            dens_JZ_array = np.zeros(gridsize*3, dtype=float)
            self.dens_JZ = np.zeros((self.nrow, self.ncol, 3), dtype=float)
        if (not check[0]) and (not check[1]) and (not check[2]):
            print('Error: Invalid Entry. Only the m, j, and J characters are supported')
            sys.exit(1)

        # Gathers the data
        pointer = 0
        skip = 6 + self.naf

        if check[0]:
            pointer = 3
            k = 0
            for i in range(0, lines):
                for j in range(0, len(data[i+pointer].split())):
                    dens_m_array[k] = data[i+pointer].split()[j]
                    k += 1
            pointer = pointer + lines + skip
            for i in range(0, lines):
                for j in range(0, len(data[i+pointer].split())):
                    dens_m_array[k] = data[i+pointer].split()[j]
                    k += 1
            pointer = pointer + lines + skip
            for i in range(0, lines):
                for j in range(0, len(data[i+pointer].split())):
                    dens_m_array[k] = data[i+pointer].split()[j]
                    k += 1
            pointer = pointer + lines + skip

            k = 0
            for dir in range(0, 3):
                for j in range(0, self.ncol):
                    for i in range(self.nrow - 1, -1, -1):
                        self.dens_m[i, j, dir] = dens_m_array[k]
                        k += 1
            del dens_m_array

        if check[1]:
            if pointer == 0:
                pointer = 3
            k = 0
            for i in range(0, lines):
                for j in range(0, len(data[i+pointer].split())):
                    dens_j_array[k] = data[i+pointer].split()[j]
                    k += 1
            pointer = pointer + lines + skip
            for i in range(0, lines):
                for j in range(0, len(data[i+pointer].split())):
                    dens_j_array[k] = data[i+pointer].split()[j]
                    k += 1
            pointer = pointer + lines + skip
            for i in range(0, lines):
                for j in range(0, len(data[i+pointer].split())):
                    dens_j_array[k] = data[i+pointer].split()[j]
                    k += 1
            pointer = pointer + lines + skip

            k = 0
            for dir in range(0, 3):
                for j in range(0, self.ncol):
                    for i in range(self.nrow - 1, -1, -1):
                        self.dens_j[i, j, dir] = dens_j_array[k]
                        k += 1
            del dens_j_array

        if check[2]:
            if pointer == 0:
                pointer = 3
            k = 0
            for i in range(0, lines):
                for j in range(0, len(data[i+pointer].split())):
                    dens_JX_array[k] = data[i+pointer].split()[j]
                    k += 1
            pointer = pointer + lines + skip
            for i in range(0, lines):
                for j in range(0, len(data[i+pointer].split())):
                    dens_JX_array[k] = data[i+pointer].split()[j]
                    k += 1
            pointer = pointer + lines + skip
            for i in range(0, lines):
                for j in range(0, len(data[i+pointer].split())):
                    dens_JX_array[k] = data[i+pointer].split()[j]
                    k += 1
            pointer = pointer + lines + skip

            k = 0
            for dir in range(0, 3):
                for j in range(0, self.ncol):
                    for i in range(self.nrow - 1, -1, -1):
                        self.dens_JX[i, j, dir] = dens_JX_array[k]
                        k += 1
            del dens_JX_array

            k = 0
            for i in range(0, lines):
                for j in range(0, len(data[i+pointer].split())):
                    dens_JY_array[k] = data[i+pointer].split()[j]
                    k += 1
            pointer = pointer + lines + skip
            for i in range(0, lines):
                for j in range(0, len(data[i+pointer].split())):
                    dens_JY_array[k] = data[i+pointer].split()[j]
                    k += 1
            pointer = pointer + lines + skip
            for i in range(0, lines):
                for j in range(0, len(data[i+pointer].split())):
                    dens_JY_array[k] = data[i+pointer].split()[j]
                    k += 1
            pointer = pointer + lines + skip

            k = 0
            for dir in range(0, 3):
                for j in range(0, self.ncol):
                    for i in range(self.nrow - 1, -1, -1):
                        self.dens_JY[i, j, dir] = dens_JY_array[k]
                        k += 1
            del dens_JY_array

            k = 0
            for i in range(0, lines):
                for j in range(0, len(data[i+pointer].split())):
                    dens_JZ_array[k] = data[i+pointer].split()[j]
                    k += 1
            pointer = pointer + lines + skip
            for i in range(0, lines):
                for j in range(0, len(data[i+pointer].split())):
                    dens_JZ_array[k] = data[i+pointer].split()[j]
                    k += 1
            pointer = pointer + lines + skip
            for i in range(0, lines):
                for j in range(0, len(data[i+pointer].split())):
                    dens_JZ_array[k] = data[i+pointer].split()[j]
                    k += 1

            k = 0
            for dir in range(0, 3):
                for j in range(0, self.ncol):
                    for i in range(self.nrow - 1, -1, -1):
                        self.dens_JZ[i, j, dir] = dens_JZ_array[k]
                        k += 1
            del dens_JZ_array

        return self

    def read_cry_ECHG(self, properties_output):
        """
        Read density profile data from a file.

        Args:
            properties_output (str): Path to the fort.25 file.
        Returns:
            self: Modified object with extracted ECHG data.
        """
        import numpy as np

        self.read_file(properties_output)
        data = self.data

        self.nrow = int(data[0].split()[1])
        self.ncol = int(data[0].split()[2])
        self.cosxy = float(data[0][42:54])
        self.a = np.array([data[1][0:12],
                           data[1][12:24],
                           data[1][24:36]], dtype=float)
        self.b = np.array([data[1][36:48],
                           data[1][48:60],
                           data[1][60:72]], dtype=float)
        self.c = np.array([data[2][0:12],
                           data[2][12:24],
                           data[2][24:36]], dtype=float)

        self.density_map = np.zeros((self.nrow, self.ncol), dtype=float)

        number_points = self.nrow * self.ncol
        lines = int(number_points / 6)
        if number_points % 6 != 0:
            lines = lines + 1

        density_temp = np.zeros(number_points, dtype=float)
        j = 0
        for i in range(0, lines):
            for k in range(0, len(data[3 + i].split())):
                density_temp[j] = data[3 + i].split()[k]
                j += 1

        k = 0
        for j in range(0, self.ncol):
            for i in range(self.nrow - 1, -1, -1):
                self.density_map[i, j] = density_temp[k]
                k += 1

        return self

    def read_cry_ECHG_delta(self, properties_output1, properties_output2):
        """
        Read density profile data from two files and plots the difference.
        It is important that the header stays the same for the two output files.

        Args:
            properties_output1 (str): Path to first fort.25 file.
            properties_output2 (str): Path to second fort.25 file.
        Returns:
            self: Modified object with extracted ECHG data.
        """

        out1 = Properties_output().read_cry_ECHG(properties_output1)
        out2 = Properties_output().read_cry_ECHG(properties_output2)

        self.a = out1.a
        self.b = out1.b
        self.c = out1.c
        self.cosxy = out1.cosxy

        self.density_map = out1.density_map - out2.density_map

        return self

    # TOPOND related functions
    def read_topond_trho(self, properties_output):
        """Read the TOPOND TRHO run output file to create associated DataFrames:
         * Coordinates for the nuclei at the unitcell (nuclei_df)
         * Data of the critical points (CP) discoverd (topo_df)

        Args:
            properties_output (str): Properties output file.
        Returns:
            Properties_output: Updated Properties_output object.
        """
        import re

        import numpy as np
        import pandas as pd

        self.read_file(properties_output)

        self.topo_filename = str(properties_output).split('.')[0]

        properties = ["type", "coord", "rho", "grho", "lap",
                      "kinetic_g", "kinetic_k", "virial", "elf",
                      "eigenval", "eigenvec", "ellip", "atom_a_id",
                      "atom_a_cell" ,"atom_a_z", "atom_b_id",
                      "atom_b_cell" ,"atom_b_z", "bp_length",
                      "distance_ab", "bp/dist"]

        self.topo_df = pd.DataFrame(columns=properties)
        self.nuclei_df = pd.DataFrame(columns=['z','coord'])

        atom_a_read = False
        bp_step = False

        file_lines = len(self.data)
        i = 0
        while  i < file_lines:
            line = self.data[i]

            if re.match(r' DIRECT LATTICE VECTOR COMPONENTS .ANGSTROM.',line) != None:
                unitcell_mat = [self.data[i+1].strip().split(),
                                self.data[i+2].strip().split(),
                                self.data[i+3].strip().split()]
                self.unitcell_mat = np.array(unitcell_mat, dtype=float)
                i += 3

            if re.match(r' N. OF ATOMS PER CELL', line) != None:
                num_atoms = int(line.strip().split()[5])

            if re.match(r'   ATOM N.AT.  SHELL', line) != None:
                i += 1
                for n_atom in range(1,num_atoms+1):
                    atom_data = self.data[i+n_atom].strip().split()
                    self.nuclei_df.loc[n_atom,'z'] = int(atom_data[1])
                    self.nuclei_df.loc[n_atom,'coord'] = np.array(atom_data[4:7],dtype=float)
                i += num_atoms

            if re.match(r' SEARCH OF BOND PATH',line) != None:
                bp_step = True

            if re.match(r' CP N.',line) != None:
                if bp_step == False:
                    crit_point_line = i
                    crit_point_number = line.strip().split()[2]

                    if crit_point_number == "X(ANG)":
                        break

                    cp_type = self.data[crit_point_line+3].strip().split()[3]
                    cp_coord = self.data[crit_point_line+4].strip().split()[5:]
                    if self.data[crit_point_line+5].strip().split()[1] == 'FRACT.':
                        crit_point_line += 1
                        i += 1
                    cp_prop = self.data[crit_point_line+5].strip().split()[3:]
                    cp_kener = self.data[crit_point_line+6].strip().split()[5:]
                    cp_virial = self.data[crit_point_line+7].strip().split()[3]
                    cp_elf = self.data[crit_point_line+8].strip().split()[2]
                    eigenval = self.data[crit_point_line+12].strip().split()[5:]
                    eigenvec = []
                    eigenvec.append(self.data[crit_point_line+13].strip().split()[2:])
                    eigenvec.append(self.data[crit_point_line+14].strip().split())
                    eigenvec.append(self.data[crit_point_line+15].strip().split())

                    # To skip already filtered lines
                    i+= 14

                    # Only bond cp has ellipticity calculated
                    if cp_type == "(3,-1)":
                        ellip = self.data[crit_point_line+17].strip().split()[2]
                        i += 1
                    else:
                        ellip = 0

                    # Populating the dataframe
                    self.topo_df.loc[crit_point_number,'type'] = cp_type
                    self.topo_df.loc[crit_point_number,'coord'] = np.array(cp_coord, dtype=float)
                    self.topo_df.loc[crit_point_number,'rho'] = float(cp_prop[0])
                    self.topo_df.loc[crit_point_number,'grho'] = float(cp_prop[1])
                    self.topo_df.loc[crit_point_number,'lap'] = float(cp_prop[2])
                    self.topo_df.loc[crit_point_number,'kinetic_g'] = float(cp_kener[0])
                    self.topo_df.loc[crit_point_number,'kinetic_k'] = float(cp_kener[1])
                    self.topo_df.loc[crit_point_number,'virial'] = float(cp_virial)
                    self.topo_df.loc[crit_point_number,'elf'] = float(cp_elf)
                    self.topo_df.loc[crit_point_number,'eigenval'] = np.array(eigenval, dtype=float)
                    self.topo_df.loc[crit_point_number,'eigenvec'] = np.array(eigenvec, dtype=float)
                    self.topo_df.loc[crit_point_number,'ellip'] = float(ellip)

                else:
                    # ecrit stands for extra critical point, li
                    ecrit_point_line = i
                    ecrit_point_number = line.strip().split()[2]

                    if ecrit_point_number == "X(ANG)":
                        break

                    cp_type = self.data[ecrit_point_line+3].strip().split()[4]
                    cp_coord = self.data[ecrit_point_line+4].strip().split()[5:]
                    if self.data[ecrit_point_line+5].strip().split()[1] == 'FRACT.':
                        ecrit_point_line += 1
                        i += 1
                    cp_prop = self.data[ecrit_point_line+5].strip().split()[3:]

                    # Populating DF
                    self.topo_df.loc[ecrit_point_number,'type'] = cp_type
                    self.topo_df.loc[ecrit_point_number,'coord'] = np.array(cp_coord, dtype=float)
                    self.topo_df.loc[ecrit_point_number,'rho'] = float(cp_prop[0])
                    self.topo_df.loc[ecrit_point_number,'grho'] = float(cp_prop[1])
                    self.topo_df.loc[ecrit_point_number,'lap'] = float(cp_prop[2])

            # Data of first attractor associated to BCP
            if ((re.match(r' CLUSTER OF ATOMS AROUND THE TERMINUS ',line) != None) and
                    bp_step and
                    (not atom_a_read)):
                atom_a_read = True
                atom_a_start = i
                atom_a_data = self.data[atom_a_start+2].strip().split()
                i += 2
                atom_a_id = atom_a_data[1]

                if len(atom_a_data) > 7:
                    atom_a_cell = atom_a_data[2:5]
                    atom_a_z = atom_a_data[5]
                else:
                    atom_a_cell = []
                    atom_a_z = atom_a_data[2]

            # Data of second attractor associated to BCP
            elif (re.match(r' CLUSTER OF ATOMS AROUND THE TERMINUS ',line) != None and
                    bp_step and
                    atom_a_read):
                atom_a_read = False
                atom_b_start = i
                atom_b_data = self.data[atom_b_start+2].strip().split()
                i += 2
                atom_b_id = atom_b_data[1]

                if len(atom_b_data) > 7:
                    atom_b_cell = atom_b_data[2:5]
                    atom_b_z = atom_b_data[5]
                else:
                    atom_b_cell = []
                    atom_b_z = atom_b_data[2]

            if re.match(r' BPL',line) != None:
                bp_step = False
                bp_length = line.strip().split()[7]
                distance_ab = line.strip().split()[8]
                bp_dist_ratio = line.strip().split()[9]

                # Populating info for atoms involved in bond CPs
                self.topo_df.loc[crit_point_number,'atom_a_id'] = int(atom_a_id)
                self.topo_df.loc[crit_point_number,'atom_a_cell'] = np.array(atom_a_cell, dtype=float)
                self.topo_df.loc[crit_point_number,'atom_a_z'] = int(atom_a_z)
                self.topo_df.loc[crit_point_number,'atom_b_id'] = int(atom_b_id)
                self.topo_df.loc[crit_point_number,'atom_b_cell'] = np.array(atom_b_cell, dtype=float)
                self.topo_df.loc[crit_point_number,'atom_b_z'] = int(atom_b_z)
                self.topo_df.loc[crit_point_number,'bp_length'] = float(bp_length)
                self.topo_df.loc[crit_point_number,'distance_ab'] = float(distance_ab)
                self.topo_df.loc[crit_point_number,'bp/dist'] = float(bp_dist_ratio)

            # Next line iterator
            i += 1

        return self

    def topond_viz_file(self, cp_type='ALL', add_atoms=False, file_type="xyz"):
        """Given the existence of TOPOND TRHO DataFrames, write files
         with the critical points (CP) coordinates for further
         visualization.

         Args:
             cp_type (str): Generate file for a given type CP.
                    Options: "ALL" to include all CP in the system,
                             "BCP" for Bond CPs,
                             "RCP" for Ring CPs,
                             "CCP" for Cage CPs,
                             "NNA" for Nuclear/Non-Nuclear Atractors.
             add_atoms (bool): Include the Nuclei atoms read from the
                    TOPOND output in the same file.
             file_type (str): Type of output file to be generated.

         Returns:
             ase_obj: ASE object.
        """
        import sys
        from ase import Atoms

        if not hasattr(self,'topo_df'):
            print('ERROR: You need first to succesfully run read_topond_trho() to use this.')
            sys.exit(1)

         # An ASE object is used as intermediate to generate different types of outputs
        if cp_type == 'ALL':
            if hasattr(self,'unitcell_mat'):
                ase_obj = Atoms(
                    symbols = ["X"] * len(self.topo_df['type']),
                    positions = list(self.topo_df['coord']*0.529177), # factor for Bohr2Angs
                    cell = self.unitcell_mat,
                    pbc=[True,True,True]
                )
            else:
                ase_obj = Atoms(
                    symbols = ["X"] * len(self.topo_df['type']),
                    positions = list(self.topo_df['coord']*0.529177)
                )
        elif cp_type == 'BCP':
            df_bcp =  self.topo_df[self.topo_df['type']=='(3,-1)']
            if hasattr(self,'unitcell_mat'):
                ase_obj = Atoms(
                    symbols = ["X"] * len(df_bcp['type']),
                    positions = list(df_bcp['coord']*0.529177),
                    cell = self.unitcell_mat,
                    pbc=[True,True,True]
                )
            else:
                ase_obj = Atoms(
                    symbols = ["X"] * len(df_bcp['type']),
                    positions = list(df_bcp['coord']*0.529177)
                )
        elif cp_type == 'RCP':
            df_rcp =  self.topo_df[self.topo_df['type']=='(3,1)']
            if hasattr(self,'unitcell_mat'):
                ase_obj = Atoms(
                    symbols = ["X"] * len(df_rcp['type']),
                    positions = list(df_rcp['coord']*0.529177),
                    cell = self.unitcell_mat,
                    pbc=[True,True,True]
                )
            else:
                ase_obj = Atoms(
                    symbols = ["X"] * len(df_rcp['type']),
                    positions = list(df_rcp['coord']*0.529177)
                )
        elif cp_type == 'CCP':
            df_ccp =  self.topo_df[self.topo_df['type']=='(3,3)']
            if hasattr(self,'unitcell_mat'):
                ase_obj = Atoms(
                    symbols = ["X"] * len(df_ccp['type']),
                    positions = list(df_ccp['coord'] * 0.529177),
                    cell = self.unitcell_mat,
                    pbc=[True,True,True]
                )
            else:
                ase_obj = Atoms(
                    symbols = ["X"] * len(df_ccp['type']),
                    positions = list(df_ccp['coord'] * 0.529177)
                )
        elif cp_type == 'NNA':
            df_nna =  self.topo_df[self.topo_df['type']=='(3,-3)']
            if hasattr(self,'unitcell_mat'):
                ase_obj = Atoms(
                    symbols = ["X"] * len(df_nna['type']),
                    positions = list(df_nna['coord'] * 0.529177), # factor for Bohr2Angs
                    cell = self.unitcell_mat,
                    pbc=[True,True,True]
                )
            else:
                ase_obj = Atoms(
                    symbols = ["X"] * len(df_nna['type']),
                    positions = list(df_nna['coord'] * 0.529177)
                )

        if add_atoms:
            for i in range(len(self.nuclei_df['z'])):
                ase_obj.append(self.nuclei_df.loc[i + 1, 'z'])
                ase_obj.positions[-1] = self.nuclei_df.loc[i + 1, 'coord']

        ase_obj.write(self.topo_filename + '_' + cp_type + '.' +file_type)
        print(self.topo_filename + '_' + cp_type + '.' + file_type + " generated.")

        return ase_obj
    
    def read_topond_tlap(self,properties_output):
        """Read the TOPOND TLAP run output file to create associated DataFrames:
         * Coordinates for the nuclei at the unitcell (nuclei_df)
         * Data of the laplacian critical points (CP) discoverd (tlap_df)

        Args:
            properties_output (str): Properties output file.
        Returns:
            Properties_output: Updated Properties_output object.
        """
        import re

        import numpy as np
        import pandas as pd

        self.read_file(properties_output)

        self.tlap_filename = str(properties_output).split('.')[0]

        properties = ["central_atom_id","central_atom_type","type", "coord", "rho", "grho", "-lap",
                    "kinetic_g", "kinetic_k", "virial", "elf",
                    "eigenval", "eigenvec"] 
                  # To be included. Attr stands for attractor.
                  # "attr_a_id","attr_a_cell" ,"attr_a_z", 
                  # "attr_b_id", "attr_b_cell" ,"attr_b_z", 
                  # "bp_length", "distance_ab", "bp/dist"]

        self.tlap_df = pd.DataFrame(columns=properties)
        self.tlap_df['eigenval'] = self.tlap_df['eigenval'].astype(object)
        self.tlap_df['eigenvec'] = self.tlap_df['eigenvec'].astype(object)
    
        # To be included
        # self.topo_df['attr_a_cell'] = self.topo_df['atom_a_cell'].astype(object)
        # self.topo_df['attr_b_cell'] = self.topo_df['atom_b_cell'].astype(object)

        self.nuclei_df = pd.DataFrame(columns=['z','coord'])

        atom_a_read = False
        bp_steps = 0
        outdiag_props = False
        global_cp_id = 0
        bp_cp_id = global_cp_id

        file_lines = len(self.data)
        i = 0
        while  i < file_lines:
            line = self.data[i]

            if re.match(r' DIRECT LATTICE VECTOR COMPONENTS .ANGSTROM.',line) != None:
                unitcell_mat = [self.data[i+1].strip().split(),
                                self.data[i+2].strip().split(),
                                self.data[i+3].strip().split()]
                self.unitcell_mat = np.array(unitcell_mat, dtype=float)
                i += 3
        
            if re.match(r' N. OF ATOMS PER CELL', line) != None:
                num_atoms = int(line.strip().split()[5])

            if re.match(r'   ATOM N.AT.  SHELL', line) != None:
                i += 1
                for n_atom in range(1,num_atoms+1):
                    atom_data = self.data[i+n_atom].strip().split()
                    self.nuclei_df.loc[n_atom,'z'] = int(atom_data[1])
                    self.nuclei_df.loc[n_atom,'coord'] = np.array(atom_data[4:7],dtype=float)
                i += num_atoms
            
            if re.match(r' ATTRACTORS OF THE UNIQUE PAIRS',line) != None:
                bp_steps = 2

            if re.match(r' \*\*\*\*\*\* CP SEARCH FOR NON',line) != None:
                c_atom = line.strip().split()[6:8]

            if re.match(r' CP N.',line) != None:
                global_cp_id += 1
                if not bp_steps:
                    bp_cp_id +=1
                    crit_point_line = i
                    crit_point_number = line.strip().split()[2]

                    if crit_point_number == "X(ANG)":
                        break
            
                    cp_type = self.data[crit_point_line+3].strip().split()[3]
                    cp_coord = self.data[crit_point_line+4].strip().split()[5:]
                    if self.data[crit_point_line+5].strip().split()[1] == 'FRACT.':
                        crit_point_line += 1
                        i += 1
                    cp_prop = self.data[crit_point_line+5].strip().split()[3:]
                    self.tlap_df.loc[global_cp_id,'central_atom_id'] = c_atom[0]
                    self.tlap_df.loc[global_cp_id,'central_atom_type'] = c_atom[1]
                    self.tlap_df.loc[global_cp_id,'type'] = cp_type
                    self.tlap_df.at[global_cp_id,'coord'] = np.array(cp_coord, dtype=float)
                    self.tlap_df.loc[global_cp_id,'rho'] = float(cp_prop[2])
                    self.tlap_df.loc[global_cp_id,'grho'] = float(cp_prop[1])
                    self.tlap_df.loc[global_cp_id,'-lap'] = float(cp_prop[0])
                    
                    if not self.data[crit_point_line+6].strip():
                        outdiag_props = False
                    elif self.data[crit_point_line+6].strip().split()[0] == 'KINETIC':
                        outdiag_props = True
                        cp_kener = self.data[crit_point_line+6].strip().split()[5:]
                        cp_virial = self.data[crit_point_line+7].strip().split()[3]
                        cp_elf = self.data[crit_point_line+8].strip().split()[2]
                        self.tlap_df.loc[global_cp_id,'kinetic_g'] = float(cp_kener[0])
                        self.tlap_df.loc[global_cp_id,'kinetic_k'] = float(cp_kener[1])
                        self.tlap_df.loc[global_cp_id,'virial'] = float(cp_virial)
                        self.tlap_df.loc[global_cp_id,'elf'] = float(cp_elf)
                        crit_point_line += 3
                        i += 3
                    else:
                        print("Non considered case, check output and code.")

                    eigenval = self.data[crit_point_line+9].strip().split()[5:]
                    eigenvec = []
                    eigenvec.append(self.data[crit_point_line+10].strip().split()[2:])
                    eigenvec.append(self.data[crit_point_line+11].strip().split())
                    eigenvec.append(self.data[crit_point_line+12].strip().split())
                    self.tlap_df.at[global_cp_id,'eigenval'] = np.array(eigenval, dtype=float)
                    self.tlap_df.at[global_cp_id,'eigenvec'] = np.array(eigenvec, dtype=float)
                    i += 12
                    

                elif bp_steps:
                    #global_cp_id += 1
                    # ecrit stands for extra critical point, li
                    ecrit_point_line = i
                    ecrit_point_number = line.strip().split()[2]

                    if ecrit_point_number == "X(ANG)":
                        break

                    cp_type = self.data[ecrit_point_line+3].strip().split()[4]
                    cp_coord = self.data[ecrit_point_line+4].strip().split()[5:]
                    if self.data[ecrit_point_line+5].strip().split()[1] == 'FRACT.':
                        ecrit_point_line += 1
                        i += 1
                    cp_prop = self.data[ecrit_point_line+5].strip().split()[3:]

                    # Populating DF
                    self.tlap_df.loc[global_cp_id,'central_atom_id'] = c_atom[0]
                    self.tlap_df.loc[global_cp_id,'central_atom_type'] = c_atom[1]
                    self.tlap_df.loc[global_cp_id,'type'] = cp_type
                    self.tlap_df.at[global_cp_id,'coord'] = np.array(cp_coord, dtype=float)
                    self.tlap_df.loc[global_cp_id,'rho'] = float(cp_prop[2])
                    self.tlap_df.loc[global_cp_id,'grho'] = float(cp_prop[1])
                    self.tlap_df.loc[global_cp_id,'-lap'] = float(cp_prop[0])

            if re.match(r' MAX. NUMBER OF STEP EXCEEDED',line) != None:
                bp_steps -= 1

            if ((re.match(r' CLUSTER OF ATOMS AROUND THE TERMINUS ',line) != None) and
                    bp_steps == 2):
                bp_steps -= 1
                # Following should be fixed, as it takes info of closest atom and not
                # of the corresponding attractor.
                # atom_a_start = i
                # atom_a_data = self.data[atom_a_start+2].strip().split()
                # i += 2
                # atom_a_id = atom_a_data[1]

                # if len(atom_a_data) > 7:
                #     atom_a_cell = atom_a_data[2:5]
                #     atom_a_z = atom_a_data[5]
                # else:
                #     atom_a_cell = []
                #     atom_a_z = atom_a_data[2]

                # self.topo_df.loc[bp_cp_id,'atom_a_id'] = int(atom_a_id)
                # self.topo_df.at[bp_cp_id,'atom_a_cell'] = np.array(atom_a_cell, dtype=float)
                # self.topo_df.loc[bp_cp_id,'atom_a_z'] = int(atom_a_z)

            # Data of second attractor associated to BCP
            elif (re.match(r' CLUSTER OF ATOMS AROUND THE TERMINUS ',line) != None and
                    bp_steps == 1):
                bp_steps -= 1
                # Following should be fixed, as it takes info of closest atom and not
                # of the corresponding attractor.
                # atom_b_start = i
                # atom_b_data = self.data[atom_b_start+2].strip().split()
                # i += 2
                # atom_b_id = atom_b_data[1]

                # if len(atom_b_data) > 7:
                #     atom_b_cell = atom_b_data[2:5]
                #     atom_b_z = atom_b_data[5]
                # else:
                #     atom_b_cell = []
                #     atom_b_z = atom_b_data[2]
                
                # self.topo_df.loc[bp_cp_id,'atom_b_id'] = int(atom_b_id)
                # self.topo_df.at[bp_cp_id,'atom_b_cell'] = np.array(atom_b_cell, dtype=float)
                # self.topo_df.loc[bp_cp_id,'atom_b_z'] = int(atom_b_z)


            if re.match(r' TTRAJ',line) != None:
                bp_steps = 0
                # This works, but first the attractors info should be fixed.
                # bp_length = line.strip().split()[6]
                # distance_ab = line.strip().split()[7]
                # bp_dist_ratio = line.strip().split()[8]

                # # Populating info for atoms involved in bond CPs
                
                # self.topo_df.loc[bp_cp_id,'bp_length'] = float(bp_length)
                # self.topo_df.loc[bp_cp_id,'distance_ab'] = float(distance_ab)
                # self.topo_df.loc[bp_cp_id,'bp/dist'] = float(bp_dist_ratio)
                bp_cp_id = global_cp_id

            # Next line iterator
            i += 1

        return self

    def topond_viz_tlap_file(self, cp_type='+3', add_atoms=False, file_type="xyz"):
        """Given the existence of TOPOND TLAP DataFrames, write files
         with the critical points (CP) coordinates for further
         visualization.

         Args:
             cp_type (str): Generate file for a given type CP.
                    Options: 'ALL' to include all CP in the system,
                             '+3' for CPs with +3 signature,
                             '+1' for CPs with +1 signature,
                             '-1' for CPs with -1 signature,
                             '+1' for CPs with -3 signature.
             add_atoms (bool): Include the Nuclei atoms read from the
                    TOPOND output in the same file.
             file_type (str): Type of output file to be generated.

         Returns:
             ase_obj: The generated ASE object.
        """
        import sys
        from ase import Atoms

        if not hasattr(self,'tlap_df'):
            print('ERROR: You need first to succesfully run read_topond_tlap() to use this.')
            sys.exit(1)

         # An ASE object is used as intermediate to generate different types of outputs
        if cp_type == 'ALL':
            if hasattr(self,'unitcell_mat'):
                ase_obj = Atoms(
                    symbols = ["X"] * len(self.tlap_df['type']),
                    positions = list(self.tlap_df['coord']*0.529177), # factor for Bohr2Angs
                    cell = self.unitcell_mat,
                    pbc=[True,True,True]
                )
            else:
                ase_obj = Atoms(
                    symbols = ["X"] * len(self.tlap_df['type']),
                    positions = list(self.tlap_df['coord']*0.529177)
                )
        elif cp_type == '-1':
            cp_type = 'minus1'
            df_bcp =  self.tlap_df[self.tlap_df['type']=='(3,-1)']
            if hasattr(self,'unitcell_mat'):
                ase_obj = Atoms(
                    symbols = ["X"] * len(df_bcp['type']),
                    positions = list(df_bcp['coord']*0.529177),
                    cell = self.unitcell_mat,
                    pbc=[True,True,True]
                )
            else:
                ase_obj = Atoms(
                    symbols = ["X"] * len(df_bcp['type']),
                    positions = list(df_bcp['coord']*0.529177)
                )
        elif cp_type == '+1':
            cp_type = 'plus1'
            df_rcp =  self.tlap_df[self.tlap_df['type']=='(3,+1)']
            if hasattr(self,'unitcell_mat'):
                ase_obj = Atoms(
                    symbols = ["X"] * len(df_rcp['type']),
                    positions = list(df_rcp['coord']*0.529177),
                    cell = self.unitcell_mat,
                    pbc=[True,True,True]
                )
            else:
                ase_obj = Atoms(
                    symbols = ["X"] * len(df_rcp['type']),
                    positions = list(df_rcp['coord']*0.529177)
                )
        elif cp_type == '+3':
            cp_type = 'plus3'
            df_ccp =  self.tlap_df[self.tlap_df['type']=='(3,+3)']
            if hasattr(self,'unitcell_mat'):
                ase_obj = Atoms(
                    symbols = ["X"] * len(df_ccp['type']),
                    positions = list(df_ccp['coord'] * 0.529177),
                    cell = self.unitcell_mat,
                    pbc=[True,True,True]
                )
            else:
                ase_obj = Atoms(
                    symbols = ["X"] * len(df_ccp['type']),
                    positions = list(df_ccp['coord'] * 0.529177)
                )
        elif cp_type == '-3':
            cp_type = 'minus3'
            df_nna =  self.tlap_df[self.tlap_df['type']=='(3,-3)']
            if hasattr(self,'unitcell_mat'):
                ase_obj = Atoms(
                    symbols = ["X"] * len(df_nna['type']),
                    positions = list(df_nna['coord'] * 0.529177), # factor for Bohr2Angs
                    cell = self.unitcell_mat,
                    pbc=[True,True,True]
                )
            else:
                ase_obj = Atoms(
                    symbols = ["X"] * len(df_nna['type']),
                    positions = list(df_nna['coord'] * 0.529177)
                )

        if add_atoms:
            for i in range(len(self.nuclei_df['z'])):
                ase_obj.append(self.nuclei_df.loc[i + 1, 'z'])
                ase_obj.positions[-1] = self.nuclei_df.loc[i + 1, 'coord']

        ase_obj.write(self.tlap_filename + '_' + cp_type + '.' +file_type)
        print(self.tlap_filename + '_' + cp_type + '.' + file_type + " generated.")

        return ase_obj


class Crystal_gui:
    """
    This class can read a CRYSTAL gui file into an object or substrate
    information of the object to generate a gui file.
    """

    def __init__(self):
        pass

    def read_gui(self, gui_file):
        """
        This is used mainly to convert the object into an ASE or pymatgen
        object.

        Args:
            gui_file (str): CRYSTAL structure (gui) file
        """
        import numpy as np

        try:
            file = open(gui_file, 'r')
            data = file.readlines()
            file.close()
        except:
            raise FileNotFoundError(
                'A CRYSTAL geometry file needs to be specified: .gui, fort.34 or opt* files.')

        self.dimensionality = int(data[0].split()[0])
        self.lattice = []
        self.symmops = []
        for i in range(1, 4):
            self.lattice.append([float(x) for x in data[i].split()])
        self.n_symmops = int(data[4].split()[0])
        for i in range(5, 5+self.n_symmops*4):
            self.symmops.append(data[i].split())
        self.symmops = np.reshape(np.array(self.symmops, dtype=float),
                                  [self.n_symmops, 4, 3])
        self.n_atoms = int(data[5+self.n_symmops*4].split()[0])
        self.atom_number = []
        self.atom_positions = []
        for i in range(6+self.n_symmops*4, 6+self.n_symmops*4+self.n_atoms):
            atom_line = data[i].split()
            self.atom_number.append(int(atom_line[0]))
            self.atom_positions.append([float(x) for x in atom_line[1:]])
        self.space_group = int(data[-1].split()[0])

        return self

    def write_gui(self, gui_file, symm=True, pseudo_atoms=[]):
        """
        Write a CRYSTAL gui file (to file)

        Args:
            gui_file (str): Name of the gui that is going to be written (including .gui).
            symm (bool): Whether to include symmetry operations (default is True).
            pseudo_atoms (list[int]): A list of atoms whose core is described
                by a pseudopotential
        """
        import numpy as np

        if symm == False:
            self.n_symmops = 1
            self.symmops = np.vstack([np.eye(3), [0.0, 0.0, 0.0]])

        file = open(gui_file, 'w')
        # First line
        file.writelines('%4s   1   1\n' % self.dimensionality)
        # Cell vectors
        for vector in self.lattice:
            file.writelines('{}\n'.format(
                ''.join(['{0: 20.12E}'.format(np.round(n, 12))
                        for n in vector])
            ))
        # N symm ops
        file.writelines('{:5d}\n'.format(self.n_symmops))

        # symm ops
        sym_list = np.reshape(self.symmops, [self.n_symmops*4, 3])
        for symmops in sym_list:
            file.writelines('{}\n'.format(
                ''.join(['{0: 20.12f}'.format(np.round(n, 12))
                        for n in symmops])
            ))
        # N atoms
        file.writelines('{:5d}\n'.format(self.n_atoms))

        # atom number (including pseudopotentials) + coordinates cart
        for i in range(self.n_atoms):
            if self.atom_number[i] in pseudo_atoms:
                file.writelines('{:5d}{}\n'.format(
                    int(self.atom_number[i])+200,
                    ''.join(['{0: 20.12E}'.format(np.round(x, 12))
                            for x in self.atom_positions[i]])
                ))
            else:
                file.writelines('{:5d}{}\n'.format(
                    int(self.atom_number[i]),
                    ''.join(['{0: 20.12E}'.format(np.round(x, 12))
                            for x in self.atom_positions[i]])
                ))

        # space group + n symm ops
        if symm == True:
            file.writelines('{:5d}{:5d}\n'.format(
                self.space_group, self.n_symmops
            ))
        else:
            file.writelines('{:5d}{:5d}\n'.format(1, 1))

        file.close()


class Crystal_density():
    """
    Read density data from a .f98 file.
    """

    def __init__(self):

        pass

    def read_cry_irr_density(self, fort98_unit):
        """
        Read density profile data from a CRYSTAL .f98 file.
        Args:
            fort98_unit (str): File containing the formatted density matrix.
        Returns:
            None
        Note:
        This is a work in progress. If you are interested in this functionality,
        please open an Issue on GitHub.
        """
        self.file_name = fort98_unit
        self.is_irr = True

        import re
        import sys

        import numpy as np

        try:
            file = open(self.file_name, 'r')
            data = file.readlines()
            file.close()
        except:
            print('EXITING: a CRYSTAL .f98 file needs to be specified')
            sys.exit(1)

        self.all_file = data

        # the keyword BASATO appears twice, this is a check to see which basato
        # is being read
        basato1 = False

        for i, line in enumerate(data):

            if re.match(r'^LIMINF LIMTOL LIMPAR', line):
                inf_vec_len, tol_vec_len, par_vec_len = [
                    int(x) for x in data[i+1].split()]

            elif re.match(r'^INF', line):
                self.inf_vec = []
                inf_n_line = int(np.ceil(inf_vec_len/8))
                for j in range(inf_n_line):
                    self.inf_vec.extend([int(x) for x in data[i+1+j].split()])
                n_symmops = self.inf_vec[0]
                n_atoms = self.inf_vec[23]
                n_shells = self.inf_vec[19]
                '''if self.inf_vec[26] == 0:
                    self.spin_pol == False
                elif self.inf_vec[26] == 1:
                    self.spin_pol == True'''
                n_prim_gto = self.inf_vec[74]
                f_irr_len = (self.inf_vec[63]+1)*self.inf_vec[227]
                p_irr_len = (self.inf_vec[63]+1)*self.inf_vec[18]
                nnnc_len = self.inf_vec[190]
                la3_len = self.inf_vec[55]
                # n_symmops_noinv = inf_vec[1]

            elif re.match(r'^TOL', line):
                self.tol_vec = []
                tol_n_line = int(np.ceil(tol_vec_len/8))
                for j in range(tol_n_line):
                    self.tol_vec.extend([int(x) for x in data[i+1+j].split()])

            elif re.match(r'^PAR', line):
                self.par_vec = []
                par_n_line = int(np.ceil(par_vec_len/4))
                for j in range(par_n_line):
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # The line below fixes that issue
                    for item in range(0, int(len(data[i+1+j])/20)):
                        self.par_vec.append(
                            float(data[i+1+j][(item)*20:(item+1)*20]))

            elif re.match(r'^XYVGVE', line):
                # This vector contains the rotations, translation,
                # lattice vectors and transformation matrix from primitive to
                # crystallographic cell
                # Read all of it first and separate later
                xyvgve_n_line = int(np.ceil((n_symmops*12+18)/4))
                xyvgve_vec = []
                for j in range(xyvgve_n_line):
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # The line below fixes that issue
                    for item in range(0, int(len(data[i+1+j])/20)):
                        xyvgve_vec.append(
                            float(data[i+1+j][(item)*20:(item+1)*20]))
                # Now let's split the xyvgve_vec
                self.rotations_vec = xyvgve_vec[0:n_symmops*9]
                self.translations_vec = xyvgve_vec[n_symmops *
                                                   9:n_symmops*9+n_symmops*3]
                self.direct_lattice_vec = xyvgve_vec[n_symmops *
                                                     12:n_symmops*12+9]
                self.transf_matrix = xyvgve_vec[-9:]

            elif re.match(r'^BASATO', line):
                if basato1 == False:
                    basato_n_line = int(
                        np.ceil((n_atoms*4+n_shells*5+n_prim_gto*7)/4))
                    basato_vec = []
                    for j in range(basato_n_line):
                        # The negative elements appear connected to the previous one
                        # eg:  0.0000000000000E+00-1.0000000000000E+00
                        # The line below fixes that issue
                        for item in range(0, int(len(data[i+1+j])/20)):
                            basato_vec.append(
                                float(data[i+1+j][(item)*20:(item+1)*20]))
                    # Extract the iformation we need from basato

                    # Atom coordinates
                    self.atom_coord = []
                    for j in range(0, 3*n_atoms, 3):
                        self.atom_coord.append(
                            basato_vec[(n_atoms+j):(n_atoms+j+3)])
                    # self.atom_coord = np.array(self.atom_coord)

                    # Assign the shell to the atom
                    self.shell_coord = []
                    # The beginning of the part of BASATO I need here
                    init = 4*n_atoms + 2*n_shells
                    for j in range(0, 3*n_shells, 3):
                        self.shell_coord.append(
                            basato_vec[(init+j):(init+j+3)])
                    # self.shell_coord = np.array(self.shell_coord)

                    # Array that defines which atom a shell belongs to
                    self.shell_to_atom = []
                    for coord in self.shell_coord:
                        self.shell_to_atom.append(self.atom_coord.index(coord))
                    basato1 = True

                elif basato1 == True:
                    self.basato2 = []
                    j = i + 1
                    while 'SPINOR' not in data[j].split()[0]:
                        # The negative elements appear connected to the previous one
                        # eg:  0.0000000000000E+00-1.0000000000000E+00
                        # As opposite to the loops above where the float read was 20
                        # characters long, this ones are 21
                        # The line below fixes that issue
                        self.basato2.extend([int(x) for x in data[j].split()])
                        j += 1
                    self.atom_shell = self.basato2[-n_shells:]

            elif re.match(r'^SPINOR', line):
                self.f_irr = []
                self.charges = []
                # self.spin = [0]*(2*n_atoms) #tmp
                self.spin = []
                self.ghost = []
                n_ghost = 0
                n_spin_line = int(np.ceil((n_atoms*2)/8))
                n_basold = 0
                if 'BASOLD' in data[i+n_spin_line+1]:
                    n_basold = 9 + 3 * \
                        self.inf_vec[1] + n_shells + 3*n_atoms+3 * \
                        n_shells+self.inf_vec[4]+1+3*self.inf_vec[78]
                n_basold_line = int(np.ceil((n_basold)/4))
                n_charge_line = int(np.ceil(n_atoms/4))
                skip = n_spin_line + n_charge_line + 1 + n_basold_line
                for j in range(n_spin_line):
                    self.spin.extend([int(x) for x in data[i + j + 1].split()])
                if 'IGHOST' in data[i+n_spin_line+1]:
                    n_ghost = int(np.ceil((n_atoms)/8)) + 1
                    skip = skip + n_ghost
                    for j in range(n_ghost-1):
                        self.ghost.extend(
                            [float(x) for x in data[i + j + n_spin_line + 2].split()])
                f_irr_n_line = int(np.ceil(f_irr_len/4))
                for j in range(n_charge_line):
                    self.charges.extend(
                        [float(x) for x in data[i+j+n_spin_line+n_basold+n_ghost+1].split()])
                for j in range(f_irr_n_line):
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # As opposite to the loops above where the float read was 20
                    # characters long, this ones are 21
                    # The line below fixes that issue
                    for item in range(0, int(len(data[i+skip+j])/21)):
                        self.f_irr.append(
                            float(data[i+skip+j][(item)*21:(item+1)*21]))
                self.p_irr = []
                p_irr_n_line = int(np.ceil(p_irr_len/4))
                skip += 1
                for k in range(i+skip+j, i+skip+j+p_irr_n_line):
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # As opposite to the loops above where the float read was 20
                    # characters long, this ones are 21
                    # The line below fixes that issue
                    for item in range(0, int(len(data[k])/21)):
                        self.p_irr.append(
                            float(data[k][(item)*21:(item+1)*21]))

            elif re.match(r'^   NCF', line):
                # The ncf vector contains the pointers to the symmetry irerducible
                # shell couples la3, la4
                self.ncf = []
                j = i+1
                while 'NSTATG' not in data[j].split()[0]:
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # As opposite to the loops above where the float read was 20
                    # characters long, this ones are 21
                    # The line below fixes that issue
                    self.ncf.extend([int(x) for x in data[j].split()])
                    j += 1

            elif re.match(r'^NSTATG', line):
                # The nstatg vector contains the pointers to the starting point
                # of each couple set in P_irr and F_irr
                self.nstatg = []
                j = i+1
                while 'NSTAFG' not in data[j].split()[0]:
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # As opposite to the loops above where the float read was 20
                    # characters long, this ones are 21
                    # The line below fixes that issue
                    self.nstatg.extend([int(x) for x in data[j].split()])
                    j += 1

            elif re.match(r'^  NNNC', line):
                # The nnnc points the starting position in the P matrix for
                # each couple, and its size corresponds to the total number
                # of shell couple in the shell couple sets
                self.nnnc = []
                nnnc_n_line = int(np.ceil(nnnc_len/8))
                for j in range(nnnc_n_line):
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # As opposite to the loops above where the float read was 20
                    # characters long, this ones are 21
                    # The line below fixes that issue
                    for item in range(0, int(len(data[i+1+j])/10)):
                        self.nnnc.append(
                            int(data[i+1+j][(item)*10:(item+1)*10]))
            elif re.match(r'^   LA3', line):
                # The nnnc points the starting position in the P matrix for
                # each couple, and its size corresponds to the total number
                # of shell couple in the shell couple sets
                self.la3 = []
                # nnnc_n_line = int(np.ceil(nnnc_len/8))
                j = i+1
                while 'LA4' not in data[j].split()[0]:
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # As opposite to the loops above where the float read was 20
                    # characters long, this ones are 21
                    # The line below fixes that issue
                    self.la3.extend([int(x) for x in data[j].split()])
                    j += 1
            elif re.match(r'^   LA4', line):
                # The nnnc points the starting position in the P matrix for
                # each couple, and its size corresponds to the total number
                # of shell couple in the shell couple sets
                self.la4 = []
                # nnnc_n_line = int(np.ceil(nnnc_len/8))
                j = i+1
                while 'IROF' not in data[j].split()[0]:
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # As opposite to the loops above where the float read was 20
                    # characters long, this ones are 21
                    # The line below fixes that issue
                    self.la4.extend([int(x) for x in data[j].split()])
                    j += 1


def cry_combine_density(density1, density2, density3, new_density='new_density.f98', spin_pol=False):
    """
    Combine density matrix files.

    Args:
        density1 (str): First density matrix file.
        density2 (str): Second density matrix file.
        density3 (str): Density matrix file for the whole system.
        new_density (str, optional): Name of the new density matrix (default is "new_density.f98").
        spin_pol (bool, optional): Specifies if the system is spin polarized (default is False).
    Returns:
        None
    Note:
        This is a work in progress. If you are interested in this functionality,
        please open an Issue on GitHub.
    """

    import sys

    import numpy as np

    try:
        density1_data = Crystal_density(density1)  # substrate
        density2_data = Crystal_density(density2)  # adsorbate

        density3_data_obj = Crystal_density(density3)
        density3_data = density3_data_obj.all_file
    except:
        print('EXITING: a CRYSTAL .f98 file needs to be specified')
        sys.exit(1)

    # Find P_irr <-> atom correspondence
    fragment_1 = []
    fragment_2 = []

    for i, j in enumerate(density1_data.ncf):
        # density1_data.la3[j] is the shell number
        # density1_data.atom_shell[density1_data.la3[j]] is the atom position number (1-6)
        # density1_data.ghost[density1_data.atom_shell[density1_data.la3[j]]] is either 0 or atomic number depending on ghost or not
        # This tells me if the shell belongs to this fragment
        n_elements = density1_data.nstatg[i] - density1_data.nstatg[i - 1]

        if density1_data.ghost[density1_data.atom_shell[density1_data.la3[j-1]-1]-1] == 0 and \
                density1_data.ghost[density1_data.atom_shell[density1_data.la4[j-1]-1]-1] == 0:
            fragment_1.extend([True]*n_elements)
        else:
            fragment_1.extend([False]*n_elements)

        if density1_data.ghost[density1_data.atom_shell[density1_data.la3[j-1]-1]-1] != 0 and \
                density1_data.ghost[density1_data.atom_shell[density1_data.la4[j-1]-1]-1] != 0:
            fragment_2.extend([True]*n_elements)
        else:
            fragment_2.extend([False]*n_elements)

    if spin_pol == True:
        spin_p1 = fragment_1.copy()
        spin_p2 = fragment_2.copy()
        fragment_1.extend(spin_p1)
        fragment_2.extend(spin_p2)

    beginning = density3_data.index('SPINOR\n')
    end = density3_data.index('   NCF\n')
    sum_density = (np.array(density1_data.p_irr) +
                   np.array(density2_data.p_irr))/2
    sum_fock = np.array(density1_data.f_irr)+np.array(density2_data.f_irr)
    sum_charges = np.array(density1_data.charges) + \
        np.array(density2_data.charges)
    spinor = ['SPINOR\n']
    charges = []
    fock = []
    density = []
    new_fock = sum_fock  # TMP
    new_fock = [0] * len(density3_data_obj.f_irr)
    new_p = []

    for i in range(len(fragment_1)):
        if fragment_1[i] == True and fragment_2[i] == False:
            # new_fock.append(density1_data.f_irr[i])
            new_p.append(density1_data.p_irr[i])
        elif fragment_1[i] == False and fragment_2[i] == True:
            # new_fock.append(density2_data.f_irr[i])
            new_p.append(density2_data.p_irr[i])
        elif fragment_1[i] == False and fragment_2[i] == False:
            # new_fock.append(0.)
            new_p.append(sum_density[i])
            # new_p.append(0.)
            # new_p.append(density3_data_obj.p_irr[i])

    for i in range(0, len(density3_data_obj.spin), 8):
        spinor.append(' '.join([str(x)
                                for x in density3_data_obj.spin[i:i+8]])+'\n')
    for i in range(0, len(sum_charges), 4):
        charges.append(' '.join(["{:.13e}".format(x)
                                 for x in sum_charges[i:i+4]])+'\n')
    for i in range(0, len(new_fock), 4):
        fock.append(' '.join(["{:.13e}".format(x)
                              for x in new_fock[i:i+4]])+'\n')
    for i in range(0, len(new_p), 4):
        density.append(' '.join(["{:.13e}".format(x)
                                 for x in new_p[i:i+4]])+'\n')

    final_fort98 = density3_data[0:beginning] + \
        spinor+charges+fock+density+density3_data[end:]
    with open(new_density, 'w') as file:
        for line in final_fort98:
            file.writelines(line)


def write_cry_density(fort98_name, new_p, new_fort98):
    """
    Write the formatted density matrix.

    Args:
        fort98_name (str): Name of the previous density matrix file.
        new_p (list): New density matrix.
        new_fort98 (str): Name of the new density matrix file.

        Returns:
        None
    Note:
        This is a work in progress. If you are interested in this functionality,
        please open an Issue on GitHub.
    """

    import numpy as np

    file = open(fort98_name, 'r')
    data = file.readlines()
    file.close()

    density = Crystal_density(fort98_name)

    n_spin_line = int(np.ceil((density.inf_vec[23] * 2) / 8))
    n_charges_line = int(np.ceil((density.inf_vec[23]) / 4))
    beginning = data.index('SPINOR\n') + n_spin_line + n_charges_line + 1
    end = data.index('   NCF\n')

    new_fock_vect = [0] * len(density.f_irr)

    new_fock = []
    for i in range(0, len(new_fock_vect), 4):
        new_fock.append(' '.join(["{:.13e}".format(x)
                                  for x in new_fock_vect[i:i + 4]]) + '\n')

    new_density = []
    for i in range(0, len(new_p), 4):
        new_density.append(' '.join(["{:.13e}".format(x)
                                     for x in new_p[i:i+4]])+'\n')

    final_fort98 = data[0:beginning]+new_fock+new_density+data[end:]
    with open(new_fort98, 'w') as file:
        for line in final_fort98:
            file.writelines(line)


class External_unit:
    # WORK IN PROGRESS
    # This class will generate an object from the CRYSTAL external units like: IRSPEC.DAT
    # RAMSPEC.DAT tha can be plotted through the corresponding function in plot.py

    def __init__(self):

        pass

    def read_external_unit(self, external_unit):
        import os

        self.file_name = external_unit
        try:
            file = open(external_unit, 'r')
            self.data = file.readlines()
            file.close()

            # directory
            dir_name = os.path.split(external_unit)[0]
            self.abspath = os.path.join(dir_name)

            # title (named "title" only to distinguish from "file_name" which means another thing)
            self.title = os.path.split(external_unit)[1]

        except:
            raise FileNotFoundError(
                'EXITING: a CRYSTAL generated .DAT unit file needs to be specified')

    def read_cry_irspec(self, external_unit):
        """
        Reads the IRSPEC.DAT unit produced by CRYSTAL (IRSPEC keyword).

        Args:
            external_unit(str): Path of the IRSPEC.DAT file.

        Returns:
            self: Modified object to be read by plot_cry_irspec().
        """
        import numpy as np

        self.read_external_unit(external_unit)

        data = self.data

        for index, line in enumerate(data):
            data[index] = line.split()

        columns = len(data[0])
        no_points = len(data)

        if columns > 3:
            self.calculation = 'solid'
        else:
            self.calculation = 'molecule'

        irspec = np.zeros((no_points, columns))

        for i, line in enumerate(data):
            for j, element in enumerate(line):
                irspec[i, j] = float(element)

        self.irspec = irspec

        return self

    def read_cry_ramspec(self, external_unit):
        """
        Reads the RAMSPEC.DAT unit produced by CRYSTAL (RAMSPEC keyword).

        Args:
            external_unit(str): Path of the RAMSPEC.DAT file.

        Returns:
            self : Modified object to be read by plot_cry_ramspec().
        """
        import numpy as np

        self.read_external_unit(external_unit)

        data = self.data

        for index, line in enumerate(data):
            data[index] = line.split()

        columns = len(data[0])
        no_points = len(data)

        ramspec = np.zeros((no_points, columns))

        for i, line in enumerate(data):
            for j, element in enumerate(line):
                ramspec[i, j] = float(element)

        self.ramspec = ramspec

        return self

    def read_phonon_band(self, band_file, output=None):
        """
        Generate bands object from CRYSTAL PHONBAND.DAT or fort.25 file.
        Energy unit: eV. E Fermi is aligned to 0.

        Args:
            band_file (str): Name of PHONBAND.DAT or fort.25 file
            output (str): Properties output file (.outp or .out). For 3D, k
                coordinates and geometry information.

        Returns:
            self.bands (BandsBASE): A Bands base object.
        """
        from CRYSTALClear.base.propout import BandsBASE, OutBASE

        self.read_external_unit(band_file)
        if '-%-' in self.data[0]:  # fort.25 file format
            self.bands = BandsBASE.f25_parser(self.data)
        else:  # BAND.DAT file format
            self.bands = BandsBASE.BAND_parser(self.data)

        if output != None:
            self.bands.geometry = OutBASE.get_geometry(output)
            self.bands.tick_pos3d, self.bands.k_point_pos3d = OutBASE.get_3dkcoord(
                output)
            self.bands.reciprocal_latt = self.bands.geometry.lattice.reciprocal_lattice.matrix

        return self.bands

    def read_phonon_dos(self, external_unit):
        """
        Generate doss object from CRYSTAL PHONONDOSS.DAT or fort.25 file.
        Energy unit: eV. E Fermi is aligned to 0.

        Args:
            properties_output (str): File name.

        Returns:
            self.doss (DOSBASE): A DOS base object.
        """
        from CRYSTALClear.base.propout import DOSBASE

        self.read_external_unit(external_unit)
        if '-%-' in self.data[0]:  # fort.25 file format
            self.doss = DOSBASE.f25_parser(self.data)
        else:  # DOSS.DAT file format
            self.doss = DOSBASE.DOSS_parser(self.data)

        return self.doss
