#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and methods to phrase 'crystal' output file.
"""


class GeomBASE():
    """
    A container of basic methods for SCF geometry.
    """
    @classmethod
    def read_geom(cls, data, countline):
        """
        Read lattice from 'A B C ALPHA BETA GAMMA' block, periodic boundary
        condition and atom positions from 'ATOMS IN THE ASYMMETRIC UNIT'
        block. It terminates at the first empty line after that block.

        Args:
            data (list[str]): output file read by readlines()
            countline (int): The starting line number

        Returns:
            countline (int): Line number of output file.
            struc (Pymatgen Structure | Molecule)
        """
        import re

        import numpy as np
        from pymatgen.core.lattice import Lattice
        from pymatgen.core.structure import Molecule, Structure

        pbc = {0: (False, False, False),
               1: (True, False, False),
               2: (True, True, False),
               3: (True, True, True)}

        species = []
        coords = []
        latt_mx = np.eye(3)
        while countline < len(data):
            line = data[countline]
            if re.match(r'^\s+A\s+B\s+C\s+ALPHA\s+BETA\s+GAMMA\s+$', line):
                countline += 1
                line_data = data[countline].strip().split()
                ndimen = 3 - len(re.findall('\(ANGSTROM\)', data[countline+3]))
                latt = Lattice.from_parameters(
                    a=float(line_data[0]), b=float(line_data[1]),
                    c=float(line_data[2]), alpha=float(line_data[3]),
                    beta=float(line_data[4]), gamma=float(line_data[5]),
                    pbc=pbc[ndimen]
                )
                latt_mx = latt.matrix
                countline += 1
            # Atom coordinates
            if re.match(r'^\s*ATOMS IN THE ASYMMETRIC UNIT', line):
                countline += 1
                # For molecules - no PRIMITIVE CELL line
                ndimen = 3 - len(re.findall('\(ANGSTROM\)', data[countline]))

                countline += 2
                line = data[countline].strip()
                coords = []
                species = []
                while line != '':
                    line_data = line.split()
                    coords.append(line_data[-3:])
                    species.append(line_data[3].capitalize())
                    countline += 1
                    line = data[countline].strip()
                countline += 1
                coords = np.array(coords, dtype=float)
                if ndimen != 0:
                    # to cartesian coords
                    coords[:, 0:ndimen] = coords[:,
                                                 0:ndimen] @ latt_mx[0:ndimen, 0:ndimen]
                    struc = Structure(lattice=latt, species=species,
                                      coords=coords, coords_are_cartesian=True)
                else:
                    struc = Molecule(species=species, coords=coords)
                break
            else:
                countline += 1

        if len(species) == 0 or len(coords) == 0:
            raise Exception('Geometry information not found.')

        return countline, struc

    @classmethod
    def read_conv_z(cls, data, countline):
        """
        Read conventional atom numbers.

        Args:
            data (list[str]): output file read by readlines()
            countline (int): The starting line number

        Returns:
            countline (int): Line number of output file.
            conv_z (array): Int array of conventional atom numbers
        """
        import re

        import numpy as np

        conv_z = []
        while countline < len(data):
            line = data[countline]
            if re.match(r'^\s*ATOMS IN THE ASYMMETRIC UNIT', line):
                countline += 3
                line = data[countline].strip()
                while line != '':
                    line_data = line.split()
                    conv_z.append(line_data[2])
                    countline += 1
                    line = data[countline].strip()
                countline += 1
                conv_z = np.array(conv_z, dtype=int)
                break
            else:
                countline += 1

        if len(conv_z) == 0:
            raise Exception('Geometry information not found.')

        return countline, conv_z


class SCFBASE():
    """
    A container of basic methods for SCF loop.
    """
    @classmethod
    def read_convergence(cls, data, countline):
        """
        Read SCF convergence.

        Returns:
            countline (int)
            ncyc (int): Number of cycles
            endflag (str): 'terminated', 'converged', 'too many cycles' and 'unknown'
            e (array): nCYC\*1 array of SCF energy. Unit: eV
            de (array): nCYC\*1 array of SCF energy difference. Unit: eV
        """
        import re
        import warnings

        import numpy as np

        from CRYSTALClear.units import H_to_eV

        e = []
        de = []
        endflag = 'terminated'
        while countline < len(data):
            line = data[countline]
            if re.match(r'^\s*CYC', line):
                line_data = line.strip().split()
                e.append(line_data[3])
                de.append(line_data[5])
                countline += 1
            elif re.match(r'^\s*== SCF ENDED', line):
                if re.search('CONVERGENCE', line):
                    endflag = 'converged'
                elif re.search('TOO MANY CYCLES', line):
                    endflag = 'too many cycles'
                else:
                    warnings.warn('Unknown termination. Check line {}.'.format(countline + 1),
                                  stacklevel=3)
                    endflag = 'unknown'
                break
            else:
                countline += 1

        if endflag != 'converged':
            warnings.warn(
                'SCF convergence not achieved or missing.', stacklevel=3)

        ncyc = len(e)
        e = np.array(e, dtype=float)
        de = np.array(de, dtype=float)

        return countline, ncyc, endflag, H_to_eV(e), H_to_eV(de)

    @classmethod
    def read_fermi_energy(cls, data, countline, history=False):
        """
        Read Fermi energy.

        Args:
            history (bool): Whether to read e fermi of all steps

        Returns:
            countline (int)
            spin (bool): Whether the system is spin-polarised.
            efermi (float | array): Fermi energy. Unit: eV
        """
        import re
        import warnings

        import numpy as np

        from CRYSTALClear.units import H_to_eV

        if history == True:
            efermi = []
        else:
            efermi = None

        spin = False
        tmp_efermi = []
        while countline >= 0:
            line = data[countline]
            # Metal spin or no spin
            if re.match(r'^ POSSIBLY CONDUCTING STATE - EFERMI', line):
                line_data = line.strip().split()
                if history == False:
                    efermi = float(line_data[5])
                    break
                else:
                    efermi.append(line_data[5])
                    countline -= 1  # Note the reversed sequence here
            # spin flag
            elif re.match(r'^\s*SUMMED SPIN DENSITY', line):
                spin = True
                tmp_efermi = []
                countline -= 1
            # insulator, use top of valence bands
            elif re.match(r'^\s*TOP OF VALENCE BANDS', line):
                line_data = line.strip().split()
                tmp_efermi.append(float(line_data[10]))
                if spin == True and len(tmp_efermi) == 2:
                    if history == False:
                        # Note the reversed sequence here
                        efermi = [tmp_efermi[1], tmp_efermi[0]]
                        break
                    else:
                        efermi.append([tmp_efermi[1], tmp_efermi[0]])
                        tmp_efermi = []
                elif spin == False:
                    if history == False:
                        efermi = tmp_efermi[0]
                        break
                    else:
                        efermi.append(tmp_efermi[0])
                        tmp_efermi = []

                countline -= 1
            # Before the SCF block
            elif re.match(r'^\s*A+$', line):
                break
            else:
                countline -= 1

        if spin == True or history == True:
            efermi = np.array(efermi, dtype=float)
        if history == True:
            # Note the reversed sequence here
            efermi = efermi[::-1]

        return countline, spin, H_to_eV(efermi)

    @classmethod
    def read_band_gap(cls, data, countline, history=False):
        """
        Read band gap.

        Args:
            history (bool): Whether to read band gap of all steps

        Returns:
            countline (int)
            spin (bool): Whether the system is spin-polarised.
            gap (float | array): Band gap. Unit: eV
        """
        import re
        import warnings

        import numpy as np

        if history == True:
            gap = []
        else:
            gap = None

        spin = False
        tmp_gap = []
        while countline >= 0:
            line = data[countline]
            # Metal spin or no spin
            if re.match(r'^ POSSIBLY CONDUCTING STATE', line):
                warnings.warn('Conducting state is identified.', stacklevel=3)
                if history == False:
                    gap = 0.
                    break
                else:
                    gap.append(0.)
                    countline -= 1
            # spin flag
            elif re.match(r'^\s*SUMMED SPIN DENSITY', line):
                spin = True
                tmp_gap = []
                countline -= 1
            # insulator
            elif re.match(r'^\s*.*DIRECT ENERGY BAND GAP', line):
                line_data = line.strip().split()
                tmp_gap.append(float(line_data[4]))
                if spin == True and len(tmp_gap) == 2:
                    if history == False:
                        # Note the reversed sequence here
                        gap = [tmp_gap[1], tmp_gap[0]]
                        break
                    else:
                        gap.append([tmp_gap[1], tmp_gap[0]])
                        tmp_gap = []
                elif spin == False:
                    if history == False:
                        gap = tmp_gap[0]
                        break
                    else:
                        gap.append(tmp_gap[0])
                        tmp_gap = []

                countline -= 1
            # Before the SCF block
            elif re.match(r'^\s*A+$', line):
                break
            else:
                countline -= 1

        if spin == True or history == True:
            gap = np.array(gap, dtype=float)
        if history == True:
            # Note the reversed sequence here
            gap = gap[::-1]

        return countline, spin, gap


class OptBASE():
    """
    A container of basic methods for Opt loop.
    """
    @classmethod
    def read_optblock(cls, data, countline):
        """
        Read optimisation blocks.

        Returns:
            countline (int): Line number of output file.
            ncyc (int): Number of cycles
            endflag (str): 'terminated', 'converged', 'failed' and 'unknown'
            e (array): nCYC\*1 array of total energy. Unit: eV
            de (array): nCYC\*1 array of total energy difference. Unit: eV
            struc (list[Structure]): nCYC\*1 list of pymatgen structures
            maxg (array): nCYC\*1 array of max energy gradient convergence.
                Unit: Hartree / Bohr
            rmsg (array): nCYC\*1 array of RMS energy gradient convergence.
                Unit: Hartree / Bohr
            maxd (array): nCYC\*1 array of max displacement convergence.
                Unit: Bohr
            rmsd (array): nCYC\*1 array of RMS displacement convergence.
                Unit: Bohr
        """
        import re
        import warnings

        import numpy as np

        from CRYSTALClear.base.crysout import GeomBASE
        from CRYSTALClear.units import H_to_eV

        e = []
        de = []
        struc = []
        maxg = []
        rmsg = []
        # Displacement: the initial step is 0
        maxd = []
        rmsd = []

        endflag = 'terminated'
        while countline < len(data):
            line = data[countline]
            if re.match(r'^\s+TOTAL ENERGY\(DFT\)\(AU\)\(.+DE \(AU\)', line):
                line_data = line.strip().split()
                e.append(line_data[3])
                de.append(line_data[6])
                countline += 1
            # Use atom coords to read molecule geometries. Go 4 lines up for periodic systems
            elif re.match(r'^\s*ATOMS IN THE ASYMMETRIC UNIT', line):
                output = GeomBASE.read_geom(data, countline - 4)
                countline = output[0]
                struc.append(output[1])
            elif re.match(r'^\s+MAX GRADIENT', line):
                line_data = line.strip().split()
                maxg.append(line_data[2])
                countline += 1
            elif re.match(r'^\s+RMS GRADIENT', line):
                line_data = line.strip().split()
                rmsg.append(line_data[2])
                countline += 1
            elif re.match(r'^\s+MAX DISPLAC\.', line):
                line_data = line.strip().split()
                maxd.append(line_data[2])
                countline += 1
            elif re.match(r'^\s+RMS DISPLAC\.', line):
                line_data = line.strip().split()
                rmsd.append(line_data[2])
                countline += 1
            elif re.match(r'\s*\*\s*OPT END', line):
                if re.search('CONVERGED', line):
                    endflag = 'converged'
                elif re.search('FAILED', line):
                    endflag = 'failed'
                else:
                    warnings.warn('Unknown termination. Check line {}.'.format(countline + 1),
                                  stacklevel=3)
                    endflag = 'unknown'
                break
            else:
                countline += 1

        if endflag != 'converged':
            warnings.warn(
                'Optimisation convergence not achieved.', stacklevel=3)

        e = np.array(e, dtype=float)
        de = np.array(de, dtype=float)
        maxg = np.array(maxg, dtype=float)
        rmsg = np.array(rmsg, dtype=float)
        maxd = np.array(maxd, dtype=float)
        rmsd = np.array(rmsd, dtype=float)

        # In case that a calculation is restarted and stopped at the initial step
        ncyc = len(rmsg)

        return countline, ncyc, endflag, H_to_eV(e), H_to_eV(de), struc, maxg, rmsg, maxd, rmsd


class PhononBASE():
    """
    A container of basic methods for phonon information.
    """
    @classmethod
    def readmode_basic(cls, data, countline):
        """
        Read basic frequency information.

        Returns:
            countline (int): Line number of output file.
            frequency (array[float]): nmode \* 1
            intens (array[float]): nmode \* 1
            IR (array[bool]): nmode \* 1
            Raman (array[bool]): nmode \* 1
        """
        import re

        import numpy as np

        frequency = []
        intens = []
        IR = []
        Raman = []

        has_spec = False
        while countline < len(data):
            line = data[countline]
            if re.match(r'^\s+MODES\s+EIGV\s+FREQUENCIES\s+IRREP', line):
                countline += 2
                continue
            elif re.match(r'^\s+\(HARTREE\*\*2\)\s+\(CM\*\*\-1\)', line):
                countline += 1
                continue
            elif line.strip() == '':
                countline += 1
                break
            line_data = line.split()
            nm_a = int(line_data[0].strip('-'))
            nm_b = int(line_data[1])
            freq = float(line_data[4])
            # IR/Raman analysis, closed by default in dispersion calcs
            if 'A' in line or 'I' in line:
                has_spec = True
                intensity = float(line_data[-2].strip(')')) / (nm_b - nm_a + 1)
                ir = line_data[-4] == 'A'
                ram = line_data[-1] == 'A'

            for mode in range(nm_a, nm_b + 1):
                frequency.append(freq)
                if has_spec:
                    intens.append(intensity)
                    IR.append(ir)
                    Raman.append(ram)
            countline += 1

        return countline, frequency, intens, IR, Raman

    @classmethod
    def readmode_eigenvector(cls, data, countline):
        """
        Get mode eigenvectors.

        Returns:
            countline (int): Line number of output file.
            eigvt (array[float]): nmode\*natom\*3 array.
        """
        import re

        import numpy as np

        # Read the eigenvector region as its original shape
        total_data = []
        nmode = 0
        while countline < len(data):
            line = data[countline]
            # Found a block
            if re.match(r'^\s+FREQ\(CM\*\*\-1\)', line):
                countline += 2
                block_data = []
                while data[countline].strip() != '':
                    # Trim annotation part (12 characters)
                    line_data = re.findall(r'\-*[\d\.]+[E\d\-\+]*',
                                           data[countline][13:])
                    if line_data:
                        block_data.append(line_data)
                    countline += 1

                nmode += len(line_data)
                total_data.append(block_data)
                countline += 1
                continue
            # Empty lines
            elif line.strip() == '':
                countline += 1
                continue
            # Other lines
            else:
                break

        eigvt = np.hstack([i for i in total_data])  # (3*natom) * nmode
        eigvt = np.array(eigvt, dtype=float)
        # 1st dimension, nmode
        eigvt = np.transpose(eigvt)  # nmode * (3*natom)
        # 2nd dimension, natom
        natom = int(nmode / 3)
        eigvt = np.reshape(eigvt, [nmode, natom, 3], order='C')

        return countline, eigvt

    @classmethod
    def classical_amplitude(cls, struc, freq):
        """
        Get classical amplitude of phonon modes

        .. math::
            x = \sqrt{\frac{\hbar}{\mu\omega}}

        Args:
            struc (Structure): Pymatgen structure
            freq (float | array): Frequency. Unit: THz
        Returns:
            classic_a (array): nfreq\*3natom\*3natom array, or 3natom\*3natom
                if ``freq`` is float. The diagonal matrix of classical amplitude.
        """
        import numpy as np

        from CRYSTALClear.units import amu_to_me, thz_to_hartree

        if type(freq) == float:
            freq = np.array([freq])

        natom = struc.num_sites
        nmode = int(3*natom)
        mass_rev = np.zeros([nmode, nmode])  # In AU
        for i in range(0, nmode, 3):
            atid = i // 3
            atmass = amu_to_me(float(struc.species[atomid].atomic_mass))
            mass_rev[i, i] = 1 / np.sqrt(atmass)
            mass_rev[i+1, i+1] = 1 / np.sqrt(atmass)
            mass_rev[i+2, i+2] = 1 / np.sqrt(atmass)

        freq_rev = 1 / np.sqrt(thz_to_hartree(freq))

        classic_a = []
        for f in freq_rev:
            classic_a.append(f * mass_rev)
        classic_a = np.array(classic_a)
        if len(freq) == 1:
            classic_a = classic_a[0]

        return classic_a

    @classmethod
    def normalize_eigenvector(cls, eigvt, amplitude=1.):
        """
        Normalize the mode of eigenvectors.

        Args:
            eigvt (array[complex]): nmode\*natom\*3 array.
            amplitude (float): Amplitude of normalization
        Returns:
            eigvt (array[complex]): Normalized eigenvector.
        """
        import numpy as np

        for idx_m, m in enumerate(eigvt):
            eigvt[idx_m] = eigvt[idx_m] / np.linalg.norm(m) * amplitude

        return eigvt

    @classmethod
    def clean_q_overlap(cls, crysout, threshold):
        """
        Remove the repeated q points at both ends of line segment when
        dispersion is read. The weight of q points will be updated here.

        Args:
            crysout (Crystal_output): :code:`CRYSTALClear.crystal_io.Crystal_output` object
            threshold (float): The q point overlap threshold.
        """
        import warnings

        import numpy as np

        if crysout.nqpoint <= 1:
            return crysout

        overlap = []
        for idx_q1 in range(crysout.nqpoint - 1):
            qvec1 = crysout.qpoint[idx_q1][0]
            for idx_q2 in range(idx_q1 + 1, crysout.nqpoint):  # q1 < q2
                qvec2 = crysout.qpoint[idx_q2][0]
                if np.linalg.norm(qvec1 - qvec2) < threshold:
                    warnings.warn('Overlap of q points is detected between q points {:3d} and {:3d}'.format(idx_q1, idx_q2),
                                  stacklevel=2)
                    overlap.append([idx_q1, idx_q2])

        overlap = np.array(overlap, dtype=int)
        weight = np.array([i[1] for i in crysout.qpoint])
        for idx_q in range(crysout.nqpoint, 1, -1):
            if len(overlap) == 0:  # No overlap
                break
            for idx_o in np.where(overlap[:, 1] == idx_q)[0]:
                crysout.nqpoint -= 1
                del crysout.qpoint[overlap[idx_o, 1]]
                weight = np.delete(weight, overlap[idx_o, 1])
                crysout.nmode = np.delete(crysout.nmode, overlap[idx_o, 1])
                crysout.frequency = np.delete(
                    crysout.frequency, overlap[idx_o, 1], axis=0)
                if len(crysout.intens) != 0:
                    crysout.intens = np.delete(
                        crysout.intens, overlap[idx_o, 1], axis=0)
                    crysout.IR = np.delete(
                        crysout.IR, overlap[idx_o, 1], axis=0)
                    crysout.Raman = np.delete(
                        crysout.Raman, overlap[idx_o, 1], axis=0)
                if len(crysout.eigenvector) != 0:
                    crysout.eigenvector = np.delete(
                        crysout.eigenvector, overlap[idx_o, 1], axis=0)
                break  # Avoid repeatly deleting

        # Update weight
        weight = weight * 1 / np.sum(weight)
        for i in range(crysout.nqpoint):
            crysout.qpoint[i][1] = weight[i]

        return crysout

    @classmethod
    def clean_imaginary(cls, crysout, threshold):
        """
        Substitute imaginary modes and corresponding eigenvectors with numpy
        NaN format and print warning message.

        Args:
            crysout (Crystal_output): :code:`CRYSTALClear.crystal_io.Crystal_output` object
            threshold (float): The threshold to identify a phonon mode as negative.
        """
        import warnings

        import numpy as np

        for q, freq in enumerate(crysout.frequency):
            if np.isnan(freq[0]) or freq[0] > threshold:
                continue

            warnings.warn('Negative frequencies detected.\n  Calculated thermodynamics might be inaccurate. Negative frequencies will be substituted by NaN.',
                          stacklevel=2)

            neg_rank = np.where(freq <= threshold)[0]
            crysout.frequency[q, neg_rank] = np.nan

            if len(crysout.eigenvector) != 0:
                natom = int(crysout.nmode[q] / 3)
                nan_eigvt = np.full([natom, 3], np.nan)
                crysout.eigenvector[q, neg_rank] = nan_eigvt

            if len(crysout.intens) != 0:
                crysout.intens[q, neg_rank] = np.nan
                for n in neg_rank:
                    crysout.IR[q][n] = False
                    crysout.Raman[q][n] = False

        return crysout
