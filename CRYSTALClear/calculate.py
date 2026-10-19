#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def cry_ads_energy(e_full_system, e_substrate, e_adsorbate):
    """
    Calculate the adsorption energy of a system.

    Args:
        e_full_system (float): Total energy of the full system.
        e_substrate (float): Energy of the substrate.
        e_adsorbate (float): Energy of the adsorbate.

    Returns:
        float: Adsorption energy calculated as the difference between the
               total energy of the full system and the sum of the energies
               of the substrate and adsorbate.
    """
    return e_full_system-(e_substrate+e_adsorbate)


def cry_shrink(structure, spacing=0.2):
    """
    Determine the number of unit cells needed to achieve a desired spacing.

    Args:
        structure (pymatgen.core.structure.Structure): The input structure.
        spacing (float): The desired spacing between unit cells. Default is 0.2.

    Returns:
        int: The number of unit cells (rounded up) needed to achieve the desired spacing.
    """

    import numpy as np

    vector = np.average(structure.lattice.reciprocal_lattice.lengths)

    return int(np.ceil(vector/spacing))


def cry_track_modes(cos, labels=None):
    """
    Follow the normal modes along a series of calculations of the same
    system, e.g. at increasing pressure, where modes change their number as
    their frequencies cross.

    Modes of consecutive calculations are paired one to one so that the sum
    of the overlaps |e_i . e_j*|^2 of their eigenvectors, as normalised by
    get_phonon(), is the largest, and the pairs
    are chained from the first calculation. The series should be fine enough
    for neighbouring calculations to be alike: a low overlap flags a mode
    that mixes with another one, near a crossing. Degenerate modes are paired
    too, but their overlaps are only meaningful summed over the degenerate set.

    Args:
        cos (list[crystal_io.Crystal_output]): Crystal output objects in the order of the series, on each of which get_phonon(read_eigvt=True) has already been called.
        labels (list, optional): Label of each calculation, e.g. the pressures. Default is None, i.e. 0, 1, 2...

    Returns:
        pandas.DataFrame: Mode number in each calculation (columns) of every mode of the first one (index).
        pandas.DataFrame: Overlap of the step that reached each calculation, 1 for the first one.
    """
    import numpy as np
    import pandas as pd
    from scipy.optimize import linear_sum_assignment

    if labels is None:
        labels = list(range(len(cos)))
    if len(labels) != len(cos):
        raise ValueError('One label per calculation is needed.')

    # Gamma eigenvectors, one row per mode, blanked modes set to zero -->
    eigvts = []
    for co in cos:
        if len(getattr(co, 'eigenvector', [])) == 0:
            raise AttributeError(
                "No eigenvectors found: call Crystal_output.get_phonon(read_eigvt=True) first.")
        eigvt = np.array(co.eigenvector[0])
        eigvts.append(np.nan_to_num(eigvt.reshape(eigvt.shape[0], -1)))
    if any(eigvt.shape != eigvts[0].shape for eigvt in eigvts):
        raise ValueError(
            'The calculations do not have the same modes and atoms.')
    # <--

    # Pairing of consecutive calculations, chained from the first one -->
    nmode = eigvts[0].shape[0]
    chain = np.arange(nmode)
    modes = [chain + 1]
    overlaps = [np.ones(nmode)]
    for previous, current in zip(eigvts, eigvts[1:]):
        overlap = np.abs(previous @ current.conj().T) ** 2
        _, match = linear_sum_assignment(-overlap)
        overlaps.append(overlap[chain, match[chain]])
        chain = match[chain]
        modes.append(chain + 1)
    # <--

    index = pd.Index(np.arange(1, nmode + 1), name='mode')

    return (pd.DataFrame(np.array(modes).T, index=index, columns=labels),
            pd.DataFrame(np.array(overlaps).T, index=index, columns=labels))
