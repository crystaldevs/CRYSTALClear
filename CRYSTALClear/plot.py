#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions to visualize physical properties computed with CRYSTAL .
"""

##############################################################################
#                                                                            #
#                       ELECTRONIC STRUCTURE AND PHONONS                     #
#                                                                            #
##############################################################################

# ------------------------------------ECHG------------------------------------#


def plot_dens_ECHG(obj_echg, levels=150, xticks=5,
                   yticks=5, cmap_max=None, cmap_min=None):
    """
    Plots the 2D ECHG density map from a fort.25 file.

    Args: 
        obj_echg (crystal_io.Properties_output): Properties output object.
        levels (int or array-like, optional): Determines the number and positions of the contour lines/regions. Default is 150.
        xticks (int, optional): Number of ticks in the x direction. Default is 5.
        yticks (int, optional): Number of ticks in the y direction. Default is 5.
        cmap_max(float, optional): Maximum value used for the colormap. Default is None.
        cmap_min(float, optional): Minimun value used for the colormap. Default is None.

   Returns:
        matplotlib.figure.Figure
        matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    vector_ab = obj_echg.a - obj_echg.b
    lenght_ab = np.sqrt(vector_ab[0]**2 + vector_ab[1]**2 + vector_ab[2]**2)
    vector_cb = obj_echg.c - obj_echg.b
    lenght_cb = np.sqrt(vector_cb[0]**2 + vector_cb[1]**2 + vector_cb[2]**2)
    points_ab = len(obj_echg.density_map[:, 0])
    points_cb = len(obj_echg.density_map[0, :])

    mesh_x = np.zeros((points_ab, points_cb), dtype=float)
    mesh_y = np.zeros((points_ab, points_cb), dtype=float)
    for i in range(0, points_ab):
        for j in range(0, points_cb):
            mesh_y[i, j] = (((lenght_ab / points_ab) * i) *
                            np.sqrt(1 - (obj_echg.cosxy**2)))
            mesh_x[i, j] = ((lenght_cb / points_cb) * j) + \
                (((lenght_ab / points_ab) * i) * obj_echg.cosxy)

    dens = obj_echg.density_map * (1.88973**2)  # Bohr to Angstrom conversion

    if cmap_max is None:
        max_data = np.amax(dens)
    else:
        max_data = cmap_max

    if cmap_min is None:
        min_data = np.amin(dens)
    else:
        min_data = cmap_min

    fig, ax = plt.subplots()
    im = ax.contourf(mesh_x, mesh_y, dens, levels, cmap='gnuplot')
    divider = make_axes_locatable(ax)
    im.set_clim(vmin=min_data, vmax=max_data)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')
    ax.set_xlabel('$\AA$')
    ax.set_xticks(np.linspace(0, lenght_cb, xticks).tolist())
    ax.set_yticks(np.linspace(0, lenght_ab, yticks).tolist())
    ax.set_ylabel('$\AA$')
    ax.set_aspect(1.0)
    ax.set_xlim(np.amin(mesh_x), np.amax(mesh_x))
    ax.set_ylim(0, np.amax(mesh_y) * np.sqrt(1 - (obj_echg.cosxy**2)))

    return fig, ax

# ----------------------------------SPIN CURRENTS------------------------------#


def plot_vecfield2D_m(header, dens, quivscale, levels=150):
    """
    Plots the 2D magnetization vector field.

    Args:
        header (list): List containing information about the fort.25 header.
        dens (numpy.ndarray): Array containing the vector field data.
        quivscale (float): Scale factor for the quiver plot.
        levels (int or array-like, optional): Determines the number and positions of the contour lines/regions.

    Returns:
        matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    import numpy as np

    for i in range(2, 20):
        if (header[0] % i) == 0:
            nrow_split = int(header[0]/i)

    for i in range(2, 20):
        if (header[1] % i) == 0:
            ncol_split = int(header[1]/i)

    # initializes the arrays
    projx_m = np.zeros((nrow_split, ncol_split), dtype=float)
    projy_m = np.zeros((nrow_split, ncol_split), dtype=float)
    mod_m = np.zeros((header[0], header[1]), dtype=float)

    # Generates the meshgrid
    mesh_x = np.zeros((header[0], header[1]), dtype=float)
    mesh_y = np.zeros((header[0], header[1]), dtype=float)
    mesh_projx = np.zeros((nrow_split, ncol_split), dtype=float)
    mesh_projy = np.zeros((nrow_split, ncol_split), dtype=float)

    T = np.array([[np.sqrt(1 - header[4]**2)*header[2], 0],
                  [header[4]*header[2], header[3]]])  # Change of basis matrix

    for i in range(0, header[0]):
        for j in range(0, header[1]):
            mesh_x[i, j] = np.dot(T, np.array([i, j]))[0]
            mesh_y[i, j] = np.dot(T, np.array([i, j]))[1]

    mesh_x = mesh_x * 0.529177249  # Bohr to Angstrom conversion
    mesh_y = mesh_y * 0.529177249

    r = 0
    s = 0
    step_nrow = int(header[0]/nrow_split)
    step_ncol = int(header[1]/ncol_split)
    for i in range(0, nrow_split):
        for j in range(0, ncol_split):
            mesh_projx[i, j] = mesh_x[r, s]
            mesh_projy[i, j] = mesh_y[r, s]
            s += step_ncol
        s = 0
        r += step_nrow

    # Creates the orthogonal vectorial basis for the projections
    CB = header[7] - header[6]
    BA = header[6] - header[5]
    if header[4] == 0:
        mod_BA = np.sqrt(BA[0]**2 + BA[1]**2 + BA[2]**2)
        mod_CB = np.sqrt(CB[0]**2 + CB[1]**2 + CB[2]**2)
        v1 = BA/mod_BA
        v2 = CB/mod_CB
    else:
        BA = BA - ((np.dot(BA, CB))/(np.dot(BA, BA)))*CB
        mod_BA = np.sqrt(BA[0]**2 + BA[1]**2 + BA[2]**2)
        mod_CB = np.sqrt(CB[0]**2 + CB[1]**2 + CB[2]**2)
        v1 = BA/mod_BA
        v2 = CB/mod_CB

    # Calculates the modulus of the vectorial field
    for i in range(0, header[0]):
        for j in range(0, header[1]):
            mod_m[i, j] = np.sqrt((dens[i, j, 0]**2) +
                                  (dens[i, j, 1]**2)+(dens[i, j, 2]**2))

    # Calculates the projections of the vectors in the ABC plane
    ABC_normal = np.cross(BA, CB)
    mod_normal = np.sqrt(ABC_normal[0]**2 +
                         ABC_normal[1]**2 + ABC_normal[2]**2)

    for i in range(0, nrow_split):
        for j in range(0, ncol_split):
            projx_m[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                                                np.cross(np.array([dens[int(i*step_nrow), int(j*step_ncol), 0], dens[int(i*step_nrow), int(j*step_ncol), 1],
                                                                                   dens[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v2)
            projy_m[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                                                np.cross(np.array([dens[int(i*step_nrow), int(j*step_ncol, 0)], dens[int(i*step_nrow), int(j*step_ncol), 1],
                                                                                   dens[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v1)

    # Plotting

    m = plt.figure()
    m = plt.contourf(mesh_x, mesh_y, mod_m, levels, cmap='cool')
    m = plt.colorbar(mappable=m)
    m = plt.quiver(mesh_projx, mesh_projy, projx_m, projy_m, scale=quivscale)
    m = plt.xlabel('$\AA$')
    m = plt.ylabel('$\AA$')

    return m


def plot_vecfield2D_j(header, dens, quivscale, levels=150):
    """
    Plots the 2D vector field of the spin current.

    Args:
        header (list): List containing information about the fort.25 header.
        dens (numpy.ndarray): Array representing the vector field.
        quivscale (float): Scale factor for the quiver plot.
        levels (int or array-like, optional): Determines the number and positions of the contour lines/regions.

    Returns:
        matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    import numpy as np

    for i in range(2, 20):
        if (header[0] % i) == 0:
            nrow_split = int(header[0]/i)

    for i in range(2, 20):
        if (header[1] % i) == 0:
            ncol_split = int(header[1]/i)

    # initializes the arrays
    projx_j = np.zeros((nrow_split, ncol_split), dtype=float)
    projy_j = np.zeros((nrow_split, ncol_split), dtype=float)
    mod_j = np.zeros((header[0], header[1]), dtype=float)

    # Generates the meshgrid
    mesh_x = np.zeros((header[0], header[1]), dtype=float)
    mesh_y = np.zeros((header[0], header[1]), dtype=float)
    mesh_projx = np.zeros((nrow_split, ncol_split), dtype=float)
    mesh_projy = np.zeros((nrow_split, ncol_split), dtype=float)

    T = np.array([[np.sqrt(1 - header[4]**2)*header[2], 0],
                  [header[4]*header[2], header[3]]])  # Change of basis matrix

    for i in range(0, header[0]):
        for j in range(0, header[1]):
            mesh_x[i, j] = np.dot(T, np.array([i, j]))[0]
            mesh_y[i, j] = np.dot(T, np.array([i, j]))[1]

    mesh_x = mesh_x * 0.529177249  # Bohr to Angstrom conversion
    mesh_y = mesh_y * 0.529177249

    r = 0
    s = 0
    step_nrow = int(header[0]/nrow_split)
    step_ncol = int(header[1]/ncol_split)
    for i in range(0, nrow_split):
        for j in range(0, ncol_split):
            mesh_projx[i, j] = mesh_x[r, s]
            mesh_projy[i, j] = mesh_y[r, s]
            s += step_ncol
        s = 0
        r += step_nrow

    # Creates the orthogonal vectorial basis for the projections
    CB = header[7] - header[6]
    BA = header[6] - header[5]
    if header[4] == 0:
        mod_BA = np.sqrt(BA[0]**2 + BA[1]**2 + BA[2]**2)
        mod_CB = np.sqrt(CB[0]**2 + CB[1]**2 + CB[2]**2)
        v1 = BA/mod_BA
        v2 = CB/mod_CB
    else:
        BA = BA - ((np.dot(BA, CB))/(np.dot(BA, BA)))*CB
        mod_BA = np.sqrt(BA[0]**2 + BA[1]**2 + BA[2]**2)
        mod_CB = np.sqrt(CB[0]**2 + CB[1]**2 + CB[2]**2)
        v1 = BA/mod_BA
        v2 = CB/mod_CB

    # Calculates the modulus of the vectorial field
    for i in range(0, header[0]):
        for j in range(0, header[1]):
            mod_j[i, j] = np.sqrt((dens[i, j, 0]**2) +
                                  (dens[i, j, 1]**2)+(dens[i, j, 2]**2))

    # Calculates the projections of the vectors in the ABC plane
    ABC_normal = np.cross(BA, CB)
    mod_normal = np.sqrt(ABC_normal[0]**2 +
                         ABC_normal[1]**2 + ABC_normal[2]**2)

    for i in range(0, nrow_split):
        for j in range(0, ncol_split):
            projx_j[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                                                np.cross(np.array([dens[int(i*step_nrow), int(j*step_ncol), 0], dens[int(i*step_nrow), int(j*step_ncol), 1],
                                                                                   dens[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v2)
            projy_j[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                                                np.cross(np.array([dens[int(i*step_nrow), int(j*step_ncol), 0], dens[int(i*step_nrow), int(j*step_ncol), 1],
                                                                                   dens[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v1)

    j = plt.figure()
    j = plt.contourf(mesh_x, mesh_y, mod_j, levels, cmap='winter')
    j = plt.colorbar(mappable=j)
    j = plt.quiver(mesh_projx, mesh_projy, projx_j, projy_j, scale=quivscale)
    j = plt.xlabel('$\AA$')
    j = plt.ylabel('$\AA$')

    return j


def plot_vecfield2D_J(header, dens_JX, dens_JY, dens_JZ, quivscale, levels=150):
    """
    Plots the 2D spin current density vector fields.

    Args:
        header (list): List containing information about the fort.25 header.
        dens_JX (numpy.ndarray): Array representing the X-component of the spin current density.
        dens_JY (numpy.ndarray): Array representing the Y-component of the spin current density.
        dens_JZ (numpy.ndarray): Array representing the Z-component of the spin current density.
        quivscale: Scale factor for the quiver plot.
        levels (int or array-like, optional): Determines the number and positions of the contour lines/regions.

    Returns:
        JX (matplotlib.figure.Figure)
        JY (matplotlib.figure.Figure)
        JZ (matplotlib.figure.Figure)
    """
    import matplotlib.pyplot as plt
    import numpy as np

    for i in range(2, 20):
        if (header[0] % i) == 0:
            nrow_split = int(header[0]/i)

    for i in range(2, 20):
        if (header[1] % i) == 0:
            ncol_split = int(header[1]/i)

    # initializes the arrays
    projx_JX = np.zeros((nrow_split, ncol_split), dtype=float)
    projy_JX = np.zeros((nrow_split, ncol_split), dtype=float)
    mod_JX = np.zeros((header[0], header[1]), dtype=float)
    projx_JY = np.zeros((nrow_split, ncol_split), dtype=float)
    projy_JY = np.zeros((nrow_split, ncol_split), dtype=float)
    mod_JY = np.zeros((header[0], header[1]), dtype=float)
    projx_JZ = np.zeros((nrow_split, ncol_split), dtype=float)
    projy_JZ = np.zeros((nrow_split, ncol_split), dtype=float)
    mod_JZ = np.zeros((header[0], header[1]), dtype=float)

    # Generates the meshgrid
    mesh_x = np.zeros((header[0], header[1]), dtype=float)
    mesh_y = np.zeros((header[0], header[1]), dtype=float)
    mesh_projx = np.zeros((nrow_split, ncol_split), dtype=float)
    mesh_projy = np.zeros((nrow_split, ncol_split), dtype=float)

    T = np.array([[np.sqrt(1 - header[4]**2)*header[2], 0],
                  [header[4]*header[2], header[3]]])  # Change of basis matrix

    for i in range(0, header[0]):
        for j in range(0, header[1]):
            mesh_x[i, j] = np.dot(T, np.array([i, j]))[0]
            mesh_y[i, j] = np.dot(T, np.array([i, j]))[1]

    mesh_x = mesh_x * 0.529177249  # Bohr to Angstrom conversion
    mesh_y = mesh_y * 0.529177249

    r = 0
    s = 0
    step_nrow = int(header[0]/nrow_split)
    step_ncol = int(header[1]/ncol_split)
    for i in range(0, nrow_split):
        for j in range(0, ncol_split):
            mesh_projx[i, j] = mesh_x[r, s]
            mesh_projy[i, j] = mesh_y[r, s]
            s += step_ncol
        s = 0
        r += step_nrow

    # Creates the orthogonal vectorial basis for the projections
    CB = header[7] - header[6]
    BA = header[6] - header[5]
    if header[4] == 0:
        mod_BA = np.sqrt(BA[0]**2 + BA[1]**2 + BA[2]**2)
        mod_CB = np.sqrt(CB[0]**2 + CB[1]**2 + CB[2]**2)
        v1 = BA/mod_BA
        v2 = CB/mod_CB
    else:
        BA = BA - ((np.dot(BA, CB))/(np.dot(BA, BA)))*CB
        mod_BA = np.sqrt(BA[0]**2 + BA[1]**2 + BA[2]**2)
        mod_CB = np.sqrt(CB[0]**2 + CB[1]**2 + CB[2]**2)
        v1 = BA/mod_BA
        v2 = CB/mod_CB

    # Calculates the modulus of the vectorial field
    for i in range(0, header[0]):
        for j in range(0, header[1]):
            mod_JX[i, j] = np.sqrt(
                (dens_JX[i, j, 0]**2)+(dens_JX[i, j, 1]**2)+(dens_JX[i, j, 2]**2))
            mod_JY[i, j] = np.sqrt(
                (dens_JY[i, j, 0]**2)+(dens_JY[i, j, 1]**2)+(dens_JY[i, j, 2]**2))
            mod_JZ[i, j] = np.sqrt(
                (dens_JZ[i, j, 0]**2)+(dens_JZ[i, j, 1]**2)+(dens_JZ[i, j, 2]**2))

    # Calculates the projections of the vectors in the ABC plane
    ABC_normal = np.cross(BA, CB)
    mod_normal = np.sqrt(ABC_normal[0]**2 +
                         ABC_normal[1]**2 + ABC_normal[2]**2)

    for i in range(0, nrow_split):
        for j in range(0, ncol_split):
            projx_JX[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                    np.cross(np.array([dens_JX[int(i*step_nrow), int(j*step_ncol), 0], dens_JX[int(i*step_nrow), int(j*step_ncol), 1],
                                                       dens_JX[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v2)
            projy_JX[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                    np.cross(np.array([dens_JX[int(i*step_nrow), int(j*step_ncol), 0], dens_JX[int(i*step_nrow), int(j*step_ncol), 1],
                                                       dens_JX[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v1)
            projx_JY[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                    np.cross(np.array([dens_JY[int(i*step_nrow), int(j*step_ncol), 0], dens_JY[int(i*step_nrow), int(j*step_ncol), 1],
                                                       dens_JY[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v2)
            projy_JY[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                    np.cross(np.array([dens_JY[int(i*step_nrow), int(j*step_ncol), 0], dens_JY[int(i*step_nrow), int(j*step_ncol), 1],
                                                       dens_JY[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v1)
            projx_JZ[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                    np.cross(np.array([dens_JZ[int(i*step_nrow), int(j*step_ncol), 0], dens_JZ[int(i*step_nrow), int(j*step_ncol), 1],
                                                       dens_JZ[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v2)
            projy_JZ[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                    np.cross(np.array([dens_JZ[int(i*step_nrow), int(j*step_ncol), 0], dens_JZ[int(i*step_nrow), int(j*step_ncol), 1],
                                                       dens_JZ[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v1)

    # Plotting

    JX = plt.figure()
    JX = plt.contourf(mesh_x, mesh_y, mod_JX, levels, cmap='summer')
    JX = plt.colorbar(mappable=JX)
    JX = plt.quiver(mesh_projx, mesh_projy, projx_JX,
                    projy_JX, scale=quivscale)
    JX = plt.xlabel('$\AA$')
    JX = plt.ylabel('$\AA$')

    JY = plt.figure()
    JY = plt.contourf(mesh_x, mesh_y, mod_JY, levels, cmap='summer')
    JY = plt.colorbar(mappable=JY)
    JY = plt.quiver(mesh_projx, mesh_projy, projx_JY,
                    projy_JY, scale=quivscale)
    JY = plt.xlabel('$\AA$')
    JY = plt.ylabel('$\AA$')

    JZ = plt.figure()
    JZ = plt.contourf(mesh_x, mesh_y, mod_JZ, levels, cmap='summer')
    JZ = plt.colorbar(mappable=JZ)
    JZ = plt.quiver(mesh_projx, mesh_projy, projx_JZ,
                    projy_JZ, scale=quivscale)
    JZ = plt.xlabel('$\AA$')
    JZ = plt.ylabel('$\AA$')

    return JX, JY, JZ


# --------------------------------BAND STRUCTURES------------------------------#

def plot_phonon_band(bands, unit='cm-1', k_labels=None, mode='single',
                     not_scaled=False, freq_range=None, k_range=None,
                     color='blue', labels=None, linestl='-', linewidth=1,
                     line_freq0=None, title=None, figsize=None,
                     scheme=None, sharex=True, sharey=True, fontsize=12):
    """
    A wrapper of plot_cry_bands for phonon band structure.

    Args:
        bands (BandsBASE|list): Bands object generated by `CRYSTALClear.crystal_io.Crystal_output.read_pband` or
            a list of BandsBASE objects.
        unit (str): The unit of frequency. Can be 'cm-1' or 'THz'.
        k_labels (list): A list of high-symmetric k point labels. Greek alphabets should be, for example, 'Gamma'.
        mode (str): The plotting mode. Possible values are 'single', 'multi', and 'compare'.
        not_scaled (bool): Whether to scale the x-axis for different volumes.
        energy_range (array): A 2x1 array specifying the energy range.
        k_range (array): A 2x1 array specifying the k-range.
        color (str|list): Color of plot lines. Should be consistent with bands.
        labels (str|list): Plot legend. Should be consistent with bands.
        linestl (str|list): Linestyle string. Should be consistent with bands.
        linewidth (float): The width of the plot lines.
        line_freq0 (str): The color of the frequency=0 line.
        title (str): The title of the plot.
        figsize (list): The figure size specified as [width, height].
        scheme (list|tuple): The layout of subplots.
        sharex (bool): Whether to share the x-axis among subplots.
        sharey (bool): Whether to share the y-axis among subplots.
        dpi (int): Dots per inch resolution of the saved file.
        fontsize (int): Fontsize of the axis labels.            
        transparency(bool): Background transparency of the saved file,

    Returns:
        Matplotlib object

    Raises:
        ValueError: If the specified unit is unknown.

    """
    import re

    import matplotlib.pyplot as plt

    from CRYSTALClear.base.plotbase import plot_cry_bands
    from CRYSTALClear.units import cm_to_thz, thz_to_cm

    if re.match(r'^cm\-1$', unit, re.IGNORECASE):
        unit = 'cm-1'
        is_thz = False
    elif re.match(r'^THz$', unit, re.IGNORECASE):
        unit = 'THz'
        is_thz = True
    else:
        raise ValueError('Unknown unit.')

    if not (isinstance(bands, list) or isinstance(bands, tuple)):
        bands = [bands]

    if line_freq0 == None:
        line_freq0 = (1., 0., 0., 0.)  # Transparent

    for b in bands:
        if unit != b.unit:
            if unit == 'cm-1':
                b.bands[:, :, :] = thz_to_cm(b.bands[:, :, :])
            else:
                b.bands[:, :, :] = cm_to_thz(b.bands[:, :, :])
            b.unit = unit
    if len(bands) == 1:
        bands = bands[0]

    energy_range = freq_range

    fig, ax = plot_cry_bands(bands, k_labels=k_labels, energy_range=energy_range, title=title,
                             not_scaled=not_scaled, mode=mode, linestl=linestl, linewidth=linewidth,
                             color=color, fermi=line_freq0, k_range=k_range, labels=labels,
                             figsize=figsize, scheme=scheme, sharex=sharex, sharey=sharey, fermialpha=1, fermiwidth=0)
    if is_thz == True:
        fig.supylabel('Frequency (THz)', fontsize=fontsize)
    else:
        fig.supylabel('Frequency (cm$^{-1}$)', fontsize=fontsize)

    return fig, ax


def plot_electron_band(bands, unit='eV', k_labels=None, mode='single',
                       not_scaled=False, energy_range=None, k_range=None,
                       color='blue', labels=None, linestl='-', linewidth=1,
                       fermi='forestgreen', fermiwidth=1.5, fermialpha=1, title=None, figsize=None,
                       scheme=None, sharex=True, sharey=True, fontsize=12):
    """
    A wrapper of plot_cry_bands for electron band structure.

    Args:
        bands (BandsBASE|list): Bands object generated by `CRYSTALClear.crystal_io.Properties_output.read_bands` or
            a list of BandsBASE objects.
        unit (str): The unit of energy. Can be 'eV' or 'Hartree'.
        k_labels (list): A list of high-symmetric k point labels. Greek alphabets should be, for example, 'Gamma'.
        mode (str): The plotting mode. Possible values are 'single', 'multi', and 'compare'.
        not_scaled (bool): Whether to scale the x-axis for different volumes.
        energy_range (array): A 2x1 array specifying the energy range.
        k_range (array): A 2x1 array specifying the k-range.
        color (str|list): Color of plot lines. Should be consistent with bands.
        labels (str|list): Plot legend. Should be consistent with bands.
        linestl (str|list): Linestyle string. Should be consistent with bands.
        linewidth (float): The width of the plot lines.
        fermi (str): The color of the Fermi level line.
        fermiwidth (float): The width of the fermi line.
        fermialpha (float): Opacity of the fermi level 0-1.
        title (str): The title of the plot.
        figsize (list): The figure size specified as [width, height].
        scheme (list|tuple): The layout of subplots.
        sharex (bool): Whether to share the x-axis among subplots.
        sharey (bool): Whether to share the y-axis among subplots.
        dpi (int): Dots per inch resolution of the saved file.
        fontsize (int): Fontsize of the axis labels 
        transparency: Background Transparency of the saved file.

    Returns:
        Matplolib object

    Raises:
        ValueError: If the specified unit is unknown.

    """
    import re

    import matplotlib.pyplot as plt

    from CRYSTALClear.base.plotbase import plot_cry_bands
    from CRYSTALClear.units import H_to_eV, eV_to_H

    if re.match(r'^eV$', unit, re.IGNORECASE):
        unit = 'eV'
        is_ev = True
    elif re.match(r'^Hartree$', unit, re.IGNORECASE):
        unit = 'Hartree'
        is_ev = False
    else:
        raise ValueError('Unknown unit.')

    if not (isinstance(bands, list) or isinstance(bands, tuple)):
        bands = [bands]

    for b in bands:
        if unit != b.unit:
            if unit == 'eV':
                b.bands[:, :, :] = H_to_eV(b.bands[:, :, :])
            else:
                b.bands[:, :, :] = eV_to_H(b.bands[:, :, :])
            b.unit = unit
    if len(bands) == 1:
        bands = bands[0]

    fig, ax = plot_cry_bands(bands, k_labels=k_labels, energy_range=energy_range, title=title,
                             not_scaled=not_scaled, mode=mode, linestl=linestl, linewidth=linewidth,
                             color=color, fermi=fermi, fermiwidth=fermiwidth, fermialpha=fermialpha, k_range=k_range, labels=labels,
                             figsize=figsize, scheme=scheme, sharex=sharex, sharey=sharey)
    if is_ev == True:
        fig.supylabel('$E-E_{F}$ (eV)', fontsize=fontsize)
    else:
        fig.supylabel('$E-E_{F}$ (Hartree)', fontsize=fontsize)

    return fig, ax


# -------------------------------DENSITY OF STATES-----------------------------#


def plot_electron_dos(doss, unit='eV', beta='up', overlap=False, prj=None,
                      energy_range=None, dos_range=None, color='blue',
                      labels=None, linestl=None, linewidth=1, fermi='forestgreen',
                      title=None, figsize=None):
    """
    A wrapper of plot_cry_doss for electron density of states.

    Args:
        doss (DOSBASE): DOS obect generated by code:`CRYSTALClear.crystal_io.Properties_output.read_doss`.
            Or a list of DOSBASE objects.
        unit (str): 'eV' or 'Hartree'
        beta (str): Plot spin-down state 'up' or 'down'
        overlap (bool): Plotting multiple lines into the same figure
        prj (list): Index of selected projection. Consistent with the
            index of the 2nd dimension of :code:`doss.doss`
        energy_range (list[float]): 2*1 list of energy range
        dos_range (list[float]): 2*1 list of DOS range
        color (str | list[str]): Color of plot lines. *Should be
            consistent with number of projections.*
        labels (str | list[str]): Plot legend. *Should be consistent with
            number of projections.*
        linestl (str | list[str]): linestyle string. *Should be consistent
            with number of projections.*
        linewidth (float)
        fermi (str): Color of Fermi level line.
        title (str)
        figsize (list[float])

    Returns:
        Matplotlib object
    """
    import re

    import matplotlib.pyplot as plt

    from CRYSTALClear.base.plotbase import plot_cry_doss
    from CRYSTALClear.units import H_to_eV, eV_to_H

    if re.match(r'^ev$', unit, re.IGNORECASE):
        unit = 'eV'
        is_ev = True
    elif re.match(r'^hartree$', unit, re.IGNORECASE):
        unit = 'Hartree'
        is_ev = False
    else:
        raise ValueError('Unknown unit.')

    if not (isinstance(doss, list) or isinstance(doss, tuple)):
        doss = [doss]

    for d in doss:
        if unit != d.unit:
            if unit == 'eV':
                d.doss[:, 0, :] = H_to_eV(d.doss[:, 0, :])
                d.doss[:, 1:, :] = eV_to_H(d.doss[:, 1:, :])
            else:
                d.doss[:, 0, :] = eV_to_H(d.doss[:, 0, :])
                d.doss[:, 1:, :] = H_to_eV(d.doss[:, 1:, :])
            d.unit = unit
    if len(doss) == 1:
        doss = doss[0]

    fig, ax = plot_cry_doss(doss, color=color, fermi=fermi, overlap=overlap,
                            labels=labels, figsize=figsize, linestl=linestl,
                            linewidth=linewidth, title=title, beta=beta,
                            energy_range=energy_range, dos_range=dos_range, prj=prj)
    if is_ev == True:
        fig.supylabel('DOS (states/eV)')
        fig.supxlabel('Energy (eV)')
    else:
        fig.supylabel('DOS (states/Hartree)')
        fig.supxlabel('Energy (Hartree)')

    return fig, ax


def plot_phonon_dos(doss, unit='cm-1', overlap=False, prj=None,
                    freq_range=None, dos_range=None, color='blue',
                    labels=None, linestl=None, linewidth=1, line_freq0=None,
                    title=None, figsize=None):
    """
    A wrapper of plot_cry_doss for electron density of states.

    Args:
        doss (DOSBASE): DOS obect generated by code:`CRYSTALClear.crystal_io.Crystal_output.read_pdos`.
            Or a list of DOSBASE objects.
        unit (str): 'cm-1' or 'THz'
        overlap (bool): Plotting multiple lines into the same figure
        prj (list): Index of selected projection. Consistent with the
            index of the 2nd dimension of :code:`doss.doss`
        freq_range (list[float]): 2*1 list of frequency range
        dos_range (list[float]): 2*1 list of DOS range
        color (str | list[str]): Color of plot lines. *Should be
            consistent with number of projections.*
        labels (str | list[str]): Plot legend. *Should be consistent with
            number of projections.*
        linestl (str | list[str]): linestyle string. *Should be consistent
            with number of projections.*
        linewidth (float)
        line_freq0 (str): Color of frequency = 0 line.
        title (str)
        figsize (list[float])

    Returns:
        Matplotlib object
    """
    import re

    import matplotlib.pyplot as plt

    from CRYSTALClear.base.plotbase import plot_cry_doss
    from CRYSTALClear.units import cm_to_thz, thz_to_cm

    if re.match(r'^cm\-1$', unit, re.IGNORECASE):
        unit = 'cm-1'
        is_thz = False
    elif re.match(r'^thz$', unit, re.IGNORECASE):
        unit = 'THz'
        is_thz = True
    else:
        raise ValueError('Unknown unit.')

    if not (isinstance(doss, list) or isinstance(doss, tuple)):
        doss = [doss]

    for d in doss:
        if unit != d.unit:
            if unit == 'cm-1':
                d.doss[:, 0, :] = thz_to_cm(d.doss[:, 0, :])
                d.doss[:, 1:, :] = cm_to_thz(d.doss[:, 1:, :])
            else:
                d.doss[:, 0, :] = cm_to_thz(d.doss[:, 0, :])
                d.doss[:, 1:, :] = thz_to_cm(d.doss[:, 1:, :])
            d.unit = unit

    if line_freq0 == None:
        line_freq0 = (1., 0., 0., 0.)  # Transparent
    if len(doss) == 1:
        doss = doss[0]

    fig, ax = plot_cry_doss(doss, color=color, fermi=line_freq0, overlap=overlap,
                            labels=labels, figsize=figsize, linestl=linestl,
                            linewidth=linewidth, title=title, beta='up',
                            energy_range=freq_range, dos_range=dos_range, prj=prj)

    if is_thz == True:
        fig.supylabel('DOS (states/THz)')
        fig.supxlabel('Frequency (THz)')
    else:
        fig.supylabel('DOS (states/cm$^{-1}$)')
        fig.supxlabel('Frequency (cm$^{-1}$)')

    return fig, ax


# -----------------------------BAND + DENSITY OF STATES------------------------#


def plot_electron_banddos(bands, doss, unit='eV', k_labels=None, dos_beta='down',
                          dos_prj=None, energy_range=None, dos_range=None,
                          color_band='blue', color_dos='blue', labels=None, linestl_band='-',
                          linestl_dos=None, linewidth=1, fermi='forestgreen',
                          title=None, figsize=None, legend=False):
    """
    A wrapper of plot_cry_es for electron band structure + dos. For spin-polarized cases, beta state.

    Args:
        bands (BandsBASE|list): Bands object generated by CRYSTALClear.crystal_io.Properties_output.read_bands
            or a list of BandsBASE objects.
        doss (DOSBASE): DOS object generated by CRYSTALClear.crystal_io.Properties_output.read_doss
            or a list of DOSBASE objects.
        unit (str): Unit of energy. Valid options are 'eV' or 'Hartree'.
        k_labels (list): A list of high-symmetric k point labels. Greek alphabets should be
            represented as strings, for example, 'Gamma'.
        dos_beta (str): Spin state to plot. Valid options are 'Up' or 'down'. If 'down', the beta
            state will be plotted on the same side as the alpha state, otherwise on the other side.
        dos_prj (list): Index of selected projection. Consistent with the index of the 2nd dimension
            of doss.doss.
        energy_range (list): A list of two values representing the energy range to be plotted.
        dos_range (list): DOS range for the y-axis.
        color_band (str): Color of the electron bands in the plot.
        color_dos (str): Color of the density of states (DOS) in the plot.
        labels (list): A list of labels for the plot legend.
        linestl_band (str): Linestyle of the electron bands.
        linestl_dos (str): Linestyle of the density of states (DOS).
        linewidth (float): Width of the lines in the plot.
        fermi (str): Color of the Fermi level line.
        title (str): Title of the plot.
        figsize (list[float]): Size of the figure in inches (width, height).
        legend (bool): Enables or disables the legend of the density of states (DOS).

    Returns:
        Matplotlib object

    Raises:
        ValueError: If the unit parameter is unknown.

    """
    import re

    import matplotlib.pyplot as plt

    from CRYSTALClear.base.plotbase import plot_cry_es
    from CRYSTALClear.units import H_to_eV, eV_to_H

    if re.match(r'^ev$', unit, re.IGNORECASE):
        unit = 'eV'
        is_ev = True
    elif re.match(r'^hartree$', unit, re.IGNORECASE):
        unit = 'Hartree'
        is_ev = False
    else:
        raise ValueError('Unknown unit.')

    if unit != doss.unit:
        if unit == 'eV':
            doss.doss[:, 0, :] = H_to_eV(doss.doss[:, 0, :])
            doss.doss[:, 1:, :] = eV_to_H(doss.doss[:, 1:, :])
        else:
            doss.doss[:, 0, :] = eV_to_H(doss.doss[:, 0, :])
            doss.doss[:, 1:, :] = H_to_eV(doss.doss[:, 1:, :])
        doss.unit = unit
    if unit != bands.unit:
        if unit == 'eV':
            bands.bands[:, :, :] = H_to_eV(bands.bands[:, :, :])
        else:
            bands.bands[:, :, :] = eV_to_H(bands.bands[:, :, :])
        bands.unit = unit

    fig, ax = plot_cry_es(bands=bands, doss=doss, k_labels=k_labels, color_bd=color_band,
                          color_doss=color_dos, fermi=fermi, energy_range=energy_range,
                          linestl_bd=linestl_band, linestl_doss=linestl_dos,
                          linewidth=linewidth, prj=dos_prj, figsize=figsize, labels=labels,
                          dos_range=dos_range, title=title, dos_beta=dos_beta, legend=legend)
    if is_ev == True:
        fig.supylabel('Energy (eV)')
    else:
        fig.supylabel('Energy (Hartree)')

    return fig, ax


def plot_phonon_banddos(bands, doss, unit='cm-1', k_labels=None, dos_prj=None,
                        freq_range=None, dos_max_range=None, color_band='blue',
                        color_dos='blue', labels=None, linestl_band='-',
                        linestl_dos=None, linewidth=1, freq0_line=None,
                        title=None, figsize=None):
    """
    A wrapper of plot_cry_es for phonon band structure + dos. Only one pair is permitted.

    Args:
        bands (BandsBASE|list): Bands object generated by CRYSTALClear.crystal_io.Properties_output.read_bands
            or a list of BandsBASE objects.
        doss (DOSBASE): DOS object generated by CRYSTALClear.crystal_io.Properties_output.read_doss
            or a list of DOSBASE objects.
        unit (str): Unit of frequency. Valid options are 'cm-1' or 'THz'.
        k_labels (list): A list of high-symmetric k point labels. Greek alphabets should be
            represented as strings, for example, 'Gamma'.
        dos_prj (list): Index of selected projection. Consistent with the index of the 2nd dimension
            of doss.doss.
        freq_range (list): A list of two values representing the frequency range to be plotted.
        dos_max_range (float): Maximum DOS range for the y-axis.
        color_band (str): Color of the phonon bands in the plot.
        color_dos (str): Color of the density of states (DOS) in the plot.
        labels (list): A list of labels for the plot legend.
        linestl_band (str): Linestyle of the phonon bands.
        linestl_dos (str): Linestyle of the density of states (DOS).
        linewidth (float): Width of the lines in the plot.
        freq0_line (str): Color of the frequency=0 line.
        title (str): Title of the plot.
        figsize (list[float]): Size of the figure in inches (width, height).

    Returns:
        Matplotlib object

    Raises:
        ValueError: If the unit parameter is unknown.

    """
    import re

    import matplotlib.pyplot as plt

    from CRYSTALClear.base.plotbase import plot_cry_es
    from CRYSTALClear.units import cm_to_thz, thz_to_cm

    if re.match(r'^cm\-1$', unit, re.IGNORECASE):
        unit = 'cm-1'
        is_thz = False
    elif re.match(r'^thz$', unit, re.IGNORECASE):
        unit = 'THz'
        is_thz = True
    else:
        raise ValueError('Unknown unit.')

    if unit != doss.unit:
        if unit == 'cm-1':
            doss.doss[:, 0, :] = thz_to_cm(doss.doss[:, 0, :])
            doss.doss[:, 1:, :] = cm_to_thz(doss.doss[:, 1:, :])
        else:
            doss.doss[:, 0, :] = cm_to_thz(doss.doss[:, 0, :])
            doss.doss[:, 1:, :] = thz_to_cm(doss.doss[:, 1:, :])
        doss.unit = unit
    if unit != bands.unit:
        if unit == 'cm-1':
            bands.bands[:, :, :] = thz_to_cm(bands.bands[:, :, :])
        else:
            bands.bands[:, :, :] = cm_to_thz(bands.bands[:, :, :])
        bands.unit = unit

    if line_freq0 == None:
        line_freq0 = (1., 0., 0., 0.)  # Transparent

    fig, ax = plot_cry_es(bands=bands, doss=doss, k_labels=k_labels, color_bd=color_band,
                          color_doss=color_dos, fermi=line_freq0, energy_range=energy_range,
                          linestl_bd=linestl_band, linestl_doss=linestl_dos,
                          linewidth=linewidth, prj=dos_prj, figsize=figsize, labels=labels,
                          dos_max_range=dos_max_range, title=title, dos_beta='up')
    if is_thz == True:
        fig.supylabel('Frequency (THz)')
    else:
        fig.supylabel('Frequency (cm$^{-1}$)')

    return fig, ax


##############################################################################
#                                                                            #
#                                     QTAIM                                  #
#                                                                            #
##############################################################################

# ----------------------------------CONTOUR PLOT-------------------------------#


def plot_cry_contour(contour_obj):
    """
    Plot a contour plot.

    Args:
        contour_obj (object): Contour object representing the contour plot.

    Returns:
        None

    Notes:
        - Plots a contour plot based on the data in the contour object.
        - Retrieves the data from the contour object and converts it to a 2D list.
        - Sets the figure size based on x_graph_param and y_graph_param attributes of the contour object.
        - Sets the x-axis and y-axis labels.
        - Creates a meshgrid using the x_points and y_points attributes of the contour object.
        - Defines contour levels, colors, linestyles, and fmt.
        - Plots the contour plot.
        - Saves the plot to a file named 'figure_TIPO_YYYY-MM-DD_HHMMSS.jpg' in the current directory.

    """
    import os
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    df = contour_obj.df
    n_punti_x = contour_obj.npx

    for i in range(0, 8):
        df[i] = df[i].astype(float)

    flat_list = [item for sublist in df.values for item in sublist]

    cleaned_list = [x for x in flat_list if ~np.isnan(x)]

    l = [cleaned_list[x:x+n_punti_x]
         for x in range(0, len(cleaned_list), n_punti_x)]

    c = contour_obj.x_graph_param
    d = contour_obj.y_graph_param

    plt.rcParams["figure.figsize"] = [c, d]

    plt.xlabel(r'$\AA$', fontsize=18)
    plt.ylabel(r'$\AA$', fontsize=18)

    X, Y = np.meshgrid(contour_obj.x_points, contour_obj.y_points)

    levels = contour_obj.levels
    colors = contour_obj.colors
    linestyles = contour_obj.linestyles
    fmt = contour_obj.fmt

    # Change here to have or not the isovalues on the plot
    iso = True
    # iso = False

    if (iso == True):
        L = plt.contour(X, Y, l, levels=levels, colors=colors, linestyles=linestyles, linewidths=0.7,
                        alpha=1)
        plt.clabel(L, inline=1, fontsize=7, fmt=fmt)
    elif (iso == False):
        L = plt.contour(X, Y, l, levels=levels, colors=colors, linestyles=linestyles, linewidths=0.7,
                        alpha=1)

    path = os.path.join('./'+'figure_' + contour_obj.tipo +
                        '_' + time.strftime("%Y-%m-%d_%H%M%S") + '.jpg')
    plt.savefig(path, bbox_inches='tight', dpi=600)
    print('\nThe image has been saved in the current directory')

    plt.show()


def plot_cry_contour_differences(contour_obj, contour_obj_ref):
    """
    Plot the differences between two contour plots.

    Args:
        contour_obj (object): Contour object representing the original contour plot.
        contour_obj_ref (object): Contour object representing the reference contour plot.

    Returns:
        None

    Notes:
        - Plots the differences between two contour plots.
        - Requires the contour objects to have a tipo attribute with values 'SURFLAPP', 'SURFLAPM', 'SURFRHOO', or 'SURFELFB'.
        - Calculates the difference between the dataframes of the two contour objects.
        - Sets the figure size based on x_graph_param and y_graph_param attributes of the contour object.
        - Sets the x-axis and y-axis labels.
        - Creates a meshgrid using the x_points and y_points attributes of the contour object.
        - Defines contour levels, colors, and linestyles.
        - Plots the contour differences.
        - Saves the plot to a file named 'figure_diff_TIPO_YYYY-MM-DD_HHMMSS.jpg' in the current directory.

    """
    import os
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    if (contour_obj.tipo == 'SURFLAPP') or (contour_obj.tipo == 'SURFLAPM') or (contour_obj.tipo == 'SURFRHOO') or (contour_obj.tipo == 'SURFELFB'):
        pass
    else:
        sys.exit(
            'Difference option only allowed for SURFLAPP, SURFLAPM, SURFRHOO and SURFELFB file')

    n_punti_x = contour_obj.npx

    df = contour_obj.df
    for i in range(0, 8):
        df[i] = df[i].astype(float)

    df_ref = contour_obj_ref.df
    for i in range(0, 8):
        df_ref[i] = df_ref[i].astype(float)

    df_diff = df - df_ref

    flat_list = [item for sublist in df_diff.values for item in sublist]

    cleaned_list = [x for x in flat_list if ~np.isnan(x)]

    l = [cleaned_list[x:x+n_punti_x]
         for x in range(0, len(cleaned_list), n_punti_x)]

    c = contour_obj.x_graph_param
    d = contour_obj.y_graph_param

    plt.rcParams["figure.figsize"] = [c, d]

    plt.xlabel(r'$\AA$', fontsize=18)
    plt.ylabel(r'$\AA$', fontsize=18)

    X, Y = np.meshgrid(contour_obj.x_points, contour_obj.y_points)

    ctr1dif = np.array([-8, -4, -2, -0.8, -0.4, -0.2, -0.08, -0.04, -0.02, -0.008, -0.004, -0.002, -0.0008, -0.0004, -0.0002, 0,
                       0.0002, 0.0004, 0.0008, 0.002, 0.004, 0.008, 0.02, 0.04, 0.08, 0.2, 0.4, 0.8, 2, 4, 8])
    colors1dif = ['b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'k', 'r', 'r', 'r', 'r', 'r', 'r', 'r',
                  'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
    ls1dif = ['--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', 'dotted', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']

    levels = ctr1dif
    colors = colors1dif
    linestyles = ls1dif
    fmt = '%1.4f'

    # Change here to have or not the isovalues on the plot
    iso = True
    # iso = False

    if (iso == True):
        L = plt.contour(X, Y, l, levels=levels, colors=colors, linestyles=linestyles, linewidths=0.7,
                        alpha=1)
        plt.clabel(L, inline=1, fontsize=7, fmt=fmt)
    elif (iso == False):
        L = plt.contour(X, Y, l, levels=levels, colors=colors, linestyles=linestyles, linewidths=0.7,
                        alpha=1)

    path = os.path.join('./'+'figure_diff_' + contour_obj.tipo +
                        '_' + time.strftime("%Y-%m-%d_%H%M%S") + '.jpg')
    plt.savefig(path, bbox_inches='tight', dpi=600)
    print('\nThe image has been saved in the current directory')

    plt.show()

# --------------------------------------XRD------------------------------------#


def plot_cry_xrd(xrd_obj):
    """
    Plot the X-ray diffraction pattern.

    Args:
        xrd_obj (object): XRD object containing the data for the X-ray diffraction pattern.

    Returns:
        None

    Notes:
        - Plots the X-ray diffraction pattern.
        - Sets the figure size to [16, 9].
        - Sets the x-axis limit to (0, 30).
        - Saves the plot to a file named 'figure_XRD_YYYY-MM-DD_HHMMSS.jpg' in the current directory.

    """
    import os
    import time

    import matplotlib.pyplot as plt

    plt.rcParams["figure.figsize"] = [16, 9]

    plt.plot(xrd_obj.x, xrd_obj.y)

    plt.xlim((0, 30))

    path = os.path.join('./'+'figure_'+'XRD_' +
                        time.strftime("%Y-%m-%d_%H%M%S") + '.jpg')
    plt.title(xrd_obj.title, fontsize=20)
    plt.savefig(path, bbox_inches='tight', dpi=600)

    plt.show()

# -------------------------------------RHOLINE---------------------------------#


def plot_cry_rholine(rholine_obj):
    """
    Plot the resistivity as a function of distance.

    Args:
        rholine_obj (object): Rholine object containing the data for the resistivity.

    Returns:
        None

    Notes:
        - Plots the resistivity as a function of distance.
        - Sets the x-axis label as 'd  [$\AA$]' and the y-axis label as r'$\rho$  [$\frac{e}{\AA^3}$]'.
        - Saves the plot to a file named 'figure_rholine_YYYY-MM-DD_HHMMSS.jpg' in the current directory.
    """
    import os
    import time

    import matplotlib.pyplot as plt

    plt.plot(rholine_obj.x, rholine_obj.y)

    plt.xlabel('d  [$\AA$]', fontsize=14)
    plt.ylabel(r'$\rho$  [$\frac{e}{\AA^3}$]', fontsize=16)

    path = os.path.join('./'+'figure_'+'rholine_' +
                        time.strftime("%Y-%m-%d_%H%M%S") + '.jpg')
    plt.title(rholine_obj.title, fontsize=15)
    plt.savefig(path, bbox_inches='tight', dpi=600)

    plt.show()

# -----------------------------------LAPLACIAN---------------------------------#


def plot_cry_lapl_profile(lapl_obj):
    """
    Plot the Laplacian profile of a crystal.

    Args:
        lapl_obj (object): Laplacian object containing the data for the Laplacian profile.

    Returns:
        None

    Notes:
        - Plots the Laplacian profile using the data from the Laplacian object.
        - The x-axis represents the distance in angstroms.
        - The y-axis represents the Laplacian in electrons per cubic angstrom to the fifth power (e/A^5).
        - The area under the curve where the Laplacian is negative is filled with a light blue color.
        - The area under the curve where the Laplacian is positive is filled with a light coral color.
    """
    import time

    import matplotlib.pyplot as plt

    plt.plot(lapl_obj.datax, lapl_obj.datay)

    plt.fill_between(lapl_obj.datax, lapl_obj.datay, where=(
        lapl_obj.datay < 0), color='lightblue', interpolate=True)
    plt.fill_between(lapl_obj.datax, lapl_obj.datay, where=(
        lapl_obj.datay > 0), color='lightcoral', interpolate=True)

    # plt.xlim(-0.5,0.5)
    # plt.ylim(-200,200)

    plt.xlabel('Distance [A]')
    plt.ylabel('Laplacian [e/A^5]')

    plt.show()

# -----------------------------DENSITY PROFILE---------------------------------#


def plot_cry_density_profile(lapl_obj):
    """
    Plot the density profile of a crystal.

    Args:
        lapl_obj (object): Laplacian object containing the data for the density profile.

    Returns:
        None

    Notes:
        - Plots the density profile using the data from the Laplacian object.
        - The x-axis represents the distance in angstroms.
        - The y-axis represents the density in electrons per cubic angstrom (e/A^3).
    """
    import time

    import matplotlib.pyplot as plt

    plt.plot(lapl_obj.datax, lapl_obj.datay)

    plt.xlabel('Distance [A]')
    plt.ylabel('Density [e/A^3]')

    plt.show()

##############################################################################
#                                                                            #
#                             TRANSPORT PROPERTIES                           #
#                                                                            #
##############################################################################

# -------------------------------------SEEBACK---------------------------------#


def plot_cry_seebeck_potential(seebeck_obj, direction, temperature):
    """
    Plot the Seebeck coefficient as a function of chemical potential.

    Args:
        seebeck_obj (object): Seebeck object containing the data for the Seebeck coefficient.
        direction (str): choose the direction to plot among 'S_xx', 'S_xy', 'S_xz', 'S_yx', 'S_yy', 'S_yz', 'S_yz', 'S_zx', 'S_zy', 'S_zz'.
        temperature (value/str): choose the temperature to be considered or 'all' to consider them all together

    Returns:
        Figure object

    Notes:
        - Plots the Seebeck coefficient as a function of chemical potential for each temperature.
        - Distinguishes between n-type and p-type conduction with dashed and solid lines, respectively.
    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    direction = direction.lower().replace('_', '')

    if direction.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid direction')

    if direction == 'sxx':
        col = 3
    elif direction == 'sxy':
        col = 4
    elif direction == 'sxz':
        col = 5
    elif direction == 'syx':
        col = 6
    elif direction == 'syy':
        col = 7
    elif direction == 'syz':
        col = 8
    elif direction == 'szx':
        col = 9
    elif direction == 'szy':
        col = 10
    elif direction == 'szz':
        col = 11

    else:
        sys.exit('please, choose a valid direction')

    listtemp = []

    for k in range(0, len(seebeck_obj.all_data)):
        listtemp.append(seebeck_obj.temp[k])

    if temperature != 'all' and temperature not in listtemp:
        sys.exit('Please, choose a valid temperature.')

    vol = seebeck_obj.volume

    x = []
    for k in range(0, len(seebeck_obj.all_data)):
        x.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: float(x.split()[0]))))

    carrier = []
    for k in range(0, len(seebeck_obj.all_data)):
        carrier.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: (float(x.split()[2])/vol))))

    y = []
    for k in range(0, len(seebeck_obj.all_data)):
        y.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: float(x.split()[col])*1000000)))

    yneg = []
    ypos = []
    xpos = []
    xneg = []
    yposfin = []
    xposfin = []
    ynegfin = []
    xnegfin = []

    for k in range(0, len(seebeck_obj.all_data)):
        for j in range(0, len(carrier[k])):
            if carrier[k][j] >= 0:
                xpos.append(x[k][j])
                ypos.append(y[k][j])
            else:
                xneg.append(x[k][j])
                yneg.append(y[k][j])
        yposfin.append(ypos)
        ynegfin.append(yneg)
        xposfin.append(xpos)
        xnegfin.append(xneg)
        xpos = []
        ypos = []
        xneg = []
        yneg = []

    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    colours = ['royalblue', 'orange', 'green', 'red',
               'purple', 'brown', 'pink', 'grey', 'olive', 'cyan']

    endx = []
    endy = []

    if temperature != 'all':
        k = listtemp.index(temperature)
        endx = [xposfin[k][-1], xnegfin[k][0]]
        endy = [yposfin[k][-1], ynegfin[k][0]]
        plt.figure()
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xposfin[k], yposfin[k], color=colours[k],
                 label=str(seebeck_obj.temp[k])+' K')
        plt.plot(xnegfin[k], ynegfin[k], '--', color=colours[k])
        plt.xlabel('Chemical Potential (eV)', fontsize=12)
        plt.ylabel('Seebeck Coefficient ($\mu$V/K)', fontsize=12)
        plt.axhline(0, color='k')
        plt.title('Seebeck at ' + str(seebeck_obj.temp[k]) + ' K')
        plt.legend(loc='upper left', fontsize=12)
        fig = plt.gcf()
        # plt.savefig('seebeck_potential_at_' + str(seebeck_obj.temp[k]) + 'K___' + time.strftime(
        # "%Y-%m-%d_%H%M%S") + '.jpg', format='jpg', dpi=600, bbox_inches='tight')
        return fig
    else:

        from matplotlib.pyplot import figure
        figure(figsize=(7, 7))
        for k in range(0, len(seebeck_obj.all_data)):
            endx = [xposfin[k][-1], xnegfin[k][0]]
            endy = [yposfin[k][-1], ynegfin[k][0]]
            plt.plot(endx, endy, color=colours[k])
            plt.plot(xposfin[k], yposfin[k], color=colours[k],
                     label=str(seebeck_obj.temp[k])+' K')
            plt.plot(xnegfin[k], ynegfin[k], '--', color=colours[k])
            plt.xlabel('Chemical Potential (eV)', fontsize=12)
            plt.axhline(0, color='k')
            plt.title('Seebeck at different T')
            fig = plt.gcf()
        # plt.savefig('seebeck_potential_different_T_' + time.strftime("%Y-%m-%d_%H%M%S") +
            # '.jpg', format='jpg', dpi=100, bbox_inches='tight')
        return fig


def plot_cry_seebeck_carrier(seebeck_obj, direction, temperature):
    """
    Plot the Seebeck coefficient as a function of charge carrier concentration.

    Args:
        seebeck_obj: Seebeck object containing the data for the Seebeck coefficient.
        direction (str): choose the direction to plot among 'S_xx', 'S_xy', 'S_xz', 'S_yx', 'S_yy', 'S_yz', 'S_yz', 'S_zx', 'S_zy', 'S_zz'.
        temperature (value/str): choose the temperature to be considered or 'all' to consider them all together

    Returns:
        Figure object

    Notes:
        - Plots the Seebeck coefficient as a function of charge carrier concentration for each temperature, distinguishing between n-type and p-type conduction.
    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    direction = direction.lower().replace('_', '')

    if direction.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid direction')

    if direction == 'sxx':
        col = 3
    elif direction == 'sxy':
        col = 4
    elif direction == 'sxz':
        col = 5
    elif direction == 'syx':
        col = 6
    elif direction == 'syy':
        col = 7
    elif direction == 'syz':
        col = 8
    elif direction == 'szx':
        col = 9
    elif direction == 'szy':
        col = 10
    elif direction == 'szz':
        col = 11
    else:
        sys.exit('please, choose a valid direction')

    listtemp = []

    for k in range(0, len(seebeck_obj.all_data)):
        listtemp.append(seebeck_obj.temp[k])

    if temperature != 'all' and temperature not in listtemp:
        sys.exit('Please, choose a valid temperature.')

    vol = seebeck_obj.volume

    x = []
    for k in range(0, len(seebeck_obj.all_data)):
        x.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: (float(x.split()[2])/vol))))
    y = []

    for k in range(0, len(seebeck_obj.all_data)):
        y.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: float(x.split()[col])*1000000)))

    yneg = []
    ypos = []
    xpos = []
    xneg = []
    yposfin = []
    xposfin = []
    ynegfin = []
    xnegfin = []

    for k in range(0, len(seebeck_obj.all_data)):
        for j in range(0, len(x[k])):
            if x[k][j] >= 0:
                xpos.append(x[k][j])
                ypos.append(y[k][j])
            else:
                xneg.append(x[k][j])
                yneg.append(y[k][j])
        yposfin.append(ypos)
        ynegfin.append(yneg)
        xposfin.append(xpos)
        xnegfin.append(xneg)
        xpos = []
        ypos = []
        xneg = []
        yneg = []

    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    colours = []
    colours = ['royalblue', 'orange', 'green', 'red',
               'purple', 'brown', 'pink', 'grey', 'olive', 'cyan']

    if temperature != 'all':
        k = listtemp.index(temperature)
        endx = [xposfin[k][-1], xnegfin[k][0]]
        endy = [yposfin[k][-1], ynegfin[k][0]]
        plt.figure()
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xposfin[k], yposfin[k], color=colours[k],
                 label=str(seebeck_obj.temp[k])+' K')
        plt.plot(abs(np.array(xnegfin[k])), ynegfin[k], '--', color=colours[k])
        plt.xlabel('Charge Carrier Concentration (cm$^{-3}$)', fontsize=12)
        plt.ylabel('Seebeck Coefficient ($\mu$V/K)', fontsize=12)
        plt.axhline(0, color='k')
        plt.title('Seebeck at ' + str(seebeck_obj.temp[k]) + ' K')
        plt.legend(loc='upper left', fontsize=12)
        plt.xscale('log')
        plt.show()
        fig = plt.gcf()
        return fig

    else:

        from matplotlib.pyplot import figure

        figure(figsize=(7, 7))
        for k in range(0, len(seebeck_obj.all_data)):
            endx = [xposfin[k][-1], xnegfin[k][0]]
            endy = [yposfin[k][-1], ynegfin[k][0]]
            plt.plot(endx, endy, color=colours[k])
            plt.plot(xposfin[k], yposfin[k], color=colours[k],
                     label=str(seebeck_obj.temp[k])+' K')
            plt.plot(abs(np.array(xnegfin[k])),
                     ynegfin[k], '--', color=colours[k])
            plt.xlabel('Charge Carrier Concentration (cm$^{-3}$)', fontsize=12)
            plt.ylabel('Seebeck Coefficient ($\mu$V/K)', fontsize=12)
            plt.axhline(0, color='k')
            plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
            plt.xscale('log')
            plt.title('Seebeck at different T')
        # plt.savefig('seebeck_carrier_different_T_' + time.strftime("%Y-%m-%d_%H%M%S") +
         #           '.jpg', format='jpg', dpi=100, bbox_inches='tight')
        plt.show()
        fig = plt.gcf()
        return fig


def plot_cry_multiseebeck(direction, temperature, minpot, maxpot, *seebeck):
    """
    Plot the seebeck coefficient from different files.

    Args:
        direction (str): choose the direction to plot among 'S_xx', 'S_xy', 'S_xz', 'S_yx', 'S_yy', 'S_yz', 'S_yz', 'S_zx', 'S_zy', 'S_zz'.
        temperature (value): choose the temperature to be considered 
        minpot (value): lower value of chemical potential you want to plot in eV  
        maxpot (value): higher value of chemical potential you want to plot in eV 
        *seebeck (obj): Variable number of seebeck objects containing the data for the Seebeck coefficient.

    Returns:
        Figure object

    Notes:
        - Plots the seebeck coefficient for each seebeck object.
        - Differentiates transport coefficients due to n-type or p-type conduction using dashed and solid lines.

    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    # k = int(input(
    #   'Insert the index of temperature you want to plot \n(i.e. if your temperature are [T1, T2, T3] indexes are [0, 1, 2])'))
    # minpot = float(
    #   input('Insert the lower value of chemical potential you want to plot in eV'))
    # maxpot = float(
    #   input('Inser the higher value of chemical potential you want to plot in eV'))

    direction = direction.lower().replace('_', '')

    if direction.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid direction')

    if direction == 'sxx':
        col = 3
    elif direction == 'sxy':
        col = 4
    elif direction == 'sxz':
        col = 5
    elif direction == 'syx':
        col = 6
    elif direction == 'syy':
        col = 7
    elif direction == 'syz':
        col = 8
    elif direction == 'szx':
        col = 9
    elif direction == 'szy':
        col = 10
    elif direction == 'szz':
        col = 11

    else:
        sys.exit('please, choose a valid direction')

    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    i = 0
    for n in seebeck:
        vol = n.volume

        listtemp = []
        for k in range(0, len(n.all_data)):
            listtemp.append(n.temp[k])
        if temperature != 'all' and temperature not in listtemp:
            sys.exit('Please, choose a valid temperature.')

        x = []
        for kq in range(0, len(n.all_data)):
            x.append(np.array(n.all_data[kq].apply(
                lambda x: float(x.split()[0]))))

        carrier = []
        for kq in range(0, len(n.all_data)):
            carrier.append(np.array(n.all_data[kq].apply(
                lambda x: (float(x.split()[2])/vol))))

        y = []
        for kq in range(0, len(n.all_data)):
            y.append(np.array(n.all_data[kq].apply(
                lambda x: float(x.split()[col])*1000000)))

        yneg = []
        ypos = []
        xpos = []
        xneg = []
        yposfin = []
        xposfin = []
        ynegfin = []
        xnegfin = []

        for kq in range(0, len(n.all_data)):
            for j in range(0, len(carrier[kq])):
                if carrier[kq][j] >= 0:
                    xpos.append(x[kq][j])
                    ypos.append(y[kq][j])
                else:
                    xneg.append(x[kq][j])
                    yneg.append(y[kq][j])
            yposfin.append(ypos)
            ynegfin.append(yneg)
            xposfin.append(xpos)
            xnegfin.append(xneg)
            xpos = []
            ypos = []
            xneg = []
            yneg = []

        colours = ['royalblue', 'orange', 'green', 'red',
                   'purple', 'brown', 'pink', 'grey', 'olive', 'cyan']

        endx = []
        endy = []

        k = listtemp.index(temperature)
        endx = [xposfin[k][-1], xnegfin[k][0]]
        endy = [yposfin[k][-1], ynegfin[k][0]]
        plt.plot(endx, endy, color=colours[i])
        plt.plot(xposfin[k], yposfin[k], color=colours[i], label=str(n.title))
        plt.plot(xnegfin[k], ynegfin[k], '--', color=colours[i])
        plt.xlabel('Chemical Potential (eV)', fontsize=12)
        plt.ylabel('Seebeck Coefficient ($\mu$V/K)', fontsize=12)
        plt.xlim(minpot, maxpot)
        plt.axhline(0, color='k')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
        i = 1+i
    plt.title('MultiSeebeck ' + str(n.temp[k]) + ' K')
    # plt.savefig('multiseebeck' + time.strftime("%Y-%m-%d_%H%M%S") +
    # '.jpg', format='jpg', dpi=100, bbox_inches='tight')
    plt.show()
    fig = plt.gcf()
    return fig


# -------------------------------------SIGMA-----------------------------------#

def plot_cry_sigma_potential(sigma_obj, direction, temperature):
    """
    Plot the electrical conductivity as a function of chemical potential.

    Args:
        sigma_obj (object): Sigma object containing the data for electrical conductivity.
        direction (str): choose the direction to plot among 'S_xx', 'S_xy', 'S_xz', 'S_yx', 'S_yy', 'S_yz', 'S_yz', 'S_zx', 'S_zy', 'S_zz'.
        temperature (value/str): choose the temperature to be considered or 'all' to consider them all together

    Returns:
        Returns figure object

    Notes:
        - Plots the electrical conductivity as a function of chemical potential the selected temperature.
        - Distinguishes between n-type and p-type conduction with dashed and solid lines, respectively.

    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    direction = direction.lower().replace('_', '')

    if direction.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid chioce')

    if direction == 'sxx':
        col = 3
    elif direction == 'sxy':
        col = 4
    elif direction == 'sxz':
        col = 5
    elif direction == 'syy':
        col = 6
    elif direction == 'syz':
        col = 7
    elif direction == 'szz':
        col = 8
    else:
        sys.exit('please, choose a valid chioce')

    listtemp = []

    for k in range(0, len(sigma_obj.all_data)):
        listtemp.append(sigma_obj.temp[k])

    if temperature != 'all' and temperature not in listtemp:
        sys.exit('Please, choose a valid temperature.')

    vol = sigma_obj.volume

    x = []
    for k in range(0, len(sigma_obj.all_data)):
        x.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: float(x.split()[0]))))

    carrier = []
    for k in range(0, len(sigma_obj.all_data)):
        carrier.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: (float(x.split()[2])/vol))))

    y = []
    for k in range(0, len(sigma_obj.all_data)):
        y.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: float(x.split()[col]))))

    yneg = []
    ypos = []
    xpos = []
    xneg = []
    yposfin = []
    xposfin = []
    ynegfin = []
    xnegfin = []

    for k in range(0, len(sigma_obj.all_data)):
        for j in range(0, len(x[k])):
            if carrier[k][j] >= 0:
                xpos.append(x[k][j])
                ypos.append(y[k][j])
            else:
                xneg.append(x[k][j])
                yneg.append(y[k][j])
        yposfin.append(ypos)
        ynegfin.append(yneg)
        xposfin.append(xpos)
        xnegfin.append(xneg)
        xpos = []
        ypos = []
        xneg = []
        yneg = []

    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')
    colours = []
    colours = ['royalblue', 'orange', 'green', 'red',
               'purple', 'brown', 'pink', 'grey', 'olive', 'cyan']
    endx = []
    endy = []

    if temperature != 'all':
        k = listtemp.index(temperature)
        endx = [xposfin[k][-1], xnegfin[k][0]]
        endy = [yposfin[k][-1], ynegfin[k][0]]
        plt.figure()
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xposfin[k], yposfin[k], color=colours[k],
                 label=str(sigma_obj.temp[k])+' K')
        plt.plot(xnegfin[k], ynegfin[k], '--', color=colours[k])
        plt.xlabel('Chemical Potential (eV)', fontsize=12)
        plt.ylabel('Electrical Conductivity (S/m)', fontsize=12)
        plt.axhline(0, color='k')
        plt.title('Sigma at '+str(sigma_obj.temp[k]) + 'K')
        plt.legend(loc='upper left', fontsize=12)
        # plt.savefig('sigma_potential_at_' + str(sigma_obj.temp[k]) + 'K___' + time.strftime(
        #   "%Y-%m-%d_%H%M%S") + '.jpg', format='jpg', dpi=600, bbox_inches='tight')
        plt.show()
        fig = plt.gcf()
        return fig

    else:
        from matplotlib.pyplot import figure

        figure(figsize=(7, 7))

        for i in range(0, len(sigma_obj.all_data)):
            endx = [xposfin[i][-1], xnegfin[i][0]]
            endy = [yposfin[i][-1], ynegfin[i][0]]
            plt.plot(endx, endy, color=colours[i])
            plt.plot(xposfin[i], yposfin[i], color=colours[i],
                     label=str(sigma_obj.temp[i])+' K')
            plt.plot(xnegfin[i], ynegfin[i], '--', color=colours[i])
            plt.xlabel('Chemical Potential (eV)', fontsize=12)
            plt.ylabel('Electrical Conductivity (S/m)', fontsize=12)
            plt.title('Sigma at different T')
            plt.axhline(0, color='k')
            plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
        # plt.savefig('sigma_potential_different_T_' + time.strftime("%Y-%m-%d_%H%M%S") +
         #           '.jpg', format='jpg', dpi=100, bbox_inches='tight')
        plt.show()
        fig = plt.gcf()
        return fig


def plot_cry_sigma_carrier(sigma_obj, direction, temperature):
    """
    Plot the electrical conductivity as a function of charge carrier concentration.

   Args:
        sigma_obj (object): Sigma object containing the data for electrical conductivity.
        direction (str): choose the direction to plot among 'S_xx', 'S_xy', 'S_xz', 'S_yx', 'S_yy', 'S_yz', 'S_yz', 'S_zx', 'S_zy', 'S_zz'.
        temperature (value/str): choose the temperature to be considered or 'all' to consider them all together


    Returns:
        Returns Figure object

    Notes:
        - Plots the electrical conductivity as a function of charge carrier concentration the selected temperature.
        - Distinguishes between n-type and p-type conduction with dashed and solid lines, respectively.
    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    direction = direction.lower().replace('_', '')

    if direction.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid chioce')

    if direction == 'sxx':
        col = 3
    elif direction == 'sxy':
        col = 4
    elif direction == 'sxz':
        col = 5
    elif direction == 'syy':
        col = 6
    elif direction == 'syz':
        col = 7
    elif direction == 'szz':
        col = 8
    else:
        sys.exit('please, choose a valid chioce')

    listtemp = []

    for k in range(0, len(sigma_obj.all_data)):
        listtemp.append(sigma_obj.temp[k])

    if temperature != 'all' and temperature not in listtemp:
        sys.exit('Please, choose a valid temperature.')

    vol = sigma_obj.volume

    x = []
    for k in range(0, len(sigma_obj.all_data)):
        x.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: (float(x.split()[2])/vol))))

    y = []
    for k in range(0, len(sigma_obj.all_data)):
        y.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: float(x.split()[col]))))

    yneg = []
    ypos = []
    xpos = []
    xneg = []
    yposfin = []
    xposfin = []
    ynegfin = []
    xnegfin = []

    for k in range(0, len(sigma_obj.all_data)):
        for j in range(0, len(x[k])):
            if x[k][j] >= 0:
                xpos.append(x[k][j])
                ypos.append(y[k][j])
            else:
                xneg.append(x[k][j])
                yneg.append(y[k][j])
        yposfin.append(ypos)
        ynegfin.append(yneg)
        xposfin.append(xpos)
        xnegfin.append(xneg)
        xpos = []
        ypos = []
        xneg = []
        yneg = []

    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    colours = []
    colours = ['royalblue', 'orange', 'green', 'red',
               'purple', 'brown', 'pink', 'grey', 'olive', 'cyan']

    if temperature != 'all':
        k = listtemp.index(temperature)
        endx = [xposfin[k][-1], xnegfin[k][0]]
        endy = [yposfin[k][-1], ynegfin[k][0]]
        plt.figure()
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xposfin[k], yposfin[k], color=colours[k],
                 label=str(sigma_obj.temp[k])+' K')
        plt.plot(abs(np.array(xnegfin[k])), ynegfin[k], '--', color=colours[k])
        plt.xlabel('Charge Carrier Concentration (cm$^{-3}$)', fontsize=12)
        plt.ylabel('Electrical Conductivity (S/m)', fontsize=12)
        plt.axhline(0, color='k')
        plt.title('Sigma at ' + str(sigma_obj.temp[k]) + 'K')
        plt.legend(loc='upper left', fontsize=12)
        plt.xscale('log')
        plt.savefig('sigma_carrier_at_' + str(sigma_obj.temp[k]) + 'K___' + time.strftime(
            "%Y-%m-%d_%H%M%S") + '.jpg', format='jpg', dpi=600, bbox_inches='tight')
        plt.show()
        fig = plt.gcf()
        return fig
    else:

        for i in range(0, len(sigma_obj.all_data)):
            endx = [xposfin[i][-1], xnegfin[i][0]]
            endy = [yposfin[i][-1], ynegfin[i][0]]
            plt.plot(endx, endy, color=colours[k])
            plt.plot(xposfin[i], yposfin[i], color=colours[i],
                     label=str(sigma_obj.temp[i])+' K')
            plt.plot(abs(np.array(xnegfin[i])),
                     ynegfin[i], '--', color=colours[i])
            plt.xlabel('Charge Carrier Concentration (cm$^{-3}$)', fontsize=12)
            plt.ylabel('Electrical Conductivity (S/m)', fontsize=12)
            plt.title('Sigma at different T')
            plt.axhline(0, color='k')
            plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
            plt.xscale('log')
        # plt.savefig('sigma_carrier_different_T_' + time.strftime("%Y-%m-%d_%H%M%S") +
        #           '.jpg', format='jpg', dpi=100, bbox_inches='tight')
        plt.show()
        fig = plt.gcf()
        return fig


def plot_cry_multisigma(direction, temperature, minpot, maxpot, *sigma):
    """
    Plot the electron cinductivity from different files.

    Args:
        direction (str): choose the direction to plot among 'S_xx', 'S_xy', 'S_xz', 'S_yx', 'S_yy', 'S_yz', 'S_yz', 'S_zx', 'S_zy', 'S_zz'.
        temperature (value): choose the temperature to be considered 
        minpot (value): lower value of chemical potential you want to plot in eV  
        maxpot (value): higher value of chemical potential you want to plot in eV 
        *seebeck (obj): Variable number of seebeck objects containing the data for the electron conductivity (sigma).

    Returns:
        Figure object

    Notes:
        - Plots the electron conductivity for each sigma object.
        - Differentiates transport coefficients due to n-type or p-type conduction using dashed and solid lines.
    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    direction = direction.lower().replace('_', '')

    if direction.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid direction')

    if direction == 'sxx':
        col = 3
    elif direction == 'sxy':
        col = 4
    elif direction == 'sxz':
        col = 5
    elif direction == 'syy':
        col = 6
    elif direction == 'syz':
        col = 7
    elif direction == 'szz':
        col = 8

    else:
        sys.exit('please, choose a valid direction')

    i = 0
    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')
    for n in sigma:
        vol = n.volume
        listtemp = []
        for k in range(0, len(n.all_data)):
            listtemp.append(n.temp[k])
        if temperature != 'all' and temperature not in listtemp:
            sys.exit('Please, choose a valid temperature.')

        x = []
        for kq in range(0, len(n.all_data)):
            x.append(np.array(n.all_data[kq].apply(
                lambda x: float(x.split()[0]))))

        carrier = []
        for kq in range(0, len(n.all_data)):
            carrier.append(np.array(n.all_data[kq].apply(
                lambda x: (float(x.split()[2])/vol))))

        y = []
        for kq in range(0, len(n.all_data)):
            y.append(np.array(n.all_data[kq].apply(
                lambda x: float(x.split()[col]))))

        yneg = []
        ypos = []
        xpos = []
        xneg = []
        yposfin = []
        xposfin = []
        ynegfin = []
        xnegfin = []

        for kq in range(0, len(n.all_data)):
            for j in range(0, len(x[kq])):
                if carrier[kq][j] >= 0:
                    xpos.append(x[kq][j])
                    ypos.append(y[kq][j])
                else:
                    xneg.append(x[kq][j])
                    yneg.append(y[kq][j])
            yposfin.append(ypos)
            ynegfin.append(yneg)
            xposfin.append(xpos)
            xnegfin.append(xneg)
            xpos = []
            ypos = []
            xneg = []
            yneg = []

        colours = []
        colours = ['royalblue', 'orange', 'green', 'red',
                   'purple', 'brown', 'pink', 'grey', 'olive', 'cyan']
        endx = []
        endy = []
        k = listtemp.index(temperature)

        endx = [xposfin[k][-1], xnegfin[k][0]]
        endy = [yposfin[k][-1], ynegfin[k][0]]
        plt.plot(endx, endy, color=colours[i])
        plt.plot(xposfin[k], yposfin[k], color=colours[i], label=str(n.title))
        plt.plot(xnegfin[k], ynegfin[k], '--', color=colours[i])
        plt.xlabel('Chemical Potential (eV)', fontsize=12)
        plt.ylabel('Electrical Conductivity (S/m)', fontsize=12)
        plt.axhline(0, color='k')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
        plt.xlim(minpot, maxpot)
        i = 1+i
    plt.title('MultiSigma ' + str(sigma[0].temp[k]) + ' K')
   # plt.savefig('multisigma' + time.strftime("%Y-%m-%d_%H%M%S") +
    # '.jpg', format='jpg', dpi=100, bbox_inches='tight')
    plt.show()
    fig = plt.gcf()
    return fig


# --------------------------------POWERFACTOR----------------------------------#

def plot_cry_powerfactor_potential(seebeck_obj, sigma_obj, direction, temperature):
    """
    Plot the power factor for different potentials.

    Args:
        seebeck_obj (obj): Seebeck object containing the data for the Seebeck coefficient.
        sigma_obj (obj): Sigma object containing the data for the electrical conductivity.
        direction (str): choose the direction to plot among 'PF_xx', 'PF_xy', 'PF_xz', 'PF_yx', 'PF_yy', 'PF_yz', 'PF_yz', 'PF_zx', 'PF_zy', 'PF_zz'.
        temperature (value/str): choose the temperature to be considered or 'all' to consider them all together

    Returns:
        Figure object

    Notes:
        - Calculates the power factor using the Seebeck coefficient and electrical conductivity data for each temperature.
        - Plots the power factor for each temperature as a function of the chemical potential, distinguishing between n-type and p-type conduction.
    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    direction = direction.lower().replace('_', '')

    if direction.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid direction')

    if direction == 'pfxx':
        col = 3
    elif direction == 'pfxy':
        col = 4
    elif direction == 'pfxz':
        col = 5
    elif direction == 'pfyx':
        col = 6
    elif direction == 'pfyy':
        col = 7
    elif direction == 'pfyz':
        col = 8
    elif direction == 'pfzx':
        col = 9
    elif direction == 'pfzy':
        col = 10
    elif direction == 'pfzz':
        col = 11
    else:
        sys.exit('please, choose a valid direction')

    if direction == 'pfxx':
        cols = 3
    elif direction == 'pfxy':
        cols = 4
    elif direction == 'pfxz':
        cols = 5
    elif direction == 'pfyx':
        cols = 4
    elif direction == 'pfyy':
        cols = 6
    elif direction == 'pfyz':
        cols = 7
    elif direction == 'pfzx':
        cols = 5
    elif direction == 'pfzy':
        cols = 7
    elif direction == 'pfzz':
        cols = 8
    else:
        sys.exit('please, choose a valid direction')

    listtemp = []
    for k in range(0, len(seebeck_obj.all_data)):
        listtemp.append(seebeck_obj.temp[k])
    if temperature != 'all' and temperature not in listtemp:
        sys.exit('Please, choose a valid temperature.')
    listtemps = []
    for k in range(0, len(sigma_obj.all_data)):
        listtemps.append(sigma_obj.temp[k])
    if temperature != 'all' and temperature not in listtemps:
        sys.exit('Please, choose a valid temperature.')

    x = []
    for k in range(0, len(seebeck_obj.all_data)):
        x.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: float(x.split()[0]))))

    vol = sigma_obj.volume
    carrier = []

    yse = []
    for k in range(0, len(seebeck_obj.all_data)):
        yse.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: float(x.split()[col]))))

    ysi = []
    for k in range(0, len(sigma_obj.all_data)):
        ysi.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: float(x.split()[cols]))))

    carrier = []
    for k in range(0, len(seebeck_obj.all_data)):
        carrier.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: (float(x.split()[2])/vol))))

    ysineg = []
    ysipos = []
    xsipos = []
    xsineg = []
    ysiposfin = []
    xsiposfin = []
    ysinegfin = []
    xsinegfin = []

    yseneg = []
    ysepos = []
    xsepos = []
    xseneg = []
    yseposfin = []
    xseposfin = []
    ysenegfin = []
    xsenegfin = []

    for k in range(0, len(seebeck_obj.all_data)):
        for j in range(0, len(carrier[k])):
            if carrier[k][j] >= 0:
                xsipos.append(x[k][j])
                ysipos.append(ysi[k][j])
            else:
                xsineg.append(x[k][j])
                ysineg.append(ysi[k][j])
        ysiposfin.append(np.array(ysipos))
        ysinegfin.append(np.array(ysineg))
        xsiposfin.append(np.array(xsipos))
        xsinegfin.append(np.array(xsineg))
        xsipos = []
        ysipos = []
        xsineg = []
        ysineg = []

    for k in range(0, len(seebeck_obj.all_data)):
        for j in range(0, len(carrier[k])):
            if carrier[k][j] >= 0:
                xsepos.append(x[k][j])
                ysepos.append(yse[k][j])
            else:
                xseneg.append(x[k][j])
                yseneg.append(yse[k][j])

        yseposfin.append(np.array(ysepos))
        ysenegfin.append(np.array(yseneg))
        xseposfin.append(np.array(xsepos))
        xsenegfin.append(np.array(xseneg))
        xsepos = []
        ysepos = []
        xseneg = []
        yseneg = []

    pf_meta_pos = []
    for i in range(0, len(yseposfin)):
        pf_meta_pos.append(yseposfin[i]*yseposfin[i])

    pf_pos = []
    for i in range(0, len(pf_meta_pos)):
        pf_pos.append(pf_meta_pos[i] * ysiposfin[i])

    pf_meta_neg = []
    for i in range(0, len(ysenegfin)):
        pf_meta_neg.append(ysenegfin[i] * ysenegfin[i])

    pf_neg = []
    for i in range(0, len(pf_meta_neg)):
        pf_neg.append(pf_meta_neg[i] * ysinegfin[i])

    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    colours = []
    colours = ['royalblue', 'orange', 'green', 'red',
               'purple', 'brown', 'pink', 'grey', 'olive', 'cyan']

    if temperature != 'all':
        k = listtemp.index(temperature)
        endx = [xsiposfin[k][-1], xsinegfin[k][0]]
        endy = [pf_pos[k][-1], pf_neg[k][0]]
        plt.figure()
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xsiposfin[k], pf_pos[k], color=colours[k],
                 label=str(seebeck_obj.temp[k])+' K')
        plt.plot(xsinegfin[k], pf_neg[k], '--', color=colours[k])
        plt.xlabel('Chemical Potential (eV)', fontsize=12)
        plt.ylabel('Power Factor (10$^{-12}$WK$^{-2}$m$^{-1}$)', fontsize=12)
        plt.axhline(0, color='k')
        plt.title('Power Factor at ' + str(seebeck_obj.temp[k]) + ' K')
        plt.legend(loc='upper left', fontsize=12)
        plt.show()
        fig = plt.gcf
        return fig

    else:

        from matplotlib.pyplot import figure

        figure(figsize=(7, 7))
        for k in range(0, len(seebeck_obj.all_data)):
            endx = [xsiposfin[k][-1], xsinegfin[k][0]]
            endy = [pf_pos[k][-1], pf_neg[k][0]]
            plt.plot(endx, endy, color=colours[k])
            plt.plot(xsiposfin[k], pf_pos[k], color=colours[k],
                     label=str(seebeck_obj.temp[k])+' K')
            plt.plot(xsinegfin[k], pf_neg[k], '--', color=colours[k])
            plt.xlabel('Chemical Potential (eV)', fontsize=12)
            plt.ylabel(
                'Power Factor (10$^{-12}$WK$^{-2}$m$^{-1}$)', fontsize=12)
            plt.axhline(0, color='k')
            plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
            plt.title('Power Factor at different T')

        # plt.savefig('powerfactor_potential_different_T_' + time.strftime(
         #   "%Y-%m-%d_%H%M%S") + '.jpg', format='jpg', dpi=100, bbox_inches='tight')
        plt.show()
        fig = plt.gcf
        return fig


def plot_cry_powerfactor_carrier(seebeck_obj, sigma_obj, direction, temperature):
    """
    Plot the power factor for different charge carrier concentrations.

    Args:
        seebeck_obj (obj): Seebeck object containing the data for the Seebeck coefficient.
        sigma_obj (obj): Sigma object containing the data for the electrical conductivity.
        direction (str): choose the direction to plot among 'PF_xx', 'PF_xy', 'PF_xz', 'PF_yx', 'PF_yy', 'PF_yz', 'PF_yz', 'PF_zx', 'PF_zy', 'PF_zz'.
        temperature (value/str): choose the temperature to be considered or 'all' to consider them all together

    Returns:
        Figure object

    Notes:
        - Calculates the power factor using the Seebeck coefficient and electrical conductivity data for each temperature.
        - Plots the power factor for each temperature as a function of the charge carrier concentration, distinguishing between n-type and p-type conduction.
    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    direction = direction.lower().replace('_', '')

    if direction.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid direction')

    if direction == 'pfxx':
        col = 3
    elif direction == 'pfxy':
        col = 4
    elif direction == 'pfxz':
        col = 5
    elif direction == 'pfyx':
        col = 6
    elif direction == 'pfyy':
        col = 7
    elif direction == 'pfyz':
        col = 8
    elif direction == 'pfzx':
        col = 9
    elif direction == 'pfzy':
        col = 10
    elif direction == 'pfzz':
        col = 11
    else:
        sys.exit('please, choose a valid direction')

    if direction == 'pfxx':
        cols = 3
    elif direction == 'pfxy':
        cols = 4
    elif direction == 'pfxz':
        cols = 5
    elif direction == 'pfyx':
        cols = 4
    elif direction == 'pfyy':
        cols = 6
    elif direction == 'pfyz':
        cols = 7
    elif direction == 'pfzx':
        cols = 5
    elif direction == 'pfzy':
        cols = 7
    elif direction == 'pfzz':
        cols = 8
    else:
        sys.exit('please, choose a valid direction')

    vol = sigma_obj.volume

    listtemp = []
    for k in range(0, len(seebeck_obj.all_data)):
        listtemp.append(seebeck_obj.temp[k])
    if temperature != 'all' and temperature not in listtemp:
        sys.exit('Please, choose a valid temperature.')
    listtemps = []
    for k in range(0, len(sigma_obj.all_data)):
        listtemps.append(sigma_obj.temp[k])
    if temperature != 'all' and temperature not in listtemps:
        sys.exit('Please, choose a valid temperature.')

    x = []
    for k in range(0, len(sigma_obj.all_data)):
        x.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: (float(x.split()[2])/vol))))

    yse = []
    for k in range(0, len(seebeck_obj.all_data)):
        yse.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: float(x.split()[col]))))

    ysi = []
    for k in range(0, len(sigma_obj.all_data)):
        ysi.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: float(x.split()[cols]))))

    pf_meta = []
    for i in range(0, len(yse)):
        pf_meta.append(yse[i] * yse[i])

    pf = []
    for i in range(0, len(pf_meta)):
        pf.append(pf_meta[i] * ysi[i])

    ysineg = []
    ysipos = []
    xsipos = []
    xsineg = []
    ysiposfin = []
    xsiposfin = []
    ysinegfin = []
    xsinegfin = []

    yseneg = []
    ysepos = []
    xsepos = []
    xseneg = []
    yseposfin = []
    xseposfin = []
    ysenegfin = []
    xsenegfin = []

    for k in range(0, len(seebeck_obj.all_data)):
        for j in range(0, len(x[k])):
            if x[k][j] >= 0:
                xsipos.append(x[k][j])
                ysipos.append(ysi[k][j])
            else:
                xsineg.append(x[k][j])
                ysineg.append(ysi[k][j])
        ysiposfin.append(np.array(ysipos))
        ysinegfin.append(np.array(ysineg))
        xsiposfin.append(np.array(xsipos))
        xsinegfin.append(np.array(xsineg))
        xsipos = []
        ysipos = []
        xsineg = []
        ysineg = []

    for k in range(0, len(seebeck_obj.all_data)):
        for j in range(0, len(x[k])):
            if x[k][j] >= 0:
                xsepos.append(x[k][j])
                ysepos.append(yse[k][j])
            else:
                xseneg.append(x[k][j])
                yseneg.append(yse[k][j])
        yseposfin.append(np.array(ysepos))
        ysenegfin.append(np.array(yseneg))
        xseposfin.append(np.array(xsepos))
        xsenegfin.append(np.array(xseneg))
        xsepos = []
        ysepos = []
        xseneg = []
        yseneg = []

    pf_meta_pos = []
    for i in range(0, len(yseposfin)):
        pf_meta_pos.append(yseposfin[i]*yseposfin[i])

    pf_pos = []
    for i in range(0, len(pf_meta_pos)):
        pf_pos.append(pf_meta_pos[i] * ysiposfin[i])

    pf_meta_neg = []
    for i in range(0, len(ysenegfin)):
        pf_meta_neg.append(ysenegfin[i] * ysenegfin[i])

    pf_neg = []
    for i in range(0, len(pf_meta_neg)):
        pf_neg.append(pf_meta_neg[i] * ysinegfin[i])

    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    colours = []
    colours = ['royalblue', 'orange', 'green', 'red',
               'purple', 'brown', 'pink', 'grey', 'olive', 'cyan']

    from matplotlib.pyplot import figure

    if temperature != 'all':
        k = listtemp.index(temperature)
        endx = [xsiposfin[k][-1], xsinegfin[k][0]]
        endy = [pf_pos[k][-1], pf_neg[k][0]]
        plt.figure()
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xsiposfin[k], pf_pos[k], color=colours[k],
                 label=str(seebeck_obj.temp[k])+' K')
        plt.plot(abs(xsinegfin[k]), pf_neg[k], '--', color=colours[k])
        plt.xlabel('Charge Carrier Concentration (cm$^{-3}$)', fontsize=12)
        plt.ylabel('Power Factor (10$^{-12}$WK$^{-2}$m$^{-1}$)', fontsize=12)
        plt.axhline(0, color='k')
        plt.title('Power Factor at ' + str(seebeck_obj.temp[k]) + ' K')
        plt.legend(loc='upper left', fontsize=12)
        plt.xscale('log')
        # plt.savefig('powerfactor_carrier_at_' + str(seebeck_obj.temp[k]) + 'K___' + time.strftime(
        #   "%Y-%m-%d_%H%M%S") + '.jpg', format='jpg', dpi=600, bbox_inches='tight')
        plt.show()
        fig = plt.gcf
        return fig
    else:

        from matplotlib.pyplot import figure
        figure(figsize=(7, 7))
        for k in range(0, len(seebeck_obj.all_data)):
            endx = [xsiposfin[k][-1], xsinegfin[k][0]]
            endy = [pf_pos[k][-1], pf_neg[k][0]]
            plt.plot(endx, endy, color=colours[k])
            plt.plot(xsiposfin[k], pf_pos[k], color=colours[k],
                     label=str(seebeck_obj.temp[k])+' K')
            plt.plot(abs(xsinegfin[k]), pf_neg[k], '--', color=colours[k])
            plt.xlabel('Charge Carrier Concentration (cm$^{-3}$)', fontsize=12)
            plt.ylabel(
                'Power Factor (10$^{-12}$WK$^{-2}$m$^{-1}$)', fontsize=12)
            plt.title('Power Factor at different T')
            plt.axhline(0, color='k')
            plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
            plt.xscale('log')
        # plt.savefig('powerfactor_carrier_different_T_' + time.strftime("%Y-%m-%d_%H%M%S") +
         #           '.jpg', format='jpg', dpi=100, bbox_inches='tight')
        plt.show()
        fig = plt.gcf
        return fig


# -------------------------------------ZT--------------------------------------#


def plot_cry_zt(seebeck_obj, sigma_obj, direction, temperature, ktot):
    """
    Plot the ZT value for different temperatures.

    Args:
        seebeck_obj (obj): Seebeck object containing the data for the Seebeck coefficient.
        sigma_obj (obj): Sigma object containing the data for the electrical conductivity.
        direction (str): choose the direction to plot among 'ZT_xx', 'ZT_xy', 'ZT_xz', 'ZT_yx', 'ZT_yy', 'ZT_yz', 'ZT_yz', 'ZT_zx', 'ZT_zy', 'ZT_zz'.
        temperature (value/str): choose the temperature to be considered or 'all' to consider them all together
        ktot (value): alue of the total thermal conductivity (ktot) in W-1K-1m-1

    Returns:
        Figure object

    Notes:
        - Calculates the ZT value using the Seebeck coefficient and electrical conductivity data.
        - Plots the ZT value for each temperature as a function of the chemical potential.
    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    direction = direction.lower().replace('_', '')

    if direction.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid direction')

    if direction == 'ztxx':
        col = 3
    elif direction == 'ztxy':
        col = 4
    elif direction == 'ztxz':
        col = 5
    elif direction == 'ztyx':
        col = 6
    elif direction == 'ztyy':
        col = 7
    elif direction == 'ztyz':
        col = 8
    elif direction == 'ztzx':
        col = 9
    elif direction == 'ztzy':
        col = 10
    elif direction == 'ztzz':
        col = 11
    else:
        sys.exit('please, choose a valid direction')

    if direction == 'ztxx':
        cols = 3
    elif direction == 'ztxy':
        cols = 4
    elif direction == 'ztxz':
        cols = 5
    elif direction == 'ztyx':
        cols = 4
    elif direction == 'ztyy':
        cols = 6
    elif direction == 'ztyz':
        cols = 7
    elif direction == 'ztzx':
        cols = 5
    elif direction == 'ztzy':
        cols = 7
    elif direction == 'ztzz':
        cols = 8
    else:
        sys.exit('please, choose a valid direction')

    listtemp = []
    for k in range(0, len(seebeck_obj.all_data)):
        listtemp.append(seebeck_obj.temp[k])
    if temperature != 'all' and temperature not in listtemp:
        sys.exit('Please, choose a valid temperature.')
    listtemps = []
    for k in range(0, len(sigma_obj.all_data)):
        listtemps.append(sigma_obj.temp[k])
    if temperature != 'all' and temperature not in listtemps:
        sys.exit('Please, choose a valid temperature.')

    x = []
    for k in range(0, len(seebeck_obj.all_data)):
        x.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: float(x.split()[0]))))

    yse = []
    for k in range(0, len(seebeck_obj.all_data)):
        yse.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: float(x.split()[col]))))

    ysi = []
    for k in range(0, len(sigma_obj.all_data)):
        ysi.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: float(x.split()[cols]))))

    pf_meta = []
    for i in range(0, len(yse)):
        pf_meta.append(yse[i] * yse[i])

    pf = []
    for i in range(0, len(pf_meta)):
        pf.append(pf_meta[i] * ysi[i])

    zt = []
    for i in range(0, len(pf_meta)):
        zt.append((pf[i] * seebeck_obj.temp[i])/ktot)

    if temperature != 'all':
        k = listtemp.index(temperature)
        plt.figure()
        plt.plot(x[k], pf[k], label=str(seebeck_obj.temp[k])+' K')
        plt.xlabel('Chemical Potential (eV)', fontsize=12)
        plt.ylabel('ZT', fontsize=12)
        plt.axhline(0, color='k')
        plt.title('ZT at ' + str(seebeck_obj.temp[k]) + ' K')
        plt.legend(loc='upper left', fontsize=12)
        # plt.savefig('zt_at_' + str(seebeck_obj.temp[k]) + 'K___' + time.strftime(
        #   "%Y-%m-%d_%H%M%S") + '.jpg', format='jpg', dpi=600, bbox_inches='tight')
        plt.show()
        fig = plt.gcf
        return fig
    else:
        for k in range(0, len(seebeck_obj.all_data)):
            plt.plot(x[k], pf[k], label=str(seebeck_obj.temp[k])+' K')
            plt.xlabel('Chemical Potential (eV)', fontsize=12)
            plt.ylabel('ZT', fontsize=12)
            plt.title('ZT at different T')
            plt.axhline(0, color='k')
            plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
        # plt.savefig('zt_different_T_' + time.strftime("%Y-%m-%d_%H%M%S") +
            # '.jpg', format='jpg', dpi=100, bbox_inches='tight')
        plt.show()
        fig = plt.gcf
        return fig


##############################################################################
#                                                                            #
#                             ELASTIC PROPERTIES                             #
#                                                                            #
##############################################################################

# --------------------------------YOUNG MODULUS--------------------------------#

def plot_cry_young(theta, phi, S):
    """
    Compute Young's modulus for each direction of the space (i.e., each pair
    of theta and phi angles).

    Args:
        theta (float): Theta value.
        phi (float): Phi value.
        S (numpy.ndarray): Compliance matrix.

    Returns:
        float: Young's modulus values.

    Notes:
        - This function is intended to be called by cry_ela_plot
    """
    import numpy as np

    # C2V = Matrix to refer the Cartesian into Voigt's notation
    # Observe that the matrix should be written as is shown below
    # C2V = np.array([[1,6,5],[6,2,4],[5,4,3]])
    # Since python start counting from zero all numbers must be subtracted by 1
    C2V = np.array(
        [
            [0, 5, 4],
            [5, 1, 3],
            [4, 3, 2]
        ]
    )
    # print("The Matrix to convert Cartesian into Voigs Notation: \n", C2V)
    # creating the 1x3 vector "a"
    a = np.array(
        [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ]
    )
    # e is a pseudo Young modulus value folowing the relation 1/e
    e = 0.0
    # i,j,k,l are updatable indices refering to cartesian notation that
    # will be converted by C2V into Voigt's
    for i in range(3):
        for j in range(3):
            v = C2V[i, j]
            for k in range(3):
                for l in range(3):
                    u = C2V[k, l]
                    # rf is a factor that must be multipled by the compliance element if
                    # certain conditions are satisfied
                    rf = 1
                    if v >= 3 and u >= 3:
                        rf = 4
                    if v >= 3 and u < 3:
                        rf = 2
                    if u >= 3 and v < 3:
                        rf = 2

                    rtmp = a[i] * a[j] * a[k] * a[l] * (S[v, u] / rf)
                    e = e + rtmp
    E_tmp = 1 / e  # is the Young Modulus of each cycle
    return E_tmp

# ----------------------------COMPRESSION PROPERTIES---------------------------#


def plot_cry_comp(theta, phi, S):
    """
    Compute linear compressibility for each direction of the space (i.e., each
    pair of theta and phi angles).

    Args:
        theta (float): Theta value.
        phi (float): Phi value.
        S (numpy.ndarray): Compliance matrix.

    Returns:
        float: Linear compressibility values.

    Notes:
        - This function is intended to be called by cry_ela_plot
    """
    import numpy as np

    C2V = np.array(
        [
            [0, 5, 4],
            [5, 1, 3],
            [4, 3, 2]
        ]
    )

    # Paper u vector
    a = np.array(
        [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ]
    )
    B = 0.0

    for i in range(3):
        for j in range(3):
            v = C2V[i, j]
            for k in range(3):
                u = C2V[k, k]
                rf = 1
                if v >= 3 and u >= 3:
                    rf = 4
                if v >= 3 and u < 3:
                    rf = 2
                if u >= 3 and v < 3:
                    rf = 2

                # beta(theta, phi) = u_i * u_j * S_{ij}
                rtmp = a[i] * a[j] * (S[v, u] / rf)
                B = B + rtmp
    return B


# --------------------------------SHEAR MODULUS--------------------------------#

def plot_cry_shear(theta_1D, phi_1D, S, ndeg, shear_choice):
    """
    For each direction of the space (i.e., for each pair
    of theta and phi angles) the shear modulus is computed for the third angle
    chi and the average, maximum and minimum values are stored.

    Args:
        theta_1D (numpy.ndarray): One-dimensional array of theta values.
        phi_1D (numpy.ndarray): One-dimensional array of phi values.
        S (numpy.ndarray): Compliance matrix.
        ndeg (int): Number of degrees for discretization.
        shear_choice (str): Type of shear property to plot. Options: "avg", "min", "max".

    Returns:
        numpy.ndarray: Shear property array.

    Notes:
        - This function is intended to be called by cry_ela_plot
    """
    import numpy as np

    C2V = np.array(
        [
            [0, 5, 4],
            [5, 1, 3],
            [4, 3, 2],
        ]
    )
    shear_chi = np.zeros(ndeg)
    shear_min = np.zeros((ndeg, ndeg))
    shear_max = np.zeros((ndeg, ndeg))
    shear_avg = np.zeros((ndeg, ndeg))
    chi_1D = np.linspace(0, 2 * np.pi, ndeg)

    for phi_idx in range(ndeg):
        phi = phi_1D[phi_idx]
        for theta_idx in range(ndeg):
            theta = theta_1D[theta_idx]
            for chi_idx in range(ndeg):
                chi = chi_1D[chi_idx]
                a = np.array(
                    [
                        np.sin(theta) * np.cos(phi),
                        np.sin(theta) * np.sin(phi),
                        np.cos(theta),
                    ]
                )
                b = np.array(
                    [
                        np.cos(theta) * np.cos(phi) * np.cos(chi)
                        - np.sin(phi) * np.sin(chi),
                        np.cos(theta) * np.sin(phi) * np.cos(chi)
                        + np.cos(phi) * np.sin(chi),
                        -np.sin(theta) * np.cos(chi),
                    ]
                )
                shear_tmp = 0
                for i in range(3):
                    for j in range(3):
                        v = C2V[i, j]
                        for k in range(3):
                            for l in range(3):
                                u = C2V[k, l]
                                rf = 1
                                if v >= 3 and u >= 3:
                                    rf = 4
                                if v >= 3 and u < 3:
                                    rf = 2
                                if u >= 3 and v < 3:
                                    rf = 2
                                rtmp = a[i] * b[j] * a[k] * \
                                    b[l] * (S[v, u] / rf)
                                shear_tmp = shear_tmp + rtmp
                shear_chi[chi_idx] = 1 / (4 * shear_tmp)
            shear_min[phi_idx, theta_idx] = np.amin(shear_chi)
            shear_max[phi_idx, theta_idx] = np.amax(shear_chi)
            shear_avg[phi_idx, theta_idx] = np.mean(shear_chi)

    if shear_choice == "avg":
        return shear_avg
    if shear_choice == "min":
        return shear_min
    if shear_choice == "max":
        return shear_max


# ------------------------------------POISSON RATIO----------------------------#

def plot_cry_poisson(theta_1D, phi_1D, S, ndeg, poisson_choice):
    """
    For each direction of the space (i.e., for each pair
    of theta and phi angles) the Poisson ratio is computed for the third angle
    chi and the average, maximum and minimum values are stored.

    Args:
        theta_1D (numpy.ndarray): One-dimensional array of theta values.
        phi_1D (numpy.ndarray): One-dimensional array of phi values.
        S (numpy.ndarray): Compliance matrix.
        ndeg (int): Number of degrees for discretization.
        poisson_choice (str): Type of Poisson's ratio to plot. Options: "avg", "min", "max".

    Returns:
        numpy.ndarray: Poisson's ratio array.

    Notes:
        - This function is intended to be called by cry_ela_plot
    """
    import numpy as np

    C2V = np.array(
        [
            [0, 5, 4],
            [5, 1, 3],
            [4, 3, 2],
        ]
    )
    poisson_chi = np.zeros(ndeg)
    poisson_min = np.zeros((ndeg, ndeg))
    poisson_max = np.zeros((ndeg, ndeg))
    poisson_avg = np.zeros((ndeg, ndeg))
    chi_1D = np.linspace(0, 2 * np.pi, ndeg)

    for phi_idx in range(ndeg):
        phi = phi_1D[phi_idx]
        for theta_idx in range(ndeg):
            theta = theta_1D[theta_idx]
            for chi_idx in range(ndeg):
                chi = chi_1D[chi_idx]
                a = np.array(
                    [
                        np.sin(theta) * np.cos(phi),
                        np.sin(theta) * np.sin(phi),
                        np.cos(theta),
                    ]
                )
                b = np.array(
                    [
                        np.cos(theta) * np.cos(phi) * np.cos(chi)
                        - np.sin(phi) * np.sin(chi),
                        np.cos(theta) * np.sin(phi) * np.cos(chi)
                        + np.cos(phi) * np.sin(chi),
                        -np.sin(theta) * np.cos(chi),
                    ]
                )
                poisson_num = 0
                poisson_den = 0
                for i in range(3):
                    for j in range(3):
                        v = C2V[i, j]
                        for k in range(3):
                            for l in range(3):
                                u = C2V[k, l]
                                rf = 1
                                if v >= 3 and u >= 3:
                                    rf = 4
                                if v >= 3 and u < 3:
                                    rf = 2
                                if u >= 3 and v < 3:
                                    rf = 2
                                num = (a[i] * a[j] * b[k] * b[l] * S[v, u])/rf
                                den = (a[i] * a[j] * a[k] * a[l] * S[v, u])/rf
                                poisson_num = poisson_num + num
                                poisson_den = poisson_den + den
                poisson_chi[chi_idx] = - poisson_num / poisson_den
            poisson_min[phi_idx, theta_idx] = np.amin(poisson_chi)
            poisson_max[phi_idx, theta_idx] = np.amax(poisson_chi)
            poisson_avg[phi_idx, theta_idx] = np.mean(poisson_chi)

    if poisson_choice == "avg":
        return poisson_avg
    if poisson_choice == "min":
        return poisson_min
    if poisson_choice == "max":
        return poisson_max


# ----------------------------------ELASTIC------------------------------------#

def plot_cry_ela(choose, ndeg,  *args, twoD=False):
    """
    Plot crystal elastic properties on the basis of the elastic tensor. A
    variable number of elastic tensors can be provided in order to get
    multiple plots in one shot, establishing a fixed color scale among them.

    Args:
        choose (str): Property to plot. Options: "young", "comp", "shear avg", "shear min", "shear max", "poisson avg", "poisson min", "poisson max".
        ndeg (int): Number of degrees for discretization.
        *args: Variable number of elastic tensors.

    Returns:
        Tuple of lists:
        - fig_list : list of matplotlib.figure.Figure
            A list containing matplotlib Figure objects for each plot.
        - ax_list : list of matplotlib.axes._axes.Axes
            A list containing the Axes objects associated with each plot.
        - plt_list : list of matplotlib.pyplot
            A list of the pyplot objects for each plot, representing the actual plot.
    """
    import math

    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import animation, cm, colors
    from mpl_toolkits.mplot3d import Axes3D, axes3d

    i = 0
    R = [None] * len(args)
    tmin = []
    tmax = []

    # Compute elastic properties for each tensor -->
    for C in args:

        # Inverse of the matrix C in GPa (Compliance)
        S = np.linalg.inv(C)

        # One dimentional array of theta from 0 to pi
        theta_1D = np.linspace(0, np.pi, ndeg)
        # One dimentional array of phi from 0 to 2pi
        phi_1D = np.linspace(0, 2 * np.pi, ndeg)
        # Make a 2D array for theta and phi
        theta_2D, phi_2D = np.meshgrid(theta_1D, phi_1D)

        # Call to function
        if choose == "young":
            R[i] = plot_cry_young(theta_2D, phi_2D, S)
        elif choose == "comp":
            R[i] = plot_cry_comp(theta_2D, phi_2D, S)
        elif choose == "shear avg":
            R[i] = plot_cry_shear(theta_1D, phi_1D, S, ndeg, "avg")
        elif choose == "shear min":
            R[i] = plot_cry_shear(theta_1D, phi_1D, S, ndeg, "min")
        elif choose == "shear max":
            R[i] = plot_cry_shear(theta_1D, phi_1D, S, ndeg, "max")
        elif choose == "poisson avg":
            R[i] = plot_cry_poisson(theta_1D, phi_1D, S, ndeg, "avg")
        elif choose == "poisson min":
            R[i] = plot_cry_poisson(theta_1D, phi_1D, S, ndeg, "min")
        elif choose == "poisson max":
            R[i] = plot_cry_poisson(theta_1D, phi_1D, S, ndeg, "max")

        i += 1
    # <--

    # Find highest and lowest values -->
    for k in range(i):
        tmin.append(np.min(R[k]))
        tmax.append(np.max(R[k]))
    vmin = min(tmin)
    vmax = max(tmax)
    # <--

    # Create plot for each tensor -->
    plt_list = []
    ax_list = []
    fig_list = []
    for k in range(i):
        X = R[k] * np.sin(theta_2D) * np.cos(phi_2D)
        Y = R[k] * np.sin(theta_2D) * np.sin(phi_2D)
        Z = R[k] * np.cos(theta_2D)

        norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=False)

        # if twoD == True:
        #     fig, ax = plt.subplots()
        #     ax.plot(X, Y)
        #     ax.set_xlabel("X")
        #     ax.set_ylabel("Y")
        #     fig_list.append(fig)
        #     ax_list.append(ax)
        #     plt_list.append(plt)
        #     fig, ax = plt.subplots()
        #     ax.plot(X, Z)
        #     ax.set_xlabel("X")
        #     ax.set_ylabel("Z")
        #     fig_list.append(fig)
        #     ax_list.append(ax)
        #     plt_list.append(plt)
        #     fig, ax = plt.subplots()
        #     ax.plot(Y, Z)
        #     ax.set_xlabel("Y")
        #     ax.set_ylabel("Z")
        #     fig_list.append(fig)
        #     ax_list.append(ax)
        #     plt_list.append(plt)

        fig, ax = plt.subplots(subplot_kw=dict(projection="3d"))

        ax.plot_surface(
            X,
            Y,
            Z,
            rstride=1,
            cstride=1,
            facecolors=cm.jet(norm(R[k])),
            antialiased=True,
            alpha=0.75,
        )

        m = cm.ScalarMappable(cmap=cm.jet, norm=norm)
        m.set_array(R[k])
        fig.colorbar(m, shrink=0.7, location="left", ax=ax)

        # Make the planes transparent
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        # Make the grid lines transparent
        #  ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        #  ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        #  ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        # Fixing limits
        ax.set_xlim(-1 * np.max(R), np.max(R))
        ax.set_ylim(-1 * np.max(R), np.max(R))
        ax.set_zlim3d(-1 * np.max(R), np.max(R))
        ax.locator_params(nbins=5)  # tight=True,
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")

        ax.set_box_aspect(aspect=(1, 1, 1))  # Fix aspect ratio

        fig_list.append(fig)
        ax_list.append(ax)
        plt_list.append(plt)

    return fig_list, ax_list, plt_list

    # <--


##############################################################################
#                                                                            #
#                             VIBRATIONAL PROPERTIES                         #
#                                                                            #
##############################################################################

# ------------------------------------HARMONIC---------------------------------#

def plot_cry_irspec(irspec, x_unit='cm-1', y_mode='LG', figsize=None, linestyle='-',
                    linewidth=1.5, color='tab:blue', freq_range=None, int_range=None,
                    label=None, dpi=100, offset=0):
    """Generates the IR spectra for the IRSPEC.DAT file produced by an IRSPEC calculation

    Args:
        irspec (External_unit object): Object (or a list of) generated by the read_cry_irspec function necessary for the plot
        x_unit (str, optional): Unit measure of the x axes. Avalilable: 'cm-1' and 'nm'. Defaults to 'cm-1'.
        y_mode (str, optional): Peak broadening modality in absorbance and reflectance. 
                                Available: 'LG'(Lorentzian-Gaussian broadening), 'V' (Voight broadening), 'RS' (Rayleigh spherical particles), 'RE' (Rayleigh with elipsoid particles), 'REFL' (Reflectance)
                                Defaults to 'LG'.
        figsize (tuple, optional): Image dimensions correspondig to matplotlib figsize. Defaults to None.
        linestyle (str/list[str], optional): linestyle corresponding to the matplotlib one it can be a list for a multiplot. Defaults to '-'.
        linewidth (float/list[float], optional): linewidth corresponding to the matplotlib one it can be a list for a multiplot. Defaults to 1.5.
        color (str/list[str], optional): Color of the spectra it can accept all matplotlib colors it can be a list for multiplots. Defaults to 'tab:blue'.
        freq_range (list, optional): Two element list [min, max], that allows to visualize the spectra in a given frequency window. Defaults to None.
        int_range (list, optional): Two element list [min, max], that allows to visualize the spectra in a given intensity window. Defaults to None.
        label (list[str], optional): List of labels for the legend of a multiplot. Defaults to None.
        dpi (int, optional): Resolution of the saved file. Defaults to 100
        offset(float, optional): Allows the user to define an offset between different spectra in a multi plot


    Returns:
        Matplotlib object
    Raises:
        ValueError: The function raises an error when the object to be plotted does not have the required y_mode  
    """

    import sys
    import warnings

    import matplotlib.pyplot as plt
    import numpy as np

    modes = ['single', 'multi']
    accepted_y = ['LG', 'V', 'RS', 'RE', 'REFL']

    if isinstance(irspec, list):
        mode = modes[1]

        if not isinstance(linestyle, list):
            style = linestyle
            linestyle = []
            for i in enumerate(irspec):
                linestyle.append(style)

        if not isinstance(linewidth, list):
            width = linewidth
            linewidth = []
            for i in enumerate(irspec):
                linewidth.append(width)

        if not isinstance(color, list):
            color = ['dimgrey', 'blue', 'indigo', 'slateblue',
                     'thistle', 'purple', 'orchid', 'crimson']

        for file in irspec:
            if (file.calculation == 'molecule') and (y_mode != accepted_y[0]):
                raise ValueError(
                    'This spectra does not contain the y_mode requested: available y_mode'+accepted_y[0])

    else:
        mode = modes[0]

        if (irspec.calculation == 'molecule') and (y_mode != accepted_y[0]):
            raise ValueError(
                'This spectra does not contain the y_mode requested: available y_mode'+accepted_y[0])

    if figsize is not None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    else:
        fig, ax = plt.subplots(dpi=dpi)

    if mode == modes[0]:

        # selection of the x axis unit
        if x_unit == 'cm-1':
            x = irspec.irspec[:, 0]

        elif x_unit == 'nm':
            x = irspec.irspec[:, 1]

        # selection of the intensities mode
        if y_mode == accepted_y[0]:
            y = irspec.irspec[:, 2]

        elif y_mode == accepted_y[1]:
            y = irspec.irspec[:, 5]

        elif y_mode == accepted_y[2]:
            y = irspec.irspec[:, 6]

        elif y_mode == accepted_y[3]:
            y = irspec.irspec[:, 7]

        elif y_mode == accepted_y[4]:
            y = irspec.irspec[:, 8]

        xmin = min(x)
        xmax = max(x)
        ymin = min(y)-1
        ymax = max(y)+10

        ax.plot(x, y, linestyle=linestyle, linewidth=linewidth, color=color)

    if mode == modes[1]:

        xmin = []
        xmax = []
        ymin = []
        ymax = []
        add_offset = 0

        for index, file in enumerate(irspec):
            # selection of the x axis unit
            if x_unit == 'cm-1':
                x = file.irspec[:, 0]

            elif x_unit == 'nm':
                x = file.irspec[:, 1]

            # selection of the intensities mode
            if y_mode == accepted_y[0]:
                y = file.irspec[:, 2] + add_offset

            elif y_mode == accepted_y[1]:
                y = file.irspec[:, 5] + add_offset

            elif y_mode == accepted_y[2]:
                y = file.irspec[:, 6] + add_offset

            elif y_mode == accepted_y[3]:
                y = file.irspec[:, 7] + add_offset

            elif y_mode == accepted_y[4]:
                y = file.irspec[:, 8] + add_offset

            xmin.append(min(x))
            xmax.append(max(x))
            ymin.append(min(y)-1)
            ymax.append(max(y)+10)

            if label is not None:
                ax.plot(x, y, linestyle=linestyle[index], linewidth=linewidth[index],
                        color=color[index], label=label[index])
                plt.legend()
            else:
                ax.plot(
                    x, y, linestyle=linestyle[index], linewidth=linewidth[index], color=color[index])

            add_offset += offset

        xmin = min(xmin)
        xmax = max(xmax)
        ymin = min(ymin)
        ymax = max(ymax)

    if freq_range is not None:
        xmin = freq_range[0]
        xmax = freq_range[1]

    if int_range is not None:
        ymin = int_range[0]
        ymax = int_range[1]

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    if x_unit == 'cm-1':
        plt.xlabel('Wavenumber (cm$^{-1}$)')
    elif x_unit == 'nm':
        plt.xlabel('Wavelength (nm)')

    if y_mode != accepted_y[4]:
        plt.ylabel('Absorbance (A.U.)')
    else:
        plt.ylabel('Reflectance (A.U.)')

    return fig, ax


def plot_cry_ramspec(ramspec,  y_mode='total', figsize=None, linestyle='-',
                     linewidth=1.5, color='tab:blue', freq_range=None, int_range=None,
                     label=None, dpi=100, offset=0):
    """Generates the RAMAN spectra for the RAMSPEC.DAT file produced by an RAMSPEC calculation

    Args:
        ramspec (External_unit object): Object (or a list of) generated by the read_cry_ramspec function necessary for the plot
        y_mode (str, optional): Polarization of the spectra for the simulated compound
                                Available: 'total', 'parallel', 'perpendicular' (for powders), 'xx', 'xy', 'xz', 'yy', 'yz', 'zz' (for single crystals)
                                Defaults to 'LG'.
        figsize (tuple, optional): Image dimensions correspondig to matplotlib figsize. Defaults to None.
        linestyle (str/list[str], optional): linestyle corresponding to the matplotlib one it can be a list for a multiplot. Defaults to '-'.
        linewidth (float/list[float], optional): linewidth corresponding to the matplotlib one it can be a list for a multiplot. Defaults to 1.5.
        color (str/list[str], optional): Color of the spectra it can accept all matplotlib colors it can be a list for multiplots. Defaults to 'tab:blue'.
        freq_range (list, optional): Two element list [min, max], that allows to visualize the spectra in a given frequency window. Defaults to None.
        int_range (list, optional): Two element list [min, max], that allows to visualize the spectra in a given intensity window. Defaults to None.
        label (list[str], optional): List of labels for the legend of a multiplot. Defaults to None.
        dpi (int, optional): Resolution of the saved file. Defaults to 300.
        offset (float): Allows the user to define an offset between different spectra in a multi plot
    Returns:
        Matplotlib object
    """

    import sys
    import warnings

    import matplotlib.pyplot as plt
    import numpy as np

    modes = ['single', 'multi']
    accepted_y = ['total', 'parallel', 'perpendicular',
                  'xx', 'xy', 'xz', 'yy', 'yz', 'zz']

    if isinstance(ramspec, list):
        mode = modes[1]
        if not isinstance(linestyle, list):
            style = linestyle
            linestyle = []
            for i in enumerate(ramspec):
                linestyle.append(style)

        if not isinstance(linewidth, list):
            width = linewidth
            linewidth = []
            for i in enumerate(ramspec):
                linewidth.append(width)

        if not isinstance(color, list):
            color = ['dimgrey', 'blue', 'indigo', 'slateblue',
                     'thistle', 'purple', 'orchid', 'crimson']

    else:
        mode = modes[0]

    if figsize is not None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    else:
        fig, ax = plt.subplots()

    if mode == modes[0]:

        x = ramspec.ramspec[:, 0]

        # selection of the intensities mode
        if y_mode == accepted_y[0]:
            y = ramspec.ramspec[:, 1]

        elif y_mode == accepted_y[1]:
            y = ramspec.ramspec[:, 2]

        elif y_mode == accepted_y[2]:
            y = ramspec.ramspec[:, 3]

        elif y_mode == accepted_y[3]:
            y = ramspec.ramspec[:, 4]

        elif y_mode == accepted_y[4]:
            y = ramspec.ramspec[:, 5]

        elif y_mode == accepted_y[5]:
            y = ramspec.ramspec[:, 6]

        elif y_mode == accepted_y[6]:
            y = ramspec.ramspec[:, 7]

        elif y_mode == accepted_y[7]:
            y = ramspec.ramspec[:, 8]

        elif y_mode == accepted_y[8]:
            y = ramspec.ramspec[:, 9]

        xmin = min(x)
        xmax = max(x)
        ymin = min(y)-1
        ymax = max(y)+10

        ax.plot(x, y, linestyle=linestyle,
                linewidth=linewidth, color=color)

    if mode == modes[1]:
        xmin = []
        xmax = []
        ymin = []
        ymax = []
        add_offset = 0

        for index, file in enumerate(ramspec):
            x = file.ramspec[:, 0]

            # selection of the intensities mode
            if y_mode == accepted_y[0]:
                y = file.ramspec[:, 1] + add_offset

            elif y_mode == accepted_y[1]:
                y = file.ramspec[:, 2] + add_offset

            elif y_mode == accepted_y[2]:
                y = file.ramspec[:, 3] + add_offset

            elif y_mode == accepted_y[3]:
                y = file.ramspec[:, 4] + add_offset

            elif y_mode == accepted_y[4]:
                y = file.ramspec[:, 5] + add_offset

            elif y_mode == accepted_y[5]:
                y = file.ramspec[:, 6] + add_offset

            elif y_mode == accepted_y[6]:
                y = file.ramspec[:, 7] + add_offset

            elif y_mode == accepted_y[7]:
                y = file.ramspec[:, 8] + add_offset

            elif y_mode == accepted_y[8]:
                y = file.ramspec[:, 9] + add_offset

            xmin.append(min(x))
            xmax.append(max(x))
            ymin.append(min(y)-1)
            ymax.append(max(y)+10)

            if label is not None:
                ax.plot(x, y, linestyle=linestyle[index], linewidth=linewidth[index],
                        color=color[index], label=label[index])
                plt.legend()
            else:
                ax.plot(
                    x, y, linestyle=linestyle[index], linewidth=linewidth[index], color=color[index])

            add_offset += offset

        xmin = min(xmin)
        xmax = max(xmax)
        ymin = min(ymin)
        ymax = max(ymax)

    if freq_range is not None:
        xmin = freq_range[0]
        xmax = freq_range[1]

    if int_range is not None:
        ymin = int_range[0]
        ymax = int_range[1]

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    plt.xlabel('Wavenumber (cm$^{-1}$)')

    if y_mode != accepted_y[4]:
        plt.ylabel('Absorbance (A.U.)')
    else:
        plt.ylabel('Reflectance (A.U.)')

    return fig, ax


# -----------------------------------ANHARMONIC--------------------------------#

def plot_cry_spec(transitions, typeS="lorentz", components=False, bwidth=5,
                  stdev=3, eta=0.5, fmin=None, fmax=None, ylim=None,
                  exp_spec=None, sep=";", export_csv=False, label=None,
                  xlabel='Wavenumber [cm$^{-1}$]', ylabel='Intensity [arb. u.]',
                  linewidth=2.0, padd=100, style=None, compstyle=None,
                  figsize=(16, 6), fig=None, ax=None, offset=0):
    """
    This function enables the simulation of vibrational spectra based on a 2D 
    NumPy array containing a list of transition frequencies and the 
    corresponding intensities. The code allows users to model spectral
    broadening according to various profiles (Gaussian, Lorentzian, 
    pseudo-Voigt), or zero broadening (Dirac deltas-like lines). Please, note
    that by turning the optional argument 'component' to `True` you can
    additionally plot contributions arising from each transition.

    Args:
        transitions (numpy.ndarray): 2D array containing transition frequencies (axis=0) and corresponding intensities (axis=1).
        typeS (str): String specifying the spectral profile: "bars", "lorentz", "gauss", "pvoigt" (default is "lorentz"). 
        components (bool, optional): Whether to plot contributions arising from each transition (default is `False`).  
        bwidth (float, optional): Half-width at half-maximum of the Lorentzian profile (default is 5).
        stdev (float, optional): Standard deviation of the Gaussian profile (default is 5).
        eta (float, optional): Fraction of Lorentzian character in pseudo-Voigt profile (default is 0.5).
        fmin (float, optional): Minimum frequency.
        fmax(float, optional): Maximum frequency.
        ylim (float, optional): Maximum intensity.
        export_csv (bool, optional): Whether to save plot in csv format (default is `False`).
        xlabel (str, optional): x-axis label (default is "Wavenumber [cm$^{-1}$]").
        ylabel (str, optional): y-axis label (default is "Intensity [arb. u.]").
        linewidth (float, optional): Linewidth (default is 2.0).
        padd (float, optional): left- and right- hand side padding expressed in the same unit of the quantity reported in x-axis (default is 100).
        style (str, optional): String specifying Matplotlib style. 
        compstyle (list[string], optional): List containing Matplotlib styles to plot each component. 
        figsize (list[real], optional): List of two numbers specifying the aspect ratio of the figure (default is [16, 6]).

    Returns:
        matplotlib.figure.Figure
        matplotlib.axes.Axes
    """

    import math
    import time
    from copy import deepcopy

    import matplotlib.pyplot as plt
    import numpy as np
    from numpy import genfromtxt

    if (ax is None):
        fig, ax = plt.subplots(figsize=figsize)

    if (ylim is not None):
        ax.set_ylim(0, ylim)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    bars = False
    lorentz = False
    gauss = False
    pseudo_voigt = False

    if typeS == 'bars':
        bars = True

    if typeS == 'lorentz':
        lorentz = True

    if typeS == 'gauss':
        gauss = True

    if typeS == 'pvoigt':
        pseudo_voigt = True

    n = 20000

    if fmin is None:
        fmin = min(transitions[:, 0] - padd)
    if fmax is None:
        fmax = max(transitions[:, 0] + padd)

    x = np.linspace(fmin, fmax, num=n)
    y = np.zeros(n)

    spec_data = np.block([[x], [y]]).T
    sbuff = np.block([[x], [y]]).T

    if bars:
        spec_data = np.concatenate((spec_data, transitions), axis=0)
        spec_data = spec_data[spec_data[:, 0].argsort()]
    elif lorentz:
        iL = 0
        L = []
        for i in range(len(transitions)):
            if transitions[i, 1] == 0:
                continue
            for j, f in enumerate(spec_data[:, 0]):
                lorentz = (1/math.pi)*bwidth / \
                    ((f-transitions[i, 0])**2+bwidth**2)*transitions[i, 1]
                sbuff[j, 1] = lorentz
            L.append(deepcopy(sbuff))
            iL = iL + 1
        if (not components):
            for i in range(len(L)):
                spec_data[:, 1] = spec_data[:, 1] + L[i][:, 1]
        else:
            for i in range(len(L)):
                ax.plot(spec_data[:, 0], L[i][:, 1], linewidth=linewidth)
            for i in range(len(L)):
                spec_data[:, 1] = spec_data[:, 1] + L[i][:, 1]

    elif gauss:
        G = []
        for i in range(len(transitions)):
            if transitions[i, 1] == 0:
                continue
            for j, f in enumerate(spec_data[:, 0]):
                gauss = (1/(stdev*math.sqrt(2*math.pi))) * \
                    math.exp(-((f-transitions[i, 0])**2) /
                             (2*stdev**2))*transitions[i, 1]
                sbuff[j, 1] = gauss
            G.append(deepcopy(sbuff))
        if (not components):
            for i in range(len(G)):
                spec_data[:, 1] = spec_data[:, 1] + G[i][:, 1]
        else:
            for i in range(len(G)):
                ax.plot(spec_data[:, 0], G[i][:, 1], linewidth=linewidth)
            for i in range(len(G)):
                spec_data[:, 1] = spec_data[:, 1] + G[i][:, 1]

    elif pseudo_voigt:
        V = []
        for i in range(len(transitions)):
            if transitions[i, 1] == 0:
                continue
            for j, f in enumerate(spec_data[:, 0]):
                gauss = (1/(stdev*math.sqrt(2*math.pi))) * \
                    math.exp(-((f-transitions[i, 0])**2) /
                             (2*stdev**2))*transitions[i, 1]
                lorentz = (1/math.pi)*bwidth / \
                    ((f-transitions[i, 0])**2+bwidth**2)*transitions[i, 1]
                sbuff[j, 1] = eta*lorentz + (1-eta)*gauss
            V.append(deepcopy(sbuff))
        if (not components):
            for i in range(len(V)):
                spec_data[:, 1] = spec_data[:, 1] + V[i][:, 1]
        else:
            for i in range(len(V)):
                if (compstyle is not None):
                    ax.plot(spec_data[:, 0], V[i][:, 1], compstyle[i],
                            linewidth=linewidth)
                else:
                    ax.plot(spec_data[:, 0], V[i][:, 1], linewidth=linewidth)
            for i in range(len(V)):
                spec_data[:, 1] = spec_data[:, 1] + V[i][:, 1]

    if (exp_spec is not None):
        exp_data = genfromtxt(exp_spec, delimiter=sep)
        area_spec_data = np.trapz(spec_data[:, 1], spec_data[:, 0])
        area_exp_data = np.trapz(exp_data[:, 1], exp_data[:, 0])
        norm_fac = area_spec_data / area_exp_data
        baseline = 0.2
        exp_data[:, 1] = exp_data[:, 1] * norm_fac - baseline  # * 0.5
        ax.plot(exp_data[:, 0], exp_data[:, 1], 'r-', linewidth=linewidth)

    if ((label is not None) and (style is None)):
        ax.plot(spec_data[:, 0], spec_data[:, 1] + offset, linewidth=linewidth,
                label=label)
    elif ((label is None) and (style is not None)):
        ax.plot(spec_data[:, 0], spec_data[:, 1] +
                offset, style, linewidth=linewidth)
    elif ((label is not None) and (style is not None)):
        ax.plot(spec_data[:, 0], spec_data[:, 1] + offset, style, linewidth=linewidth,
                label=label)
    else:
        ax.plot(spec_data[:, 0], spec_data[:, 1] + offset, linewidth=linewidth)

    if (export_csv):
        np.savetxt(typeS + time.strftime("%Y-%m-%d_%H%M%S.") + 'csv',
                   spec_data, delimiter=';')
    return fig, ax


def plot_cry_spec_multi(files, typeS="lorentz", components=False, bwidth=5,
                        stdev=3, eta=0.5, fmin=None, fmax=None, ylim=None,
                        label=None, xlabel='Wavenumber [cm$^{-1}$]',
                        ylabel='Intensity [arb. u.]', linewidth=2.0, padd=100,
                        style=None, figsize=(16, 6), exp_spec=None, norm_fac=1,
                        sep=';', offset=0):
    """
    This function is a wrapper for `plot_spec` function, enablng the simulation 
    of many vibrational spectra coming from a list of NumPy array.  

    Args:
        transitions (list[numpy.ndarray]): List of 2D arrays containing transition frequencies (axis=0) and corresponding intensities (axis=1).
        typeS (str): String specifying the spectral profile: "bars", "lorentz", "gauss", "pvoigt". 
        components (bool, optional): Whether to plot contributions arising from each transition (default is `False`).  
        bwidth (float, optional): Half-width at half-maximum of the Lorentzian profile (default is 5).
        stdev (float, optional): Standard deviation of the Gaussian profile (default is 5).
        eta (float, optional): Fraction of Lorentzian character in pseudo-Voigt profile (default is 0.5).
        fmin (float, optional): Minimum frequency.
        fmax(float, optional): Maximum frequency
        ylim (float, optional): Maximum intensity.
        xlabel (str, optional): x-axis label (default is "Wavenumber [cm$^{-1}$]").
        ylabel (str, optional): y-axis label (default is "Intensity [arb. u.]").
        linewidth (float, optional): Linewidth (default is 2.0).
        padd (float, optional): left- and right- hand side padding expressed in the same unit of the quantity reported in x-axis (default is 100).
        style (str, optional): String specifying Matplotlib style. 
        figsize (list[float], optional): List of two numbers specifying the aspect ratio of the figure (default is [16, 6]).
        offset (float, optional) : Offset along the y axis (default is 0). 

    Returns:
        matplotlib.figure.Figure
        matplotlib.axes.Axes
    """

    import matplotlib.pyplot as plt
    from numpy import genfromtxt

    fig, ax = plt.subplots(figsize=figsize)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if (exp_spec is not None):
        exp_data = genfromtxt(exp_spec, delimiter=sep)
        # area_spec_data = np.trapz(spec_data[:, 1], spec_data[:, 0])
        # area_exp_data = np.trapz(exp_data[:, 1], exp_data[:, 0])
        # norm_fac = area_spec_data / area_exp_data
        baseline = 0.2
        exp_data[:, 1] = exp_data[:, 1] * norm_fac - baseline
        ax.plot(exp_data[:, 0], exp_data[:, 1], 'r-', linewidth=linewidth)

    compstyle = []
    if (style is not None):
        for i in range(len(style)):
            compstyle.append([style[i]] * 100)

    for i, transitions in enumerate(files):
        if ((label is not None) and (style is None)):
            plot_cry_spec(transitions, typeS, components, bwidth, stdev, eta,
                          fmin, fmax, ylim, label=label[i], linewidth=linewidth,
                          padd=padd, xlabel=xlabel, ylabel=ylabel, fig=fig,
                          ax=ax, offset=offset*i)
        elif ((style is not None) and (label is None)):
            plot_cry_spec(transitions, typeS, components, bwidth, stdev, eta,
                          fmin, fmax, ylim, linewidth=linewidth, padd=padd,
                          style=style[i], xlabel=xlabel, ylabel=ylabel,
                          compstyle=compstyle[i], fig=fig, ax=ax,
                          offset=offset*i)
        elif ((style is not None) and (label is not None)):
            plot_cry_spec(transitions, typeS, components, bwidth, stdev, eta,
                          fmin, fmax, ylim, linewidth=linewidth, padd=padd,
                          style=style[i], label=label[i], xlabel=xlabel,
                          ylabel=ylabel, compstyle=compstyle[i], fig=fig,
                          ax=ax, offset=offset*i)
        else:
            plot_cry_spec(transitions, typeS, components, bwidth, stdev, eta, fmin,
                          fmax, ylim, linewidth=linewidth, padd=padd,
                          xlabel=xlabel, ylabel=ylabel, fig=fig, ax=ax,
                          offset=offset*i)

    if (label is not None):
        ax.legend(loc='upper left')

    return fig, ax


def plot_cry_anscan(co, scale_wf=None, scale_prob=None, harmpot=False,
                    scanpot=True, figsize=[10, 10]):
    """
    This function provides a plotting tool for the ANSCAN keyword.  

    Args:
        co (crystal_io.Crystal_output): Crystal output object.
        scale_wf (float, optional): Scaling factor for wavefunctions plot. By 
        default, wavefunctions are not plotted. Providing a value for this 
        argument enables the wavefunction plot and scales it accordingly.  
        scale_prob (float, optional): Scaling factor for probability density 
        plot. By default, probability densities are not plotted. Providing a 
        value for this argument enables the probability density plot and 
        scales it accordingly.  
        harmpot (bool, optional): A logical flag to activate the plotting of 
        the harmonic potential (default is False). 
        scanpot (bool, optional): A logical flag to activate the plotting of 
        the scan potential provided by ANSCAN (default is True). 
        figsize (list[float])

    Notes:
        - This is a work in progress.
    """

    import math

    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import special

    # Unpack co
    harm_freq = co.harm_freq
    force_const = co.force_const
    energy = co.energy
    wf = co.wf
    scan_energy = co.anhpot
    alpha = co.alpha
    rangescan = co.rangescan

    # debug
    print(co.alpha)
    # debug

    # npts
    npts = 10000

    # Matplotlib aspect ratio
    plt.figure(figsize=figsize)

    # Define harmonic freq
    amu_me = 1822.88848
    Ha2wn = 219473.5152
    lambda_AU = abs(harm_freq / Ha2wn) * amu_me

    # debug
    print(harm_freq)
    print(lambda_AU)
    # debug

    # Define coordinates (basis set)
    x = np.linspace(-1000, 1000, npts)
    xi = x * alpha

    # Plot levels
    Nlevel = 10
    for i in range(Nlevel):
        plt.hlines(y=energy[i], xmin=rangescan[0], xmax=rangescan[1],
                   colors='k', linewidth=0.3)

    # Set number of basis functions and wf
    N = len(wf)
    Nwf = N

    # Compute Gaussian functions
    G = np.exp(-(xi)**2/2)

    # Compute harmonic wf
    wfHO = np.zeros([len(xi), N])
    for m in range(N):
        norm = np.sqrt((alpha) / ((np.sqrt(math.pi))
                       * (2**m) * math.factorial(m)))
        Herm = special.hermite(m, monic=False)
        wfHO[:, m] = norm * Herm(xi) * G

    # Build anharmonic wf
    wfANH = np.zeros([len(wfHO), N])

    # Define coordinates (wf)
    x = np.linspace(-100, 100, npts)
    xi = x * lambda_AU**0.25

    # for s in range(N):
    for s in range(Nlevel):
        for i in range(N):
            wfANH[:, s] = wfANH[:, s] + wf[i, s]*wfHO[:, i]

    # Wavefunctions -->
    if (scale_wf is not None):

        # Plot wf
        for i in range(Nwf):
            yp = wfANH[:, i]*scale_wf + energy[i]
            plt.plot(xi, yp, "m-", linewidth=1)
    # <-- Wavefunctions

    # Probability density  -->
    if (scale_prob is not None):
        for s in range(Nlevel):
            prob = wfANH[:, s]**2*scale_prob**2 + energy[s]
            lower_bound = energy[s] + 0*xi
            plt.fill_between(xi, lower_bound, prob, color='c', alpha=0.3)

    # Plot anharmonic potential
    anhpot = 1/2 * force_const[2] * xi**2 \
        + 1/6 * force_const[3] * xi**3 \
        + 1/24 * force_const[4] * xi**4
    plt.plot(xi, anhpot, 'b', linewidth=2)

    # Plot harmonic potential
    if (harmpot):
        HOpot = 1/2 * harm_freq * xi**2
        plt.plot(xi, HOpot, 'r--', linewidth=2)

    # Plot scan potential
    if (scanpot):
        # Define scaled_x
        scaled_x = np.linspace(rangescan[0], rangescan[1], len(scan_energy))
        plt.plot(scaled_x, scan_energy, 'bo')

    plt.ylabel('$\Delta E$ [cm$^{-1}$]')
    plt.xlabel(r'$\xi$')
    plt.xlim([rangescan[0], rangescan[1]])
    plt.ylim([harm_freq, abs(harm_freq)*5])

    return plt

# --------------------------------------EOS-----------------------------------#


def plot_cry_EOS(eos, formula_unit=None, plot='VvsE', color='tab:blue', figsize=(5, 5), dpi=72, marker='o', fontsize=12, legend=None):
    """
    This function provides a tool to plot quantities extracted from the EOS module

    Args:
        eos[obj/list]: CRYSTAL object or list of two CRYSTAL objects containing quantities extracted 
                       from the EOS module generated by the get_EOS() function
        formula_unit[list]: List containing the formula unit associated with each CRYSTAL object,
                            the order should be the same as in the eos list (default:None)
        plot[str]: String identifing the Equation of State/Computed data the user wants to plot (default:VvsE)
        color[str/list]: String or list of defining the colors of the plotted series of data requested (default:'tab:blue')
        figsize[tuple]: Tuple of floats defining the size of the figure (default:(5,5))
        dpi[int]: Integer defining the resolution of the figure (default:72)
        fontsize[int]: Integer defining the fontsize of the text (default:12)
        legend[list]: List of strings that will be used as legend in the plot (default:None)

    Returns:
        Matplotlib object
    """
    import matplotlib.pyplot as plt
    import numpy as np

    plot_type = ['VvsE', 'Murnaghan',
                 'Birch-Murnaghan', 'Poirier-Tarantola', 'Vinet']

    plt.rcParams.update({'font.size': fontsize})

    # Check the type of plot required by the user -->
    if type(eos) is list:
        phase_transition = True
        if type(color) is not list:
            color = ['tab:blue', 'tab:orange']
        if type(marker) is not list:
            marker = ['o', 'o']
        if legend == None:
            legend = ['Phase 1', 'Phase 2']
        if formula_unit is None:
            formula_unit = [1, 1]
    else:
        phase_transition = False
    # <--

    # Plotting -->
    if phase_transition:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        # Plot of the Volumes against the Energy -->
        if plot == plot_type[0]:
            for index, phase in enumerate(eos):
                ax.plot(phase.VvsE[:, 0]/formula_unit[index], phase.VvsE[:, 1]/formula_unit[index],
                        color=color[index], marker=marker[index], label=legend[index])

            plt.xlabel('V/Z ($\AA$)')
            plt.ylabel('E/Z (a.u)')
        # <--

        # Murnaghan EOS plot of delta of Gibbs free energy between the two phases-->
        if plot == plot_type[1]:
            maxp = []
            minp = []
            for _, phase in enumerate(eos):
                maxp.append(max(phase.murnaghan[:, 1]))
                minp.append(min(phase.murnaghan[:, 1]))

            maxp = max(maxp)
            minp = min(minp)

            referencex = np.linspace(minp-10, maxp+10, 2)
            referencey = np.zeros(2)

            ax.plot(referencex, referencey, color='black', label=legend[0])

            deltag = eos[1].murnaghan[:, 3]/formula_unit[1] - \
                eos[0].murnaghan[:, 3]/formula_unit[0]

            ax.plot(eos[1].murnaghan[:, 1], deltag,
                    color=color[1], label=legend[1])
        # <--

        # Birch-Murnaghan EOS plot of delta of Gibbs free energy between the two phases-->
        if plot == plot_type[2]:
            maxp = []
            minp = []
            for _, phase in enumerate(eos):
                maxp.append(max(phase.bmurnaghan[:, 1]))
                minp.append(min(phase.bmurnaghan[:, 1]))

            maxp = max(maxp)
            minp = min(minp)

            referencex = np.linspace(minp-10, maxp+10, 2)
            referencey = np.zeros(2)

            ax.plot(referencex, referencey, color='black', label=legend[0])

            deltag = eos[1].bmurnaghan[:, 3]/formula_unit[1] - \
                eos[0].bmurnaghan[:, 3]/formula_unit[0]

            ax.plot(eos[1].bmurnaghan[:, 1], deltag,
                    color=color[1], label=legend[1])
        # <--

        # Poirier-Tarantola EOS plot of delta of Gibbs free energy between the two phases-->
        if plot == plot_type[3]:
            maxp = []
            minp = []
            for _, phase in enumerate(eos):
                maxp.append(max(phase.pt[:, 1]))
                minp.append(min(phase.pt[:, 1]))

            maxp = max(maxp)
            minp = min(minp)

            referencex = np.linspace(minp-10, maxp+10, 2)
            referencey = np.zeros(2)

            ax.plot(referencex, referencey, color='black', label=legend[0])

            deltag = eos[1].pt[:, 3]/formula_unit[1] - \
                eos[0].pt[:, 3]/formula_unit[0]

            ax.plot(eos[1].pt[:, 1], deltag,
                    color=color[1], label=legend[1])
        # <--

        # Vinet EOS plot of delta of Gibbs free energy between the two phases-->
        if plot == plot_type[4]:
            maxp = []
            minp = []
            for _, phase in enumerate(eos):
                maxp.append(max(phase.vinet[:, 1]))
                minp.append(min(phase.vinet[:, 1]))

            maxp = max(maxp)
            minp = min(minp)

            referencex = np.linspace(minp-10, maxp+10, 2)
            referencey = np.zeros(2)

            ax.plot(referencex, referencey, color='black', label=legend[0])

            deltag = eos[1].vinet[:, 3]/formula_unit[1] - \
                eos[0].vinet[:, 3]/formula_unit[0]

            ax.plot(eos[1].vinet[:, 1], deltag,
                    color=color[1], label=legend[1])
        # <--
        if plot != plot_type[0]:
            plt.xlabel('P (GPa)')
            plt.ylabel('$\Delta$H/Z (a.u)')

        if legend != ['Phase 1', 'Phase 2']:
            plt.legend()

    else:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        # Plot of the Volume against the Energy -->
        if plot == plot_type[0]:
            ax.plot(eos.VvsE[:, 0], eos.VvsE[:, 1],
                    color=color, marker=marker)
            plt.xlabel('V/Z ($\AA$)')
            plt.ylabel('E/Z (a.u)')
        # <--

        # Murnaghan Pressure vs Gibbs Free Energy -->
        elif plot == plot_type[1]:
            ax.plot(eos.murnaghan[:, 1], eos.murnaghan[:, 3],
                    color=color, marker=marker)
        # <--

        # Birch-Murnaghan Pressure vs Gibbs Free Energy -->
        elif plot == plot_type[2]:
            ax.plot(eos.bmurnaghan[:, 1], eos.bmurnaghan[:, 3],
                    color=color, marker=marker)
        # <--

        # Poirier-Tarantola Pressure vs Gibbs Free Energy -->
        elif plot == plot_type[3]:
            ax.plot(eos.pt[:, 1], eos.pt[:, 3], color=color, marker=marker)
        # <--

        # Vinet Pressure vs Gibbs Free Energy -->
        elif plot == plot_type[4]:
            ax.plot(eos.vinet[:, 1], eos.vinet[:, 3],
                    color=color, marker=marker)
        # <--

        if ploi != plot_type[0]:
            plt.xlabel('P (GPa)')
            plt.ylabel('H (a. u.)')

    # <--

    return fig, ax
