from CRYSTALClear.crystal_io import Properties_output
from CRYSTALClear.units import au_to_angstrom, angstrom_to_au
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import math

class ScalarField(Properties_output):
    def __init__(self, fort25_path, units="Angstrom"):
        super().read_cry_ECHG(fort25_path)

        self.density_map = angstrom_to_au(angstrom_to_au(self.density_map))
        self.__get_minmaxvalue()
        self.__generate_mesh()
        self.units = "Angstrom"
        self.change_units(units)

        self.__scaled = False

    def __generate_mesh(self):
        ab = np.linalg.norm(self.a - self.b)
        cb = np.linalg.norm(self.c - self.b)
        self.__meshx = np.zeros((self.nrow, self.ncol), dtype=float)
        self.__meshy = np.zeros((self.nrow, self.ncol), dtype=float)
        for i in range(0, self.nrow):
            for j in range(0, self.ncol):
                self.__meshy[i, j] = ((ab/self.nrow)*i) * np.sqrt(1 - self.cosxy**2)
                self.__meshx[i, j] = (((cb/self.ncol)*j) * np.sqrt(1 - self.cosxy**2)) + \
                        (((ab/self.nrow)*i) * self.cosxy)
        self.__meshx_max = np.amax(self.__meshx)
        self.__meshy_max = np.amax(self.__meshy)
        self.__meshx_min = np.amin(self.__meshx)
        self.__meshy_min = np.amin(self.__meshy)

    def __get_minmaxvalue(self):
        self.dens_min_value = np.amin(self.density_map)
        self.dens_max_value = np.amax(self.density_map)
 
    def change_units(self, units: str):
        if (self.units == "Angstrom" and units == "Bohr"):
            self.units = "Bohr"
            self.a = angstrom_to_au(self.a)
            self.b = angstrom_to_au(self.b)
            self.c = angstrom_to_au(self.c)
            self.density_map = au_to_angstrom(au_to_angstrom(self.density_map))
            self.__get_minmaxvalue()
            self.__generate_mesh()
        elif (self.units == "Bohr" and units == "Angstrom"):
            self.units = "Angstrom"
            self.a = au_to_angstrom(self.a)
            self.b = au_to_angstrom(self.b)
            self.c = au_to_angstrom(self.c)
            self.density_map = angstrom_to_au(angstrom_to_au(self.density_map))
            self.__get_minmaxvalue()
            self.__generate_mesh()
        elif (self.units == units):
            pass
        else:
            raise ValueError("Only the keywords Bohr and Angstrom are allowed")

    def add_map(self, map, sign="p"):
        if (self.units != map.units):
            raise ValueError("Maps units mismatch")
        if ((self.ncol != map.ncol) or (self.nrow != map.nrow)):
            raise ValueError("Maps have different sizes")
        if (self.__scaled ^ map.__scaled):
            raise ValueError("Summing maps with different scaling")

        if (sign == "p"):
            self.density_map = self.density_map + map.density_map
        elif (sign == "m"):
            self.density_map = self.density_map - map.density_map
        else:
            raise ValueError("Error: add or subtract with p or m")

        self.__get_minmaxvalue()

    def log_scale(self):
        self.__scaled = True

        sign = 1
        for i in range(0, self.nrow):
            for j in range(0, self.ncol):
                sign = np.sign(self.density_map[i, j])
                self.density_map[i, j] = sign * math.log10(abs(self.density_map[i, j]))
        self.__get_minmaxvalue()

    def set_cmaprange(self, cmap_range: list, units="Angstrom"):
        if (len(cmap_range) > 2):
            raise ValueError("Error: invalid range selected")

        self.dens_min_value = float(cmap_range[0])
        self.dens_max_value = float(cmap_range[1])

    def plot(self, cmap="gnuplot", levels=100, add_colorbar=True):

        fig, ax = plt.subplots()
        ax.contourf(self.__meshx, self.__meshy, self.density_map, levels=levels, cmap=cmap)
        if self.units == "Angstrom":
            ax.set_xlabel("$\AA$")
            ax.set_ylabel("$\AA$")
        elif self.units == "Bohr":
            ax.set_xlabel("a.u.")
            ax.set_ylabel("a.u.")
        if self.__scaled:
            print('Scaled Plot')
        ax.set_aspect(1.0)
        ax.set_xlim(self.__meshx_min, self.__meshx_max)
        ax.set_ylim(0, self.__meshy_max * np.sqrt(1 - self.cosxy**2))
        if add_colorbar and not self.__scaled:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            mappable = ax.contourf(self.__meshx, self.__meshy, self.density_map, levels=levels, cmap=cmap)
            mappable.set_clim(self.dens_min_value, self.dens_max_value)
            fig.colorbar(mappable=mappable, cax=cax, orientation="vertical")

        return fig, ax
