from CRYSTALClear.crystal_io import Properties_output
from CRYSTALClear.units import au_to_angstrom, angstrom_to_au
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import math

class ElectricChargeDensity(Properties_output):
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

class OrbitalMagnetizationField(Properties_output):
    def __init__(self, fort25_path, units="Angstrom", quivdens=20):
        self.which_prop='m'
        super().read_vecfield(fort25_path, which_prop=self.which_prop)
        
        self.dens_m = angstrom_to_au(angstrom_to_au(self.dens_m))
        self.__get_modulo()
        self.__quiver_density(quivdens)
        self.__generate_mesh()
        self.__generate_projections()
        self.__get_minmaxvalue()
        self.units = "Angstrom"
        self.change_units(units)

        self.__scaled = False

    def __quiver_density(self, quivdens):
        """
        Large values of quivdens reduce the number of arrows.
        """

        for i in range(1, quivdens):
            if (self.nrow % i) == 0:
                self.__nrow_split = int(self.nrow/i)
                self.__step_nrow = i

        for i in range(1, quivdens):
            if (self.ncol % i) == 0:
                self.__ncol_split = int(self.ncol/i)
                self.__step_ncol = i

    def __get_modulo(self):
        self.mod_dens = np.zeros((self.nrow, self.ncol), dtype=float)

        for i in range(0, self.nrow):
            for j in range(0, self.ncol):
                self.mod_dens[i,j] = np.sqrt(self.dens_m[i,j,0]**2 + self.dens_m[i,j,1]**2 + \
                        self.dens_m[i,j,2]**2)
    
    def __generate_mesh(self):

        # Mesh for the modulo plot 
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

        # Mesh for the projections
        self.__meshprojx = np.zeros((self.__nrow_split, self.__ncol_split), dtype=float)
        self.__meshprojy = np.zeros((self.__nrow_split, self.__ncol_split), dtype=float)
        r = 0
        s = 0
        for i in range(0, self.__nrow_split):
            for j in range(0, self.__ncol_split):
                self.__meshprojx[i, j] = self.__meshx[r, s]
                self.__meshprojy[i, j] = self.__meshy[r, s]
                s += self.__step_ncol
            s = 0
            r += self.__step_nrow

    def __generate_projections(self):
        """
        Calculates the projections of the vector V from the field on the ABC plane
        """
        
        cb = self.c - self.b
        ab = self.a - self.b
        if self.cosxy == 0:
            mod_ab = np.linalg.norm(ab)
            mod_cb = np.linalg.norm(cb)
            v1 = ab/mod_ab
            v2 = cb/mod_cb
        else:
            ab = ab - ((np.dot(ab, cb))/(np.dot(ab, ab)))*cb
            mod_ab = np.linalg.norm(ab)
            mod_cb = np.linalg.norm(cb)
            v1 = ab/mod_ab
            v2 = cb/mod_cb
        
        abc_normal = np.cross(ab, cb)
        abc_normal_mod = np.linalg.norm(abc_normal)

        self.projx = np.zeros((self.__nrow_split, self.__ncol_split), dtype=float)
        self.projy = np.zeros((self.__nrow_split, self.__ncol_split), dtype=float)
        for i in range(0, self.__nrow_split):
            for j in range(0, self.__ncol_split):
                vec = np.array([self.dens_m[i*self.__step_nrow, j*self.__step_ncol, 0],
                                self.dens_m[i*self.__step_nrow, j*self.__step_ncol, 1],
                                self.dens_m[i*self.__step_nrow, j*self.__step_ncol, 2]])
                step = (1/(abc_normal_mod**2))*np.cross(abc_normal, np.cross(vec, abc_normal))
                self.projx[i, j] = np.dot((step), v2)
                self.projy[i, j] = np.dot((step), v1)

    def __get_minmaxvalue(self):
        self.dens_min_value = np.amin(self.mod_dens)
        self.dens_max_value = np.amax(self.mod_dens)
 
    def change_units(self, units: str):
        if (self.units == "Angstrom" and units == "Bohr"):
            self.units = "Bohr"
            self.a = angstrom_to_au(self.a)
            self.b = angstrom_to_au(self.b)
            self.c = angstrom_to_au(self.c)
            self.dens_m = au_to_angstrom(au_to_angstrom(self.dens_m))
            self.__get_modulo()
            self.__generate_mesh()
            self.__generate_projections()
            self.__get_minmaxvalue()
        elif (self.units == "Bohr" and units == "Angstrom"):
            self.units = "Angstrom"
            self.a = au_to_angstrom(self.a)
            self.b = au_to_angstrom(self.b)
            self.c = au_to_angstrom(self.c)
            self.dens_m = angstrom_to_au(angstrom_to_au(self.dens_m))
            self.__get_modulo()
            self.__generate_mesh()
            self.__generate_projections()
            self.__get_minmaxvalue()
        elif (self.units == units):
            pass
        else:
            raise ValueError("Only the keywords Bohr and Angstrom are allowed")

    def add_map(self, map, sign="p"):
        if (self.units != map.units):
            raise ValueError("Maps units mismatch")
        if (self.which_prop != map.which_prop):
            raise ValueError("Trying to sum two different quantities, you silly goose!")
        if ((self.ncol != map.ncol) or (self.nrow != map.nrow)):
            raise ValueError("Maps have different sizes")
        if (self.__scaled ^ map.__scaled):
            raise ValueError("Summing maps with different scaling")

        if (sign == "p"):
            self.dens_m = self.dens_m + map.dens_m
        elif (sign == "m"):
            self.dens_m = self.dens_m - map.dens_m
        else:
            raise ValueError("Error: add or subtract with p or m")

        self.__get_modulo()
        self.__generate_projections()
        self.__get_minmaxvalue()

    def log_scale(self):
        self.__scaled = True

        sign = 1
        for i in range(0, self.nrow):
            for j in range(0, self.ncol):
                sign = np.sign(self.mod_dens[i, j])
                self.mod_dens[i, j] = sign * math.log10(abs(self.mod_dens[i, j]))

        self.__get_modulo()
        self.__generate_projections()
        self.__get_minmaxvalue()

    def set_cmaprange(self, cmap_range: list):
        if (len(cmap_range) > 2):
            raise ValueError("Error: invalid range selected")

        self.dens_min_value = float(cmap_range[0])
        self.dens_max_value = float(cmap_range[1])

    def plot(self, cmap="gnuplot", levels=100, add_colorbar=True, scale=0.8):

        fig, ax = plt.subplots()
        ax.contourf(self.__meshx, self.__meshy, self.mod_dens, levels=levels, cmap=cmap)
        ax.quiver(self.__meshprojx, self.__meshprojy, self.projx, self.projy, scale=scale)
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
            mappable = ax.contourf(self.__meshx, self.__meshy, self.mod_dens, levels=levels, cmap=cmap)
            mappable.set_clim(self.dens_min_value, self.dens_max_value)
            fig.colorbar(mappable=mappable, cax=cax, orientation="vertical")
            ax.quiver(self.__meshprojx, self.__meshprojy, self.projx, self.projy, scale=scale)

        return fig, ax

class ChargeCurrentField(Properties_output):
    def __init__(self, fort25_path, units="Angstrom", quivdens=20):
        self.which_prop='j'
        super().read_vecfield(fort25_path, which_prop=self.which_prop)
        
        self.dens_j = angstrom_to_au(angstrom_to_au(self.dens_j))
        self.__get_modulo()
        self.__quiver_density(quivdens)
        self.__generate_mesh()
        self.__generate_projections()
        self.__get_minmaxvalue()
        self.units = "Angstrom"
        self.change_units(units)

        self.__scaled = False

    def __quiver_density(self, quivdens):
        """
        Large values of quivdens reduce the number of arrows.
        """

        for i in range(1, quivdens):
            if (self.nrow % i) == 0:
                self.__nrow_split = int(self.nrow/i)
                self.__step_nrow = i

        for i in range(1, quivdens):
            if (self.ncol % i) == 0:
                self.__ncol_split = int(self.ncol/i)
                self.__step_ncol = i

    def __get_modulo(self):
        self.mod_dens = np.zeros((self.nrow, self.ncol), dtype=float)

        for i in range(0, self.nrow):
            for j in range(0, self.ncol):
                self.mod_dens[i,j] = np.sqrt(self.dens_j[i,j,0]**2 + self.dens_j[i,j,1]**2 + \
                        self.dens_j[i,j,2]**2)
    
    def __generate_mesh(self):

        # Mesh for the modulo plot 
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

        # Mesh for the projections
        self.__meshprojx = np.zeros((self.__nrow_split, self.__ncol_split), dtype=float)
        self.__meshprojy = np.zeros((self.__nrow_split, self.__ncol_split), dtype=float)
        r = 0
        s = 0
        for i in range(0, self.__nrow_split):
            for j in range(0, self.__ncol_split):
                self.__meshprojx[i, j] = self.__meshx[r, s]
                self.__meshprojy[i, j] = self.__meshy[r, s]
                s += self.__step_ncol
            s = 0
            r += self.__step_nrow

    def __generate_projections(self):
        """
        Calculates the projections of the vector V from the field on the ABC plane
        """
        
        cb = self.c - self.b
        ab = self.a - self.b
        if self.cosxy == 0:
            mod_ab = np.linalg.norm(ab)
            mod_cb = np.linalg.norm(cb)
            v1 = ab/mod_ab
            v2 = cb/mod_cb
        else:
            ab = ab - ((np.dot(ab, cb))/(np.dot(ab, ab)))*cb
            mod_ab = np.linalg.norm(ab)
            mod_cb = np.linalg.norm(cb)
            v1 = ab/mod_ab
            v2 = cb/mod_cb
        
        abc_normal = np.cross(ab, cb)
        abc_normal_mod = np.linalg.norm(abc_normal)

        self.projx = np.zeros((self.__nrow_split, self.__ncol_split), dtype=float)
        self.projy = np.zeros((self.__nrow_split, self.__ncol_split), dtype=float)
        for i in range(0, self.__nrow_split):
            for j in range(0, self.__ncol_split):
                vec = np.array([self.dens_j[i*self.__step_nrow, j*self.__step_ncol, 0],
                                self.dens_j[i*self.__step_nrow, j*self.__step_ncol, 1],
                                self.dens_j[i*self.__step_nrow, j*self.__step_ncol, 2]])
                step = (1/(abc_normal_mod**2))*np.cross(abc_normal, np.cross(vec, abc_normal))
                self.projx[i, j] = np.dot((step), v2)
                self.projy[i, j] = np.dot((step), v1)

    def __get_minmaxvalue(self):
        self.dens_min_value = np.amin(self.mod_dens)
        self.dens_max_value = np.amax(self.mod_dens)
 
    def change_units(self, units: str):
        if (self.units == "Angstrom" and units == "Bohr"):
            self.units = "Bohr"
            self.a = angstrom_to_au(self.a)
            self.b = angstrom_to_au(self.b)
            self.c = angstrom_to_au(self.c)
            self.dens_j = au_to_angstrom(au_to_angstrom(self.dens_j))
            self.__get_modulo()
            self.__generate_mesh()
            self.__generate_projections()
            self.__get_minmaxvalue()
        elif (self.units == "Bohr" and units == "Angstrom"):
            self.units = "Angstrom"
            self.a = au_to_angstrom(self.a)
            self.b = au_to_angstrom(self.b)
            self.c = au_to_angstrom(self.c)
            self.dens_j = angstrom_to_au(angstrom_to_au(self.dens_j))
            self.__get_modulo()
            self.__generate_mesh()
            self.__generate_projections()
            self.__get_minmaxvalue()
        elif (self.units == units):
            pass
        else:
            raise ValueError("Only the keywords Bohr and Angstrom are allowed")

    def add_map(self, map, sign="p"):
        if (self.units != map.units):
            raise ValueError("Maps units mismatch")
        if (self.which_prop != map.which_prop):
            raise ValueError("Trying to sum two different quantities, you silly goose!")
        if ((self.ncol != map.ncol) or (self.nrow != map.nrow)):
            raise ValueError("Maps have different sizes")
        if (self.__scaled ^ map.__scaled):
            raise ValueError("Summing maps with different scaling")

        if (sign == "p"):
            self.dens_j = self.dens_j + map.dens_j
        elif (sign == "m"):
            self.dens_j = self.dens_j - map.dens_j
        else:
            raise ValueError("Error: add or subtract with p or m")

        self.__get_modulo()
        self.__generate_projections()
        self.__get_minmaxvalue()

    def log_scale(self):
        self.__scaled = True

        sign = 1
        for i in range(0, self.nrow):
            for j in range(0, self.ncol):
                sign = np.sign(self.mod_dens[i, j])
                self.mod_dens[i, j] = sign * math.log10(abs(self.mod_dens[i, j]))

        self.__get_modulo()
        self.__generate_projections()
        self.__get_minmaxvalue()

    def set_cmaprange(self, cmap_range: list):
        if (len(cmap_range) > 2):
            raise ValueError("Error: invalid range selected")

        self.dens_min_value = float(cmap_range[0])
        self.dens_max_value = float(cmap_range[1])

    def plot(self, cmap="gnuplot", levels=100, add_colorbar=True, scale=0.8):

        fig, ax = plt.subplots()
        ax.contourf(self.__meshx, self.__meshy, self.mod_dens, levels=levels, cmap=cmap)
        ax.quiver(self.__meshprojx, self.__meshprojy, self.projx, self.projy, scale=scale)
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
            mappable = ax.contourf(self.__meshx, self.__meshy, self.mod_dens, levels=levels, cmap=cmap)
            mappable.set_clim(self.dens_min_value, self.dens_max_value)
            fig.colorbar(mappable=mappable, cax=cax, orientation="vertical")
            ax.quiver(self.__meshprojx, self.__meshprojy, self.projx, self.projy, scale=scale)

        return fig, ax

class SpinCurrentField(Properties_output):
    def __init__(self, fort25_path, units="Angstrom", quivdens=20):
        self.which_prop='J'
        super().read_vecfield(fort25_path, which_prop=self.which_prop)
        
        self.dens_JX = angstrom_to_au(angstrom_to_au(self.dens_JX))
        self.dens_JY = angstrom_to_au(angstrom_to_au(self.dens_JY))
        self.dens_JZ = angstrom_to_au(angstrom_to_au(self.dens_JZ))
        self.__get_modulo()
        self.__quiver_density(quivdens)
        self.__generate_mesh()
        self.__generate_projections()
        self.__get_minmaxvalue()
        self.units = "Angstrom"
        self.change_units(units)

        self.__scaled = False

    def __quiver_density(self, quivdens):
        """
        Large values of quivdens reduce the number of arrows.
        """

        for i in range(1, quivdens):
            if (self.nrow % i) == 0:
                self.__nrow_split = int(self.nrow/i)
                self.__step_nrow = i

        for i in range(1, quivdens):
            if (self.ncol % i) == 0:
                self.__ncol_split = int(self.ncol/i)
                self.__step_ncol = i

    def __get_modulo(self):
        self.mod_dens_JX = np.zeros((self.nrow, self.ncol), dtype=float)
        self.mod_dens_JY = np.zeros((self.nrow, self.ncol), dtype=float)
        self.mod_dens_JZ = np.zeros((self.nrow, self.ncol), dtype=float)

        for i in range(0, self.nrow):
            for j in range(0, self.ncol):
                self.mod_dens_JX[i,j] = np.sqrt(self.dens_JX[i,j,0]**2 + self.dens_JX[i,j,1]**2 + \
                        self.dens_JX[i,j,2]**2)
                self.mod_dens_JY[i,j] = np.sqrt(self.dens_JY[i,j,0]**2 + self.dens_JY[i,j,1]**2 + \
                        self.dens_JY[i,j,2]**2)
                self.mod_dens_JZ[i,j] = np.sqrt(self.dens_JZ[i,j,0]**2 + self.dens_JZ[i,j,1]**2 + \
                        self.dens_JZ[i,j,2]**2)
    
    def __generate_mesh(self):

        # Mesh for the modulo plot 
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

        # Mesh for the projections
        self.__meshprojx = np.zeros((self.__nrow_split, self.__ncol_split), dtype=float)
        self.__meshprojy = np.zeros((self.__nrow_split, self.__ncol_split), dtype=float)
        s = 0
        r = 0
        for i in range(0, self.__nrow_split):
            for j in range(0, self.__ncol_split):
                self.__meshprojx[i, j] = self.__meshx[r, s]
                self.__meshprojy[i, j] = self.__meshy[r, s]
                s += self.__step_ncol
            s = 0
            r += self.__step_nrow

    def __generate_projections(self):
        """
        Calculates the projections of the vector V from the field on the ABC plane
        """
        
        cb = self.c - self.b
        ab = self.a - self.b
        if self.cosxy == 0:
            mod_ab = np.linalg.norm(ab)
            mod_cb = np.linalg.norm(cb)
            v1 = ab/mod_ab
            v2 = cb/mod_cb
        else:
            ab = ab - ((np.dot(ab, cb))/(np.dot(ab, ab)))*cb
            mod_ab = np.linalg.norm(ab)
            mod_cb = np.linalg.norm(cb)
            v1 = ab/mod_ab
            v2 = cb/mod_cb
        
        abc_normal = np.cross(ab, cb)
        abc_normal_mod = np.linalg.norm(abc_normal)

        self.projx_JX = np.zeros((self.__nrow_split, self.__ncol_split), dtype=float)
        self.projy_JX = np.zeros((self.__nrow_split, self.__ncol_split), dtype=float)
        self.projx_JY = np.zeros((self.__nrow_split, self.__ncol_split), dtype=float)
        self.projy_JY = np.zeros((self.__nrow_split, self.__ncol_split), dtype=float)
        self.projx_JZ = np.zeros((self.__nrow_split, self.__ncol_split), dtype=float)
        self.projy_JZ = np.zeros((self.__nrow_split, self.__ncol_split), dtype=float)
        for i in range(0, self.__nrow_split):
            for j in range(0, self.__ncol_split):
                vec = np.array([self.dens_JX[i*self.__step_nrow, j*self.__step_ncol, 0],
                                self.dens_JX[i*self.__step_nrow, j*self.__step_ncol, 1],
                                self.dens_JX[i*self.__step_nrow, j*self.__step_ncol, 2]])
                step = (1/(abc_normal_mod**2))*np.cross(abc_normal, np.cross(vec, abc_normal))
                self.projx_JX[i, j] = np.dot((step), v2)
                self.projy_JX[i, j] = np.dot((step), v1)
                vec = np.array([self.dens_JY[i*self.__step_nrow, j*self.__step_ncol, 0],
                                self.dens_JY[i*self.__step_nrow, j*self.__step_ncol, 1],
                                self.dens_JY[i*self.__step_nrow, j*self.__step_ncol, 2]])
                step = (1/(abc_normal_mod**2))*np.cross(abc_normal, np.cross(vec, abc_normal))
                self.projx_JY[i, j] = np.dot((step), v2)
                self.projy_JY[i, j] = np.dot((step), v1)
                vec = np.array([self.dens_JZ[i*self.__step_nrow, j*self.__step_ncol, 0],
                                self.dens_JZ[i*self.__step_nrow, j*self.__step_ncol, 1],
                                self.dens_JZ[i*self.__step_nrow, j*self.__step_ncol, 2]])
                step = (1/(abc_normal_mod**2))*np.cross(abc_normal, np.cross(vec, abc_normal))
                self.projx_JZ[i, j] = np.dot((step), v2)
                self.projy_JZ[i, j] = np.dot((step), v1)

    def __get_minmaxvalue(self):
        self.dens_JX_min_value = np.amin(self.dens_JX)
        self.dens_JX_max_value = np.amax(self.dens_JX)
        self.dens_JY_min_value = np.amin(self.dens_JY)
        self.dens_JY_max_value = np.amax(self.dens_JY)
        self.dens_JZ_min_value = np.amin(self.dens_JZ)
        self.dens_JZ_max_value = np.amax(self.dens_JZ)

        self.dens_min_value = np.amin([self.mod_dens_JX, self.mod_dens_JY,
                                      self.mod_dens_JZ])
        self.dens_max_value = np.amax([self.mod_dens_JX, self.mod_dens_JY,
                                      self.mod_dens_JZ])
 
    def change_units(self, units: str):
        if (self.units == "Angstrom" and units == "Bohr"):
            self.units = "Bohr"
            self.a = angstrom_to_au(self.a)
            self.b = angstrom_to_au(self.b)
            self.c = angstrom_to_au(self.c)
            self.dens_JX = au_to_angstrom(au_to_angstrom(self.dens_JX))
            self.dens_JY = au_to_angstrom(au_to_angstrom(self.dens_JY))
            self.dens_JZ = au_to_angstrom(au_to_angstrom(self.dens_JZ))
            self.__get_modulo()
            self.__generate_mesh()
            self.__generate_projections()
            self.__get_minmaxvalue()
        elif (self.units == "Bohr" and units == "Angstrom"):
            self.units = "Angstrom"
            self.a = au_to_angstrom(self.a)
            self.b = au_to_angstrom(self.b)
            self.c = au_to_angstrom(self.c)
            self.dens_JX = angstrom_to_au(angstrom_to_au(self.dens_JX))
            self.dens_JY = angstrom_to_au(angstrom_to_au(self.dens_JY))
            self.dens_JZ = angstrom_to_au(angstrom_to_au(self.dens_JZ))
            self.__get_modulo()
            self.__generate_mesh()
            self.__generate_projections()
            self.__get_minmaxvalue()
        elif (self.units == units):
            pass
        else:
            raise ValueError("Only the keywords Bohr and Angstrom are allowed")

    def add_map(self, map, sign="p"):
        if (self.units != map.units):
            raise ValueError("Maps units mismatch")
        if (self.which_prop != map.which_prop):
            raise ValueError("Trying to sum two different quantities, you silly goose!")
        if ((self.ncol != map.ncol) or (self.nrow != map.nrow)):
            raise ValueError("Maps have different sizes")
        if (self.__scaled ^ map.__scaled):
            raise ValueError("Summing maps with different scaling")

        if (sign == "p"):
            self.dens_JX = self.dens_JX + map.dens_JX
            self.dens_JY = self.dens_JY + map.dens_JY
            self.dens_JZ = self.dens_JZ + map.dens_JZ
        elif (sign == "m"):
            self.dens_JX = self.dens_JX - map.dens_JX
            self.dens_JY = self.dens_JY - map.dens_JY
            self.dens_JZ = self.dens_JZ - map.dens_JZ
        else:
            raise ValueError("Error: add or subtract with p or m")

        self.__get_modulo()
        self.__generate_projections()
        self.__get_minmaxvalue()

    def log_scale(self):
        self.__scaled = True

        sign = 1
        for i in range(0, self.nrow):
            for j in range(0, self.ncol):
                sign = np.sign(self.mod_dens_JX[i, j])
                self.mod_dens_JX[i, j] = sign * math.log10(abs(self.mod_dens_JX[i, j]))
                sign = np.sign(self.mod_dens_JY[i, j])
                self.mod_dens_JY[i, j] = sign * math.log10(abs(self.mod_dens_JY[i, j]))
                sign = np.sign(self.mod_dens_JZ[i, j])
                self.mod_dens_JZ[i, j] = sign * math.log10(abs(self.mod_dens_JZ[i, j]))

        self.__get_modulo()
        self.__generate_projections()
        self.__get_minmaxvalue()

    def set_cmaprange(self, cmap_range: list):
        if (len(cmap_range) > 2):
            raise ValueError("Error: invalid range selected")

        self.dens_min_value = float(cmap_range[0])
        self.dens_max_value = float(cmap_range[1])

    def plot(self, cmap="gnuplot", levels=100, add_colorbar=True, scale=0.8):

        fig1, ax1 = plt.subplots()
        ax1.contourf(self.__meshx, self.__meshy, self.mod_dens_JX, levels=levels, cmap=cmap)
        ax1.quiver(self.__meshprojx, self.__meshprojy, self.projx_JX, self.projy_JX, scale=scale)
        if self.units == "Angstrom":
            ax1.set_xlabel("$\AA$")
            ax1.set_ylabel("$\AA$")
        elif self.units == "Bohr":
            ax1.set_xlabel("a.u.")
            ax1.set_ylabel("a.u.")
        if self.__scaled:
            print('Scaled Plot')
        ax1.set_aspect(1.0)
        ax1.set_xlim(self.__meshx_min, self.__meshx_max)
        ax1.set_ylim(0, self.__meshy_max * np.sqrt(1 - self.cosxy**2))
        if add_colorbar:  
            divider1 = make_axes_locatable(ax1)
            mappable1 = ax1.contourf(self.__meshx, self.__meshy, self.mod_dens_JX, levels=levels, cmap=cmap)
            mappable1.set_clim(self.dens_min_value, self.dens_max_value)
            cax1 = divider1.append_axes('right', size='5%', pad=0.05)
            fig1.colorbar(mappable=mappable1, cax=cax1, orientation="vertical")
            ax1.quiver(self.__meshprojx, self.__meshprojy, self.projx_JX, self.projy_JX, scale=scale)

        fig2, ax2 = plt.subplots()
        ax2.contourf(self.__meshx, self.__meshy, self.mod_dens_JY, levels=levels, cmap=cmap)
        ax2.quiver(self.__meshprojx, self.__meshprojy, self.projx_JY, self.projy_JY, scale=scale)
        if self.units == "Angstrom":
            ax2.set_xlabel("$\AA$")
            ax2.set_ylabel("$\AA$")
        elif self.units == "Bohr":
            ax2.set_xlabel("a.u.")
            ax2.set_ylabel("a.u.")
        if self.__scaled:
            print('Scaled Plot')
        ax2.set_aspect(1.0)
        ax2.set_xlim(self.__meshx_min, self.__meshx_max)
        ax2.set_ylim(0, self.__meshy_max * np.sqrt(1 - self.cosxy**2))
        if add_colorbar:  
            divider2 = make_axes_locatable(ax2)
            mappable2 = ax2.contourf(self.__meshx, self.__meshy, self.mod_dens_JY, levels=levels, cmap=cmap)
            mappable2.set_clim(self.dens_min_value, self.dens_max_value)
            cax2 = divider2.append_axes('right', size='5%', pad=0.05)
            fig2.colorbar(mappable=mappable2, cax=cax2, orientation="vertical")
            ax2.quiver(self.__meshprojx, self.__meshprojy, self.projx_JY, self.projy_JY, scale=scale)

        fig3, ax3 = plt.subplots()
        ax3.contourf(self.__meshx, self.__meshy, self.mod_dens_JZ, levels=levels, cmap=cmap)
        ax3.quiver(self.__meshprojx, self.__meshprojy, self.projx_JZ, self.projy_JZ, scale=scale)
        if self.units == "Angstrom":
            ax3.set_xlabel("$\AA$")
            ax3.set_ylabel("$\AA$")
        elif self.units == "Bohr":
            ax3.set_xlabel("a.u.")
            ax3.set_ylabel("a.u.")
        if self.__scaled:
            print('Scaled Plot')
        ax3.set_aspect(1.0)
        ax3.set_xlim(self.__meshx_min, self.__meshx_max)
        ax3.set_ylim(0, self.__meshy_max * np.sqrt(1 - self.cosxy**2))
        if add_colorbar:  
            divider3 = make_axes_locatable(ax3)
            mappable3 = ax3.contourf(self.__meshx, self.__meshy, self.mod_dens_JZ, levels=levels, cmap=cmap)
            mappable3.set_clim(self.dens_min_value, self.dens_max_value)
            cax3 = divider3.append_axes('right', size='5%', pad=0.05)
            fig3.colorbar(mappable=mappable3, cax=cax3, orientation="vertical")
            ax3.quiver(self.__meshprojx, self.__meshprojy, self.projx_JZ, self.projy_JZ, scale=scale)

        return fig1, ax1, fig2, ax2, fig3, ax3
