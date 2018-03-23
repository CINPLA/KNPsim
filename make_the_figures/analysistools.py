import h5py
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri


class Ion:
    def __init__(self, name, diffusion_constant, charge, index):
        self.name = name
        self.diffusion_constant = diffusion_constant
        self.charge = charge
        self.index = index

class AnalysisTools:
    def __init__(self, filename):
        self.F = 9.648e4 # Faradays constant, C/mol
        hdfpy = h5py.File(filename, 'r')
        ions = hdfpy['ions']
        self.N = len(ions.keys())
        self.ion_list = [0]*self.N

        for ion_name in ions.keys():
            ion_dataset = ions[ion_name]
            D = float(ion_dataset.attrs['diffusion constant'])
            z = float(ion_dataset.attrs['charge'])
            index = int(ion_dataset.attrs['index'])
            ion = Ion(ion_name, D, z, index)
            self.ion_list[index] = ion

        self.time_series = hdfpy['time'][()] # read whole dataset to numpy array
        element = hdfpy['FunctionSpace'].attrs['element']
        order = int(hdfpy['FunctionSpace'].attrs['order'])
        hdfpy.close()

        self.hdf = HDF5File(mpi_comm_world(), filename, 'r')
        self.mesh = Mesh()
        self.hdf.read(self.mesh,"geometry/mesh", False)
        # self.dim = self.mesh.geometry().dim()
        # element_types = [interval, triangle, tetrahedron]
        # element = element_types[self.dim-1]
        self.dim = self.mesh.geometry().dim()

        self.element_type_list = [interval, triangle, tetrahedron]
        self.P1 = FiniteElement("P", self.element_type_list[self.dim-1],
                                order)
        self.R0 = FiniteElement('R', self.element_type_list[self.dim-1], 0)
        self.V = FunctionSpace(self.mesh, self.P1)
        self.R = FunctionSpace(self.mesh, self.R0)

        # Create mixed element
        element_list = []
        for i in range(self.N+2):
            element_list.append(self.P1)

        for i in range(2):
            element_list.append(self.R0)

        TH = MixedElement(element_list)
        self.W = FunctionSpace(self.mesh, TH)

        self.attributes = self.hdf.attributes('attributes')
        self.vertex_to_dof_map = vertex_to_dof_map(self.V)

    def plot_initial_condition_fenics(self):
        u = Function(self.W)
        self.hdf.read(u, '/initial_state')
        plot(u.sub(self.N), title='Initial Potential [V]')
        rho = Function(self.V)

        for ion in self.ion_list:
            plot(u.sub(ion.index), title="Initial " + str(ion.name) + " Concentration [mM/L]")
            rho += project(u.sub(ion.index)*ion.charge)
        rho = self.F*rho
        plot(rho, title='Charge concentration [C/L]')

    def plot_nth_state_fenics(self, n):
        u = Function(self.W)
        self.hdf.read(u, '/solution/vector_'+str(n))
        time = self.time_series[n]
        plot(u.sub(self.N), title='Potential at t='+str(time)+' ms')
        rho = Function(self.V)
        for ion in self.ion_list:
            plot(u.sub(ion.index), title=str(ion.name)+' concentration at t='+str(time)+' ms')
            rho += project(u.sub(ion.index)*ion.charge)

        rho = self.F*rho
        plot(rho, title='Charge concentration [C/L]')

    def plot_initial_condition_matplotlib(self):
        u = Function(self.W)
        self.hdf.read(u, '/initial_state')
        coor = self.mesh.coordinates()

        if self.mesh.geometry().dim()==1:
            phi = Function(self.V)
            phi.assign(project(u.sub(self.N), self.V))

            phi_array = phi.vector().array()
            phi_array = phi_array[self.vertex_to_dof_map]
            plt.plot(coor, phi_array)
            plt.title('Initial potential')
            plt.ylabel('Potential [V]')
            plt.xlabel('Position [um]')
            plt.figure()
            plt.hold('on')
            rho = Function(self.V)
            for ion in self.ion_list:
                c = project(u.sub(ion.index), self.V)
                c_array = c.vector().array()
                c_array[:] = c_array[self.vertex_to_dof_map]
                plt.plot(coor,c_array,label=str(ion.name))
                rho += project(u.sub(ion.index)*ion.charge)

            plt.title('Initial ion concentations')
            plt.ylabel('Concentration [mM/L]')
            plt.xlabel('Position [um]')
            plt.legend()

            rho = self.F*rho
            rho_array = project(rho,self.V).vector().array()
            rho_array = rho_array[self.vertex_to_dof_map]
            plt.figure()
            plt.plot(coor, rho_array)
            plt.title('Initial charge concentation')
            plt.ylabel('Charge concentration [C/L]')
            plt.xlabel('Position [um]')

        elif self.mesh.geometry().dim()==2:
            phi = project(u.sub(self.N), self.V)
            phi_array = phi.vector().array()
            x = coor[:,0]
            y = coor[:,1]
            phi = phi_array[self.vertex_to_dof_map]

            z = phi
            fig = plt.figure()
            ax=fig.add_subplot(111)
            triang = tri.Triangulation(x,y)
            C = ax.tricontourf(triang,z, 20, cmap='RdBu', extend='both')
            cbar = fig.colorbar(C)
            cbar.set_label("Potential [V]")
            ax.tricontour(triang,z, 20, colors='k')
            ax.axis('image')
            plt.xlabel('x position [um]')
            plt.ylabel('y position [um]')
            plt.title('Initial Potential')

            for ion in self.ion_list:
                fig = plt.figure()
                ax=fig.add_subplot(111)
                ax.tricontourf(x,y,z, 20, cmap='RdBu', extend='both')
                ax.axis('image')
                plt.xlabel('x position [um]')
                plt.ylabel('y position [um]')
                plt.title('Inititial ' + str(ion.name) + ' concentration [mM]')


    def make_potential_difference_time_course(self, from_idx, to_idx):
        u = Function(self.W)
        t = []
        phi_t = []
        phi = Function(self.V)
        for n in range(from_idx, to_idx):
            if n % 100 == 0:
                print(n)
            self.hdf.read(u, '/solution/vector_'+str(n))
            phi.assign(project(u.sub(self.N), self.V))
            phi_array = phi.vector().array()
            t.append(self.time_series[n])
            phi_t.append(phi_array.max() - phi_array.min())
        return phi_t, t


    def plot_potential_difference_time_course(self, from_idx, to_idx):
        u = Function(self.W)
        t = []
        phi_t = []
        phi = Function(self.V)
        for n in range(from_idx, to_idx):
            if n % 100 == 0:
                print(n)
            self.hdf.read(u, '/solution/vector_'+str(n))
            phi.assign(project(u.sub(self.N), self.V))
            phi_array = phi.vector().array()
            t.append(self.time_series[n])
            phi_t.append(phi_array.max() - phi_array.min())
        fig = plt.figure()
        plt.plot(t, phi_t)
        plt.xlabel("Time [ms]")
        plt.ylabel("Potential difference [V]")
        return phi_t, t


    def plot_nth_state_matplotlib(self, n):
        u = Function(self.W)
        self.hdf.read(u, '/solution/vector_'+str(n))
        time = self.time_series[n]
        coor = self.mesh.coordinates()
        # plot(u.sub(self.N), title='Initial Potential')

        if self.mesh.geometry().dim()==1:
            phi = Function(self.V)
            phi.assign(project(u.sub(self.N), self.V))

            phi_array = phi.vector().array()
            phi_array = phi_array[self.vertex_to_dof_map]
            plt.plot(coor, phi_array)
            plt.title('Potential at t='+str(time)+' ms')
            plt.ylabel('Potential [V]')
            plt.xlabel('Position [um]')
            plt.figure()
            plt.hold('on')
            rho = Function(self.V)
            for ion in self.ion_list:
                # self.hdf.read(u, '/ions/'+ion['name']+'/initial_concentration')
                c = project(u.sub(ion.index), self.V)
                c_array = c.vector().array()
                c_array[:] = c_array[self.vertex_to_dof_map]
                plt.plot(coor,c_array,label=str(ion.name))
                rho += project(u.sub(ion.index)*ion.charge)

            plt.title('Ion concentrations at t='+str(time)+' ms')
            plt.ylabel('Concentration [mM/L]')
            plt.xlabel('Position [um]')
            plt.legend()

            rho = self.F*rho
            rho_array = project(rho,self.V).vector().array()
            rho_array = rho_array[self.vertex_to_dof_map]
            plt.figure()
            plt.plot(coor, rho_array)
            plt.title('Charge concentration at t='+str(time)+' ms')
            plt.ylabel('Charge concentration [C/L]')
            plt.xlabel('Position [um]')

        elif self.mesh.geometry().dim()==2:
            phi = project(u.sub(self.N), self.V)
            phi_array = phi.vector().array()
            x = coor[:,0]
            y = coor[:,1]
            phi = phi_array[self.vertex_to_dof_map]

            z = phi
            fig = plt.figure()
            ax=fig.add_subplot(111)
            # ax.tripcolor(x,y,z, cmap='RdBu')
            triang = tri.Triangulation(x,y)
            C = ax.tricontourf(triang,z, 20, cmap='RdBu', extend='both') # choose 20 contour levels, just to show how good its interpolation is
            cbar = fig.colorbar(C)
            cbar.set_label("Potential [V]")
            ax.tricontour(triang,z, 20, colors='k')
            ax.axis('image')
            plt.xlabel('x position [um]')
            plt.ylabel('y position [um]')
            plt.title('Potential at t='+str(time)+' ms')

            for ion in self.ion_list:
                fig = plt.figure()
                ax=fig.add_subplot(111)
                # ax.tripcolor(x,y,z, cmap='RdBu')
                # triang = tri.Triangulation(x,y)
                ax.tricontourf(x,y,z, 20, cmap='RdBu', extend='both') # choose 20 contour levels, just to show how good its interpolation is
                ax.axis('image')
                plt.xlabel('x position [um]')
                plt.ylabel('y position [um]')
                fig.colorbar(C)
                plt.title(str(ion.name) + ' concentrations at t='+str(time)+' ms')
