import sys
from dolfin import *
from ion import Ion
from time_solver import Time_solver



class Simulator:
    def __init__(self,geometry, T=300.):
        self.ion_list = []
        self.dt = Expression("dt", dt=0)
        self.geometry = geometry
        self.time_solver = None
        self.potential = None
        self.state_saver = None
        self.live_plotter = None
        self.T = T
        self.R = 8.314 # gas constant, J/(K*mol)
        self.F = 9.648e4 # Faradays constant, C/mol
        self.eps_0 = 8.854187e-12 # F/m
        self.psi = self.R*self.T/self.F
        self.eps_r = 80
        self.epsilon = self.eps_0*self.eps_r
        self.deltas = []
        self.N = 0 # number of ions in simulator


    def add_ion(self, ion):
        assert(isinstance(ion, Ion))
        self.ion_list.append(ion)
        ion.index = self.N
        self.N += 1

    def set_time_solver(self, time_solver):
        assert(isinstance(time_solver, Time_solver))
        if self.time_solver != None:
            print "Can only set one time solver for system! Exiting..."
            sys.exit(1)
        else:
            self.time_solver = time_solver

    def set_potential(self, potential):
        assert(isinstance(potential, Potential))
        if self.potential != None:
            print "Can only set one time solver for system! Exiting..."
            sys.exit(1)
        else:
            self.potential = potential


    def initialize_simulator(self):
        assert(self.time_solver != None)
        assert(self.potential != None)
        self.bcs = []
        # Function spaces for concentrations:
        function_space_list = [self.geometry.V]*(self.N)
        # Function spaces for potential:
        function_space_list.extend([self.geometry.V,self.geometry.V])
        # Function spaces for potential point sources:
        function_space_list.extend([self.geometry.R,self.geometry.R])

        self.geometry.W = MixedFunctionSpace(function_space_list)
        self.u = Function(self.geometry.W)
        self.u_new = Function(self.geometry.W)
        self.u_res = Function(self.geometry.W)
        self.v_list = TestFunctions(self.geometry.W)

        for i, ion in enumerate(self.ion_list):
            if ion.boundary != None:
                bc = DirichletBC(self.geometry.W.sub(i), ion.boundary_condition, ion.boundary)
                self.bcs.append(bc)
            ion.c = self.u[i]
            assign(self.u.sub(i), interpolate(ion.initial_condition, self.geometry.V))
            ion.c_new = self.u_new[i]
            assign(self.u_new.sub(i), interpolate(ion.initial_condition, self.geometry.V))
            ion.v = self.v_list[i]

        if self.potential.bc != None:
            self.bcs.append(DirichletBC(self.geometry.W.sub(self.N), self.potential.bc, "on_boundary"))

        [bc.apply(self.u.vector()) for bc in self.bcs]
        [bc.apply(self.u_new.vector()) for bc in self.bcs]
        n = len(self.ion_list)

        self.potential.phi = self.u[n]
        self.potential.phi_new = self.u_new[n]

        self.potential.dummy = self.u[n+2]
        self.potential.dummy_new = self.u_new[n+2]

        self.potential.phi_ps = self.u[n+1]
        self.potential.phi_ps_new = self.u_new[n+1]

        self.potential.dummy_ps = self.u[n+3]
        self.potential.dummy_ps_new = self.u_new[n+3]

        self.v_phi, self.d_phi = self.v_list[n], self.v_list[n+2]
        self.v_phi_ps, self.d_phi_ps = self.v_list[n+1], self.v_list[n+3]

        self.conductance = 0
        for i, ion in enumerate(self.ion_list):
            self.conductance = self.conductance + self.F*ion.D*ion.z**2*ion.c_new/self.psi
        self.set_form()

    def set_form(self):
        """
        This function is called by initialize_simulator in order to set up the variational form
        """
        self.form = 0
        psi, dt = Constant(self.psi), Constant(self.time_solver.dt)
        phi_new = self.potential.phi_new + self.potential.phi_ps_new
        for i, ion in enumerate(self.ion_list):
            c, c_new, f, D, z = ion.c, ion.c_new, ion.f, Constant(ion.D), Constant(ion.z)
            v = self.v_list[i]
            k = Constant(1/self.time_solver.dt)
            self.form += (k*(c_new - c)*v + inner(D*nabla_grad(c_new) + \
                D*c_new*z*nabla_grad(phi_new)/psi, nabla_grad(v)) - f*v)*dx

        self.form = self.potential.set_form(self.form)
        self.w = TrialFunction(self.geometry.W)
        self.Jac = derivative(self.form, self.u_new, self.w)
        self.A = assemble(self.Jac)

    def add_point_source(self,delta):
        self.deltas.append(delta)
