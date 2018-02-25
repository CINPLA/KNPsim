import sys
from dolfin import *
from .ion import Ion
from .time_solver import Time_solver


class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class DuplicationError(Error):
    """Duplication error"""
    pass


class Simulator:
    """
    This class serves as a foundation for running a simulation. It should be
    the second class initialized by the user (after Geometry). It keeps track
    of all the other parts (Geometry, Time_solver, Potential, Ions, etc), as
    well as physical constants etc.
    The method initialize_simulator() is called when all these things are in
    place to set up the variational form and prepare for the call(s) to
    time_solver.solve_for_time_step().

    Args:
        geometry (Geometry): The geometry instance that holds the mesh for
            the simulation.
        temperature (float, optional): The temperature of the system.

    Note:
        All physical constants (appart from temperature) are hard-coded in the
        initializer.
    """
    def __init__(self, geometry, temperature=300.):
        self.ion_list = []
        self.dt = Expression("dt", dt=0, degree=2)
        self.geometry = geometry
        self.time_solver = None
        self.potential = None
        self.state_saver = None
        self.live_plotter = None
        self.T = temperature                    # temperature, K
        self.R = 8.314                          # gas constant, J/(K*mol)
        self.F = 9.648e4                        # Faradays constant, C/mol
        self.eps_0 = 8.854187e-12               # vacuum permittivity, F/m
        self.psi = self.R*self.T/self.F
        self.eps_r = 80
        self.epsilon = self.eps_0*self.eps_r
        self.deltas = []  # deltas are added later by making a Delta object
        self.N = 0  # number of ion species currently in the simulator

    def add_ion(self, ion):
        """
        This function adds an ion to the simulator instance. It is called by
        the constructor in Ion, and should never be called by the user in
        normal use.

        Args:
            ion (Ion): The Ion instance.
        """
        assert(isinstance(ion, Ion))
        self.ion_list.append(ion)
        ion.index = self.N
        self.N += 1

    def set_time_solver(self, time_solver):
        """
        This function sets the time solver for the system. Will raise an error
        if a time solver already exists. This function is called by the
        constructor in Time_solver, and should never be called by the user.

        Args:
            time_solver (Time_solver): The Time_solver instance.
        """
        assert(isinstance(time_solver, Time_solver))
        if self.time_solver is not None:
            raise DuplicationError("Can only set one time solver for the \
                system!")
        else:
            self.time_solver = time_solver

    def set_potential(self, potential):
        """
        This function sets the potential type for the system. Will raise an
        error if a time solver already exists. This function is called by the
        constructor in Potential, and should never be called by the user.

        Args:
            time_solver (Time_solver): The Time_solver instance.
        """
        assert(isinstance(potential, Potential))
        if self.potential is not None:
            raise DuplicationError("Can only set one potential for \
                                    the system!")
        else:
            self.potential = potential

    def add_point_source(self, delta):
        """
        This function adds a point source to the system.

        Args:
            delta (Delta): the point source.
        """
        self.deltas.append(delta)

    def initialize_simulator(self):
        """
        This function initialized the simulator, and should be called by the
        user when the Time_solver, Potential and all Ions and Deltas have been
        set.
        """
        assert(self.time_solver is not None)
        assert(self.potential is not None)
        self.bcs = []

        # Create mixed element
        element_list = []
        for i in range(self.N+2):
            element_list.append(self.geometry.P1)

        for i in range(2):
            element_list.append(self.geometry.R0)

        TH = MixedElement(element_list)
        self.geometry.W = FunctionSpace(self.geometry.mesh, TH)

        # Set up functions, testfunctions
        self.u = Function(self.geometry.W)
        self.u_new = Function(self.geometry.W)
        self.u_res = Function(self.geometry.W)
        self.v_list = TestFunctions(self.geometry.W)

        # Initialize all the ion functions with boundary cond and init cond.
        for i, ion in enumerate(self.ion_list):
            if ion.boundary is not None:
                bc = DirichletBC(self.geometry.W.sub(i),
                                 ion.boundary_condition,
                                 ion.boundary)
                self.bcs.append(bc)
            ion.c = self.u[i]
            assign(self.u.sub(i),
                   interpolate(ion.initial_condition, self.geometry.V))
            ion.c_new = self.u_new[i]
            assign(self.u_new.sub(i),
                   interpolate(ion.initial_condition, self.geometry.V))
            ion.v = self.v_list[i]

        # This should never happen:
        if self.potential.bc is not None:
            if MPI.rank(mpi_comm_world()) == 0:
                print("Why are you setting a dirichlet BC for the potential??")
            self.bcs.append(DirichletBC(self.geometry.W.sub(self.N),
                                        self.potential.bc, "on_boundary"))

        [bc.apply(self.u.vector()) for bc in self.bcs]
        [bc.apply(self.u_new.vector()) for bc in self.bcs]
        n = len(self.ion_list)

        # declare some variables related to the potential:
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

        # define conductance
        self.conductance = 0
        for i, ion in enumerate(self.ion_list):
            self.conductance = (self.conductance +
                                self.F*ion.D*ion.z**2*ion.c_new/self.psi)

        # set up variational form:
        self.set_form()

    def set_form(self):
        """
        This function is called by initialize_simulator() in order to set up
        the variational form. Should never be called by the user.
        """
        self.form = 0
        psi, dt = Constant(self.psi), Constant(self.time_solver.dt)
        phi_new = self.potential.phi_new + self.potential.phi_ps_new
        for i, ion in enumerate(self.ion_list):
            c, c_new, f, D, = ion.c, ion.c_new, ion.f, Constant(ion.D)
            z = Constant(ion.z)
            v = self.v_list[i]
            k = Constant(1/self.time_solver.dt)
            self.form += (k*(c_new - c)*v + inner(
                                            D*nabla_grad(c_new) +
                                            D*c_new*z*nabla_grad(phi_new)/psi,
                                            nabla_grad(v)) - f*v)*dx

        self.form = self.potential.set_form(self.form)
        self.w = TrialFunction(self.geometry.W)
        self.Jac = derivative(self.form, self.u_new, self.w)
        self.A = assemble(self.Jac)
