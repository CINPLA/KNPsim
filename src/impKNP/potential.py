import sys
from dolfin import *
class Potential:
    """
    Contains the electric field used in the electro-diffusion solver
    """
    def __init__(self, simulator, bc=None):
        self.simulator = simulator
        simulator.potential = self
        self.bc = bc

    def set_form(self, form):
        print "Error! The Potential super class should be thought of as abstract!"
        sys.exit(1)

#
class KirchoffPotential(Potential):
    """
    Contains details on potential spesific for the KNP formalism
    """
    def __init__(self, simulator, bc=None):
        Potential.__init__(self, simulator, bc)

    def find_initial_potential(self):
        V = self.simulator.geometry.V
        R = self.simulator.geometry.R
        W = MixedFunctionSpace([V, R])
        v,d = TestFunctions(W)

        a = interpolate(Constant(0), V)
        b = interpolate(Constant(0), V)
        f = interpolate(Constant(0), V)

        F = self.simulator.F
        R = self.simulator.R
        T = self.simulator.T

        psi = self.simulator.R*self.simulator.T/self.simulator.F


        for ion in self.simulator.ion_list:
            a += (1./psi)*ion.z**2*ion.D*ion.c
            b += ion.z*ion.D*ion.c
            ion.f.t = self.simulator.time_solver.t
            f += ion.z*ion.f

        (phi_init, c) = TrialFunction(W)
        form = F*(inner(a*nabla_grad(phi_init) + nabla_grad(b), nabla_grad(v)) + phi_init*d + v*c - f*v)*dx
        phi_W = Function(W)
        (self.a, self.L) = system(form)
        self.A = assemble(self.a)
        self.b = assemble(self.L)


        for delta in self.simulator.deltas:
            I = 0
            for current in delta.currents:
                i_j = current.magnitude_function(self.simulator.time_solver.t)
                I += i_j
                p = PointSource(V, delta.point, I)
            p.apply(self.b)

        solve(self.A, phi_W.vector(), self.b)
        assign(self.simulator.u.sub(self.simulator.N),phi_W.sub(0))


    def set_form(self,form):
        v, d = self.simulator.v_phi, self.simulator.d_phi
        a = interpolate(Constant(0), self.simulator.geometry.V)
        b = interpolate(Constant(0), self.simulator.geometry.V)
        f = interpolate(Constant(0), self.simulator.geometry.V)
        F = self.simulator.F
        R = self.simulator.R
        T = self.simulator.T
        psi = Constant(self.simulator.R*self.simulator.T/self.simulator.F)
        theta = self.simulator.time_solver.theta


        for ion in self.simulator.ion_list:
            a += (1./psi)*ion.z**2*ion.D*ion.c_new
            b += ion.z*ion.D*ion.c_new
            f += ion.z*ion.f

        form += F*(inner(a*nabla_grad(self.phi_new) + nabla_grad(b), nabla_grad(v)) + self.phi_new*d + v*self.dummy_new - f*v)*dx
        return form



class PoissonPotential(Potential):
    """
    Contains details on potential spesific for the PNP formalism
    """
    def __init__(self, simulator, bc=None):
        Potential.__init__(self, simulator, bc)

    def find_initial_potential(self):
        V = self.simulator.geometry.V
        R = self.simulator.geometry.R
        W = MixedFunctionSpace([V, R])
        v,d = TestFunctions(W)

        F = Constant(self.simulator.F)
        eps = Constant(self.simulator.epsilon)

        rho = Function(self.simulator.geometry.V)
        for ion in self.simulator.ion_list:
            rho += F*ion.z*ion.c

        (phi_init, c) = TrialFunction(W)
        g = Constant(0)
        form = (inner(nabla_grad(phi_init), nabla_grad(v)) + c*v + phi_init*d - rho*v/eps)*dx + g*v*ds
        phi_W = Function(W)
        (self.a, self.L) = system(form)
        self.A = assemble(self.a)
        self.b = assemble(self.L)

        solve(self.A, phi_W.vector(), self.b)
        assign(self.simulator.u.sub(self.simulator.N),phi_W.sub(0))





    def set_form(self, form):

        v, d = self.simulator.v_phi, self.simulator.d_phi
        psi = self.simulator.psi
        a = interpolate(Constant(0), self.simulator.geometry.V)
        b = interpolate(Constant(0), self.simulator.geometry.V)
        for ion in self.simulator.ion_list:
            a += (1./psi)*ion.z**2*ion.D*ion.c
            b += ion.z*ion.D*ion.c

        g = inner(nabla_grad(b)/a, FacetNormal(self.simulator.geometry.mesh))
        F = Constant(self.simulator.F)
        eps = Constant(self.simulator.epsilon)

        rho = Function(self.simulator.geometry.V)
        for ion in self.simulator.ion_list:
            rho += F*ion.z*ion.c_new


        form += (inner(nabla_grad(self.phi_new), nabla_grad(v)) + self.dummy_new*v + self.phi_new*d - rho*v/eps)*dx
        form += g*v*ds
        return form


#
# class NoField(Potential):
#     """
#     Returns a zero potential
#     """
#     def __init__(self, system):
#         Potential.__init__(self, system)
#
#     def set_form(self):
#         (phi, c) = TrialFunction(self.W)
#         (v, d) = TestFunctions(self.W)
