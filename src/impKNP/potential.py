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

    def set_form(self,form):
        v, d = self.simulator.v_phi, self.simulator.d_phi
        sigma = interpolate(Constant(0), self.simulator.geometry.V)
        b = interpolate(Constant(0), self.simulator.geometry.V)
        f = interpolate(Constant(0), self.simulator.geometry.V)
        F = self.simulator.F
        R = self.simulator.R
        T = self.simulator.T
        psi = Constant(self.simulator.R*self.simulator.T/self.simulator.F)

        for ion in self.simulator.ion_list:
            sigma += (1./psi)*ion.z**2*ion.D*ion.c_new
            b += ion.z*ion.D*ion.c_new
            f += ion.z*ion.f

        form += F*(inner(sigma*nabla_grad(self.phi_new) + nabla_grad(b), nabla_grad(v)) + self.phi_new*d + v*self.dummy_new - f*v)*dx
        return form



class PoissonPotential(Potential):
    """
    Contains details on potential spesific for the PNP formalism
    """
    def __init__(self, simulator, bc=None):
        Potential.__init__(self, simulator, bc)

    def set_form(self, form):
        v, d = self.simulator.v_phi, self.simulator.d_phi
        psi = self.simulator.psi
        sigma = interpolate(Constant(0), self.simulator.geometry.V)
        b = interpolate(Constant(0), self.simulator.geometry.V)
        for ion in self.simulator.ion_list:
            sigma += (1./psi)*ion.z**2*ion.D*ion.c
            b += ion.z*ion.D*ion.c

        g = inner(nabla_grad(b)/sigma, FacetNormal(self.simulator.geometry.mesh))
        F = Constant(self.simulator.F)
        eps = Constant(self.simulator.epsilon)

        rho = Function(self.simulator.geometry.V)
        for ion in self.simulator.ion_list:
            rho += F*ion.z*ion.c_new


        form += (inner(nabla_grad(self.phi_new), nabla_grad(v)) + self.dummy_new*v + self.phi_new*d - rho*v/eps)*dx
        form += g*v*ds
        return form




class ZeroPotential(Potential):
    """
    Returns a zero potential
    """
    def __init__(self, system):
        Potential.__init__(self, system)

    def set_form(self, form):
        v, d = self.simulator.v_phi, self.simulator.d_phi
        form += (inner(nabla_grad(self.phi_new), nabla_grad(v)) + self.dummy_new*v + self.phi_new*d)*dx
        return form
