import sys
from dolfin import Function, Constant, interpolate, inner, nabla_grad, dx

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class Potential:
    """
    Contains the electric field used in the electro-diffusion solver. Should
    always be called from a subclass.
    """
    def __init__(self, simulator, bc=None):
        self.simulator = simulator
        simulator.potential = self
        self.bc = bc

    def set_form(self, form):
        """
        This function is called by simulator.set_form to set up the equations
        for the electric field. This function will
        """
        raise Error("Error! You should always use a subclass of Potential!")



class KirchoffPotential(Potential):
    """
    Contains details on potential spesific for the KNP formalism
    """
    def __init__(self, simulator, bc=None):
        Potential.__init__(self, simulator, bc)

    def set_form(self, form):
        rho = Function(self.simulator.geometry.V)
        for ion in self.simulator.ion_list:
            rho += self.simulator.F*ion.z*ion.c_new

        self.rho = rho
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

        form += Constant(F)*(inner(sigma*nabla_grad(self.phi_new) +
                                   nabla_grad(b), nabla_grad(v)) +
                             self.phi_new*d + v*self.dummy_new - f*v)*dx

        v, d = self.simulator.v_phi_ps, self.simulator.d_phi_ps
        form += F*(inner(sigma*nabla_grad(self.phi_ps_new), nabla_grad(v)) +
                   self.phi_ps_new*d + v*self.dummy_ps_new)*dx

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

        g = inner(nabla_grad(b)/sigma,
                  FacetNormal(self.simulator.geometry.mesh))
        F = Constant(self.simulator.F)
        eps = Constant(self.simulator.epsilon)

        rho = Function(self.simulator.geometry.V)
        for ion in self.simulator.ion_list:
            rho += F*ion.z*ion.c_new

        self.rho = rho
        form += (inner(nabla_grad(self.phi_new), nabla_grad(v)) +
                 self.dummy_new*v + self.phi_new*d - rho*v/eps)*dx
        form += g*v*ds

        v, d = self.simulator.v_phi_ps, self.simulator.d_phi_ps
        form += (inner(nabla_grad(self.phi_ps_new), nabla_grad(v)) +
                 self.phi_ps_new*d + v*self.dummy_ps_new)*dx

        return form


class ZeroPotential(Potential):
    """
    Returns a zero potential
    """
    def __init__(self, system):
        Potential.__init__(self, system)

    def set_form(self, form):
        # first set rho (not used in calculations for this potential)
        F = Constant(self.simulator.F)
        rho = Function(self.simulator.geometry.V)
        for ion in self.simulator.ion_list:
            rho += F*ion.z*ion.c_new
        self.rho = rho

        # set equations so that phi_diff and phi_vc are 0 everywhere
        v, d = self.simulator.v_phi, self.simulator.d_phi
        form += (inner(nabla_grad(self.phi_new), nabla_grad(v)) +
                 self.dummy_new*v + self.phi_new*d)*dx

        v, d = self.simulator.v_phi_ps, self.simulator.d_phi_ps
        form += (inner(nabla_grad(self.phi_ps_new), nabla_grad(v)) +
                 self.phi_ps_new*d + v*self.dummy_ps_new)*dx

        return form
