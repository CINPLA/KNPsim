from dolfin import *


class Geometry:
    """
    This class keeps track of mesh and the function space of the electro
    diffusion solver.
    """
    def __init__(self, mesh, space='CG', order=1):
        self.meshtype = mesh
        self.space = space
        self.order = order
        if isinstance(mesh, list):
            domain_type = [UnitIntervalMesh, UnitSquareMesh, UnitCubeMesh]
            self.dim = len(mesh)
            if len(mesh) == self.dim:
                self.mesh = domain_type[self.dim-1](*mesh)
            else:
                print("dimension mismatch in set_geometry! mesh does not match \
                    dimension")
                print(str(self.dim))
                print(str(len(mesh)))
                sys.exit()

        elif isinstance(mesh, str):
            print("interpreting mesh input as filename...")
            try:
                self.mesh = Mesh(mesh)
                self.dim = self.mesh.geometry().dim()
            except IOError:
                print("Could not find the file spesified, exiting....")
                sys.exit(1)

        elif isinstance(mesh, Mesh):
            self.mesh = mesh
            self.dim = self.mesh.geometry().dim()
        else:
            print("input not understood! Exiting...")
            sys.exit(1)

        self.V = FunctionSpace(mesh, space, order)
        self.R = FunctionSpace(mesh, "R", 0)
