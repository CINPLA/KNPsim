from dolfin import *


class Geometry:
    """
    This class keeps track of mesh and the function space of the electro-
    diffusion solver. It should be the first class initialized by the user.

    Args:
        mesh (Mesh, list, or str): The mesh to simulate on. Can be a FEniCS
            Mesh instance, a string to a file containing a mesh, or a list of
            numbers specifying the number of elements in each spatial
            direction.
        element (str, optional): A string spesifying the element to be used.
            Passed to the FEniCS constructor FiniteElement().
        order (int, optional): The order of the element. Passed to
            FiniteElement()
    """
    def __init__(self, mesh, element='P', order=1):
        self.meshtype = mesh
        self.element = element
        self.order = order
        self.element_type_list = [interval, triangle, tetrahedron]

        # Set mesh
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

        # Set elements:
        self.P1 = FiniteElement(element, self.element_type_list[self.dim-1],
                                order)
        self.R0 = FiniteElement('R', self.element_type_list[self.dim-1], 0)
        self.V = FunctionSpace(mesh, self.P1)
        self.R = FunctionSpace(mesh, "R", self.R0)
