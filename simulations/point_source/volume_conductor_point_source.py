import sys
sys.path.append("/home/andreavs/Dropbox/knpsim_backward_euler")
from impKNP import *
from dolfin import *
import time
# parameters['form_compiler']['optimize'] = True

p0 = Point(0,0)
p1 = Point(100,100)
mesh = RectangleMesh(p0, p1, 100, 100)
geometry = Geometry(mesh)
simulator = Simulator(geometry)

print("loaded mesh and made spaces!")

def boundary(x, on_boundary):
    return on_boundary

init_cond_Na = Expression('150', degree=2)
z_Na = 1
D_Na = 1.33
init_Na = init_cond_Na
c_boundary_Na = init_cond_Na
ion_Na = Ion(simulator, z_Na, D_Na, init_Na, c_boundary_Na, boundary, "Na")

init_cond_K = Expression('3', degree=2)
z_K = 1
D_K = 2.
init_K = init_cond_K
c_boundary_K = init_cond_K
ion_K = Ion(simulator, z_K, D_K, init_K, c_boundary_K, boundary, "K")


def mag_func(t):
    return 1*(t<1e2)


def neg_mag_func(t):
    return -mag_func(t)

current_1 = Current(mag_func, ion_K)
currents = [current_1]
delta = Delta(Point(20,50), currents)
simulator.add_point_source(delta)

current_2 = Current(neg_mag_func, ion_K)
currents = [current_2]
delta = Delta(Point(80,50), currents)
simulator.add_point_source(delta)

# delta = Delta(ion_K, Point(20,50), mag_func)
# simulator.add_point_source(delta)
#
# delta = Delta(ion_K, Point(80,50), neg_mag_func)
# simulator.add_point_source(delta)

init_cond_Cl = Expression('153', degree=2)
z_Cl = -1
D_Cl = 1.
init_Cl = init_cond_Cl
c_boundary_Cl = init_cond_Cl
# f_Cl = Expression("x[0]*x[0]*(100-x[0])*(100-x[0])*t", t=0)
ion_Cl = Ion(simulator, z_Cl, D_Cl, init_Cl, c_boundary_Cl, boundary, "Cl")

dt = 1e-3
time_solver = Time_solver(simulator, dt, theta=1, t_stop=1e-0)
potential = KirchoffPotential(simulator)

print("initializing")

simulator.initialize_simulator()

print("initialized!")

conductance = Function(simulator.geometry.V)
for ion in simulator.ion_list:
    conductance += ion.D*ion.z*ion.z*ion.c

conductance = conductance/simulator.psi

W = MixedFunctionSpace([simulator.geometry.V, simulator.geometry.R])

(u, c) = TrialFunction(W)
(v, d) = TestFunctions(W)

a = (conductance*inner(nabla_grad(u), nabla_grad(v)) + c*v + u*d)*dx
L = Constant(0)*v*dx

A, b = assemble_system(a,L)
for delta in simulator.deltas:
    I = 0
    for current in delta.currents:
        i_j = current.magnitude_function(simulator.time_solver.t)
        I += i_j
        # if current.ion != None:
        #     pointsources.append(PointSource(self.system.geometry.W.sub(current.ion.index), delta.point, i_j/self.simulator.F))
    p = PointSource(W.sub(0), delta.point, I)
    p.apply(b)

u = Function(W)
solve(A, u.vector(), b)

u_phi = project(u.sub(0),simulator.geometry.V)

fname = "/home/andreavs/Dropbox/knpsim_backward_euler/simulations/point_source_all_modes/vc.h5"
notes = "This simulation considers a point source in a 2d grid, with vc"
state_saver = State_saver(fname, simulator, notes)

assign(simulator.u.sub(simulator.N), u_phi)

time_solver.t = 0
t_stop = 1e2
dt = 1e-1

while time_solver.t < t_stop:
    state_saver.save_state()
    time_solver.t += dt
    time_solver.t_list.append(time_solver.t)

t_stop = 2e2
u = Function(W)
u_phi = project(u.sub(0),simulator.geometry.V)
assign(simulator.u.sub(simulator.N), u_phi)

while time_solver.t < t_stop:
    state_saver.save_state()
    time_solver.t += dt
    time_solver.t_list.append(time_solver.t)

state_saver.finalize()
# plot(u_phi)
# interactive()
