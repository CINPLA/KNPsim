from fenics import *
import numpy as np
import time


def Newton_manual(J, F, u, u_res, bcs, deltas, atol, rtol, max_iter,
                  relax, report_convergence):
    parameters['form_compiler']['optimize'] = True
    parameters['form_compiler']['cpp_optimize'] = True
    # parameters['form_compiler']['cpp_optimize_flags'] = '-O3'

    # Reset counters
    Iter = 0
    residual = 1
    rel_res = residual

    # Iterate until the residual criteria is meet, or max iterations
    while rel_res > rtol and residual > atol and Iter < max_iter:
        # Assemble system
        t0 = time.clock()
        A = assemble(J)
        b = assemble(-F)
        t1 = time.clock()
        if MPI.rank(Mesh().mpi_comm()) == 0:
            print("Assemble Jacobian took {:02.03f} seconds!".format(t1 - t0))

        # Solve linear system
        [bc.apply(A, b, u.vector()) for bc in bcs]
        [delta.apply(b) for delta in deltas]
        t0 = time.clock()

        solve(A, u_res.vector(), b)
        t1 = time.clock()
        if MPI.rank(Mesh().mpi_comm()) == 0:
            print("Linear solve took {:02.03f} seconds!".format(t1 - t0))

        # Update solution
        u.vector().axpy(relax, u_res.vector())
        [bc.apply(u.vector()) for bc in bcs]

        # Compute residual
        residual = b.norm('l2')

        if Iter == 0:
            rel_res0 = residual
            rel_res = 1
        else:
            rel_res = residual/rel_res0

        if MPI.rank(Mesh().mpi_comm()) == 0:
            if report_convergence:
                print(("Newton iteration %d: r(atol) = %.3e (tol=%.3e), r" +
                      "(rel) = %.3e (tol=%.3e)\n") % (Iter, residual, atol,
                      rel_res, rtol))
        Iter += 1
