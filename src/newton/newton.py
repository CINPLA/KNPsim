from fenics import *
import numpy as np
import time


def Newton_manual(
        J,
        F,
        u,
        u_res,
        bcs=[],
        deltas=[],
        atol=1e-12,
        rtol=1e-12,
        max_it=20,
        relax=1.0,
        report_convergence=True):
    parameters['form_compiler']['optimize'] = True
    parameters['form_compiler']['cpp_optimize'] = True
    # parameters['form_compiler']['cpp_optimize_flags'] = '-O3'
    # Reset counters
    Iter = 0
    residual = 1
    rel_res = residual
    while rel_res > rtol and residual > atol and Iter < max_it:
        # Assemble system
        t0 = time.clock()
        A = assemble(J)
        b = assemble(-F)
        t1 = time.clock()
        if MPI.rank(mpi_comm_world()) == 0:
            print("assembly took " + str(t1 - t0) + " seconds!")

        # Solve linear system
        [bc.apply(A, b, u.vector()) for bc in bcs]
        [delta.apply(b) for delta in deltas]
        t0 = time.clock()

        solve(A, u_res.vector(), b)
        t1 = time.clock()
        if MPI.rank(mpi_comm_world()) == 0:
            print("solving took " + str(t1 - t0) + " seconds!")

        # Update solution
        u.vector().axpy(relax, u_res.vector())
        [bc.apply(u.vector()) for bc in bcs]

        # Compute residual
        residual = b.norm('l2')

        if Iter == 0:
            rel_res0 = residual
            rel_res = 1
        else:
            rel_res = residual / rel_res0

        if MPI.rank(mpi_comm_world()) == 0:
            if report_convergence:
                print(
                    "Newton iteration %d: r (atol) = %.3e (tol = %.3e), \
                    r (rel) = %.3e (tol = %.3e) "
                    % (Iter, residual, atol, rel_res, rtol)
                    )
        Iter += 1
