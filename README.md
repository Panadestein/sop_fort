# SOP-FBR-FORT

The FORTRAN implementation of the objective function for the SOP-FBR optimization.
Now all function evaluations (including the loop computing the energies of a list of points)
happen at the FORTRAN level. To compile the FORTRAN module use the following F2PY command:

	f2py -m rho_sop --f90flags='-fopenmp' -lgomp -c sop_fbr.f90

As can be seen, the code has been parallelized with OPENMP. The average CPU time gain with
respect to Python is ~90 times.
