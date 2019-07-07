2019-07-07

Trevor and David,

I think I'm all caught up now. The latest version is 1.4 which is
published at the CPC Program Library, along with the paper. I've
downloaded all the source code and CG5 zip files.

Here's my view of the overall architecture of the system. The
computations are organized in a pipeline of three processes:

1. Compute the CG coefficients for SO(5) > SO(3) and save the results
in data files. This was done in Trevor's C program using the
techniques described in the 2009 Caprio+Rowe+Welsh paper. The
currently computed files support a range of Hamiltonians. More exotic
Hamiltionian will require this computation to be extended. I gather
this was a time-consuming calculation, so the results are computed
once and then re-used many times.

2. Input a Hamiltonian that models some nucleus and compute its matrix
in the ACM basis. This is done in the Maple code which reads the CG5
files, performs some computer algebra and ultimately produces a
numerical approximation to the matrix. The symbolic algebra may be
hard to translate to, say Python, although there may be suitable
libraries.

3. Call a linear algebra package to find eigenvalues, etc. of the
matrix. Maple calls CLAPACK which is a C library generated from LAPACK
(Fortran). The Python numpy package uses LAPACK. So this part should
have similar performance.

-- Arthur