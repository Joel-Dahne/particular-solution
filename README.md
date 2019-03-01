# Particular-solution

Particular-solution is a C++ program for computing rigorous enclosures
of eigenvalues of the Laplacian on spherical triangles. It was
originally written during my internship for Bruno Salvy at ENS Lyon
the summer of 2018 but has been occasionally been updated since then.

More information about the algorithms used by the program and the
performance on some examples can be found in the report I write during
the internship. At the moment the report is not available online but
feel free to contact me if you are interested in it.

# Example

The main program is `build/examples/triangles` which supports
computing enclosures of eigenvalues for some predefined spherical
triangles. Calling the program without arguments give information
about how it can be used.

```
> ./build/examples/triangles
Compute eigenvalues using particular_solution.
Usage: triangles -i n [-prec p] [-N_beg b] [...]

-i n       - compute for triangle n (0 <= n <= 9), or "-i all",
           - or "-i regular" or "-i singular"
-output m    - output type
             - 0: no output
             - 1: enclosure for eigenvalue
             - 2: enclosure for nu
             - 3: midpoint of enclosure for eigenvalue
             - 4: midpoint of enclosure for nu
             - 5: width of enclosure for eigenvalue
             - 6: width of enclosure for nu
             - 7: coefficients for the approximate eigenfunction
             - 8: plot of the approximate eigenfunction
             - 9: minimizing nu value
-prec p      - precision in bits (default p = 64)
-N_beg b     - N value to start at (default 4)
-N_end e     - N value to stop at (default 16)
-N_step s    - step to take with N each iteration (default 2)
-tol eps     - relative accuracy goal each iteration (default 1e-5)
-plot n      - determines how to enclose eigenfunction with output 8 (default 0)
             -  -1: non-rigorous plot on points on the boundary
             -   0: simple interval enclosure
             - > 0: enclosure with Taylor expansion with n terms
-final       - flag for plotting only for the last N value (default off)
Implemented triangles:
T0 = (3pi/4, pi/3, pi/2)
T1 = (2pi/3, pi/3, pi/2)
T2 = (2pi/3, pi/4, pi/2)
T3 = (2pi/3, pi/3, pi/3)
T4 = (3pi/4, pi/4, pi/3)
T5 = (2pi/3, pi/4, pi/4)
T6 = (2pi/3, 3pi/4, 3pi/4)
T7 = (2pi/3, 2pi/3, 2pi/3)
T8 = (pi/2, 2pi/3, 3pi/4)
T9 = (pi/2, 2pi/3, 2pi/3)

```

The following computes rigorous enclosures of the eigenvalue close to
21.3094 for the spherical triangle having angles 2π/3, π/3 and π/3 (T3
in the list above).

```
> ./build/examples/triangles -i 3
4 [21.3094 +/- 1.49e-5]
6 [21.3094076 +/- 4.10e-8]
8 [21.3094076302 +/- 3.42e-11]
10 [21.309407630190 +/- 5.17e-13]
12 [21.309407630190445 +/- 4.88e-16]
14 [21.309407630190445259 +/- 8.18e-19]
16 [21.30940763019044525895 +/- 6.34e-21]
```

The output show a list of tighter and tighter enclosures, the first
number on each line indicates the number of terms used in the
expansion during the computation. So this proves that there is an
eigenvalue in the interval [21.30940763019044525895 +/- 6.34e-21]. An
important note here is that this is most likely the first eigenvalue
but that is not something that the program can prove. To prove that it
is the first eigenvalue you would need to find a lower bound, greater
than [21.30940763019044525895 +/- 6.34e-21], for the second eigenvalue
using some other method.

# Dependencies, compilation and tests

The program depends on Arb (http://www.arblib.org/), MPFR
(https://www.mpfr.org/), MPFR C++
(http://www.holoborodko.com/pavel/mpfr/) and Eigen
(https://eigen.tuxfamily.org/index.php?title=Main_Page).

If all of these dependencies are installed, on Debian they are all
available as packages, compiling is as simple as

``` shell
make examples
```

The program also has a built in test that runs the computations for
all predefined triangles comparing the results to previously computed
ones. The tests can be run with

``` shell
make check
```
