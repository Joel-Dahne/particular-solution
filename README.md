# Particular-solution

Particular-solution is a C++ program for computing rigorous enclosures
of eigenvalues of the Laplacian on spherical triangles. It was written
during my internship for Bruno Salvy at ENS Lyon the summer of 2018.

In the current state the program should mainly be seen as a proof of
concept of what can be done when computing eigenvalues of the
Laplacian on spherical triangles. While the program, up to bugs, gives
rigorous results it is not made to be user friendly. It also contains
no automatic tests so is very likely to contain at least some bugs.

# Example

The following computes rigorous enclosures of the eigenvalue close to
24.4569 for the spherical triangle having angles $3\pi/4$, $\pi/4$ and
$\pi/3$. The computations are done by starting with an expansion with
2 terms and then every step increasing the number of terms by 2 up to
32.

    ./build/particular-solution -b 2 -s 2 -e 32 -o 4 -n 4.470604591 -w 0.1 -p 128 -- 3 4 1 4 1 3

The output is:

    2 [+/- 49.9]
    4 [2.4e+1 +/- 0.833]
    6 [24.46 +/- 7.12e-3]
    8 [24.46 +/- 4.36e-3]
    10 [24.4569 +/- 6.85e-5]
    12 [24.4569 +/- 2.67e-5]
    14 [24.45691 +/- 5.10e-6]
    16 [24.456914 +/- 3.57e-7]
    18 [24.4569138 +/- 3.35e-8]
    20 [24.45691380 +/- 5.52e-9]
    22 [24.456913796 +/- 9.50e-10]
    24 [24.4569137963 +/- 4.51e-11]
    26 [24.4569137963 +/- 1.48e-11]
    28 [24.45691379630 +/- 2.59e-12]
    30 [24.456913796299 +/- 4.08e-13]
    32 [24.4569137962991 +/- 6.45e-14]

So this proves that there is an eigenvalue in the interval
[24.4569137962991 +/- 6.45e-14]. An important note here is that this
is most likely the first eigenvalue but that is not something that the
program can prove. To prove that it is the first eigenvalue you would
need to find a lower bound, greater than [24.4569137962991 +/-
6.45e-14], for the second eigenvalue using some other method.

# Dependencies and compilation

The program depends on Arb (http://www.arblib.org/), MPFR
(https://www.mpfr.org/), MPFR C++
(http://www.holoborodko.com/pavel/mpfr/) and Eigen
(https://eigen.tuxfamily.org/index.php?title=Main_Page).

If all of these dependencies are installed, on Debian they are all
available as packages, compiling is as simple as

    make all
