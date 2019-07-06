<h1 align="center">Additive FFT</h1>


## Overview

These routines (additive_fft and add_fft) computes the additive fft on a computational friendly basis of a finite field extension.

The first routines mimics the approach of the fft in multiplicative mode (with a n-th primitive root of 1). The second approach is faster but requires precomputations that make it unpracticable on this implementation (there are way to make it fast but we have to get rid of sage).

There is also an algorithm to compute irreducible polynomials of degree charac^n where n is chosen and charac is the characteristic. It outpus a computational friendly basis.

## Usage

*For the moment the characteristic 2 is not supported*.

Example:
```bash
//setup an extension of GF(3) of degree 3**2
field, basis = setup_basis(3, 2)
R.<x> = PolynomialRing(field)

// defining a polynomial
f = 0
for i in range(3**5):
    f += (i%3)*x**i

// evaluate the polynomial on basis[:8] (3**8 elements)
t = additive_fft(f, basis[:8], field)
elapsed time: 4.93681716919

naive_evaluation(f, basis[:8], 3)
elapsed time: 113.064409018
```

The reference that I used is this [paper](https://core.ac.uk/download/pdf/82655328.pdf)

