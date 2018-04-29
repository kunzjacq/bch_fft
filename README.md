# BCH encoder/decoder

C libraries and demonstration programs for BCH encoding and decoding on a characteristic-2 finite field, using fast operations on polynomials:
- fast binary polynomial multiplication provided by the external library [gf2x](https://gforge.inria.fr/projects/gf2x/)
- Fast polynomial evaluation through either additive or multiplicative Fast Fourier transform. The performance of the multiplicative fast Fourier transform depends on the smoothness of 2^n-1 if the finite field used has 2^n elements.

a compiled gf2x library should be put into external/lib in order to compile the project with the help of the CMake build file provided. The gf2x header in external/include should be updated if necessary.