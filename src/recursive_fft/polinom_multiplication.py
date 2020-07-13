"""
The purpose of this file is to demonstrate the Cooley–Tukey FFT algorithm. This
is not meant to be an efficient implementation. If you wish to use an efficient
one, use sympys own fft.

The scripts gets its arguments from the command line. Each command line argument
corresponds to a coefficient in increasing order. You may use the script by
running

    $ python3 recursive_fft.py 2 4 2 -1
    $ python3 recursive_fft.py "1 + 2*I" 2 3 "4*I"

where the numbers after `recursive_fft.py` are the command line arguments. A
separate bash script is provided, which runs the script with sample values.
"""

from sympy import simplify
from recursive_fft import fft
from recursive_fft_inverse import recursive_fft_inverse as fft_inverse


def recursive_fft_polynom_multiplication(polynom_a, polynom_b):
    """Polynom multiplication using Cooley-Turkey and its inverse algorithm"""
    len_polynom = len_polynom_a = len(polynom_a)
    len_polynom_b = len(polynom_b)
    if len_polynom & (len_polynom - 1) or not len_polynom_a == len_polynom_b:
        raise ValueError("Input criteria not met.")
    fft_a = fft(polynom_a + len_polynom * [0])
    fft_b = fft(polynom_b + len_polynom * [0])
    fft_a_times_b = [a * b for a, b in zip(fft_a, fft_b)]

    return fft_inverse(fft_a_times_b)


if __name__ == "__main__":
    POLYNOM_A = [-10, 1, -1, 7]
    POLYNOM_B = [3, -6, 0, 8]
    POLYNOM_A_TIMES_B = recursive_fft_polynom_multiplication(
        POLYNOM_A, POLYNOM_B)
    print([simplify(e) for e in POLYNOM_A_TIMES_B])
