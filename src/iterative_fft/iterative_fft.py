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
import sys

from sympy import exp, pi, I
from sympy.parsing.sympy_parser import parse_expr


def bit_reverse_copy(coefficients, input_length):
    """Calculates the bit-reversal permutation of the coefficients"""
    bit_length = input_length.bit_length() - 1
    return [coefficients[int("{:0{}b}".format(i, bit_length)[::-1], 2)]
            for i in range(input_length)]


def fft(coefficients):
    """Implementation of Cooley–Tukey FFT algorithm"""
    input_length = len(coefficients)
    bit_length = input_length.bit_length() - 1
    bit_reversed_permutation = [
        coefficients[int("{:0{}b}".format(i, bit_length)[::-1], 2)]
        for i in range(input_length)
    ]
    for running_s in range(1, bit_length + 1):
        partition_size = 2 ** running_s
        principal_mth_root = exp(2*pi*I/partition_size)
        for k in range(0, input_length, partition_size):
            running_root = 1
            for j in range(int(partition_size/2)):
                sub_expr_t = running_root * \
                    bit_reversed_permutation[int(k + j + partition_size/2)]
                sub_expr_u = bit_reversed_permutation[k + j]
                bit_reversed_permutation[k + j] = sub_expr_u + sub_expr_t
                bit_reversed_permutation[int(
                    k + j + partition_size/2)] = sub_expr_u - sub_expr_t
                running_root *= principal_mth_root
    return bit_reversed_permutation


if __name__ == "__main__":
    COEFFICIENTS = sys.argv[1:]
    ARG_LENGTH = len(COEFFICIENTS)
    if ARG_LENGTH & (ARG_LENGTH - 1):
        raise ValueError("Number of coefficients is not a power of 2")
    print(fft([parse_expr(i) for i in COEFFICIENTS]))
