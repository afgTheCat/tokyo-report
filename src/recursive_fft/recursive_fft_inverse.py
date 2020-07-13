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
from sympy import nsimplify


def recursive_fft_inverse(coefficients, original_length=0):
    """Implementation of the inverse Cooley–Tukey FFT algorithm"""
    input_length = len(coefficients)
    original_length = max(input_length, original_length)
    if input_length == 1:
        return [nsimplify(coefficients[0]/original_length)]
    principal_nth_root = exp(-2*pi*I/input_length)
    running_root = 1
    even_indexed_coefficients = [coefficients[i]
                                 for i in range(0, input_length, 2)]
    odd_indexed_coefficients = [coefficients[i]
                                for i in range(1, input_length, 2)]
    even_fft = recursive_fft_inverse(
        even_indexed_coefficients, original_length=original_length)
    odd_fft = recursive_fft_inverse(
        odd_indexed_coefficients, original_length=original_length)
    left_fft, right_fft = [], []
    for k in range(0, int(input_length/2)):
        sub_expr = running_root * odd_fft[k]
        left_fft.append(even_fft[k] + sub_expr)
        right_fft.append(even_fft[k] - sub_expr)
        running_root *= principal_nth_root
    return left_fft + right_fft


if __name__ == "__main__":
    COEFFICIENTS = sys.argv[1:]
    ARG_LENGTH = len(COEFFICIENTS)
    if ARG_LENGTH & (ARG_LENGTH - 1):
        raise ValueError("Number of coefficients is not a power of 2")
    print(recursive_fft_inverse([parse_expr(i) for i in COEFFICIENTS]))
