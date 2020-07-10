"""
This script tests my own implementation of the fft against the one sympy uses.
To use this, run

    $ python3 test_recursive_fft.py

"""

import unittest

from parameterized import parameterized
from sympy.discrete.transforms import fft as sympy_fft
from sympy.parsing.sympy_parser import parse_expr
from sympy import simplify
from recursive_fft import fft


def read_coefficients(filename):
    """Reads the coefficients from the file name given"""
    with open(filename) as coefficients_file:
        set_of_raw_coefficients = [line.split() for line in coefficients_file]
    return [[i, [parse_expr(coefficient) for coefficient in coefficients]]
            for i, coefficients in enumerate(set_of_raw_coefficients)]


class TestingClass(unittest.TestCase):
    """Testing the correctness of the recursive fft"""

    @parameterized.expand(read_coefficients('coefficients.txt'))
    def test_fft(self, test_number, coefficients):
        """Testing each cases separately"""
        print("\nCoefficient test #{}. Coefficients tested: {}".format(
            test_number, coefficients))
        recursive_solution = fft(coefficients)
        sympy_solution = sympy_fft(coefficients)
        self.assertEqual(len(recursive_solution), len(sympy_solution))
        for r_coef, s_coef in zip(recursive_solution, sympy_solution):
            self.assertEqual(simplify(r_coef - s_coef), 0)


if __name__ == "__main__":
    unittest.main()
