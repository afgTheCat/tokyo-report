"""A working fft implementation"""
from sympy import exp, pi, I

def fft(terms):
    """Implementation of Cooley–Tukey FFT algorithm"""
    input_length = len(terms)
    if input_length == 1:
        return terms
    minus_nth_root_of_unity = exp(-2*pi*I/input_length)
    running_root = 1
    even_indexed_terms = terms[0::2]
    odd_indexed_terms = terms[1::2]
    even_fft = fft(even_indexed_terms)
    odd_fft = fft(odd_indexed_terms)
    left_fft, right_fft = [], []
    for k in range(0, input_length // 2):
        sub_expr = running_root * odd_fft[k]
        left_fft.append(even_fft[k] + sub_expr)
        right_fft.append(even_fft[k] - sub_expr)
        running_root *= minus_nth_root_of_unity
    return left_fft + right_fft
