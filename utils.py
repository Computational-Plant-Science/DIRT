import numpy as np
import scipy.optimize as sp
import yaml

from options import TraitOptions


def print_cli_header():
    print('------------------------------------------------------------')
    print('DIRT 1.1 - An automatic highthroughput root phenotyping platform')
    print('(c) 2014 Alexander Bucksch - bucksch@uga.edu')
    print('Web application by Abhiram Das - abhiram.das@gmail.com')
    print(' ')
    print('http://dirt.iplantcollaborative.org')
    print(' ')
    print('University of Georgia')
    print('------------------------------------------------------------')
    print(' ')


def get_traits(path: str = 'traits.yaml') -> TraitOptions:
    with open(path) as file:
        traits = yaml.safe_load(file)
        return TraitOptions.from_dict(**traits)


def compare_lists(l1, l2):
    l1 = l1[:len(l2)] if len(l1) > len(l2) else l1
    l2 = l2[:len(l1)] if len(l1) < len(l2) else l2
    half = len(l1) // 2

    if l1[half] == l2[half]:
        if half == 0:
            return 0
        elif half == len(l1) - 1:
            return len(l1) - 1
        elif l1[half + 1] == l2[half + 1]:
            split = compare_lists(l1[half:], l2[half:])
        else:
            split = 0
        return half + split
    else:
        if half == 0: return -1
        split = compare_lists(l1[:half], l2[:half])
        return split


def model_func(t, a, k, c):
    return a * np.exp(k * t) + c


def model_func_dia(t, a, k, c, dia):
    return a * dia ** (k * t) + c


def fit_exp_nonlinear(t, y):
    opt_parms, _ = sp.curve_fit(model_func, t, y, maxfev=100000)
    a, k, c = opt_parms
    fit_y = model_func(t, a, k, c)
    return fit_y, a, k, c


def fit_exp_linear(t, y, c):
    y = np.log(np.array(y) - c)
    k, a_log = np.polyfit(t, y, 1)
    a = np.exp(a_log)
    fit_y = model_func(np.array(t), a, k, c)
    return fit_y, a, k, c



