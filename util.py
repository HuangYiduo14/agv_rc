import numpy as np

def get_n_nodes(n_col, n_row):
    n_intersect = n_col * n_row
    n_half_shelf = 2 * (n_row - 1)
    n_full_shelf = (n_col - 2) * (n_row - 1)
    n_workstation = n_col - 1
    return n_intersect, n_half_shelf, n_full_shelf, n_workstation
