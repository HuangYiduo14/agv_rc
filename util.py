import numpy as np

def get_n_nodes(n_col, n_row):
    n_intersect = n_col * n_row
    n_half_shelf = 2 * (n_row - 1)
    n_full_shelf = (n_col - 2) * (n_row - 1)
    n_workstation = n_col - 1
    return n_intersect, n_half_shelf, n_full_shelf, n_workstation

def get_node_lists(n_col,n_row):
    n_intersect, n_half_shelf, n_full_shelf, n_workstation = get_n_nodes(n_col, n_row)
    intersect_list = list(range(n_intersect))
    half_shelf_list = list(range(n_intersect,n_intersect+n_half_shelf))
    full_shelf_list = list(range(n_intersect+n_half_shelf,n_intersect+n_half_shelf+n_full_shelf))
    workstaion_list = list(range(n_intersect+n_half_shelf+n_full_shelf,n_intersect+n_half_shelf+n_full_shelf+n_workstation))
    return intersect_list, half_shelf_list, full_shelf_list, workstaion_list