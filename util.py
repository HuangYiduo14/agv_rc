import numpy as np
from heapq import heappush, heappop, heapify
big_M = 99999
crossing_time = {'straight':2, 'angle':2, 'source':1}
AGV_length = 1
# network parameters
h_block_length = 6
v_block_length = 12
n_col = 4
n_row = 4
lp = 2
# simulation parameters
p = 0.005  # intensity of flow
n_t = 500 # simulation time
alpha = 0.4  # relative demand from station to station

def demand_transform(demand_list, type1: str, type2: str, network):
    valid_agv = np.where(demand_list)
    from_node_list = network.node_type_list[type1][valid_agv[0]]
    to_node_list = network.node_type_list[type2][valid_agv[1]]
    return from_node_list, to_node_list

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

def net_crossing_type(n0,n1,n2):
    return 'straight'