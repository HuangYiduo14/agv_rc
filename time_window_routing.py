from agv_class import *
from rc_simulation import warehouse_simulation

# create network and find all-point shortest path using Floyd Warshall
network1 = create_network(n_col, n_row, v_block_length, h_block_length)
network1.floyd_warshall()
network1.find_netflow_shortest_path()
# initialize time-window network
tw_network1 = TimeNetwork(network1)
# simulation experiment
demand_time = simulation_demand()
# hfs: half shelf, fs: full shelf, station: workstation
warehouse_simulation(tw_network1, routing_type='tw_rc_heuristic_gen', save_anim=False, ani_name='tw_routing')

'''
heuristic for RC
'''

