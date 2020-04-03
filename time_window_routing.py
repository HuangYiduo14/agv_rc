from tw_network import *
# create network and find all-point shortest path using Floyd Warshall
network1 = create_network(n_col, n_row, v_block_length, h_block_length)
network1.floyd_warshall()
network1.find_netflow_shortest_path()
# initialize time-window network
tw_network1 = TimeNetwork(network1)
info, travel_time, delay = tw_network1.time_window_routing(22, 23, 1)
info2, travel_time2, delay2 = tw_network1.time_window_routing(23, 22, 1)

