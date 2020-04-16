from agv_class import *
from rc_simulation import warehouse_simulation
import logging
# create network and find all-point shortest path using Floyd Warshall
network1 = create_network(n_col, n_row, v_block_length, h_block_length)
network1.floyd_warshall()
network1.find_netflow_shortest_path()
# initialize time-window network
tw_network1 = TimeNetwork(network1)
# simulation experiment
demand_time = simulation_demand()
# hfs: half shelf, fs: full shelf, station: workstation
#warehouse_simulation(tw_network1, routing_type='tw', save_anim=True, ani_name='tw_routing_find_max_max')
simulation_record = pd.DataFrame(columns=['prob','total_delay','total_travel_time','agv_number'])
#p_list = [0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02]
p_list = [0.003,0.0035, 0.004]
for i, p in enumerate(p_list):
    print('tw routing', p, '=' * 50)
    total_delay, total_travel_time, agv_index= warehouse_simulation(tw_network1, routing_type='tw', p=p)
    simulation_result = pd.DataFrame({'prob':p,'total_delay':[total_delay], 'total_travel_time':[total_travel_time],'agv_number':[agv_index]})
    simulation_record = simulation_record.append(simulation_result, ignore_index=True)
simulation_record.to_csv('expdata/exp_tw1.csv')




'''
result:
total delay 7460.0
total travel time 20510.0
total agvs 565
saving...
'''