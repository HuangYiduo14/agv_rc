from tw_network import *
# create network and find all-point shortest path using Floyd Warshall
network1 = create_network(n_col, n_row, v_block_length, h_block_length)
network1.floyd_warshall()
network1.find_netflow_shortest_path()
# initialize time-window network
tw_network1 = TimeNetwork(network1)
info, travel_time, delay = tw_network1.time_window_routing(22, 23, 1)
info2, travel_time2, delay2 = tw_network1.time_window_routing(23, 22, 1)



"""
# simulation experiment
n_intersect, n_half_shelf, n_full_shelf, n_workstation = get_n_nodes(n_col, n_row)

# hfs: half shelf, fs: full shelf, station: workstation
print('start simulation:')
print('generating random matrices...')
demand_time = {}
demand_time[('half shelf', 'workstation')] = np.random.rand(n_t, n_half_shelf, n_workstation) < (p / 2)
demand_time[('full shelf', 'workstation')] = np.random.rand(n_t, n_full_shelf, n_workstation) < p
demand_time[('workstation', 'half shelf')] = np.random.rand(n_t, n_workstation, n_half_shelf) < (p / 2)
demand_time[('workstation', 'full shelf')] = np.random.rand(n_t, n_workstation, n_full_shelf) < p
demand_time[('workstation', 'workstation')] = np.random.rand(n_t, n_workstation, n_workstation) < (p * alpha)
possible_od = [('half shelf', 'workstation'), ('full shelf', 'workstation'), ('workstation', 'half shelf'),
               ('workstation', 'full shelf'), ('workstation', 'workstation')]

print('start experiment')
total_travel_time = 0
total_delay = 0
travel_time_record = []
delay_record = []
trip_v_agv_record = []
X_record = []
Y_record = []
agv_index = 0
all_trip = pd.DataFrame(columns=['n0', 'n1', 't0', 't01', 't1', 'agv_ind', 'current_loc'])
for t in range(100, 100 + n_t):
    print('time:', t, '-' * 50)
    # generate trips
    for od in possible_od:
        from_node, to_node = demand_transform(demand_time[od][t - 100], od[0], od[1], network1)
        if from_node.shape[0] == 0:
            continue
        for i in range(from_node.shape[0]):
            start = from_node[i]
            end = to_node[i]
            if start == end:
                continue
            new_agv = AGV(start, end, network1.node_list[start].x, network1.node_list[start].y, t, agv_index)
            # here we need to do AGV routing using timewindow based algorithm
            #trip_record, travel_time, delay, traj_table = new_agv.rc_shortest_path_routing(network1, traj_table)
            total_travel_time += travel_time
            total_delay += delay
            travel_time_record.append(travel_time)
            delay_record.append(delay)
            trip_record_i = pd.DataFrame(trip_record, columns=['n0', 'n1', 't0', 't01', 't1'])
            trip_record_i['agv_ind'] = agv_index
            trip_record_i['current_loc'] = -1
            trip_record_i.loc[0, 'current_loc'] = 0
            agv_index += 1
            all_trip = all_trip.append(trip_record_i, ignore_index=True)
    # keep record of locations
    # print(all_trip)
    existing_agv_table = all_trip.loc[(all_trip['current_loc'] >= 0) & (all_trip['current_loc'] < 100)].copy()
    print(existing_agv_table)
    if existing_agv_table.shape[0]==0:
        continue
    existing_agv_table['current_x'] = existing_agv_table.apply(
        lambda x: network1.arc_list[network1.node_arc_dict[(int(x.n0), int(x.n1))]].interpolate_x(x.current_loc),axis=1)
    existing_agv_table['current_y'] = existing_agv_table.apply(
        lambda x: network1.arc_list[network1.node_arc_dict[(int(x.n0), int(x.n1))]].interpolate_y(x.current_loc),axis=1)

    X_record.append(existing_agv_table['current_x'].values)
    Y_record.append(existing_agv_table['current_y'].values)
    # move to next
    all_trip.loc[(all_trip['t0'] <= t) & (all_trip['t01'] >= t), 'current_loc'] = 0
    all_trip.loc[all_trip['t1'] <= t, 'current_loc'] = 999
    all_trip.loc[(all_trip['current_loc'] >= 0) & (all_trip['current_loc'] < 100) & (all_trip['t01'] <= t) & (
                all_trip['t1'] >= t), 'current_loc'] += 1

plt.figure()
plt.plot(delay_record)
print('total delay', total_delay)
print('total travel time', total_travel_time)
print('total agvs', agv_index)

# animation
x_max = (network1.n_col - 1) * network1.h_block_length
y_max = (network1.n_row - 1) * network1.v_block_length
fig = plt.figure()
ax = fig.add_subplot(111,autoscale_on = False, xlim=(-1,x_max+1),ylim=(-1,y_max+1))
agvs, = ax.plot([],[],'bo',ms=10)
def init():
    agvs.set_data([],[])
    return agvs
def animate(i):
    agvs.set_data(X_record[i],Y_record[i])
    return agvs
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=400, interval=1)

print('saving...')
plt.rcParams['animation.ffmpeg_path'] = 'C:\\ffmpeg\\bin\\ffmpeg.exe'
writer = animation.FFMpegWriter()
anim.save('simulation_basic.mp4',writer=writer)
"""
