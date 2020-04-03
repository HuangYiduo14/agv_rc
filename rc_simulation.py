from network import *
from agv_class import *
from matplotlib import animation
import pandas as pd

np.random.seed(10)


def rc_fleet_info(network, t_hat, lp: 1):
    max_sim = 10000
    virtual_fleet = {}
    agv_index = 0
    for t in range(max_sim):
        if t % t_hat == 0:
            # release horizontal AGVs
            for k in range(network.n_row // 2):
                start1 = k * 2
                end1 = k * 2 + (network.n_col - 1) * network.n_row
                start2 = k * 2 + 1 + (network.n_col - 1) * network.n_row
                end2 = k * 2 + 1
                for j in range(lp):
                    virtual_fleet[agv_index] = [[start1, end1], t + j]
                    agv_index += 1
                for j in range(lp):
                    virtual_fleet[agv_index] = [[start2, end2], t + j]
                    agv_index += 1
        if t / t_hat % 1 == 0.5:
            # release verticle AGVs
            for k in range(network.n_col // 2):
                start1 = k * 2 * network.n_row + network.n_row - 1
                end1 = k * 2 * network.n_row
                start2 = (k * 2 + 1) * network.n_row
                end2 = (k * 2 + 1) * network.n_row + network.n_row - 1
                for j in range(lp):
                    virtual_fleet[agv_index] = [[start1, end1], t + j]
                    agv_index += 1
                for j in range(lp):
                    virtual_fleet[agv_index] = [[start2, end2], t + j]
                    agv_index += 1
    return virtual_fleet


def validate_deterministic_rc(network: Network, t_hat: 0.0, lp: 1):
    def agv_out_bd(network: Network, agv: AGV):
        x_max = (network.n_col - 1) * network.h_block_length
        y_max = (network.n_row - 1) * network.v_block_length
        if agv.x < 0 or agv.y < 0 or agv.x > x_max or agv.y > y_max:
            return True
        else:
            return False

    max_sim = 100
    agv_dict = {}
    start_dict = {}
    agv_index = 0
    X_record = []
    Y_record = []
    dir_record = []
    t_hat = int(t_hat)
    for t in range(max_sim):
        purge_agv_list = []
        if t % t_hat == 0:
            # release horizontal AGVs
            for k in range(network.n_row // 2):
                start1 = k * 2
                end1 = k * 2 + (network.n_col - 1) * network.n_row
                start2 = k * 2 + 1 + (network.n_col - 1) * network.n_row
                end2 = k * 2 + 1
                for j in range(lp):
                    agv_dict[agv_index] = AGV(start1, end1, network.node_list[start1].x, network.node_list[start1].y,
                                              enter_time=t + j, index=agv_index, direction=(1, 0))
                    start_dict[agv_index] = j
                    agv_index += 1
                for j in range(lp):
                    agv_dict[agv_index] = AGV(start2, end2, network.node_list[start2].x, network.node_list[start2].y,
                                              enter_time=t + j, index=agv_index, direction=(-1, 0))
                    start_dict[agv_index] = j
                    agv_index += 1
        if t / t_hat % 1 == 0.5:
            # release verticle AGVs
            for k in range(network.n_col // 2):
                start1 = k * 2 * network.n_row + network.n_row - 1
                end1 = k * 2 * network.n_row
                start2 = (k * 2 + 1) * network.n_row
                end2 = (k * 2 + 1) * network.n_row + network.n_row - 1
                for j in range(lp):
                    agv_dict[agv_index] = AGV(start1, end1, network.node_list[start1].x, network.node_list[start1].y,
                                              enter_time=t + j, index=agv_index, direction=(0, -1))
                    start_dict[agv_index] = j
                    agv_index += 1
                for j in range(lp):
                    agv_dict[agv_index] = AGV(start2, end2, network.node_list[start2].x, network.node_list[start2].y,
                                              enter_time=t + j, index=agv_index, direction=(0, 1))
                    start_dict[agv_index] = j
                    agv_index += 1
        for agv_ind, agv in agv_dict.items():
            if start_dict[agv_ind] == 0:
                agv.move()
            else:
                start_dict[agv_ind] -= 1
            if agv_out_bd(network, agv):
                purge_agv_list.append(agv_ind)
        for agv_ind in purge_agv_list:
            del agv_dict[agv_ind]
        X_record.append({agv_ind: agv.x for agv_ind, agv in agv_dict.items()})
        dir_record.append({agv_ind: agv.direction for agv_ind, agv in agv_dict.items()})
        Y_record.append({agv_ind: agv.y for agv_ind, agv in agv_dict.items()})

    return X_record, Y_record, dir_record

# create network and find all-point shortest path using Floyd Warshall
network1 = create_network(n_col, n_row, v_block_length, h_block_length)
network1.floyd_warshall()
network1.find_netflow_shortest_path()
# find RC scheme for this network
t_hat = network1.rc_prototype(lp)
virtual_fleet = rc_fleet_info(network1, t_hat, lp)
v_fleet_traj = {ind: network1.deterministic_shortest_path(od[0][0], od[0][1], time0=od[1], has_time=True) for ind, od in
                virtual_fleet.items()}
# use pandas dataframe as the RC trajectory database
temp_list = []
for v_agv_ind, v_agv_traj in v_fleet_traj.items():
    for i in range(len(v_agv_traj[0]) - 1):
        temp_list.append([v_agv_traj[0][i], v_agv_traj[0][i + 1], v_agv_traj[1][i], v_agv_traj[1][i + 1], v_agv_ind, 0])
traj_table = pd.DataFrame(temp_list, columns=['node0', 'node1', 'time0', 'time1', 'agv_ind', 'occupied'])
traj_table.sort_values('time0', inplace=True)



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
            trip_record, travel_time, delay, traj_table = new_agv.rc_shortest_path_routing(network1, traj_table)
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
    if existing_agv_table.shape[0] == 0:
        continue
    existing_agv_table['current_x'] = existing_agv_table.apply(
        lambda x: network1.arc_list[network1.node_arc_dict[(int(x.n0), int(x.n1))]].interpolate_x(x.current_loc),
        axis=1)
    existing_agv_table['current_y'] = existing_agv_table.apply(
        lambda x: network1.arc_list[network1.node_arc_dict[(int(x.n0), int(x.n1))]].interpolate_y(x.current_loc),
        axis=1)

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
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-1, x_max + 1), ylim=(-1, y_max + 1))
agvs, = ax.plot([], [], 'bo', ms=10)


def init():
    agvs.set_data([], [])
    return agvs


def animate(i):
    agvs.set_data(X_record[i], Y_record[i])
    return agvs


anim = animation.FuncAnimation(fig, animate, init_func=init, frames=400, interval=1)

print('saving...')
plt.rcParams['animation.ffmpeg_path'] = 'C:\\ffmpeg\\bin\\ffmpeg.exe'
writer = animation.FFMpegWriter()
anim.save('simulation_basic.mp4', writer=writer)
