from network import *
from agv_class import *
from matplotlib import animation
import pandas as pd

def rc_fleet_info(network, t_hat, lp: 1):
    '''
    create virtual fleet
    :param network: Network
    :param t_hat: rhythme cycle
    :param lp: platoon size
    :return: {v agv index: [start node index, end node index, enter time]}
    '''
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
    # validate deterministic rc and create virtual fleet trajectories
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


def warehouse_simulation(network1, routing_type='rc', traj_table=None, save_anim=False, ani_name='simulation_basic',p=prob):
    # simulation experiment


    if routing_type=='tw_rc_heuristic_gen':
        demand_time = deterministic_demand()
    else:
        demand_time = simulation_demand(p)
    print('start experiment')
    total_travel_time = 0
    total_delay = 0
    travel_time_record = []
    delay_record = []
    trip_v_agv_record = []
    X_record = []
    Y_record = []
    agv_index = 0
    all_trip = pd.DataFrame(columns=['n0', 'n1', 't0', 't01', 't1', 'agv_ind', 'current_loc', 'arc_length'])
    for t in range(100, 100 + n_t):
        print('time:', t)
        # generate trips
        del_index = []
        from_node = np.array([])
        to_node = np.array([])
        # get all od demand at time t
        for od in possible_od:
            from_node_i, to_node_i = demand_transform(demand_time[od][t - 100], od[0], od[1], network1)
            from_node = np.append(from_node, from_node_i)
            to_node = np.append(to_node, to_node_i)
        #print('fron_node',from_node)
        #print('to_node',to_node)
        if routing_type == 'tw_rc_heuristic_gen':
            # if we are doing time window shortest path heuristic rc, we assign high priority to long travel distance
            distance = np.array([network1.floyd_warshall_dist[int(from_node[i]), int(to_node[i])] for i in range(len(from_node))])
            sorted_index = np.argsort(-distance)
            from_node = from_node[sorted_index]
            to_node = to_node[sorted_index]

        for i in range(from_node.shape[0]):
            if i in del_index:
                # if we have deleted this demand because it is covered by other trips in tw_rc_heuristic_gen
                continue
            start = int(from_node[i])
            end = int(to_node[i])
            if start == end:
                continue
            new_agv = AGV(start, end, network1.node_list[start].x, network1.node_list[start].y, t, agv_index)
            if routing_type == 'rc':
                trip_record, travel_time, delay, traj_table = new_agv.rc_shortest_path_routing(network1, traj_table)
            else:
                trip_record, travel_time, delay, traj = new_agv.tw_routing(network1)

            total_travel_time += travel_time
            total_delay += delay
            travel_time_record.append(travel_time)
            delay_record.append(delay)
            trip_record_i = pd.DataFrame(trip_record, columns=['n0', 'n1', 't0', 't01', 't1'])
            trip_record_i['agv_ind'] = agv_index
            trip_record_i['current_loc'] = -1
            trip_record_i.loc[0, 'current_loc'] = 0

            trip_record_i['arc_length'] = trip_record_i.apply(
                lambda x: network1.arc_list[network1.node_arc_dict[(int(x.n0), int(x.n1))]].length, axis=1)
            agv_index += 1
            all_trip = all_trip.append(trip_record_i, ignore_index=True)

            # if we are trying to find rc using time window based method, we need to update od demand
            if routing_type == 'tw_rc_heuristic_gen':
                covered_od_list = network1.find_covered_od_demand(traj)
                #  [node_o, node_d, type_o, type_d, index of o in node_type_list[type_o], index of d in node_type_list[type_d]]
                # update demand in the demand list at time t
                valid_cover_index = []
                for i in range(len(from_node)):
                    for j in range(len(covered_od_list)):
                        if from_node[i] == covered_od_list[j][0] and to_node[i] == covered_od_list[j][1]:
                            del_index.append(i)
                        else:
                            valid_cover_index.append(j)
                # update future demand
                for i in valid_cover_index:
                    t_cycle = (t - 100) % cycle_length
                    if t_cycle == cycle_length - 1:
                        break
                    type_o = covered_od_list[i][2]
                    type_d = covered_od_list[i][3]
                    ind_o_type = covered_od_list[i][4]
                    ind_d_type = covered_od_list[i][5]
                    for t1 in range(1, cycle_length-t_cycle):
                        if demand_time[(type_o,type_d)][t1+t][ind_o_type,ind_d_type]:
                            demand_time[(type_o,type_d)][t1+t][ind_o_type,ind_d_type]=False
                            break

        # keep record of locations
        # print(all_trip)
        if save_anim:
            existing_agv_table = all_trip.loc[(all_trip['current_loc'] >= 0) & (all_trip['current_loc'] < 100)].copy()
            #print('existing_table:')
            #print(existing_agv_table)

            if existing_agv_table.shape[0] == 0:
                continue
            existing_agv_table['current_x'] = existing_agv_table.apply(
                lambda x: network1.arc_list[network1.node_arc_dict[(int(x.n0), int(x.n1))]].interpolate_x(
                    x.current_loc),
                axis=1)
            existing_agv_table['current_y'] = existing_agv_table.apply(
                lambda x: network1.arc_list[network1.node_arc_dict[(int(x.n0), int(x.n1))]].interpolate_y(
                    x.current_loc),
                axis=1)

            X_record.append(existing_agv_table['current_x'].values)
            Y_record.append(existing_agv_table['current_y'].values)
            # move to next
            all_trip.loc[(all_trip['t0'] <= t) & (all_trip['t01'] >= t), 'current_loc'] = 0
            all_trip.loc[all_trip['t1'] <= t, 'current_loc'] = 999
            if routing_type == 'rc':
                all_trip.loc[
                    (all_trip['current_loc'] >= 0) & (all_trip['current_loc'] < 100) & (all_trip['t01'] <= t) & (
                            all_trip['t1'] >= t), 'current_loc'] += 1.
            else:
                index_of_interest = (all_trip['t01'] <= t) & (all_trip['t1'] >= t)
                all_trip.loc[index_of_interest, 'current_loc'] = 1. * (t - all_trip.loc[index_of_interest, 't01']) / (
                        all_trip.loc[index_of_interest, 't1'] - all_trip.loc[index_of_interest, 't01'])
                all_trip.loc[index_of_interest, 'current_loc'] = all_trip.loc[index_of_interest, 'current_loc'] * \
                                                                 all_trip.loc[index_of_interest, 'arc_length']

    plt.figure()
    plt.plot(delay_record)
    plt.xlabel('agv_index')
    plt.ylabel('delay')
    plt.title('{0},p={1},t={2}'.format(routing_type,p,n_t))
    plt.savefig('expfig/{0},p={1},t={2}.png'.format(routing_type,p,n_t))
    print('prob',p)
    print('total delay', total_delay)
    print('total travel time', total_travel_time)
    print('total agvs', agv_index)


    def init():
        agvs.set_data([], [])
        return agvs
    def animate(i):
        agvs.set_data(X_record[i], Y_record[i])
        return agvs
    if save_anim:
        # animation
        x_max = (network1.n_col - 1) * network1.h_block_length
        y_max = (network1.n_row - 1) * network1.v_block_length
        fig = plt.figure()
        ax = fig.add_subplot(111, autoscale_on=False, xlim=(-1, x_max + 1), ylim=(-1, y_max + 1))
        agvs, = ax.plot([], [], 'bo', ms=10)
        anim = animation.FuncAnimation(fig, animate, init_func=init, frames=n_t - 50, interval=1)
        print('saving...')
        plt.rcParams['animation.ffmpeg_path'] = 'C:\\ffmpeg\\bin\\ffmpeg.exe'
        writer = animation.FFMpegWriter()
        anim.save('{0}.mp4'.format(ani_name), writer=writer)
    return total_delay, total_travel_time, agv_index


if __name__ == '__main__':
    # create network and find all-point shortest path using Floyd Warshall
    network1 = create_network(n_col, n_row, v_block_length, h_block_length)
    network1.floyd_warshall()
    network1.find_netflow_shortest_path()
    # find RC scheme for this network
    t_hat = network1.rc_prototype(lp)
    virtual_fleet = rc_fleet_info(network1, t_hat, lp)
    v_fleet_traj = {ind: network1.deterministic_shortest_path(od[0][0], od[0][1], time0=od[1], has_time=True) for
                    ind, od in
                    virtual_fleet.items()}
    # use pandas dataframe as the RC trajectory database
    temp_list = []
    for v_agv_ind, v_agv_traj in v_fleet_traj.items():
        for i in range(len(v_agv_traj[0]) - 1):
            temp_list.append(
                [v_agv_traj[0][i], v_agv_traj[0][i + 1], v_agv_traj[1][i], v_agv_traj[1][i + 1], v_agv_ind, 0])
    traj_table = pd.DataFrame(temp_list, columns=['node0', 'node1', 'time0', 'time1', 'agv_ind', 'occupied'])
    traj_table.sort_values('time0', inplace=True)
    # do simulation
    simulation_record = pd.DataFrame(columns=['prob','total_delay','total_travel_time','agv_number'])
    #p_list = [0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02]
    #p_list = [0.013, 0.0135, 0.014, 0.0145, 0.015]
    p_list = [0.003,0.0035, 0.004]
    for i, p in enumerate(p_list):
        print('rc routing', p, '=' * 50)
        temp_table = traj_table.copy()
        total_delay, total_travel_time, agv_index= warehouse_simulation(network1, routing_type='rc',traj_table=temp_table, p=p)
        simulation_result = pd.DataFrame({'prob':p,'total_delay':[total_delay], 'total_travel_time':[total_travel_time],'agv_number':[agv_index]})
        simulation_record = simulation_record.append(simulation_result, ignore_index=True)
    simulation_record.to_csv('expdata/exp_rc2.csv')
'''
result
total delay 8701.0
total travel time 32809.0
total agvs 565
saving...
'''