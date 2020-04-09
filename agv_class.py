import pandas as pd
from tw_network import *


class AGV:
    def __init__(self, start, end, x0, y0, enter_time, index, direction=(0, 0), max_v=1, size=1):
        self.size = size
        self.index = index
        self.max_v = max_v
        self.start = start
        self.end = end
        self.enter_time = enter_time
        self.x = x0
        self.y = y0
        self.direction = direction
        self.path = []

    def rc_shortest_path_routing(self, network: Network, traj_table: pd.DataFrame):
        # try to find the shortest path using rc
        # traj_table = pd.DataFrame(temp_list, columns=['time0','time1','node0','node1','agv_ind','occupied'])
        temp_list = network.deterministic_shortest_path(self.start, self.end, time0=0, has_time=True)
        self.path = temp_list[0]
        optimal_travel_time = temp_list[1][-1]
        t = self.enter_time
        trip_v_agv = []
        trip_record = []
        for i in range(len(self.path) - 1):
            n0 = self.path[i]
            n1 = self.path[i + 1]
            # print(n0,n1)
            next_v_agv_trip = traj_table.loc[
                (traj_table['node0'] == n0) & (traj_table['node1'] == n1) & (traj_table['time0'] >= t) & (
                        traj_table['occupied'] == 0)].iloc[0]
            traj_table.loc[next_v_agv_trip.name, 'occupied'] = 1
            t0 = t
            t01 = next_v_agv_trip['time0']
            t = next_v_agv_trip['time1']
            # trip_v_agv.append(int(next_v_agv_trip['agv_ind']))
            trip_record.append([n0, n1, t0, t01, t])
        travel_time = t - self.enter_time
        delay = travel_time - optimal_travel_time
        return trip_record, travel_time, delay, traj_table

    def tw_routing(self, t_network: TimeNetwork):
        # try to find shortest path
        # if there is no tw path available, wait
        enter_time = self.enter_time
        while True:
            info, travel_time, delay = t_network.time_window_routing(self.start, self.end, enter_time)
            if type(delay) != bool:
                break
            else:
                enter_time += 1
        delay += (enter_time - self.enter_time)
        self.path = info[0]
        node0 = info[1]
        trip_record = []
        for i in range(len(self.path) - 1):
            next_node = t_network.reserved_time_window[node0[0]][node0[1]].next_tw
            n0 = node0[0]
            n1 = next_node[0]
            t0 = t_network.reserved_time_window[node0[0]][node0[1]].start_time + 1
            t01 = t0
            t1 = t_network.reserved_time_window[next_node[0]][next_node[1]].start_time + 1
            trip_record.append([n0, n1, t0, t01, t1])
            node0 = next_node
        return trip_record, travel_time, delay, info[0]

    def move(self):
        self.x += self.direction[0]
        self.y += self.direction[1]
        return

    def plot_agv(self):
        plt.scatter(self.x, self.y, marker='s', linewidths=self.size)
        x_list = [self.x - self.size / 2, self.x - self.size / 2, self.x + self.size / 2, self.x + self.size / 2,
                  self.x - self.size / 2]
        y_list = [self.y - self.size / 2, self.y + self.size / 2, self.x + self.size / 2, self.x - self.size / 2,
                  self.y - self.size / 2]
        # plt.plot(x_list, y_list)
        # plt.text(self.x, self.y, 'agv ' + str(self.index))
        return x_list, y_list
