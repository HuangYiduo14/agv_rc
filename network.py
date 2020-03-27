import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from util import *

from matplotlib import animation

big_M = 99999
crossing_time = {'straight':2, 'angle':2, 'source':1}
AGV_length = 1

class Node:
    def __init__(self, x, y, index, type, shelf_station_index=-1):
        self.x = x
        self.y = y
        self.index = index
        self.is_full_shelf = (type == 'full shelf')
        self.is_half_shelf = (type == 'half shelf')
        self.is_workstation = (type == 'workstation')
        self.shelf_index = -1
        self.station_index = -1
        self.queue_length = 0
        if type == 'full shelf' or type == 'half shelf':
            self.shelf_index = shelf_station_index
        elif type == 'workstation':
            self.station_index = shelf_station_index

    @staticmethod
    def distance(node1, node2):
        dx = node1.x - node2.x
        dy = node1.y - node2.y
        return np.sqrt(dx ** 2. + dy ** 2.)



class Arc:
    def __init__(self, node1, node2):
        self.node1 = node1.index
        self.node2 = node2.index
        self.length = Node.distance(node1, node2)
        self.start_xy = [node1.x, node1.y]
        self.end_xy = [node2.x, node2.y]
        self.x = (node2.x + node1.x) / 2.
        self.y = (node2.y + node1.y) / 2.

    def get_node(self):
        return (self.node1, self.node2)

    def cost_func(self, flow):
        return 1. * flow ** 2

    def interpolate_x(self, loc):
        x = self.start_xy[0] + (self.end_xy[0] - self.start_xy[0]) * loc / self.length
        return x

    def interpolate_y(self, loc):
        y = self.start_xy[1] + (self.end_xy[1] - self.start_xy[1]) * loc / self.length
        return y


class Network:
    def __init__(self, node_list, arc_list, n_col, n_row, v_block_length, h_block_length):
        self.node_list = node_list
        self.arc_list = arc_list
        intersect_list, half_shelf_list, full_shelf_list, workstaion_list = get_node_lists(n_col, n_row)
        self.node_type_list = {}
        self.node_type_list['intersect'] = np.array(intersect_list)
        self.node_type_list['half shelf'] = np.array(half_shelf_list)
        self.node_type_list['full shelf'] = np.array(full_shelf_list)
        self.node_type_list['workstation'] = np.array(workstaion_list)
        self.node_arc_dict = {(arc.node1, arc.node2): ind for ind, arc in enumerate(arc_list)}
        self.n_row = n_row
        self.n_col = n_col
        self.v_block_length = v_block_length
        self.h_block_length = h_block_length
        self.dist_matrix = big_M * np.ones((len(node_list), len(node_list)))
        self.gen_dist_matrix()
        self.floyd_warshall_dist = np.copy(self.dist_matrix)
        self.floyd_warshall_next = np.copy(self.dist_matrix)
        self.expected_flow = dict()

    def gen_dist_matrix(self):
        for arc in self.arc_list:
            self.dist_matrix[arc.node1, arc.node2] = arc.length

    def floyd_warshall(self):
        # find all points shortest path
        n_node = len(self.node_list)
        dist = np.copy(self.dist_matrix)
        next = 99999 * np.ones((n_node, n_node), dtype=int)
        for arc in self.arc_list:
            next[arc.node1, arc.node2] = arc.node2
        for vertex in self.node_list:
            next[vertex.index, vertex.index] = 0
            next[vertex.index, vertex.index] = vertex.index
        for k in range(n_node):
            for i in range(n_node):
                for j in range(n_node):
                    if dist[i, j] > dist[i, k] + dist[k, j]:
                        dist[i, j] = dist[i, k] + dist[k, j]
                        next[i, j] = next[i, k]
        self.floyd_warshall_dist = dist
        self.floyd_warshall_next = next
        return dist, next

    def find_netflow_shortest_path(self):
        # calculate the expected network flow
        flow_expectation = {(arc.node1, arc.node2): 0 for arc in self.arc_list}
        n_intersect, n_half_shelf, n_full_shelf, n_workstation = get_n_nodes(self.n_col, self.n_row)
        # half shelf to workstation
        for i in range(n_half_shelf):
            for j in range(n_workstation):
                o_node = n_intersect + i
                d_node = n_intersect + n_half_shelf + n_full_shelf + j
                amount = 1
                self.add_path_flow(flow_expectation, o_node, d_node, amount)
                self.add_path_flow(flow_expectation, d_node, o_node, amount)
        # full shelf to workstation
        for i in range(n_full_shelf):
            for j in range(n_workstation):
                o_node = n_intersect + n_half_shelf + i
                d_node = n_intersect + n_half_shelf + n_full_shelf + j
                amount = 2
                self.add_path_flow(flow_expectation, o_node, d_node, amount)
                self.add_path_flow(flow_expectation, d_node, o_node, amount)
        # workstation to workstation
        for i in range(n_workstation):
            for j in range(n_workstation):
                if i == j:
                    continue
                o_node = n_intersect + n_half_shelf + n_full_shelf + i
                d_node = n_intersect + n_half_shelf + n_full_shelf + j
                amount = 1
                self.add_path_flow(flow_expectation, o_node, d_node, amount)
                self.add_path_flow(flow_expectation, d_node, o_node, amount)
        plt.figure()
        for arc in self.arc_list:
            x1 = self.node_list[arc.node1].x
            x2 = self.node_list[arc.node2].x
            y1 = self.node_list[arc.node1].y
            y2 = self.node_list[arc.node2].y
            plt.arrow(x1, y1, x2 - x1, y2 - y1, color='blue', width=0.1, length_includes_head=True)
            plt.plot([x1, x2], [y1, y2], linewidth=flow_expectation[(arc.node1, arc.node2)], color='red', alpha=0.3)
            plt.text(arc.x, arc.y, str(flow_expectation[(arc.node1, arc.node2)]))
        plt.show()
        self.expected_flow = flow_expectation
        return flow_expectation

    def add_path_flow(self, flow_expectation, o_node, d_node, amount):
        # assign shortest path flow to arcs from o_node to d_node
        u = o_node
        while u != d_node:
            v = self.floyd_warshall_next[u, d_node]
            flow_expectation[(u, v)] += amount
            pre_u = u
            u = v
        return flow_expectation

    def deterministic_shortest_path(self, o_node, d_node, time0=0, has_time=False):
        path = [o_node]
        time_list = [time0]
        u = o_node
        t = time0
        while u != d_node:
            v = self.floyd_warshall_next[u, d_node]
            t += self.dist_matrix[u, v]
            u = v
            path.append(u)
            time_list.append(t)
        if has_time:
            return [path, time_list]
        else:
            return path

    def rc_prototype(self, lp):
        t_hat = 2 * (lp + 1)
        if self.v_block_length % t_hat != 0 or self.h_block_length % t_hat != 0:
            print('horizontal/vertical block size cannot be divided by t_hat!')
            return 0.0
        return t_hat

class ReservedTimeWindow:
    # tw_type in {first, last, mid}
    def __init__(self, node_i, start_time, end_time, tw_index, tw_type='first', prev_tw= [-big_M,-big_M], next_tw = [-big_M,-big_M]):
        self.node_i = node_i
        self.start_time = start_time
        self.end_time = end_time
        self.tw_index = (node_i, tw_index) # index is like (node_i, p)
        self.tw_type = tw_type
        self.prev_tw = prev_tw
        self.next_tw = next_tw




class TimeNetwork(Network):
    def __init__(self, network: Network, time_span = 1000, bidirect = True):
        Network.__init__(network.node_list, network.arc_list, network.n_col, network.n_row, network.v_block_length,
                         network.h_block_length)
        # we use free time window here
        self.free_time_window = [[[0,time_span]] for node in network.node_list]
        # reserve time window
        self.reserved_time_window = [[ReservedTimeWindow()] for node in network.node_list]
        # update the distance matrix if the network is bidirectional
        if bidirect:
            valid_ind = np.where(self.dist_matrix < big_M//2)
            self.dist_matrix[valid_ind[1], valid_ind[0]] = self.dist_matrix[valid_ind[0], valid_ind[1]]
            self.arc_list += [Arc(self.node_list[arc.node2], self.node_list[arc.node1]) for arc in network.arc_list]
            self.gen_dist_matrix()
            self.node_arc_dict = {(arc.node1, arc.node2): ind for ind, arc in enumerate(self.arc_list)}

    def time_calculate(self, node0, node1, label0, crossing_type, time_window0, time_window1):
        # return tau, lambda, alpha
        def reverse_test(n0, n1, tw0, tw1):
            for tw in self.reserved_time_window[n0]:
                if tw.start_time > tw0[1] and tw.prev_tw[0]== n1:
                    if self.reserved_time_window[tw.prev_tw[0]][tw.prev_tw[1]].end_time< tw1[0]:
                        return True
            return False
        arc_length = self.dist_matrix[node0, node1] - AGV_length
        # check time feasibility
        if label0 + crossing_time[crossing_type]> time_window0[1]:
            return big_M, -big_M, big_M
        alpha_ji = label0 + crossing_time[crossing_type] + arc_length
        # reachable test
        if node1 != node0:
            tau_ji = max(time_window1[0], alpha_ji)
            if arc_length > big_M//2:
                return big_M, (-big_M,-big_M), big_M
            # revese test
            if reverse_test(node0, node1, time_window0, time_window1):
                return big_M, (-big_M,-big_M), big_M
            if reverse_test(node1, node0, time_window1, time_window0):
                return big_M, (-big_M,-big_M), big_M
            lambda_ji = (node0, node1)
            return tau_ji, lambda_ji, alpha_ji
        else:
            tau_ji = time_window1[0]
            available_j = np.where(self.dist_matrix[node0]<big_M//2)[0]
            if available_j.shape[0]==0:
                return big_M, (-big_M,-big_M), big_M
            available_j = set(available_j)
            for tw in self.reserved_time_window[node0]:
                if tw.start_time > time_window0[1] and tw.end_time<time_window1[0]:
                    n0 = tw.prev_tw[0]
                    n1 = tw.next_tw[0]
                    available_j = available_j- {n0,n1}
            if len(available_j)==0:
                return big_M, (-big_M,-big_M), big_M
            lambda_ji = (node0, available_j[0])
            return tau_ji, lambda_ji, alpha_ji


    def time_window_routing(self, o, d, enter_time):
        # initialize
        label_list = [[big_M for f_tw in f_tw_node] for f_tw_node in self.free_time_window]
        prev_tw_list = [[(-big_M, -big_M) for f_tw in f_tw_node] for f_tw_node in self.free_time_window]
        prev_lane_list = [[(-big_M, -big_M) for f_tw in f_tw_node] for f_tw_node in self.free_time_window]
        set_F = set()
        dict_U = dict()
        set_T = set()
        for i in range(len(self.node_list)):
            for j in range(len(self.free_time_window[i])):
                set_F.add((i,j))

        for f_tw_ind, f_tw in enumerate(self.free_time_window[o]):
            if f_tw[0] <= enter_time and f_tw[1]>=enter_time:
                label_list[o][f_tw_ind] = enter_time
                prev_lane_list[o][f_tw_ind] = (o, o)
                prev_tw_list[o][f_tw_ind] = (o, f_tw_ind)
                set_U.add((o, f_tw_ind))
                set_F -= set_T



        return trip_record, travel_time, delay



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


def create_network(n_col, n_row, v_block_length, h_block_length=3):
    """
    :param n_col: number of vertical streets
    :param n_row: number of horizontal streets
    :param v_block_length: length of vertical blocks
    :param h_block_length: length of horizontal blocks
    :return: Network object
    """
    index_node = 0
    index_shelf = 0
    index_workstation = 0
    node_list = []
    arc_list = []
    n_intersect, n_half_shelf, n_full_shelf, n_workstation = get_n_nodes(n_col, n_row)

    # step 1. create intersection nodes
    x = h_block_length * np.arange(n_col)
    y = v_block_length * np.arange(n_row)
    X, Y = np.meshgrid(x, y, sparse=False, indexing='ij')
    for i in range(n_col):
        for j in range(n_row):
            node_list.append(Node(X[i, j], Y[i, j], index_node, 'intersection'))
            plt.scatter(X[i, j], Y[i, j])
            plt.text(X[i, j], Y[i, j], str(index_node))
            index_node += 1

    # step 2. create half shelf nodes
    x = h_block_length * np.array([0., n_col - 1.])
    y = v_block_length * (np.arange(n_row - 1.) + 0.5)
    X, Y = np.meshgrid(x, y, sparse=False, indexing='ij')
    for i in range(len(x)):
        for j in range(len(y)):
            node_list.append(Node(X[i, j], Y[i, j], index_node, 'half shelf', index_shelf))
            plt.scatter(X[i, j], Y[i, j], marker='*')
            plt.text(X[i, j], Y[i, j], index_node)
            index_node += 1
            index_shelf += 1

    # step 3. create full shelf nodes
    x = h_block_length * np.arange(1, n_row - 1)
    y = v_block_length * (np.arange(n_row - 1.) + 0.5)
    X, Y = np.meshgrid(x, y, sparse=False, indexing='ij')
    for i in range(len(x)):
        for j in range(len(y)):
            node_list.append(Node(X[i, j], Y[i, j], index_node, 'full shelf', index_shelf))
            plt.scatter(X[i, j], Y[i, j], marker='*')
            plt.text(X[i, j], Y[i, j], index_node)
            index_node += 1
            index_shelf += 1

    # step 4. create workstation nodes
    y = 0.
    x = h_block_length * (np.arange(n_col - 1) + 0.5)
    for i in range(len(x)):
        node_list.append(Node(x[i], y, index_node, 'workstation', index_workstation))
        plt.scatter(x[i], y, marker='^')
        plt.text(x[i], y, index_node)
        index_node += 1
        index_workstation += 1

    # step 5. create horizontal arcs
    # for the bottom
    for i in range(n_col - 1):
        mid_index = n_intersect + n_half_shelf + n_full_shelf + i
        arc_list.append(Arc(node_list[i * n_row], node_list[mid_index]))
        plt.plot([node_list[i * n_row].x, node_list[mid_index].x], [node_list[i * n_row].y, node_list[mid_index].y])
        arc_list.append(Arc(node_list[mid_index], node_list[(i + 1) * n_row]))
        plt.plot([node_list[mid_index].x, node_list[(i + 1) * n_row].x],
                 [node_list[mid_index].y, node_list[(i + 1) * n_row].y])
    # other arcs
    for i in range(n_col - 1):
        for j in range(n_row - 1):
            direction = (j + 1) % 2
            ind1 = i * n_row + j + 1
            ind2 = (i + 1) * n_row + j + 1
            if direction == 0:
                arc_list.append(Arc(node_list[ind1], node_list[ind2]))
            else:
                arc_list.append(Arc(node_list[ind2], node_list[ind1]))
            plt.plot([node_list[ind1].x, node_list[ind2].x], [node_list[ind1].y, node_list[ind2].y])

    # step 6. create vertical arcs
    for i in range(n_col):
        for j in range(n_row - 1):
            direction = i % 2
            ind1 = i * n_row + j
            ind2 = ind1 + 1
            if i == 0:
                mid_index = n_intersect + j
            elif i == n_col - 1:
                mid_index = n_intersect + n_row - 1 + j
            else:
                mid_index = n_intersect + n_half_shelf + (i - 1) * (n_row - 1) + j
            if direction == 0:
                arc_list.append(Arc(node_list[ind2], node_list[mid_index]))
                arc_list.append(Arc(node_list[mid_index], node_list[ind1]))
            else:
                arc_list.append(Arc(node_list[ind1], node_list[mid_index]))
                arc_list.append(Arc(node_list[mid_index], node_list[ind2]))
            plt.plot([node_list[ind1].x, node_list[mid_index].x], [node_list[ind1].y, node_list[mid_index].y])
            plt.plot([node_list[mid_index].x, node_list[ind2].x], [node_list[mid_index].y, node_list[ind2].y])
    plt.show()

    network = Network(node_list, arc_list, n_col, n_row, v_block_length, h_block_length)
    return network
