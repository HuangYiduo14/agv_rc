import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import animation
from util import get_n_nodes, get_node_lists


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
        self.x = (node2.x + node1.x) / 2.
        self.y = (node2.y + node1.y) / 2.

    def get_node(self):
        return (self.node1, self.node2)

    def cost_func(self, flow):
        return 1. * flow ** 2


class Network:
    def __init__(self, node_list, arc_list, n_col, n_row, v_block_length, h_block_length):
        self.node_list = node_list
        self.arc_list = arc_list
        intersect_list, half_shelf_list, full_shelf_list, workstaion_list = get_node_lists(n_col,n_row)
        self.node_type_list ={}
        self.node_type_list['intersect'] = np.array(intersect_list)
        self.node_type_list['half shelf'] = np.array(half_shelf_list)
        self.node_type_list['full shelf'] = np.array(full_shelf_list)
        self.node_type_list['workstation'] = np.array(workstaion_list)
        self.n_row = n_row
        self.n_col = n_col
        self.v_block_length = v_block_length
        self.h_block_length = h_block_length
        self.dist_matrix = 99999 * np.ones((len(node_list), len(node_list)))
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

    def deterministic_shortest_path(self, o_node,d_node, time0 = 0, has_time=False):
        path=[o_node]
        time_list = [time0]
        u = o_node
        t = time0
        while u!=d_node:
            v = self.floyd_warshall_next[u,d_node]
            t += self.dist_matrix[u,v]
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
        for i in range(len(self.path)-1):
            n0 = self.path[i]
            n1 = self.path[i+1]
            print(n0,n1)
            next_v_agv_trip = traj_table.loc[(traj_table['node0']==n0)&(traj_table['node1']==n1)&(traj_table['time0']>=t)&(traj_table['occupied']==0)].iloc[0]
            traj_table.loc[next_v_agv_trip.name,'occupied']=1
            t = next_v_agv_trip['time1']
            trip_v_agv.append(int(next_v_agv_trip['agv_ind']))
        travel_time = t - self.enter_time
        delay = travel_time - optimal_travel_time
        return trip_v_agv, travel_time, delay, traj_table

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





