from network import *


class ReservedTimeWindow:
    # tw_type in {first, last, mid}
    def __init__(self, node_i, start_time, end_time, tw_type='source', prev_tw=(-big_M, -big_M),
                 next_tw=(-big_M, -big_M)):
        self.node_i = node_i
        self.start_time = start_time
        self.end_time = end_time
        self.tw_type = tw_type
        self.prev_tw = prev_tw
        self.next_tw = next_tw


class TimeNetwork(Network):
    def __init__(self, network: Network, time_span=1000, bidirect=True):
        Network.__init__(self, network.node_list, network.arc_list, network.n_col, network.n_row,
                         network.v_block_length,
                         network.h_block_length)
        # we use free time window here
        self.free_time_window = [[[0, time_span]] for node in network.node_list]
        # reserve time window
        self.reserved_time_window = [[] for node in network.node_list]
        # update the distance matrix if the network is bidirectional
        if bidirect:
            valid_ind = np.where(self.dist_matrix < big_M // 2)
            self.dist_matrix[valid_ind[1], valid_ind[0]] = self.dist_matrix[valid_ind[0], valid_ind[1]]
            self.arc_list += [Arc(self.node_list[arc.node2], self.node_list[arc.node1]) for arc in network.arc_list]
            self.gen_dist_matrix()
            self.node_arc_dict = {(arc.node1, arc.node2): ind for ind, arc in enumerate(self.arc_list)}
        self.floyd_warshall()

    def time_calculate(self, node0, node1, label0, crossing_type, time_window0, time_window1):
        # return tau, lambda, alpha
        # lambda is the next lane, tau is the new label, alpha is the earliest arrival time

        def reverse_test(n0, n1, tw0, tw1):
            for tw in self.reserved_time_window[n0]:
                if tw.start_time >= tw0[1] and tw.prev_tw[0] == n1:
                    if self.reserved_time_window[tw.prev_tw[0]][tw.prev_tw[1]].end_time <= tw1[0]:
                        return True
            return False

        # check time feasibility
        if label0 + crossing_time[crossing_type] > time_window0[1]:
            return big_M, -big_M, big_M
        arc_length = self.dist_matrix[node0, node1] - AGV_length
        alpha_ji = label0 + crossing_time[crossing_type] + arc_length - 1
        # reachable test
        if node1 != node0:
            tau_ji = max(time_window1[0], alpha_ji)
            if arc_length > big_M // 2:
                return big_M, (-big_M, -big_M), big_M
            # revese test
            if reverse_test(node0, node1, time_window0, time_window1):
                return big_M, (-big_M, -big_M), big_M
            if reverse_test(node1, node0, time_window1, time_window0):
                return big_M, (-big_M, -big_M), big_M
            lambda_ji = (node0, node1)
            return tau_ji, lambda_ji, alpha_ji
        else:
            tau_ji = time_window1[0]
            available_j = np.where(self.dist_matrix[node0] < big_M // 2)[0]
            if available_j.shape[0] == 0:
                return big_M, (-big_M, -big_M), big_M
            available_j = set(available_j)
            for tw in self.reserved_time_window[node0]:
                if tw.start_time > time_window0[1] and tw.end_time < time_window1[0]:
                    n0 = tw.prev_tw[0]
                    n1 = tw.next_tw[0]
                    available_j = available_j - {n0, n1}
            if len(available_j) == 0:
                return big_M, (-big_M, -big_M), big_M
            lambda_ji = (node0, available_j[0])

            return tau_ji, lambda_ji, alpha_ji

    def update_tw(self, prev_tw_list, prev_lane_list, dist_list, o: (1, 1), d: (1, 1), enter_time):
        # t0: start time of this trip, t01: when agv start moving, t: when agv arrives
        # trip_record.append([n0, n1, t0, t01, t])
        # return trip_record, travel_time, delay
        # index of tw: (ind_node,ind_tw), e.g. o[0] is the node index of o, o[1] is the tw index of o
        # reconstruct the path
        u = d
        path = [u]
        while (u[0] != o[0]) or (u[1] != o[1]):
            v = prev_tw_list[u[0]][u[1]]
            u = v
            path.append(u)
        path.reverse()
        travel_time = dist_list[d[0]][d[1]] - enter_time
        # calculate delay
        temp_list = self.deterministic_shortest_path(o[0], d[0], time0=0, has_time=True)
        optimal_travel_time = temp_list[1][-1]
        delay = travel_time - optimal_travel_time
        traj = [u[0] for u in path]
        # update time windows
        prev_r_tw_ind = (-big_M, -big_M)
        for i in range(len(path)):
            u = path[i]  # here u is the free tw node index (node_ind, f_tw_ind)
            start_time = dist_list[u[0]][u[1]]
            if i == len(path) - 1:
                crossing_type = 'end'
            elif i == 0:
                crossing_type = 'source'
            else:
                # determine crossing type
                prev_lane = prev_lane_list[u[0]][u[1]]
                prev_node = self.get_prev_node(prev_lane, u[0])
                crossing_type = get_crossing_type(prev_node, u[0], path[i + 1][0])
            end_time = start_time + crossing_time[crossing_type]
            # update reserved tw list
            self.reserved_time_window[u[0]].append(
                ReservedTimeWindow(u[0], start_time, end_time, tw_type=crossing_type, prev_tw=prev_r_tw_ind))
            this_r_tw_ind = (u[0], len(self.reserved_time_window[u[0]]) - 1)  # keep record of the tw index
            # print(this_r_tw_ind)
            if i >= 1:
                self.reserved_time_window[prev_r_tw_ind[0]][prev_r_tw_ind[1]].next_tw = this_r_tw_ind
            prev_r_tw_ind = this_r_tw_ind
            # update free tw list
            old_ftw = self.free_time_window[u[0]][u[1]]
            if old_ftw[1] - end_time > 1 and start_time - old_ftw[0] > 1:
                self.free_time_window[u[0]][u[1]] = [end_time, old_ftw[1]]
                self.free_time_window[u[0]].insert(u[1], [old_ftw[0], start_time])
            elif old_ftw[1] - end_time > 1 and start_time - old_ftw[0] <= 1:
                self.free_time_window[u[0]][u[1]] = [end_time, old_ftw[1]]
            elif old_ftw[1] - end_time > 1 and start_time - old_ftw[0] <= 1:
                self.free_time_window[u[0]][u[1]] = [old_ftw[0], start_time]
            else:
                self.free_time_window[u[0]].pop(u[1])
        return traj, travel_time, delay

    def get_prev_node(self, prev_lane, node_ind):
        # find previous node
        # this function can help us find crossing type
        if prev_lane[0] == node_ind:
            prev_node = prev_lane[1]
        else:
            prev_node = prev_lane[0]
        return prev_node

    def time_window_routing(self, o, d, enter_time):
        # Shortest path algorithm for time window network
        # sort free time window
        for node_ind in range(len(self.node_list)):
            self.free_time_window[node_ind].sort(key=lambda x: x[0])
        # Dijkstra's algorithm initialization the first node
        # find an available time window for the first node
        o_tw_ind = -big_M
        o_tw_t0 = big_M
        for f_tw_ind, f_tw in enumerate(self.free_time_window[o]):
            if f_tw[1] >= enter_time + crossing_time['source']:  # if the enter time is valid
                o_tw_ind = f_tw_ind
                o_tw_t0 = f_tw[0]
                break
        if o_tw_t0 >= big_M:  # if there is no time window possible, then reject this demand
            return False, False, False
        actual_enter_time = max(enter_time, o_tw_t0)
        # initialize list
        prev_tw_list = [[(-big_M, -big_M) for f_tw in f_tw_node] for f_tw_node in
                        self.free_time_window]  # (index_node, index_time_window)
        prev_lane_list = [[(-big_M, -big_M) for f_tw in f_tw_node] for f_tw_node in self.free_time_window]
        dist_list = [[big_M for f_tw in f_tw_node] for f_tw_node in self.free_time_window]
        pq_q = Heapdict()
        for node_ind, f_tw_node in enumerate(self.free_time_window):
            for f_tw_ind, f_tw in enumerate(f_tw_node):
                if f_tw[1] > actual_enter_time:
                    pq_q[(node_ind, f_tw_ind)] = big_M
        # initialize for the first node
        prev_lane_list[o][o_tw_ind] = (o, o)
        prev_tw_list[o][o_tw_ind] = (o, o_tw_ind)
        dist_list[o][o_tw_ind] = actual_enter_time
        pq_q[(o, o_tw_ind)] = actual_enter_time
        # apply Dijkstra algorithm
        for loop in range(10000):
            if len(pq_q) == 0:
                return False, False, False
            u, label_u = pq_q.popitem()
            node_ind = u[0]
            tw_ind = u[1]
            tw0 = self.free_time_window[node_ind][tw_ind]
            prev_lane = prev_lane_list[node_ind][tw_ind]
            prev_node = self.get_prev_node(prev_lane, node_ind)
            # if we found the target node
            if u[0] == d:
                traj, travel_time, delay = self.update_tw(prev_tw_list, prev_lane_list, dist_list, (o, o_tw_ind), u,
                                                          enter_time)
                return traj, travel_time, delay
            neighbor_u_node = np.where(self.dist_matrix[node_ind] < big_M // 2)
            for node_ind1 in neighbor_u_node[0]:
                for tw_ind1, tw1 in enumerate(self.free_time_window[node_ind1]):
                    crossing_type = get_crossing_type(prev_node, node_ind, node_ind1)
                    tau_ji, lambda_ji, alpha_ji = self.time_calculate(node_ind, node_ind1, label_u, crossing_type, tw0,
                                                                      tw1)
                    if tau_ji < dist_list[node_ind1][tw_ind1] and tau_ji < big_M:
                        dist_list[node_ind1][tw_ind1] = tau_ji
                        prev_lane_list[node_ind1][tw_ind1] = lambda_ji
                        prev_tw_list[node_ind1][tw_ind1] = (node_ind, tw_ind)
                        pq_q[(node_ind1, tw_ind1)] = tau_ji
        return False, False, False
