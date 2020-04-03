import numpy as np


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
