import numpy as np

# network parameters
h_block_length = 6
v_block_length = 12
n_col = 4
n_row = 4
lp = 2
# simulation parameters
prob = 0.01  # intensity of flow
n_t = 1000  # simulation time
alpha = 1./4.  # relative demand from station to station
cycle_length = 10 # length of cycle
possible_od = [('half shelf', 'workstation'), ('full shelf', 'workstation'), ('workstation', 'half shelf'),
               ('workstation', 'full shelf'), ('workstation', 'workstation')]

def deterministic_demand():
    # we assume cycle length > 4
    n_cycles = n_t // cycle_length
    n_intersect, n_half_shelf, n_full_shelf, n_workstation = get_n_nodes(n_col, n_row)
    print('setting task list')
    half_ws = np.zeros((cycle_length, n_half_shelf, n_workstation))
    half_ws[:2, :, :] = 1
    full_ws = np.zeros((cycle_length, n_full_shelf, n_workstation))
    full_ws[:4, :, :] = 1
    ws_half = np.zeros((cycle_length, n_workstation, n_half_shelf))
    ws_half[:2, :, :] = 1
    ws_full = np.zeros((cycle_length, n_workstation, n_full_shelf))
    ws_full[:4, :, :] = 1
    ws_ws = np.zeros((cycle_length, n_workstation, n_workstation))
    ws_ws[:1, :, :] = 1
    demand_time = dict()
    demand_time[('half shelf', 'workstation')] = np.tile(half_ws, (n_cycles, 1, 1)) > 0.5
    demand_time[('full shelf', 'workstation')] = np.tile(full_ws, (n_cycles, 1, 1)) > 0.5
    demand_time[('workstation', 'half shelf')] = np.tile(ws_half, (n_cycles, 1, 1)) > 0.5
    demand_time[('workstation', 'full shelf')] = np.tile(ws_full, (n_cycles, 1, 1)) > 0.5
    demand_time[('workstation', 'workstation')] = np.tile(ws_ws, (n_cycles, 1, 1)) > 0.5
    return demand_time


def simulation_demand(p=prob):
    '''
    pre generate demand
    :return: demand_time[(possible od type)][time=t][o,d] = True/False
    '''
    np.random.seed(10)
    n_intersect, n_half_shelf, n_full_shelf, n_workstation = get_n_nodes(n_col, n_row)
    print('start simulation:')
    print('generating random matrices...')
    demand_time = dict()
    demand_time[('half shelf', 'workstation')] = np.random.rand(n_t, n_half_shelf, n_workstation) < (p / 2)
    demand_time[('full shelf', 'workstation')] = np.random.rand(n_t, n_full_shelf, n_workstation) < p
    demand_time[('workstation', 'half shelf')] = np.random.rand(n_t, n_workstation, n_half_shelf) < (p / 2)
    demand_time[('workstation', 'full shelf')] = np.random.rand(n_t, n_workstation, n_full_shelf) < p
    demand_time[('workstation', 'workstation')] = np.random.rand(n_t, n_workstation, n_workstation) < (p * alpha)
    return demand_time


def demand_transform(demand_list, type1: str, type2: str, network):
    '''
    :param demand_list: od matrix
    :param type1: type of node, e.g. half shelf
    :param type2: same
    :param network: Network
    :return: list of o and d
    '''
    valid_agv = np.where(demand_list)
    from_node_list = network.node_type_list[type1][valid_agv[0]]
    to_node_list = network.node_type_list[type2][valid_agv[1]]
    return np.array(from_node_list), np.array(to_node_list)


def get_n_nodes(n_col, n_row):
    # get the number of different types of nodes
    n_intersect = n_col * n_row
    n_half_shelf = 2 * (n_row - 1)
    n_full_shelf = (n_col - 2) * (n_row - 1)
    n_workstation = n_col - 1
    return n_intersect, n_half_shelf, n_full_shelf, n_workstation


def get_node_lists(n_col, n_row):
    # get the index of different types of nodes
    n_intersect, n_half_shelf, n_full_shelf, n_workstation = get_n_nodes(n_col, n_row)
    intersect_list = list(range(n_intersect))
    half_shelf_list = list(range(n_intersect, n_intersect + n_half_shelf))
    full_shelf_list = list(range(n_intersect + n_half_shelf, n_intersect + n_half_shelf + n_full_shelf))
    workstaion_list = list(
        range(n_intersect + n_half_shelf + n_full_shelf, n_intersect + n_half_shelf + n_full_shelf + n_workstation))
    return intersect_list, half_shelf_list, full_shelf_list, workstaion_list


def get_crossing_type(n0, n1, n2):
    return 'straight'


# the following code is from heapdict package
try:
    from collections.abc import MutableMapping
except ImportError:
    from collections import MutableMapping


def doc(s):
    if hasattr(s, '__call__'):
        s = s.__doc__

    def f(g):
        g.__doc__ = s
        return g

    return f


class Heapdict(MutableMapping):
    __marker = object()

    def __init__(self, *args, **kw):
        self.heap = []
        self.d = {}
        self.update(*args, **kw)

    @doc(dict.clear)
    def clear(self):
        del self.heap[:]
        self.d.clear()

    @doc(dict.__setitem__)
    def __setitem__(self, key, value):
        if key in self.d:
            self.pop(key)
        wrapper = [value, key, len(self)]
        self.d[key] = wrapper
        self.heap.append(wrapper)
        self._decrease_key(len(self.heap) - 1)

    def _min_heapify(self, i):
        n = len(self.heap)
        h = self.heap
        while True:
            # calculate the offset of the left child
            l = (i << 1) + 1
            # calculate the offset of the right child
            r = (i + 1) << 1
            if l < n and h[l][0] < h[i][0]:
                low = l
            else:
                low = i
            if r < n and h[r][0] < h[low][0]:
                low = r

            if low == i:
                break

            self._swap(i, low)
            i = low

    def _decrease_key(self, i):
        while i:
            # calculate the offset of the parent
            parent = (i - 1) >> 1
            if self.heap[parent][0] < self.heap[i][0]:
                break
            self._swap(i, parent)
            i = parent

    def _swap(self, i, j):
        h = self.heap
        h[i], h[j] = h[j], h[i]
        h[i][2] = i
        h[j][2] = j

    @doc(dict.__delitem__)
    def __delitem__(self, key):
        wrapper = self.d[key]
        while wrapper[2]:
            # calculate the offset of the parent
            parentpos = (wrapper[2] - 1) >> 1
            parent = self.heap[parentpos]
            self._swap(wrapper[2], parent[2])
        self.popitem()

    @doc(dict.__getitem__)
    def __getitem__(self, key):
        return self.d[key][0]

    @doc(dict.__iter__)
    def __iter__(self):
        return iter(self.d)

    def popitem(self):
        """D.popitem() -> (k, v), remove and return the (key, value) pair with lowest\nvalue; but raise KeyError if D is empty."""
        wrapper = self.heap[0]
        if len(self.heap) == 1:
            self.heap.pop()
        else:
            self.heap[0] = self.heap.pop()
            self.heap[0][2] = 0
            self._min_heapify(0)
        del self.d[wrapper[1]]
        return wrapper[1], wrapper[0]

    @doc(dict.__len__)
    def __len__(self):
        return len(self.d)

    def peekitem(self):
        """D.peekitem() -> (k, v), return the (key, value) pair with lowest value;\n but raise KeyError if D is empty."""
        return (self.heap[0][1], self.heap[0][0])


