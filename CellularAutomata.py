"""
Author: Justin Deterding
Created: Tue Mar 10 15:30:20 2020
Description:
    CellularAuntonima provides a simple class for building a CA on top of a 
    networkX graph. This file also served as a location to store functions and
    usefull in a number of the scripts used to generate plots.
"""
# %% IMPORTS
import numpy as np
import networkx as nx
from networkx import Graph


# %% FUNCTION DEFINITIONS
def deterministic_SIR_rule(neighboor_states, current_state, thresh=1):
    if current_state == 'I' or current_state == 'R':
        return 'R'
    elif neighboor_states.count('I') >= thresh:
        return 'I'
    else:
        return 'S'
    

def probablistic_SIR_rule(neighboor_states, current_state, 
                          transmission_prob=0.5, mutation_prob=0.0):
        
    if current_state == 'I':
        return 'R'
    elif current_state == 'R':
        if np.random.rand() < mutation_prob:
            return 'S'
        return 'R'
    else: 
        for i in range(neighboor_states.count('I')):
            if np.random.rand() < transmission_prob:
                return 'I'
        return 'S'


def graph_with_moore_neighboorhood(m, n, periodic=False):
    g = nx.grid_2d_graph(m, n, periodic=periodic)
    
    # Add right diagonal connections
    for r in range(m - 1):
        for c in range(n - 1):
            g.add_edge((r, c), (r + 1, c + 1))
    
    # Add left diagonal
    for r in range(0, m - 1):
        for c in range(1, n):
            g.add_edge((r, c), (r + 1, c - 1))

    return g


def grid_layout(g):
    pos = {}
    for key in g.nodes:
        pos[key] = np.array(key)
    return pos


def gridToRBG(g):
    n, m = max(g.nodes)
    n += 1
    m += 1
    rgb_img = np.zeros(shape=(m, n, 3), dtype=np.int)
    for node_inx in g.nodes:
        if g.nodes[node_inx]['state'] == 'S':
            rgb_img[node_inx][1] = 255
        elif g.nodes[node_inx]['state'] == 'I':
            rgb_img[node_inx][0] = 255
        elif g.nodes[node_inx]['state'] == 'R':
            rgb_img[node_inx][2] = 255
        else:
            raise ValueError()
    return rgb_img


# %% CELLULAR AUTOMATA CLASS 
class CellularAutomata:

    def __init__(self, rule=None, graph=Graph(), node_states=['on', 'off']):
        
        self.rule = rule
        self.graph = graph
        self.node_states = node_states
        self.generation = 0
        
    def setNodeStates(self, state_init='random'):
        
        state_count = [0] * len(self.node_states)
        
        if state_init == 'random':
            for node_inx in self.graph.nodes:
                rand_inx = np.random.randint(0, len(self.node_states))
                state_count[rand_inx] += 1
                self.graph.nodes[node_inx]['state'] = self.node_states[rand_inx]
                
        if isinstance(state_init, dict):
            for node_inx in state_init:
                self.graph.nodes[node_inx]['state'] = state_init[node_inx]
                count_inx = self.node_states.index(state_init[node_inx])
                state_count[count_inx] += 1
        
        self.state_count = [tuple(state_count)]
                
        return [tuple(state_count)]
    
    def getNodeColorsByState(self, colorDict):
        node_color = []
        for node_inx in self.graph.nodes:
            node_color.append(colorDict[self.graph.nodes[node_inx]['state']])
        return node_color
    
    def getNeighboorNodesIndx(self, node_inx, size=1):
        adj_nodes_set = set(self.graph.adj[node_inx])
        if size > 1:
            neighboor_adj_nodes = set()
            for neighboor_inx in adj_nodes_set:        
                neighboor_adj_nodes.update(set(self.getNeighboorNodesIndx(neighboor_inx, 
                                                                          size=size - 1)))
            adj_nodes_set.update(neighboor_adj_nodes)
            adj_nodes_set.discard(node_inx)
        return list(adj_nodes_set)
    
    def getNeighboorNodesStates(self, node_inx, size=1):
        neighboor_inxs = self.getNeighboorNodesIndx(node_inx, size=size)
        states = [self.graph.nodes[inx]['state'] for inx in neighboor_inxs]
        return states
        
    def step(self, n=1):
        state_counts = []
        for _ in range(n):
            new_states = []
            for node_inx in self.graph.nodes:
                current_state = self.graph.nodes[node_inx]['state']
                neighbood_states = self.getNeighboorNodesStates(node_inx)
                new_states.append(self.rule(neighbood_states, 
                                            current_state)) 
            
            over_all_state_change = False
            for node_inx, new_state in zip(self.graph.nodes, new_states):
                if not self.graph.nodes[node_inx]['state'] == new_state:
                    self.graph.nodes[node_inx]['state'] = new_state
                    over_all_state_change = True
            
            if not over_all_state_change:
                break
            state_counts.append(tuple([new_states.count(state) for state in self.node_states]))
        
        self.state_count += state_counts
        
        return state_counts
                        
    def maxDegreeNode(self):
        degree_dict = dict(nx.degree(self.graph))
        max_degree = max(degree_dict.values())
        max_degree_nodes = [node for node in degree_dict if degree_dict[node] == max_degree]
        return max_degree_nodes

    def minDegreeNode(self):
        degree_dict = dict(nx.degree(self.graph))
        min_degree = min(degree_dict.values())
        min_degree_nodes = [node for node in degree_dict if degree_dict[node] == min_degree]
        return min_degree_nodes
    
    def removeSingleNodes(self):
        nodes_for_removal = [ node for node in self.graph.nodes if len(self.graph.adj[node]) == 0]
        self.graph.remove_nodes_from(nodes_for_removal)