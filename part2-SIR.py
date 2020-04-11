"""
Author: Justin Deterding
Created: Tue Mar 17 15:09:24 2020
Description:
    This script generates a figure with 3 subfigures. The first Is the initial 
    configuration of the CA, the second is the state of the CA some number of 
    generations into the CA and finnaly a plot of the SIR populations for each 
    generation.
    This is done for both the probablistic and deterministic models.
"""
# %% Imports 
import CellularAutomata as ca
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import matplotlib 

matplotlib.rcParams.update({'font.size': 12})

# %% CA Initilization
  
prob_rule = lambda n, c: ca.probablistic_SIR_rule(n, c, 0.4, 0.4)

ca_prob = ca.CellularAutomata(rule=prob_rule,
                              graph=ca.graph_with_moore_neighboorhood(40, 40),
                              node_states=['S', 'I', 'R'])

deter_rule = lambda n, c: ca.deterministic_SIR_rule(n, c, 1.0)

ca_deter = ca.CellularAutomata(rule=deter_rule,
                               graph=ca.graph_with_moore_neighboorhood(40, 40),
                               node_states=['S', 'I', 'R'])

initial_state = {node_inx: 'S' for node_inx in ca_prob.graph.nodes}
initial_state[(19, 19)] = 'I'

ca_prob.setNodeStates(initial_state)
ca_deter.setNodeStates(initial_state)

# %% Plot initilization

# So long as they are all initilized with the same size grid this works
kwargs = {'pos': ca.grid_layout(ca_prob.graph),
          'node_size': 10}

for ca_type in [ca_prob, ca_deter]:
    
    fig, axarr = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))
    
    plt.axes(axarr[0])
    axarr[0].set_title('Initial Configuration')
    node_color = ca_type.getNodeColorsByState(colorDict={'S': 'b',
                                                         'I': 'r',
                                                         'R': 'g'})
    nx.draw(ca_type.graph,
            node_color=node_color, **kwargs)
    
    ca_type.step(n=25)
    
    plt.axes(axarr[1])
    axarr[1].set_title('25 generations')
    node_color = ca_type.getNodeColorsByState(colorDict={'S': 'b',
                                                         'I': 'r',
                                                         'R': 'g'})
    nx.draw(ca_type.graph,
            node_color=node_color, **kwargs)
    
    ca_type.step(n=100)
    
    plt.axes(axarr[2])
    
    axarr[2].set_title('SIR Dynamics')
    axarr[2].set_xlabel('generation')
    axarr[2].set_ylabel('Population Count')
    
    state_count_arr = np.array(ca_type.state_count)
    ln_SIR = plt.plot(state_count_arr[:, 0], '--bo', label='S')
    ln_SIR += plt.plot(state_count_arr[:, 2], '--go', label='R')
    plt.axes(plt.twinx(axarr[2]))
    ln_SIR += plt.plot(state_count_arr[:, 1], '--ro', label='I')
    plt.legend(ln_SIR, [ln.get_label() for ln in ln_SIR])