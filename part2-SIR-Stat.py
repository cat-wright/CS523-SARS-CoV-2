"""
Author: Justin Deterding
Created: Tue Mar 17 15:09:24 2020
Description:
    In this script multiple runs of the SIR population dydnamics are evaluated.
    From the multiple runs the average or the S, I, and R populations are 
    plotted. The 95% confidence interval of the I population is also displayed.
    
    This script only considers S->I->R, with no transition of a population from
    infected to recovered.
"""
# %% Imports 
import CellularAutomata as ca
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from random import sample
import matplotlib 

matplotlib.rcParams.update({'font.size': 12})

# %% Rule Definitions
   

def results_to_array(results):
    """
    Takes rusults from the ca_driver and converts to an array of common length.
    """
    test = results.copy()    
    max_len = max([len(pop) for pop in test])
    
    for inx, pop in enumerate(test):
        append_len = max_len - len(test[inx])
        if append_len > 0:
            test[inx] = np.append(test[inx], [test[inx][-1, :]] * append_len, 
                                  axis=0)
        
    return np.array(test)


def ca_driver(graph, n_initial=1, thresh=0.5):
    """
    Simple function for implimenting the multiple trials within a loop.
    """
    rule = lambda n, c: ca.probablistic_SIR_rule(n, c, thresh) 
        
    ca_prob = ca.CellularAutomata(rule=rule,
                                  graph=graph,
                                  node_states=['S', 'I', 'R'])

    initial_state = {node_inx: 'S' for node_inx in ca_prob.graph.nodes}
    sample_nodes = sample(list(initial_state), n_initial)
    for node in sample_nodes:
        initial_state[node] = 'I'

    ca_prob.setNodeStates(initial_state)
    grid_size = max(graph.nodes)
    ca_prob.step(n=grid_size[0] * grid_size[1] * 10)

    return np.array(ca_prob.state_count)


# %% Statistical analysis of problistic SIR
graph = ca.graph_with_moore_neighboorhood(40, 40)
print('Remaining Nodes: ', len(graph))
print(nx.info(graph))
print('Density: ', nx.density(graph))
N_samples = 20

data = []

for initial_infected in [1, 3, 5]:
    for prob_transmission in [0.2, 0.4, 0.6]:
        data.append({'graph': graph.copy(),
                     'N initial': initial_infected,
                     'P(S->T)': prob_transmission})
        
        
for i, d in enumerate(data):
    print('Starting data set: {}/{}'.format(i, len(data)))        
    data[i]['data'] = [ca_driver(d['graph'], d['N initial'], d['P(S->T)']) for i in range(N_samples)]
    data[i]['data'] = results_to_array(data[i]['data'])


# %% PLOTTING
fig, axarr = plt.subplots(nrows=3, ncols=3, figsize=(12, 12))

inx = 0
for d, ax in zip(data, axarr.ravel()):
    
    plt.axes(ax)
    
    ax.set_title('P(I)={:2.1f} N={}'.format(d['P(S->T)'], d['N initial']))
    
    if inx in [6, 7, 8]:
        ax.set_xlabel('generation')
    if inx in [0, 3, 6]:
        ax.set_ylabel('S & R Population Count')
    
    ln_SIR = plt.plot(np.average(d['data'], axis=0)[:, 0], '--bo', label='S')
    ln_SIR += plt.plot(np.average(d['data'], axis=0)[:, 2], '--go', label='R')
    
    axt = plt.twinx(ax)    
    plt.axes(axt)
    
    if inx in [2, 5, 8]:
        axt.set_ylabel('I Population Count')
    
    ln_SIR += plt.plot(np.average(d['data'], axis=0)[:, 1], '--ro', label='I')
    
    over_line  = np.average(d['data'], axis=0)[:, 1] + np.std(d['data'], 
                                                              axis=0)[:, 1]
    under_line = np.average(d['data'], axis=0)[:, 1] - np.std(d['data'], 
                                                              axis=0)[:, 1]
    plt.fill_between(np.arange(len(d['data'][0])), under_line, over_line, 
                     color='r', alpha=.1)
    ln_SIR += plt.fill(np.NaN, np.NaN, 'r', alpha=0.1, label='95% C.I.')
    inx += 1
    
plt.legend(ln_SIR, [ln.get_label() for ln in ln_SIR])
