"""
Author: Justin Deterding
Created: Tue Mar 17 15:09:24 2020
Description:
    In this script multiple runs of the SIR population dydnamics are evaluated.
    From the multiple runs the average or the S, I, and R populations are 
    plotted. The 95% confidence interval of the I population is also displayed.
    
    This script only considers S->I->R->S, with no transition of a population 
    from infected to recovered with the possibility of becoming re suseptable.
"""
# %% Imports 
import CellularAutomata as ca
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from random import sample
import matplotlib 

matplotlib.rcParams.update({'font.size': 12})
# %% Function Definitions
   

def results_to_array(results):
    test = results.copy()    
    max_len = max([len(pop) for pop in test])
    
    for inx, pop in enumerate(test):
        append_len = max_len - len(test[inx])
        if append_len > 0:
            test[inx] = np.append(test[inx], [test[inx][-1, :]] * append_len, 
                                  axis=0)
        
    return np.array(test)


def ca_driver(graph, n_initial=1, thresh=0.5, mute=0.5):
    
    rule = lambda n, c: ca.probablistic_SIR_rule(n, c, thresh, mute) 
        
    ca_prob = ca.CellularAutomata(rule=rule,
                                  graph=graph,
                                  node_states=['S', 'I', 'R'])

    initial_state = {node_inx: 'S' for node_inx in ca_prob.graph.nodes}
    sample_nodes = sample(list(initial_state), n_initial)
    for node in sample_nodes:
        initial_state[node] = 'I'

    ca_prob.setNodeStates(initial_state)
    ca_prob.step(n=80)

    return np.array(ca_prob.state_count)


# %% Statistical analysis of problistic SIR
graph = ca.graph_with_moore_neighboorhood(40, 40)
print('Remaining Nodes: ', len(graph))
print(nx.info(graph))
print('Density: ', nx.density(graph))
N_samples = 10

data = []

for prob_mutation in [0.08, 0.1, 0.12]:
    for prob_transmission in [0.4, 0.6, 0.8]:
        data.append({'graph': graph.copy(),
                     'N initial': 5,
                     'P(S->T)': prob_transmission,
                     'P(R->S)': prob_mutation})
        
        
for i, d in enumerate(data):
    print('Starting data set: {}/{}'.format(i + 1, len(data)))
    results = []
    for j in range(N_samples):
        print('\t Sample {} of {}'.format(j + 1, N_samples))
        results += [ca_driver(d['graph'], d['N initial'], 
                              d['P(S->T)'], d['P(R->S)'])]   
    data[i]['data'] = results_to_array(results)

# %%

fig, axarr = plt.subplots(nrows=3, ncols=3, figsize=(12, 12))

inx = 0
for d, ax in zip(data, axarr.ravel()):
    tot_pop = len(d['graph'].nodes)
    plt.axes(ax)
    
    ax.set_title('P(I)={:2.1f} P(S)={:3.2f}'.format(d['P(S->T)'],
                                                    d['P(R->S)']))
    
    if inx in [6, 7, 8]:
        ax.set_xlabel('generation')
    if inx in [0, 3, 6]:
        ax.set_ylabel('Percent of Population \n Recovered or Suseptable')
    
    ln_SIR  = plt.plot(np.average(d['data'] / tot_pop, axis=0)[:, 0], 
                       '--bo', label='S')
    ln_SIR += plt.plot(np.average(d['data'] / tot_pop, axis=0)
                       [:, 2], '--go', label='R')
    
    axt = plt.twinx(ax)    
    plt.axes(axt)
    
    if inx in [2, 5, 8]:
        axt.set_ylabel('Percent of Population \n Infected')
    
    ln_SIR += plt.plot(np.average(d['data'], axis=0)[:, 1] / tot_pop, 
                       '--ro', label='I')
    
    over_line = (np.average(d['data'], axis=0)[:, 1] / tot_pop 
                 + np.std(d['data'], axis=0)[:, 1] / tot_pop)
    under_line = (np.average(d['data'], axis=0)[:, 1] / tot_pop
                  - np.std(d['data'], axis=0)[:, 1] / tot_pop)
    
    plt.fill_between(np.arange(len(d['data'][0])), under_line, over_line, 
                     color='r', alpha=.1)
    ln_SIR += plt.fill(np.NaN, np.NaN, 'r', alpha=0.1, label='95% C.I.')
    inx += 1
    
plt.legend(ln_SIR, [ln.get_label() for ln in ln_SIR])
