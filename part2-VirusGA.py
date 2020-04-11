"""
Author: Justin Deterding
Created: Sat Mar 21 07:33:32 2020
Description:
"""
# Imports
# DEEP imports
from deap import base, creator
from deap import tools
# CA Imports
from CellularAutomata import CellularAutomata as CA
from CellularAutomata import graph_with_moore_neighboorhood as grid
from CellularAutomata import probablistic_SIR_rule as rule
# Other imports
import random

# Function Definitions


def mutGaussian_prob(individual, mu, sigma, indpb):
    
    for inx in range(len(individual)):
        if random.random() < indpb:
            cur_state = individual[inx]
            new_state = cur_state + random.gauss(mu, sigma)
            if new_state > 1.0:
                individual[inx] = 1.0
            elif new_state < 0.0:
                individual[inx] = 0.0
            else:
                individual[inx] = new_state
    return (individual,)


def evaluate(individual):
    p_trans = individual[0]
    p_mutat = individual[1]
    print('P(I)={:4.3f} \t P(M)={:4.3f}'.format(p_trans, p_mutat))
    
    ca = CA(rule=lambda n_states, c_state: rule(n_states, c_state,
                                                p_trans, p_mutat), 
            graph=grid(40, 40), node_states=['S', 'I', 'R'])
    
    Tot_S, Tot_I, Tot_R = 0, 0, 0
    num_samples = 10
    max_steps = 40
    
    for n in range(num_samples):
        
        inital_state = {node: 'S' for node in ca.graph}
        for node in random.sample(list(ca.graph), 5):
            inital_state[node] = 'I'
        ca.setNodeStates(inital_state)

        ca.step(max_steps)
        num_steps = len(ca.state_count) - 1
        print('\t{}:{}/{}'.format(n + 1, num_steps, max_steps), end='')
        N_S, N_I, N_R = ca.state_count[-1]
        Tot_S += N_S / (N_S + N_I + N_R)
        Tot_I += N_I / (N_S + N_I + N_R)
        Tot_R += N_R / (N_S + N_I + N_R)
    print()
    print('\tS:{:4.3}\tI:{:4.3}\tR:{:4.3}\t'.format(Tot_S / num_samples,
                                                    Tot_I / num_samples,
                                                    Tot_R / num_samples,))
    print('Score:', Tot_I / num_samples - p_mutat)
    return Tot_I / num_samples,


creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

# 2 7-bit long probabilities
IND_SIZE = 2

toolbox = base.Toolbox()

#pool = multiprocessing.Pool(4)
toolbox.register("map", map)

toolbox.register("attribute", random.random)
toolbox.register("individual", tools.initRepeat, creator.Individual,
                 toolbox.attribute, n=IND_SIZE)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", mutGaussian_prob, mu=0, sigma=0.125, indpb=0.5)
toolbox.register("select", tools.selTournament, tournsize=3)
toolbox.register("evaluate", evaluate)

# Set up history logbook
history = tools.History()



toolbox.decorate("mate", history.decorator)
toolbox.decorate("mutate", history.decorator)

logbook = tools.Logbook()

hall_of_fame = tools.HallOfFame(10)


def main():
    pop = toolbox.population(n=10)
    history.update(pop)    
    
    # CXPB: Cross Over Probability
    # MUTPB: Mutation Probability
    # NGEN: Number of Generations
    CXPB, MUTPB, NGEN = 0.5, 0.5, 50
    
    # Evaluate the entire population
    fitnesses = list(toolbox.map(toolbox.evaluate, pop))
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit
    
    all_fitnesses = [ind.fitness.values for ind in pop]

    print('Initial max fitness: {}'.format(max(all_fitnesses)))
    
    logbook.record(gen=0, 
                   min=min(all_fitnesses), 
                   max=max(all_fitnesses),
                   max_ind=pop[fitnesses.index(max(fitnesses))])
    hall_of_fame.update(pop)

    for g in range(NGEN):
        
        print("Evaluating Generation: {}".format(g + 1))
        # Select the next generation individuals
        offspring = toolbox.select(pop, len(pop))
        # Clone the selected individuals
        offspring = list(toolbox.map(toolbox.clone, offspring))

        # Apply crossover and mutation on the offspring
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < CXPB:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values

        for mutant in offspring:
            if random.random() < MUTPB:
                toolbox.mutate(mutant)
                del mutant.fitness.values

        # Evaluate the individuals with an invalid fitness
        #invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        #fitnesses = list(map(toolbox.evaluate, invalid_ind))
        
        fitnesses = list(map(toolbox.evaluate, offspring))
        
        
        for ind, fit in zip(offspring, fitnesses):
            ind.fitness.values = fit
              
        all_fitnesses = [ind.fitness.values for ind in pop]
        print('Max fitness: {}'.format(max(all_fitnesses)))
        logbook.record(gen=g + 1, 
                       min=min(all_fitnesses), 
                       max=max(all_fitnesses),
                       max_ind=pop[fitnesses.index(max(fitnesses))])
                                   
        hall_of_fame.update(pop)

        # The population is entirely replaced by the offspring
        pop[:] = offspring

    return pop


if __name__ == "__main__":
    from time import time
    
    start_time = time()
    pop = main()
    end_time = time()
    print(end_time - start_time)
    # %%
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    matplotlib.rcParams.update({'font.size': 12})
    
    gen, max_fit, max_ind = logbook.select('gen', 'max', 'max_ind')

    p_tra = [ind[0] for ind in max_ind]    
    p_mut = [ind[1] for ind in max_ind]
    
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)
    
    l1, = ax1.plot(gen, np.array(max_fit), '--or', label='Max Fitness')
    ax1.set_ylabel('Percent Infected')
    ax1.legend()  
    l2, = ax2.plot(gen, p_tra, '--ok', label='P(I)')
    l3, = ax2.plot(gen, p_mut, '--og', label='P(S)')
    ax2.set_ylabel('Virsus Attributes')
    ax2.set_xlabel('Generation')
    ax2.legend()
    #plt.legend((l1, l2, l3), ('Max Fitness', 'P(I)', 'P(S)'))