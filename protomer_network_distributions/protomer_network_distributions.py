import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import tellurium as te


class protomer_network():


  def __init__(self, protomer_dict):

    """
    Initialize protomer_network class with a protomer_dict
    This protomer_dict contains the dissociation constants of interactions between all species that may be used
    The protomer_dict can optionally include kinetics values as well, but will otherwise assume diffusion-limited on rates
    The class has a self.simulations list, which is used to record the input sub_protomer_dicts for each simulation
    The self.simulation_distributions list records the equilibrium data (dicts) of all simulations

    Attributes
    ----------
    protomer_dict: dictionary of dimer species and their dissociation constants in the form str: float, e.g. 'AB': 10 * 10**-9
                   Can optionally contain on rates in the form 'on_AB': 1 * 10**9
                   This dictionary at a minimum should contain the dissociation constants for every possible interaction between protomers that may be used
    simulations: recorded list of input subset protomer dictionaries for each simulation that is run
    simulation_distributions: recorded list of equilibrium distribution dictionaries for each simulation that is run
    
    """

    self.protomer_dict = protomer_dict
    self.simulations = []
    self.simulation_distributions = []


  def SimulateEquilibriumSpeciesDistribution(self, sub_protomer_dict, model_summary=True, plot=True, equilibrium_time=1):
    
    """
    Simulates the equilibrium species distribution for a desired protomer network
    This protomer network must consist of a subset of the protomers with which the class is initialized (or all protomers)
    The protomer network is input as a dictionary of protomer species, where each species is a key for its own dictionary
    Each protomer dictionary is either empty or contains a concentration and/or label value
    Concentrations are 10 nM by default if not provided
    Labels indicate if the protomer species has one half of a heterodimeric split reporter, denoted as booleans
    These labels are used to predict the binary output if the protomer network is used as a logic gate
    Optionally prints out a model_summary and plots the simulation over time (default=True)
    Simulation time is also tunable (default=1) to modify how long to run the simulation to reach equilibrium

    Parameters
    ----------
    sub_protomer_dict: dictionary of the form str: dict, with each protomer species in the desired network to simulate
                       Each protomer dict within this dict has its optional concentration and/or label
                       e.g. 'A': {'concentration': 4 * 10**-9, 'label' : False}, or 'B' : {}, or 'C' : {'label' : True}
    model_summary: boolean, default=True
    plot: boolean, default=True
    equilibrium_time: int, default=1
    
    Returns
    -------
    equilibrium_distribution: dictionary of all species and their concentrations at equilibrium
    
    """

    # make a model_str for the protomer network to simulate
    model_str = ''

    # first, add all reactions to the model_str
    for a_ind,protomer_A in enumerate(sub_protomer_dict):
      for b_ind,protomer_B in enumerate(sub_protomer_dict):
        # for all pairs of reactions, each pair only once
        if b_ind > a_ind:
          # add forward and reverse reaction for protomer pair binding
          model_str += f'{protomer_A} + {protomer_B} -> {protomer_A}_{protomer_B}; on_{protomer_A}_{protomer_B}*{protomer_A}*{protomer_B} \n'
          model_str += f'{protomer_A}_{protomer_B} -> {protomer_A} + {protomer_B}; off_{protomer_A}_{protomer_B}*{protomer_A}_{protomer_B} \n'
    model_str += ' \n'

    # second, add all rate constants to the model_str
    for a_ind,protomer_A in enumerate(sub_protomer_dict):
      for b_ind,protomer_B in enumerate(sub_protomer_dict):
        # for all pairs of reactions, each pair only once
        if b_ind > a_ind:
          # find dissociation constant in protomer_dict, checking both ways in case input is reversed
          if f'{protomer_A}_{protomer_B}' in self.protomer_dict:
            affinity_AB = self.protomer_dict[f'{protomer_A}_{protomer_B}']
          else:
            affinity_AB = self.protomer_dict[f'{protomer_B}_{protomer_A}']
          # if on rate is provided, use it, otherwise use diffusion-limited on rate (will not affect equilibrium distribution, only kinetics)
          if f'on_{protomer_A}_{protomer_B}' in self.protomer_dict or f'on_{protomer_B}_{protomer_A}' in self.protomer_dict:
            if f'on_{protomer_A}_{protomer_B}' in self.protomer_dict:
              on_rate_AB = self.protomer_dict[f'on_{protomer_A}_{protomer_B}']
            else:
              on_rate_AB = self.protomer_dict[f'on_{protomer_B}_{protomer_A}']
          else:
            # standard diffusion-limited on rate
            on_rate_AB = 8 * 10**9 
          # calculate off rate from on rate and dissociation constant
          off_rate_AB = on_rate_AB * affinity_AB
          model_str += f'on_{protomer_A}_{protomer_B} = {str(on_rate_AB)} \n'
          model_str += f'off_{protomer_A}_{protomer_B} = {str(off_rate_AB)} \n'
    model_str += ' \n'

    # third, add concentrations of all species to the model_str
    for protomer in sub_protomer_dict:
      if 'concentration' in sub_protomer_dict[protomer]:
        # add custom protomer concentrations if provided
        protomer_conc = sub_protomer_dict[protomer]['concentration']
      else:
        # otherwise, standard 10 nM concentrations
        protomer_conc = 10 * 10**-9
      model_str += f'{protomer} = {str(protomer_conc)} \n'

    # print model_str to check
    if model_summary:
      print('Model Summary:')
      print()
      print(model_str)
      print()

    # run simulation
    rr = te.loada(model_str)
    data = rr.simulate(0,equilibrium_time,100)
    if plot:
      rr.plot()
    equilibrium_distribution = {}
    for species in data.colnames:
      if species != 'time':
        equilibrium_distribution[species[1:-1]] = data[species][-1]

    # record inputs and outputs in class for further processing
    self.simulations.append(sub_protomer_dict)
    self.simulation_distributions.append(equilibrium_distribution)

    return equilibrium_distribution


  def PredictBinaryOutput(self, simulation_ind=-1, signal_threshold= 5 * 10**-9, return_cumulative_signal=False):

    """
    For a protomer network (including labels) and simulated equilibrium distribution, predicts binary output from split reporter
    Cumulative signal is added for dimeric species with complementary labels by adding their concentrations at equilibrium
    If cumulative signal is above the threshold, output is True (1), otherwise False (0)

    Parameters
    ----------
    simulation_ind: int, index of simulation and equilibrium distribution within class stored lists (default=-1, most recent simulation)
    signal_threshold: float, default = 5 * 10**-9
    return_cumulative_signal: boolean, default=False
    
    Returns
    -------
    binary_out: boolean, representing binary logic gate output from split reporter
    cumulative_signal: if set to return, return raw signal rather than boolean

    """

    sub_protomer_dict = self.simulations[simulation_ind]
    equilibrium_distribution = self.simulation_distributions[simulation_ind]

    cumulative_signal = 0
    for a_ind,protomer_A in enumerate(sub_protomer_dict):
      if 'label' in sub_protomer_dict[protomer_A]:
        for b_ind,protomer_B in enumerate(sub_protomer_dict):
          # for all pairs of reactions, each pair only once, where both are labelled
          if b_ind > a_ind and 'label' in sub_protomer_dict[protomer_B]:
            if (sub_protomer_dict[protomer_A]['label'] == True and sub_protomer_dict[protomer_B]['label'] == False) or \
               (sub_protomer_dict[protomer_A]['label'] == False and sub_protomer_dict[protomer_B]['label'] == True):
               if f'{protomer_A}_{protomer_B}' in equilibrium_distribution:
                 cumulative_signal += equilibrium_distribution[f'{protomer_A}_{protomer_B}']
               else:
                 cumulative_signal += equilibrium_distribution[f'{protomer_B}_{protomer_A}']

    # return binary output depending on if cumulative signal is above or below threshold
    if cumulative_signal > signal_threshold:
      binary_output = True
    else:
      binary_output = False

    # optionally return cumulative signal, otherwise return binary output
    if return_cumulative_signal:
      return cumulative_signal
    else:
      return binary_output


  def VisualizeEquilibriumDistribution(self, simulation_ind=-1):

    """
    For a given simulation, plots concentrations of all species to visualize relative distributions

    Parameters
    ----------
    simulation_ind: int, index of simulation and equilibrium distribution within class stored lists (default=-1, most recent simulation)

    Returns
    -------
    nothing, but displays a plot in notebook

    
    """

    sub_protomer_dict = self.simulations[simulation_ind]
    equilibrium_distribution = self.simulation_distributions[simulation_ind]

    # find all reporter dimers for visualizing in different color
    reporter_dimers = []
    for a_ind,protomer_A in enumerate(sub_protomer_dict):
      if 'label' in sub_protomer_dict[protomer_A]:
        for b_ind,protomer_B in enumerate(sub_protomer_dict):
          # for all pairs of reactions, each pair only once, where both are labelled
          if b_ind > a_ind and 'label' in sub_protomer_dict[protomer_B]:
            if (sub_protomer_dict[protomer_A]['label'] == True and sub_protomer_dict[protomer_B]['label'] == False) or \
               (sub_protomer_dict[protomer_A]['label'] == False and sub_protomer_dict[protomer_B]['label'] == True):
              reporter_dimers.append(f'{protomer_A}_{protomer_B}')

    species = []
    values = []
    labels = []
    # split by monomer vs. dimer, displaying monomers first on bar plot
    for name in equilibrium_distribution:
      if '_' not in name:
        species.append(name)
        values.append(equilibrium_distribution[name])
        labels.append('no signal') # no monomers will have signal
    # now add dimers
    for name in equilibrium_distribution:
      if '_' in name:
        species.append(name)
        values.append(equilibrium_distribution[name])
        if name in reporter_dimers:
          labels.append('signal')
        else:
          labels.append('no signal')
    
    # plot concentrations
    sns.barplot(x=species, y=values, hue=labels, dodge=False)
    plt.yscale('log')
    plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    plt.xlabel('species')
    plt.ylabel('concentration (M)')

    return

