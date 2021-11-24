# 'new_ms_model.py' is a script to determine the number of cells in each stress level compartment at each time
# and will also calculate the amount of protein shed into the plasma based on the cell growth model
#
# This script/equation will be used in a Multiple Sclerosis compartmental model of protein shedding kinetics

# Author: Corey Chitwood
# v0.1 - Oct 2021

# NOTE:
# oligodendrocyte is abbreviated as "L" throughout
# neuron is abbreviated as "N" throughout

from tqdm import tqdm_notebook as tqdm
import numpy as np
from ms_stress_function import *
import math


def oligodendrocyte_growth(parameters, results):
    """ inputs: parameters & results dictionaries, parameters contain: n, IC_L, maxIter, stress function, rates
     also inputting: stress and remission data

     stress_prob(), birth_rate(), death_rate() are other functions within this script nested within
     oligodendrocyte_growth() oxidative_stress_equation() and remission_equation() are 'ms_stress_function.py'
     functions nested within oligodendrocyte_growth()"""

    IC_L = parameters['IC_L']
    n = parameters['n']
    maxIter = parameters['maxIter']
    num_cycles = parameters['num_cycles']
    stress_length = parameters['stress_length']
    stresst0 = parameters['stresst0']
    remission_length = parameters['remission_length']
    completed_stress_cycles = 0
    completed_remission_cycles = 0

    # To keep track of L's in each compartment at each time:
    # number of column = time, number of row = compartment
    compartments = np.zeros((maxIter + 1, n))

    # To keep track of L's that become stressed at each time:
    L_stressed = np.zeros(maxIter + 1)
    L_initially_stressed_only = np.zeros(maxIter + 1)

    # To keep track of L's that die at each time:
    L_deaths = np.zeros(maxIter + 1)

    # To keep track of new L's at each time:
    L_births = np.zeros(maxIter + 1)

    # To keep track of L's with stress = 0 at each time:
    total_L_unstressed = np.zeros(maxIter + 1)

    # To keep track of L's with stress > 0 at each time:
    total_L_stressed = np.zeros(maxIter + 1)

    # To keep track of total number of living L's
    total_L_alive = np.zeros(maxIter + 1)
    total_L_dead = np.zeros(maxIter + 1)

    total_L_alive[0] = IC_L
    total_L_unstressed[0] = IC_L
    compartments[0, 0] = IC_L
    # also initial conditions, but already set at 0 by using np.zeros()
    # total_L_dead[0] = 0
    # total_L_stressed[0] = 0

    stress_range_1 = np.linspace(1, stress_length, num=stress_length)
    stress_1 = oxidative_stress_equation(parameters, stresst0, 0, 0)
    # stress level will always begin at 0

    for t in tqdm(stress_range_1, desc='Running Oligodendrocyte Stress, Relapse #1', total=len(stress_range_1),
                  mininterval=1):
        # newly stressed = (# unstressed at t - 1) * (stress probability)
        t = int(t)

        L_stressed[t] = stress_prob(total_L_stressed[t - 1], parameters) * total_L_unstressed[t - 1]
        L_initially_stressed_only[t] = L_stressed[t]
        total_L_stressed[t] = total_L_stressed[t - 1] + L_stressed[t]
        total_L_unstressed[t] = total_L_unstressed[t - 1] - L_stressed[t]
        compartments[t, 0] += total_L_unstressed[t]

        total_L_dead[t] = total_L_dead[t - 1]
        total_L_alive[t] = total_L_alive[t - 1]

        # new births
        L_births[t] = birth_rate(total_L_dead[t], parameters) * total_L_alive[t]
        total_L_alive[t] += L_births[t]
        total_L_unstressed[t] += L_births[t]
        compartments[t, 0] += L_births[t]

        # death and compartment loop
        for i in range(1, t + 1):
            current_stress_level = stress_1[t - i]
            previous_stress_level = stress_1[t - i - 1]
            kD = death_rate(current_stress_level, parameters)
            deaths = kD * L_stressed[i]
            if deaths > 0:
                L_stressed[i] -= deaths
                L_deaths[t] += deaths
                total_L_stressed[t] -= deaths
                total_L_alive[t] -= deaths
                total_L_dead[t] += deaths

            if current_stress_level < (100 / n):
                compartments[t, 0] += L_stressed[i]
            elif current_stress_level < (100 / n) * 2:
                compartments[t, 1] += L_stressed[i]
            elif current_stress_level < (100 / n) * 3:
                compartments[t, 2] += L_stressed[i]
            elif current_stress_level < (100 / n) * 4:
                compartments[t, 3] += L_stressed[i]
            else:
                compartments[t, 4] += L_stressed[i]

    # end of relapse 1 loop
    completed_stress_cycles += 1

    remission_range_1 = np.linspace(stress_length + 1, stress_length + remission_length, num=remission_length)
    for t in tqdm(remission_range_1, desc='Running Oligodendrocyte Stress, Remission #1',
                  total=len(remission_range_1), mininterval=1):
        t = int(t)

        total_L_dead[t] = total_L_dead[t - 1]
        total_L_alive[t] = total_L_alive[t - 1]
        total_L_stressed[t] = total_L_stressed[t - 1] + L_stressed[t]
        total_L_unstressed[t] = total_L_unstressed[t - 1] - L_stressed[t]
        compartments[t, 0] += total_L_unstressed[t]

        # new births
        L_births[t] = birth_rate(total_L_dead[t], parameters) * total_L_alive[t]
        total_L_alive[t] += L_births[t]
        total_L_unstressed[t] += L_births[t]
        compartments[t, 0] += L_births[t]

        # death and compartment loop
        for i in range(1, int(max(stress_range_1))):
            current_stress_level = remission_equation(parameters, stress_length + 1, stress_1[stress_length - i], i)[
                t - (stress_length + 1)]
            kD = death_rate(current_stress_level, parameters)
            deaths = kD * L_stressed[i]
            if deaths > 0:
                L_stressed[i] -= deaths
                L_deaths[t] += deaths
                total_L_stressed[t] -= deaths
                total_L_alive[t] -= deaths
                total_L_dead[t] += deaths

            if current_stress_level < (100 / n):
                compartments[t, 0] += L_stressed[i]
            elif current_stress_level < (100 / n) * 2:
                compartments[t, 1] += L_stressed[i]
            elif current_stress_level < (100 / n) * 3:
                compartments[t, 2] += L_stressed[i]
            elif current_stress_level < (100 / n) * 4:
                compartments[t, 3] += L_stressed[i]
            else:
                compartments[t, 4] += L_stressed[i]

    # end of remission 1 loop
    completed_remission_cycles += 1

    # update results dictionary
    results.update(compartments=compartments)
    results.update(total_L_dead=total_L_dead)
    results.update(total_L_stressed=total_L_stressed)
    results.update(total_L_alive=total_L_alive)
    results.update(total_L_unstressed=total_L_unstressed)
    results.update(L_stressed=L_stressed)
    results.update(L_initially_stressed_only=L_initially_stressed_only)
    results.update(L_births=L_births)
    results.update(L_deaths=L_deaths)
    results.update(completed_stress_cycles=completed_stress_cycles)
    results.update(completed_remission_cycles=completed_remission_cycles)

    if num_cycles == 1:
        return


def neuron_growth(parameters, results):
    IC_N = parameters['IC_N']
    n = parameters['n']
    maxIter = parameters['maxIter']
    num_cycles = parameters['num_cycles']
    stress_length = parameters['stress_length']
    remission_length = parameters['remission_length']
    stresst0 = parameters['stresst0']
    N_completed_stress_cycles = 0
    N_completed_remission_cycles = 0
    neuron_extension = math.ceil(stress_length / 3)

    # To keep track of L's in each compartment at each time:
    # number of column = time, number of row = compartment
    N_compartments = np.zeros((maxIter + 1, n))

    # To keep track of L's that become stressed at each time:
    N_stressed = np.zeros(maxIter + 1)
    N_initially_stressed_only = np.zeros(maxIter + 1)

    # To keep track of L's that die at each time:
    N_deaths = np.zeros(maxIter + 1)

    # To keep track of L's with stress = 0 at each time:
    total_N_unstressed = np.zeros(maxIter + 1)

    # To keep track of L's with stress > 0 at each time:
    total_N_stressed = np.zeros(maxIter + 1)

    # To keep track of total number of living L's
    total_N_alive = np.zeros(maxIter + 1)
    total_N_dead = np.zeros(maxIter + 1)

    total_N_alive[0] = IC_N
    total_N_unstressed[0] = IC_N
    N_compartments[0, 0] = IC_N
    # also initial conditions, but already set at 0 by using np.zeros()
    # total_N_dead[0] = 0
    # total_N_stressed[0] = 0

    stress_range_1 = np.linspace(1, stress_length + neuron_extension, num=(stress_length + neuron_extension))
    stress_1 = oxidative_stress_equation(parameters, stresst0, 0, neuron_extension)
    # stress level will always begin at 0

    for t in tqdm(stress_range_1, desc='Running Neuron Stress, Relapse #1', total=len(stress_range_1),
                  mininterval=1):
        # newly stressed = (# unstressed at t - 1) * (stress probability)
        t = int(t)

        N_stressed[t] = neuron_stress_prob(t, parameters, results) * total_N_unstressed[t - 1]
        # N_stressed[t] = (results['compartments'][t, 1] - results['compartments'][t - 1, 1]) / 5
        # if N_stressed[t] < 0:
        #     N_stressed[t] = 0
        N_initially_stressed_only[t] = N_stressed[t]
        total_N_stressed[t] = total_N_stressed[t - 1] + N_stressed[t]
        total_N_unstressed[t] = total_N_unstressed[t - 1] - N_stressed[t]
        N_compartments[t, 0] += total_N_unstressed[t]

        total_N_dead[t] = total_N_dead[t - 1]
        total_N_alive[t] = total_N_alive[t - 1]

        # death and compartment loop
        for i in range(1, t + 1):
            current_stress_level = stress_1[t - i]
            kD = death_rate(current_stress_level, parameters)
            deaths = kD * N_stressed[i]
            if deaths > 0:
                N_stressed[i] -= deaths
                N_deaths[t] += deaths
                total_N_stressed[t] -= deaths
                total_N_alive[t] -= deaths
                total_N_dead[t] += deaths

            if current_stress_level < (100 / n):
                N_compartments[t, 0] += N_stressed[i]
            elif current_stress_level < (100 / n) * 2:
                N_compartments[t, 1] += N_stressed[i]
            elif current_stress_level < (100 / n) * 3:
                N_compartments[t, 2] += N_stressed[i]
            elif current_stress_level < (100 / n) * 4:
                N_compartments[t, 3] += N_stressed[i]
            else:
                N_compartments[t, 4] += N_stressed[i]

    # end of relapse 1 loop
    N_completed_stress_cycles += 1

    remission_range_1 = np.linspace(stress_length + neuron_extension + 1, stress_length + remission_length,
                                    num=(remission_length - neuron_extension))
    for t in tqdm(remission_range_1, desc='Running Neuron Stress, Remission #1',
                  total=len(remission_range_1), mininterval=1):
        t = int(t)

        total_N_dead[t] = total_N_dead[t - 1]
        total_N_alive[t] = total_N_alive[t - 1]
        total_N_stressed[t] = total_N_stressed[t - 1] + N_stressed[t]
        total_N_unstressed[t] = total_N_unstressed[t - 1] - N_stressed[t]
        N_compartments[t, 0] += total_N_unstressed[t]

        # death and compartment loop
        for i in range(1, int(max(stress_range_1))):
            current_stress_level = remission_equation(parameters, stress_length + 1 + neuron_extension, stress_1[stress_length - i], i)[
                t - (stress_length + 1 + neuron_extension)]
            kD = death_rate(current_stress_level, parameters)
            deaths = kD * N_stressed[i]
            if deaths > 0:
                N_stressed[i] -= deaths
                N_deaths[t] += deaths
                total_N_stressed[t] -= deaths
                total_N_alive[t] -= deaths
                total_N_dead[t] += deaths

            if current_stress_level < (100 / n):
                N_compartments[t, 0] += N_stressed[i]
            elif current_stress_level < (100 / n) * 2:
                N_compartments[t, 1] += N_stressed[i]
            elif current_stress_level < (100 / n) * 3:
                N_compartments[t, 2] += N_stressed[i]
            elif current_stress_level < (100 / n) * 4:
                N_compartments[t, 3] += N_stressed[i]
            else:
                N_compartments[t, 4] += N_stressed[i]

    # end of remission 1 loop
    N_completed_remission_cycles += 1

    # update results dictionary
    results.update(N_compartments=N_compartments)
    results.update(total_N_dead=total_N_dead)
    results.update(total_N_stressed=total_N_stressed)
    results.update(total_N_alive=total_N_alive)
    results.update(total_N_unstressed=total_N_unstressed)
    results.update(N_stressed=N_stressed)
    results.update(N_initially_stressed_only=N_initially_stressed_only)
    results.update(N_deaths=N_deaths)
    results.update(N_completed_stress_cycles=N_completed_stress_cycles)
    results.update(N_completed_remission_cycles=N_completed_remission_cycles)

    if num_cycles == 1:
        return


def stress_prob(total_L_stressed, parameters):
    """# function uses: current number of L's stressed, prob_stress_0, prob_stress_max, IC_L, total_L_stressed,
    prob_type """

    prob_stress_0 = parameters['prob_stress_0']
    prob_stress_max = parameters['prob_stress_max']
    IC_L = parameters['IC_L']
    prob_type = parameters['prob_type']

    x = np.linspace(0, 100, num=101)

    if prob_type == 1:  # exponential
        prob_eqn = np.exp(0.05 * x)
        prob_eqn = ((prob_stress_max - prob_stress_0) / max(prob_eqn)) * prob_eqn
        prob_eqn = prob_eqn + prob_stress_0

    elif prob_type == 2:  # linear
        prob_eqn = ((prob_stress_max - prob_stress_0) / 100) * x + prob_stress_0

    else:
        print('Stress Probability Error')
        prob_eqn = x * 'error'

    percent_stressed = round(total_L_stressed / IC_L)
    if percent_stressed > 100:
        percent_stressed = 100

    prob = prob_eqn[percent_stressed]

    return prob


def neuron_stress_prob(t, parameters, results):
    prob_stress_max = parameters['prob_stress_max']
    stress_length = parameters['stress_length']
    compartments = results['compartments']

    comparison = stress_length / 2
    #
    for i in range(0, stress_length):
        if compartments[i, 1] > 0:
            comparison = i
            break
        else:
            comparison = comparison

    # comparison = 25

    if t < comparison:
        prob = 0
    else:
        prob = prob_stress_max / 2

    return prob


def death_rate(current_stress_level, parameters):
    """# function uses: current stress level, kDL_0, kDL_max, current_stress_level, death_rate_type"""

    death_rate_type = parameters['death_rate_type']
    kDL_0 = parameters['kDL_0']
    kDL_max = parameters['kDL_max']
    n = parameters['n']

    x = np.linspace(0, 100, num=101)

    if death_rate_type == 1:  # exponential
        death_rate_eqn = np.exp(0.05 * x)
        death_rate_eqn = (kDL_max - kDL_0) / max(death_rate_eqn) * death_rate_eqn
        kDL = death_rate_eqn[current_stress_level]

    elif death_rate_type == 2:  # linear
        death_rate_eqn = (kDL_max - kDL_0) / 100 * x
        kDL = death_rate_eqn[current_stress_level]

    elif death_rate_type == 3:  # CnDeathsOnly
        CnMin = 100 / n * (n - 1)
        if current_stress_level < CnMin:
            kDL = kDL_0
        else:
            kDL = kDL_max

    else:
        print('Death Rate Calculation Error')
        kDL = 'error'

    return kDL


def birth_rate(total_L_dead, parameters):
    """# birth rate dependent on how many L's are currently dead"""

    kBL_0 = parameters['kBL_0']
    kBL_max = parameters['kBL_max']
    birth_rate_type = parameters['birth_rate_type']
    IC_L = parameters['IC_L']

    percent_dead = total_L_dead / IC_L
    if percent_dead > 100:
        percent_dead = 100

    if birth_rate_type == 1:  # exponential
        kBL = np.exp(0.05 * percent_dead) * ((kBL_max - kBL_0) / 148)

    elif birth_rate_type == 2:  # linear
        kBL = (kBL_max - kBL_0) / 100 * percent_dead

    else:
        print('Birth Rate Calculation Error')
        kBL = 'error'

    return kBL


def protein_shedding(parameters, results):
    # inputs: cells alive/dead at each time, shedding rates, normal shedding rates
    # normal shedding rates

    # this model does not incorporate a spatial component, so all proteins shed from the cells will instantaneously be
    # "dumped" into the plasma

    num_cycles = parameters['num_cycles']
    stress_length = parameters['stress_length']
    remission_length = parameters['remission_length']

    shedding_length = (remission_length + stress_length) * num_cycles

    total_L_alive = results['total_L_alive']
    total_N_alive = results['total_N_alive']

    L_deaths = results['L_deaths']
    N_deaths = results['N_deaths']

    compartments = results['compartments']
    N_compartments = results['N_compartments']

    IC_L = parameters['IC_L']
    IC_N = parameters['IC_N']

    maxIter = parameters['maxIter']

    # protein shedding rates
    uh_ec = parameters['uh_ec']
    uh_nec = parameters['uh_nec']
    un_ec = parameters['un_ec']
    un_nec = parameters['un_nec']

    # protein degradation rates
    kE_ec = parameters['kE_ec']
    kE_nec = parameters['kE_nec']

    current_L_ec = np.zeros(shedding_length + 1)
    current_N_ec = np.zeros(shedding_length + 1)
    current_L_nec = np.zeros(shedding_length + 1)
    current_N_nec = np.zeros(shedding_length + 1)
    dumped_L_nec = np.zeros(shedding_length + 1)
    dumped_N_nec = np.zeros(shedding_length + 1)
    dumped_L_ec = np.zeros(shedding_length + 1)
    dumped_N_ec = np.zeros(shedding_length + 1)
    cleared_L_nec = np.zeros(shedding_length + 1)
    cleared_N_nec = np.zeros(shedding_length + 1)
    cleared_L_ec = np.zeros(shedding_length + 1)
    cleared_N_ec = np.zeros(shedding_length + 1)

    a = parameters['a']  # a is the ratio of oligodendrocyte to neuron shedding

    shedding_range = np.linspace(1, shedding_length, num=shedding_length)

    for t in tqdm(shedding_range, desc='Running Protein Shedding Analysis', total=(len(shedding_range)), mininterval=1):
        t = int(t)

        # shedding due to stress level
        N_compartment_1_ec = uh_ec * N_compartments[t, 0]
        N_compartment_2_ec = uh_ec * N_compartments[t, 1] * 2
        N_compartment_3_ec = uh_ec * N_compartments[t, 2] * 3
        N_compartment_4_ec = uh_ec * N_compartments[t, 3] * 4
        N_compartment_5_ec = uh_ec * N_compartments[t, 4] * 5
        N_compartment_total_ec = N_compartment_1_ec + N_compartment_2_ec + N_compartment_3_ec + N_compartment_4_ec + N_compartment_5_ec

        L_compartment_1_ec = uh_ec * compartments[t, 0]
        L_compartment_2_ec = uh_ec * compartments[t, 1] * 2
        L_compartment_3_ec = uh_ec * compartments[t, 2] * 3
        L_compartment_4_ec = uh_ec * compartments[t, 3] * 4
        L_compartment_5_ec = uh_ec * compartments[t, 4] * 5
        L_compartment_total_ec = L_compartment_1_ec + L_compartment_2_ec + L_compartment_3_ec + L_compartment_4_ec + L_compartment_5_ec

        N_compartment_1_nec = uh_nec * N_compartments[t, 0]
        N_compartment_2_nec = uh_nec * N_compartments[t, 1] * 2
        N_compartment_3_nec = uh_nec * N_compartments[t, 2] * 3
        N_compartment_4_nec = uh_nec * N_compartments[t, 3] * 4
        N_compartment_5_nec = uh_nec * N_compartments[t, 4] * 5
        N_compartment_total_nec = N_compartment_1_nec + N_compartment_2_nec + N_compartment_3_nec + N_compartment_4_nec + N_compartment_5_nec

        L_compartment_1_nec = uh_nec * compartments[t, 0]
        L_compartment_2_nec = uh_nec * compartments[t, 1] * 2
        L_compartment_3_nec = uh_nec * compartments[t, 2] * 3
        L_compartment_4_nec = uh_nec * compartments[t, 3] * 4
        L_compartment_5_nec = uh_nec * compartments[t, 4] * 5
        L_compartment_total_nec = L_compartment_1_nec + L_compartment_2_nec + L_compartment_3_nec + L_compartment_4_nec + L_compartment_5_nec

        # normal shedding (uh)
        dumped_N_ec[t] = N_compartment_total_ec * a
        dumped_L_ec[t] = L_compartment_total_ec
        dumped_N_nec[t] = N_compartment_total_nec * a
        dumped_L_nec[t] = L_compartment_total_nec

        # shedding due to necrosis (un)
        dumped_N_ec[t] += un_ec * N_deaths[t] * a
        dumped_L_ec[t] += un_ec * L_deaths[t]
        dumped_N_nec[t] += un_nec * N_deaths[t] * a
        dumped_L_nec[t] += un_nec * L_deaths[t]

        # protein elimination (kE)
        cleared_N_ec[t] = kE_ec * current_N_ec[t - 1]
        cleared_L_ec[t] = kE_ec * current_L_ec[t - 1]
        cleared_N_nec[t] = kE_nec * current_N_nec[t - 1]
        cleared_L_nec[t] = kE_nec * current_L_nec[t - 1]

        # current protein in plasma
        current_N_ec[t] = dumped_N_ec[t] - cleared_N_ec[t] + current_N_ec[t - 1]
        if current_N_ec[t] < 0:
            current_N_ec[t] = 0
        current_L_ec[t] = dumped_L_ec[t] - cleared_L_ec[t] + current_L_ec[t - 1]
        if current_L_ec[t] < 0:
            current_L_ec[t] = 0
        current_N_nec[t] = dumped_N_nec[t] - cleared_N_nec[t] + current_N_nec[t - 1]
        if current_N_nec[t] < 0:
            current_N_nec[t] = 0
        current_L_nec[t] = dumped_L_nec[t] - cleared_L_nec[t] + current_L_nec[t - 1]
        if current_L_nec[t] < 0:
            current_L_nec[t] = 0

    # update results dictionary
    results.update(current_L_nec=current_L_nec)
    results.update(current_L_ec=current_L_ec)
    results.update(dumped_L_nec=dumped_L_nec)
    results.update(dumped_L_ec=dumped_L_ec)
    results.update(cleared_L_nec=cleared_L_nec)
    results.update(cleared_L_ec=cleared_L_ec)
    results.update(current_N_ec=current_N_ec)
    results.update(current_N_nec=current_N_nec)
    results.update(dumped_N_nec=dumped_N_nec)
    results.update(dumped_N_ec=dumped_N_ec)
    results.update(cleared_N_nec=cleared_N_nec)
    results.update(cleared_N_ec=cleared_N_ec)
