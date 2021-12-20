# 'new_ms_plotter.py' is a script to create plots associated with simulation data
# This script/equation will be used in a Multiple Sclerosis compartmental model of protein shedding kinetics

# Author: Corey Chitwood
# v0.1 - Oct 2021


import matplotlib.pyplot as plt
import numpy as np
import os
import sys


def plot_stress_equation(os_eqn, n, stress_type, t, length):
    x = np.linspace(t, t + length, num=(length + 1))

    if stress_type == 1:
        stress_str = 'sigmoidal'
    elif stress_type == 2:
        stress_str = 'linear'
    elif stress_type == 3:
        stress_str = 'exponential'
    elif stress_type == 4:
        stress_str = 'logarithmic'
    else:
        print('Invalid stress type entered')
        return

    plt.plot(x, os_eqn)  # os_eqn = oxidative stress equation

    plt.xlabel('Time')
    plt.ylabel('Oxidative Stress')

    title_str = stress_str + ' Model of Oxidative Stress'
    plt.title(title_str.upper())

    # plot horizontal lines corresponding to division between compartments

    for i in range(n + 1):
        y = i * (100 / n)
        plt.axhline(y=y, color='r', linestyle='dotted')

    plt.show()


def plot_remission_equation(rem_eqn, n, remission_type, t, length):
    x = np.linspace(t + 1, t + length, num=length)

    if remission_type == 1:
        remission_str = 'sigmoidal'
    elif remission_type == 2:
        remission_str = 'linear'
    elif remission_type == 3:
        remission_str = 'exponential'
    elif remission_type == 4:
        remission_str = 'logarithmic'
    else:
        print('Invalid remission type entered')
        return


    plt.plot(x, rem_eqn)

    plt.xlabel('Time')
    plt.ylabel('Oxidative Stress')

    rem_title_str = remission_str + ' Model of Oxidative Stress Remission'
    plt.title(rem_title_str.upper())

    # plot horizontal lines corresponding to division between compartments

    for i in range(n + 1):
        y = i * (100 / n)
        plt.axhline(y=y, color='r', linestyle='dotted')

    plt.show()


def plot_overall_stress(overall_stress, stress_type, remission_type, n, t0, tmax):
    x = np.linspace(t0, tmax, num=(tmax + 1))

    if stress_type == 1:
        stress_str = 'sigmoidal'
    elif stress_type == 2:
        stress_str = 'linear'
    elif stress_type == 3:
        stress_str = 'exponential'
    elif stress_type == 4:
        stress_str = 'logarithmic'
    else:
        print('Invalid stress type entered')
        return

    if remission_type == 1:
        remission_str = 'sigmoidal'
    elif remission_type == 2:
        remission_str = 'linear'
    elif remission_type == 3:
        remission_str = 'exponential'
    elif remission_type == 4:
        remission_str = 'logarithmic'
    else:
        print('Invalid remission type entered')
        return

    plt.plot(x, overall_stress)

    plt.xlabel('Time')
    plt.ylabel('Oxidative Stress')

    if stress_str == remission_str:
        title_type = stress_str
    else:
        title_type = stress_str + ' & ' + remission_str

    overall_title_str = 'overall ' + title_type + ' Model of Oxidative Stress'
    plt.title(overall_title_str.upper())

    # plot horizontal lines corresponding to division between compartments

    for i in range(n + 1):
        y = i * (100 / n)
        plt.axhline(y=y, color='r', linestyle='dotted')
    stripped_title = overall_title_str.replace(' ', '') + '.png'
    plt.savefig(stripped_title, bbox_inches='tight')

    plt.show()


def plot_compartments_vs_time(parameters, results):
    cwd = 'C:/Users/Corey/Desktop/MS Model/'
    folder = os.path.join(cwd, 'compartments')

    try:
        os.mkdir(folder)
    except:
        pass

    # oligodendrocyte plot
    plt.figure(0)

    compartments = results['compartments']

    n = parameters['n']
    q = np.shape(compartments)[0]
    x = np.linspace(0, q - 1, num=q)
    plt.suptitle('Oligodendrocyte Compartments')

    for i in range(0, n):
        label_text = 'Compartment #' + str(i + 1)
        plt.subplot(5, 1, i+1)
        plt.plot(x, compartments[:, i], label=label_text)
        plt.legend()

    plt.savefig('compartments/oligodendrocyte_compartments_vs_time.png', bbox_inches='tight')
    plt.show()

    # neuron plot
    plt.figure(1)

    N_compartments = results['N_compartments']

    q = np.shape(N_compartments)[0]
    x = np.linspace(0, q - 1, num=q)
    plt.suptitle('Neuron Compartments')

    for i in range(0, n):
        label_text = 'Compartment #' + str(i + 1)
        plt.subplot(5, 1, i+1)
        plt.plot(x, N_compartments[:, i], label=label_text)
        plt.legend()

    plt.savefig('compartments/neuron_compartments_vs_time.png', bbox_inches='tight')
    plt.show()


def plot_alive_dead_vs_time(parameters, results):
    cwd = 'C:/Users/Corey/Desktop/MS Model/'
    folder = os.path.join(cwd, 'alive vs dead')
    try:
        os.mkdir(folder)
    except:
        pass

    n = parameters['n']

    plt.figure(0)
    total_L_alive = results['total_L_alive']
    total_L_dead = results['total_L_dead']

    if len(total_L_dead) != len(total_L_alive):
        print('Error: Olg. Total Alive and Dead Array Sizes do not match')
        return

    x = np.linspace(1, len(total_L_dead), num=(len(total_L_dead)))
    plt.plot(x, total_L_dead, label='Total Oligodendrocytes Dead')
    plt.plot(x, total_L_alive, label='Total Oligodendrocytes Alive')
    plt.legend()
    plt.title('Alive and Dead Oligodendrocytes vs Time')
    plt.savefig('alive vs dead/oligodendrocytes_vs_time.png', bbox_inches='tight')
    plt.show()

    plt.figure(1)
    total_N_alive = results['total_N_alive']
    total_N_dead = results['total_N_dead']

    if len(total_N_dead) != len(total_N_alive):
        print('Error: Neuron Total Alive and Dead Array Sizes do not match')
        return

    x = np.linspace(1, len(total_N_dead), num=(len(total_N_dead)))
    plt.plot(x, total_N_dead, label='Total Neurons Dead')
    plt.plot(x, total_N_alive, label='Total Neurons Alive')
    plt.legend()
    plt.title('Alive and Dead Neurons vs Time')
    plt.savefig('alive vs dead/neurons_vs_time.png', bbox_inches='tight')
    plt.show()

    plt.figure(2)
    IC_N = parameters['IC_N']
    IC_L = parameters['IC_L']
    L_deaths = results['L_deaths'] / IC_L
    L_births = results['L_births'] / IC_L
    N_deaths = results['N_deaths'] / IC_N

    if len(L_deaths) != len(L_births) != len(N_deaths):
        print('Error: Olg. Birth and Death Array Sizes do not match')

    x = np.linspace(1, len(L_births), num=(len(L_births)))
    plt.plot(x, L_births, label='Oligodendrocyte Births')
    plt.plot(x, L_deaths, label='Oligodendrocyte Deaths')
    plt.plot(x, N_deaths, label='Neuron Deaths')
    plt.legend()
    plt.title('Number of Births and Deaths at Each Time')
    plt.savefig('alive vs dead/births_deaths_vs_time.png', bbox_inches='tight')
    plt.show()


def plot_stressed_unstressed_vs_time(parameters, results):
    cwd = 'C:/Users/Corey/Desktop/MS Model/'
    folder = os.path.join(cwd, 'stressed unstressed')
    try:
        os.mkdir(folder)
    except:
        pass

    n = parameters['n']

    plt.figure(0)
    total_L_stressed = results['total_L_stressed']
    total_L_unstressed = results['total_L_unstressed']

    if len(total_L_stressed) != len(total_L_unstressed):
        print('Olg Stressed/Unstressed array sizes do not match')
        return

    x = np.linspace(1, len(total_L_unstressed), num=(len(total_L_unstressed)))
    plt.plot(x, total_L_unstressed, label='Oligodendrocytes Unstressed')
    plt.plot(x, total_L_stressed, label='Oligodendrocytes Stressed')
    plt.legend()
    plt.title('Total Number of Oligodendrocytes Stressed and Unstressed vs Time')
    plt.savefig('stressed unstressed/oligodendrocytes_total_stressed_unstressed_vs_time.png', bbox_inches='tight')
    plt.show()

    plt.figure(1)
    L_stressed = results['L_stressed']
    L_initially_stressed_only = results['L_initially_stressed_only']

    if len(L_stressed) != len(L_initially_stressed_only):
        print('Olg Stressed vs Initially Stressed array sizes do not match')
        return

    x = np.linspace(1, len(L_stressed), num=(len(L_stressed)))
    plt.plot(x, L_stressed, label='Oligodendrocytes Stressed Remaining')
    plt.plot(x, L_initially_stressed_only, label='Oligodendrocytes Initially Stressed at each Time')
    plt.legend()
    plt.title('Initially Stressed Olg. Compared to those Olg. Remaining Alive after Simulation')
    plt.savefig('stressed unstressed/oligodendrocytes_initially_stressed_vs_still_alive.png', bbox_inches='tight')
    plt.show()

    plt.figure(2)
    total_N_stressed = results['total_N_stressed']
    total_N_unstressed = results['total_N_unstressed']

    if len(total_N_stressed) != len(total_N_unstressed):
        print('Neuron Stressed/Unstressed array sizes do not match')
        return

    x = np.linspace(1, len(total_L_unstressed), num=(len(total_L_unstressed)))
    plt.plot(x, total_N_unstressed, label='Neurons Unstressed')
    plt.plot(x, total_N_stressed, label='Neurons Stressed')
    plt.legend()
    plt.title('Total Number of Neurons Stressed and Unstressed vs Time')
    plt.savefig('stressed unstressed/neurons_total_stressed_unstressed_vs_time.png', bbox_inches='tight')
    plt.show()

    plt.figure(3)
    N_stressed = results['N_stressed']
    N_initially_stressed_only = results['N_initially_stressed_only']

    if len(L_stressed) != len(L_initially_stressed_only):
        print('Neuron Stressed vs Initially Stressed array sizes do not match')
        return

    x = np.linspace(1, len(L_stressed), num=(len(L_stressed)))
    plt.plot(x, N_stressed, label='Neurons Stressed Remaining')
    plt.plot(x, N_initially_stressed_only, label='Neurons Initially Stressed at each Time')
    plt.legend()
    plt.title('Initially Stressed Neurons Compared to those Neurons Remaining Alive after Simulation')
    plt.savefig('stressed unstressed/Neurons_initially_stressed_vs_still_alive.png', bbox_inches='tight')
    plt.show()


def plot_protein_shedding(parameters, results):
    cwd = 'C:/Users/Corey/Desktop/MS Model/'
    folder = os.path.join(cwd, 'protein shedding analyses')
    try:
        os.mkdir(folder)
    except:
        pass

    maxIter = parameters['maxIter']

    current_L_nec = results['current_L_nec']
    current_N_nec = results['current_N_nec']
    current_L_ec = results['current_L_ec']
    current_N_ec = results['current_N_ec']
    dumped_L_nec = results['current_L_nec']
    dumped_N_nec = results['current_N_nec']
    dumped_L_ec = results['current_L_ec']
    dumped_N_ec = results['current_L_ec']
    cleared_L_nec = results['cleared_L_nec']
    cleared_N_nec = results['cleared_N_nec']
    cleared_L_ec = results['cleared_L_ec']
    cleared_N_ec = results['cleared_N_ec']

    total_nec_prot = current_N_nec + current_L_nec
    total_ec_prot = current_N_ec + current_L_ec
    total_L_prot = current_L_nec + current_L_ec
    total_N_prot = current_N_nec + current_N_ec
    total_prot = total_N_prot + total_L_prot
    total_dumped_ec = dumped_L_ec + dumped_N_ec
    total_dumped_nec = dumped_L_nec + dumped_N_nec
    total_dumped_L = dumped_L_nec + dumped_N_ec
    total_dumped_N = dumped_N_nec + dumped_N_ec
    total_dumped = total_dumped_N + total_dumped_L
    total_cleared_ec = cleared_L_ec + cleared_N_ec
    total_cleared_nec = cleared_L_nec + cleared_N_nec
    total_cleared_L = cleared_L_ec + cleared_L_nec
    total_cleared_N = cleared_N_ec + cleared_N_nec
    total_cleared = total_cleared_N + total_cleared_L



    x = np.linspace(1, maxIter+1, num=maxIter+1)

    # EC vs NEC vs total in plasma vs time
    plt.figure(0)
    plt.subplot(2, 1, 1)
    plt.plot(x, total_prot, label='Total Protein')
    plt.legend()
    plt.subplot(2, 1, 2)
    plt.plot(x, total_nec_prot, label='NEC Protein')
    plt.plot(x, total_ec_prot, label='EC Protein')
    plt.legend()
    plt.suptitle('Type of Protein in Plasma vs Time')
    plt.savefig('protein shedding analyses/NEC_EC_Total_vs_time.png', bbox_inches='tight')
    plt.show()

    results.update(total_nec_prot=total_nec_prot)
    results.update(total_ec_prot=total_ec_prot)
    results.update(total_prot=total_prot)

    # L vs N vs total in plasma vs time
    plt.figure(1)
    plt.subplot(2, 1, 1)
    plt.plot(x, total_prot, label='Total Protein')
    plt.legend()
    plt.subplot(2, 1, 2)
    plt.plot(x, total_L_prot, label='Oligodendrocyte Protein')
    plt.plot(x, total_N_prot, label='Neuron Protein')
    plt.legend()
    plt.suptitle('Source of Protein in plasma vs Time')
    plt.savefig('protein shedding analyses/L_N_Total_vs_time.png', bbox_inches='tight')
    plt.show()

    results.update(total_L_prot=total_L_prot)
    results.update(total_N_prot=total_N_prot)

    # L only
    # plot doesn't yield much valuable data

    # plt.figure(2)
    # plt.subplot(2, 1, 1)
    # plt.plot(x, total_L_prot, label='Total Oligodendrocyte Protein')
    # plt.legend()
    # plt.subplot(2, 1, 2)
    # plt.plot(x, current_L_ec, label='Oligodendrocyte EC Protein')
    # plt.plot(x, current_L_nec, label='Oligodendrocyte NEC Protein')
    # plt.legend()
    # plt.suptitle('Type of Oligodendrocyte Protein vs Time')
    # plt.savefig('protein shedding analyses/LEC_LNEC_vs_time.png', bbox_inches='tight')
    # plt.show()

    # N only
    # plot doesn't yield much valuable data

    # plt.figure(3)
    # plt.subplot(2, 1, 1)
    # plt.plot(x, total_N_prot, label='Total Neuron Protein')
    # plt.legend()
    # plt.subplot(2, 1, 2)
    # plt.plot(x, current_N_ec, label='Neuron EC Protein')
    # plt.plot(x, current_N_nec, label='Neuron NEC Protein')
    # plt.legend()
    # plt.suptitle('Type of Neuron Protein vs time')
    # plt.savefig('protein shedding analyses/NEC_NNEC_vs_time.png', bbox_inches='tight')
    # plt.show()

    # EC only
    plt.figure(4)
    plt.subplot(2, 1, 1)
    plt.plot(x, total_ec_prot, label='Total EC Protein')
    plt.legend()
    plt.subplot(2, 1, 2)
    plt.plot(x, current_N_ec, label='Neuron EC Protein')
    plt.plot(x, current_L_ec, label='Oligodendrocyte EC Protein')
    plt.legend()
    plt.suptitle('Source of EC Protein vs Time')
    plt.savefig('protein shedding analyses/NEC_LEC_vs_time.png', bbox_inches='tight')
    plt.show()

    # NEC only
    plt.figure(5)
    plt.subplot(2, 1, 1)
    plt.plot(x, total_nec_prot, label='Total NEC Protein')
    plt.legend()
    plt.subplot(2, 1, 2)
    plt.plot(x, current_N_nec, label='Neuron NEC Protein')
    plt.plot(x, current_L_nec, label='Oligodendrocyte NEC Protein')
    plt.legend()
    plt.suptitle('Source of NEC Protein vs Time')
    plt.savefig('protein shedding analyses/NNEC_LNEC_vs_time.png', bbox_inches='tight')
    plt.show()



# protein shedding plots
# *** need total protein in bloodstream/time
# *** need EC vs NEC protein/time
# comparison to normal shedding rates? Assume Uh for all cells in model and only run as EC shedding?


# do we want a sensitivity analysis for each parameter?'''

