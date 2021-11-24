# 'ms_stress_function.py' is a script to determine an equation for oxidative stress in oligodendrocytes
# This script/equation will be used in a Multiple Sclerosis compartmental model of protein shedding kinetics

# Author: Corey Chitwood
# v0.1 - Oct 2021


import numpy as np


def oxidative_stress_equation(parameters, t, stress, extension):
    # Inputs needed: parameters, stresst0, current stress level,

    length = parameters['stress_length']
    n = parameters['n']
    stress_type = parameters['stress_type']
    sigmoidal_steepness = parameters['sigmoidal_steepness']

    x = np.linspace(t, t + length + extension, num=(length + 1 + extension))

    offset = stress
    # stress curve shifted up or down based on starting stress

    if stress_type == 1:  # sigmoidal
        a = n - 1
        numerator = 100 / a
        other_term = np.zeros(a)

        os_eqn = np.zeros(length + 1 + extension)

        for i in range(len(other_term)):
            j = i + 1
            other_term[i] = -(j * (length / n) + t)

        for j in range(len(os_eqn)):
            for i in range(len(other_term)):
                os_eqn[j] += (numerator / (1 + np.exp(-1 * sigmoidal_steepness * (x[j] + other_term[i]))))

        os_eqn = os_eqn + offset

        for i in range(len(os_eqn)):
            if os_eqn[i] > 100:
                os_eqn[i] = 100  # can't go above 100 stress
            elif os_eqn[i] < 0:
                os_eqn[i] = 0  # can't go under 0 stress

    elif stress_type == 2:  # linear

        os_eqn = 100 / length * (x - t) + stress

        for i in range(len(os_eqn)):
            if os_eqn[i] < 0:
                os_eqn[i] = 0  # can't go below 0 stress

    elif stress_type == 3:  # exponential

        os_eqn = (np.exp(0.05 * x))
        os_eqn = (100 / max(os_eqn)) * os_eqn

    elif stress_type == 4:  # logarithmic
        #         os_eqn = (np.log(x + 1))
        #         os_eqn = (100 / max(os_eqn)) * os_eqn
        os_eqn = (np.log(x + 1 - t))  # + t + length +
        # the + 1 adjusts the intercept of the ln() curve such that the stress level goes to 0 at exactly t + length
        os_eqn = (100 / max(os_eqn)) * os_eqn

    else:
        print('Invalid stress type entered')
        os_eqn = 'error'

    return os_eqn


def remission_equation(parameters, t, stress, time_initially_stressed):
    # Inputs needed: parameters, remissiont0, current_stress_level, and time initially stressed
    # Thus can- plot from current stress level to 0, starting at x = current time and ending before max iter
    # Shift in

    n = parameters['n']
    stress_type = parameters['stress_type']
    stress_length = parameters['stress_length']
    remission_type = parameters['remission_type']
    length = parameters['remission_length']
    sigmoidal_steepness = parameters['sigmoidal_steepness']

    other_term = np.zeros(n - 1)
    for i in range(len(other_term)):
        j = i + 1
        other_term[i] = (j * (length / n)) + t
    # this if statement adjusts the remission function to account for groups of cells stressed of different times having
    # the same final stress simulation value in sigmoidal stress functions
    # before this fix, all cells ending with the same final stress value would have the same remission function
    # and the compartments vs time plots would have unrealistic unit step like jumps as the large groups moved between
    # compartments at the same rate
    if stress_type == 1:
        if (time_initially_stressed > 0) and (time_initially_stressed < (other_term[0] - (3 / sigmoidal_steepness))):
            minimum = int(1)
            maximum = int(other_term[0] - (3 / sigmoidal_steepness))
            diff = int(maximum - minimum + 1)
            x = np.linspace(minimum, maximum, num=diff).astype(int)
            adjustment = np.where(x == time_initially_stressed)[0]
            # adjustment = x.index(time_initially_stressed)
        elif (time_initially_stressed > (other_term[0] + 3 / sigmoidal_steepness)) and (
                    time_initially_stressed < (other_term[1] - 3 / sigmoidal_steepness)):
            minimum = int(other_term[0] + (3 / sigmoidal_steepness))
            maximum = int(other_term[1] - (3 / sigmoidal_steepness))
            diff = int(maximum - minimum + 1)
            x = np.linspace(minimum, maximum, num=diff).astype(int)
            # adjustment = x.index(time_initially_stressed)
            adjustment = np.where(x == time_initially_stressed)[0]
        elif (time_initially_stressed > (other_term[1] + 3 / sigmoidal_steepness)) and (
                    time_initially_stressed < (other_term[2] - 3 / sigmoidal_steepness)):
            minimum = int(other_term[1] + (3 / sigmoidal_steepness))
            maximum = int(other_term[2] - (3 / sigmoidal_steepness))
            diff = int(maximum - minimum + 1)
            x = np.linspace(minimum, maximum, num=diff).astype(int)
            # adjustment = x.index(time_initially_stressed)
            adjustment = np.where(x == time_initially_stressed)[0]
        elif (time_initially_stressed > (other_term[2] + 3 / sigmoidal_steepness)) and (
                    time_initially_stressed < (other_term[3] - 3 / sigmoidal_steepness)):
            minimum = int(other_term[2] + (3 / sigmoidal_steepness))
            maximum = int(other_term[3] - (3 / sigmoidal_steepness))
            diff = int(maximum - minimum + 1)
            x = np.linspace(minimum, maximum, num=diff).astype(int)
            # adjustment = x.index(time_initially_stressed)
            adjustment = np.where(x == time_initially_stressed)[0]
        elif (time_initially_stressed > (
                    other_term[3] + 3 / sigmoidal_steepness)) and time_initially_stressed < stress_length:
            minimum = int(other_term[3] + (3 / sigmoidal_steepness))
            maximum = int(100)
            diff = int(maximum - minimum + 1)
            x = np.linspace(minimum, maximum, num=diff).astype(int)
            # adjustment = x.index(time_initially_stressed)
            adjustment = np.where(x == time_initially_stressed)[0]
        else:
            adjustment = 0
    else:
        adjustment = 0
    x = np.linspace(t + 1, t + length, num=length)

    offset = 100 - stress
    # sigmoidal curve shifted down rather than such that the rate of decrease of stress is the same as the increase
    # in stress

    if remission_type == 1:  # sigmoidal
        a = n - 1
        numerator = 100 / a

        rem_eqn = np.zeros(length)

        for j in range(len(rem_eqn)):
            for i in range(len(other_term)):
                rem_eqn[j] += (numerator / (1 + np.exp(sigmoidal_steepness * (x[j] - other_term[i] + adjustment))))

        rem_eqn = rem_eqn - offset

        for i in range(len(rem_eqn)):
            if rem_eqn[i] < 0:
                rem_eqn[i] = 0  # can't go below 0 stress

    elif remission_type == 2:  # linear

        rem_eqn = -stress / length * (x - t) + stress

        for i in range(len(rem_eqn)):
            if rem_eqn[i] < 0:
                rem_eqn[i] = 0  # can't go below 0 stress

    elif remission_type == 3:  # exponential

        rem_eqn = (np.exp(-0.05 * x))
        rem_eqn = (stress / max(rem_eqn)) * rem_eqn

    elif remission_type == 4:  # logarithmic

        rem_eqn = (np.log(-x + t + length + 1))
        # the + 1 adjusts the intercept of the ln() curve such that the stress level goes to 0 at exactly t + length
        rem_eqn = (stress / max(rem_eqn)) * rem_eqn

    else:
        print('Invalid stress type entered')
        rem_eqn = 'error'

    return rem_eqn
