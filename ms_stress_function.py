# 'ms_stress_function.py' is a script to determine an equation for oxidative stress in oligodendrocytes
# This script/equation will be used in a Multiple Sclerosis compartmental model of protein shedding kinetics

# Author: Corey Chitwood
# v0.1 - Oct 2021


import numpy as np


def oxidative_stress_equation(t, length, stress, stress_type, n):
    # Inputs needed: remission type, max iter, the current stress level, and current time
    # Thus can- plot from current stress level to 0, starting at x = current time and ending before max iter

    x = np.linspace(t, t + length, num=(length + 1))

    offset = stress
    # sigmoidal curve shifted down rather than such that the rate of decrease of stress is the same as the increase
    # in stress

    if stress_type == 1: # sigmoidal
        a = n - 1
        numerator = 100 / a
        other_term = np.zeros(a)

        os_eqn = np.zeros(length + 1)

        for i in range(len(other_term)):
            j = i + 1
            other_term[i] = -(j * (length / n) + t)

        for j in range(len(os_eqn)):
            for i in range(len(other_term)):
                os_eqn[j] += (numerator / (1 + np.exp(-(x[j] + other_term[i]))))

        os_eqn = os_eqn + offset

        for i in range(len(os_eqn)):
            if os_eqn[i] > 100:
                os_eqn[i] = 100  # can't go above 100 stress
            elif os_eqn[i] < 0:
                os_eqn[i] = 0  # can't go under 0 stress

    elif stress_type == 2: # linear

        os_eqn = 100 / length * (x - t) + stress

        for i in range(len(os_eqn)):
            if os_eqn[i] < 0:
                os_eqn[i] = 0  # can't go below 0 stress

    elif stress_type == 3: # exponential

        os_eqn = (np.exp(0.05 * x))
        os_eqn = (100 / max(os_eqn)) * os_eqn

    elif stress_type == 4: # logarithmic
        #         os_eqn = (np.log(x + 1))
        #         os_eqn = (100 / max(os_eqn)) * os_eqn
        os_eqn = (np.log(x + 1 - t))  # + t + length +
        # the + 1 adjusts the intercept of the ln() curve such that the stress level goes to 0 at exactly t + length
        os_eqn = (100 / max(os_eqn)) * os_eqn

    else:
        print('Invalid stress type entered')
        os_eqn = 'error'

    return os_eqn


def remission_equation(t, length, stress, remission_type, n):
    # Inputs needed: remissiont0, remission_length, current_stress_level, remission_type, and n
    # Thus can- plot from current stress level to 0, starting at x = current time and ending before max iter

    x = np.linspace(t + 1, t + length, num=length)

    offset = 100 - stress
    # sigmoidal curve shifted down rather than such that the rate of decrease of stress is the same as the increase
    # in stress

    if remission_type == 1: # sigmoidal

        a = n - 1
        numerator = 100 / a
        other_term = np.zeros(a)

        rem_eqn = np.zeros(length)

        for i in range(len(other_term)):
            j = i + 1
            # other_term[i] = -(j * ((t + length) / (a + 1)))
            other_term[i] = -(j * (t / (a + 1))) - t

        for j in range(len(rem_eqn)):
            for i in range(len(other_term)):
                rem_eqn[j] += (numerator / (1 + np.exp((x[j] + other_term[i]))))

        rem_eqn = rem_eqn - offset

        for i in range(len(rem_eqn)):
            if rem_eqn[i] < 0:
                rem_eqn[i] = 0  # can't go below 0 stress

    elif remission_type == 2: # linear

        rem_eqn = -stress / length * (x - t) + stress

        for i in range(len(rem_eqn)):
            if rem_eqn[i] < 0:
                rem_eqn[i] = 0  # can't go below 0 stress

    elif remission_type == 3: # exponential

        rem_eqn = (np.exp(-0.05 * x))
        rem_eqn = (stress / max(rem_eqn)) * rem_eqn

    elif remission_type == 4: # logarithmic

        rem_eqn = (np.log(-x + t + length + 1))
        # the + 1 adjusts the intercept of the ln() curve such that the stress level goes to 0 at exactly t + length
        rem_eqn = (stress / max(rem_eqn)) * rem_eqn

    else:
        print('Invalid stress type entered')
        rem_eqn = 'error'

    return rem_eqn
