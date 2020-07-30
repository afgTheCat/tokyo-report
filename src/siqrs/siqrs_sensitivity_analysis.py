"""
Simulating sensitivity analysis on beta and <k> according to model described in
chapter two of Epidemic-logistics Modeling: A New Perspective on Operations
Research by Ming Liu, Jie Cao, Jing Liang and MingJun Chen
"""

import numpy as np
import matplotlib.pyplot as plt
from jitcdde import y, t, jitcdde


def return_differential_equations(delta, k):
    """Returns the differential equations"""
    initially_susceptible = 99980
    initially_infected = 20
    beta = 0.000001
    gamma = 0.0002
    # pylint: disable=invalid-name
    mu = 0.1
    d_one = 0.005
    d_two = 0.0001
    tau = 5
    list_of_odes = [
        gamma * y(4) - beta * k * y(0) * y(2),
        beta * k * y(0) * y(2) - beta * k * y(0, t - tau) * y(2, t-tau),
        beta * k * y(0, t - tau) * y(2, t-tau) - d_one * y(2) - delta * y(2),
        delta * y(2) - d_two * y(3) - mu * y(3),
        mu * y(3) - gamma * y(4)
    ]

    delayed_differential_equation = jitcdde(list_of_odes)
    delayed_differential_equation.add_past_point(
        -10.0, [initially_susceptible, 0, initially_infected, 0, 0], [0, 0, 0, 0, 0])
    delayed_differential_equation.add_past_point(
        0, [initially_susceptible, 0, initially_infected, 0, 0], [0, 0, 0, 0, 0])

    delayed_differential_equation.step_on_discontinuities()
    return delayed_differential_equation


def display_sensitivity_analysis():
    """Display how a typical endemic progresses in time"""

    list_of_k = list(range(4, 11, 2))
    list_of_delta = [0.2, 0.3, 0.4, 0.5]

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    fig = plt.figure()
    subplot_k = fig.add_subplot(1, 2, 1)
    subplot_k.set_xlabel('Days from outbreak')
    subplot_k.set_ylabel('Population')

    subplot_delta = fig.add_subplot(1, 2, 2)
    subplot_delta.set_xlabel('Days from outbreak')
    subplot_delta.set_ylabel('Population')

    lines = ['', '-', '--', '-.', ':']

    for index, k in enumerate(list_of_k):
        delayed_differential_equation = return_differential_equations(
            0.3, k)
        stoptime = 180
        numpoints = 18000
        times = delayed_differential_equation.t + \
            np.linspace(1, stoptime, numpoints)
        data = []
        for time in times:
            data.append(delayed_differential_equation.integrate(time))
        npdata = np.array(data)
        infected = npdata[:, 2]
        subplot_k.plot(
            times, infected, lines[index], label=r'\(\langle k \rangle\)' + '={}'.format(k))

    for index, delta in enumerate(list_of_delta):
        delayed_differential_equation = return_differential_equations(delta, 6)
        stoptime = 180
        numpoints = 18000
        times = delayed_differential_equation.t + \
            np.linspace(1, stoptime, numpoints)
        data = []
        for time in times:
            data.append(delayed_differential_equation.integrate(time))
        npdata = np.array(data)
        infected = npdata[:, 2]
        subplot_delta.plot(
            times, infected, lines[index], label=r'\(\delta\)' + '={}'.format(delta))

    subplot_k.legend()
    subplot_delta.legend()
    plt.savefig('sensitivity_analysis.svg', format="svg", dpi=1200)


if __name__ == "__main__":
    display_sensitivity_analysis()
