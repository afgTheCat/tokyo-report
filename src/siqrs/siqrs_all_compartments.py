"""
Simulating siqrs according to model described in chapter two of
Epidemic-logistics Modeling: A New Perspective on Operations Research by Ming
Liu, Jie Cao, Jing Liang and MingJun Chen
"""

import matplotlib.pyplot as plt
from jitcdde import y, t, jitcdde
import numpy as np


def display_endemic_progress():
    """Display how a typical endemic progresses in time"""
    beta = 0.000001
    k = 6
    gamma = 0.0002
    delta = 0.3
    # pylint: disable=invalid-name
    mu = 0.1
    d_one = 0.005
    d_two = 0.001
    tau = 5
    initially_susceptible = 99980
    initially_infected = 20

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

    stoptime = 180
    numpoints = 18000
    times = delayed_differential_equation.t + \
        np.linspace(1, stoptime, numpoints)

    data = []
    for time in times:
        data.append(delayed_differential_equation.integrate(time))

    npdata = np.array(data)
    susceptible = npdata[:, 0]
    exposed = npdata[:, 1]
    infected = npdata[:, 2]
    quarantined = npdata[:, 3]
    recovered = npdata[:, 4]

    fig = plt.figure()
    subplot = fig.add_subplot(1, 1, 1)
    subplot.plot(times, susceptible, label='Susceptible people')
    subplot.plot(times, exposed, '-', label='Exposed people')
    subplot.plot(times, infected, '--', label='Infected people')
    subplot.plot(times, quarantined, '-.', label='Quarantined people')
    subplot.plot(times, recovered, ':', label='Recovered people')
    subplot.set_xlabel('Days from outbreak')
    subplot.set_ylabel('Population')
    subplot.legend(loc='best')
    plt.show()


if __name__ == "__main__":
    display_endemic_progress()
