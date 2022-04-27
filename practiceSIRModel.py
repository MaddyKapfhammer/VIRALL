"""A practice program with SIR modeling in Python."""

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
# matplotlib inline
# import mpld3
# mpld3.enable_notebook()

def deriv(compartment, t, N, beta, gamma):
    S, I, R = compartment
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I
    dRdt = gamma * I
    
    print(dSdt)
    return dSdt, dIdt, dRdt


def specify_compartments(N, beta, D, S0, I0, R0):
    gamma = 1.0 / D
    t = np.linspace(0, 49, 50)
    y0 = S0, I0, R0
    ret = odeint(deriv, y0, t, args=(N, beta, gamma))
    S, I, R = ret.T

    # print(S, I, R)
    return t, S, I, R


def plotsir(t, S, I, R):
    f, ax = plt.subplots(1,1,figsize=(10,4))
    ax.plot(t, S, 'b', alpha=0.7, linewidth=2, label='Susceptible')
    ax.plot(t, I, 'y', alpha=0.7, linewidth=2, label='Infected')
    ax.plot(t, R, 'g', alpha=0.7, linewidth=2, label='Recovered')

    ax.set_xlabel('Time (days)')

    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which='major', c='w', lw=2, ls='-')
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)
    
    plt.show()


if __name__ == "__main__":
    N = int(input("Enter size of population: "))
    beta = int(input("How many people does an infected individual infect per day?: "))
    D = int(input("How long does infection last?: "))

    S0 = int(input("Enter initial susceptible population: "))
    I0 = int(input("Enter initial infected population: "))
    R0 = int(input("Enter initial recovered population: "))

    t, S, I, R = specify_compartments(N, beta, D, S0, I0, R0)
    
    plotsir(t, S, I, R)
