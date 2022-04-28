"""A practice program with SIR modeling in Python."""

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
# matplotlib inline
# import mpld3
# mpld3.enable_notebook()

def deriv(compartment, t, N, beta, gamma, host, vector, mu):
    S, I, R, Vect = compartment
    dS = -(S * Vect * beta)
    dI = (S * Vect * beta)
    dR = (I * gamma) - (R * mu)
    dVect = (host*vector) - (host * mu)
    return dS, dI, dR, dVect


def specify_compartments(N, beta, host, gamma, vector, mu, S0, I0, R0, Vec):
    t = np.linspace(0, 49, 50)
    y0 = S0, I0, R0, Vec
    ret = odeint(deriv, y0, t, args=(N, beta, gamma, host, vector, mu))
    S, I, R, Vec = ret.T

    return t, S, I, R, Vec


def plotsir(t, S, I, R, Vec):
    f, ax = plt.subplots(1,1,figsize=(10,4))
    ax.plot(t, S, 'blue', alpha=0.7, linewidth=2, label='Susceptible')
    ax.plot(t, I, 'green', alpha=0.7, linewidth=2, label='Infected')
    ax.plot(t, R, 'red', alpha=0.7, linewidth=2, label='Recovered')
    ax.plot(t, Vec, 'yellow', alpha=0.7, linewidth=2, label="Vector Population")

    ax.set_xlabel('Time')

    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which='major', c='w', lw=2, ls='-')
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)
    
    plt.show()


if __name__ == "__main__":
    N = 1000
    S0 = 9000
    I0 = 200
    R0 = 300
    Vec = 50

    beta = 4.29
    gamma = 0.121
    host = 190000
    vector_transmission = 0.4
    mu = 4.83

    t, S, I, R, Vector = specify_compartments(N, beta, gamma, host, vector_transmission, mu, S0, I0, R0, Vec)
    
    plotsir(t, S, I, R, Vector)
