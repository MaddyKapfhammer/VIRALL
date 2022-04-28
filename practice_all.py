from black import diff
from scipy.integrate import odeint
from textx import metamodel_from_file
from textx.export import metamodel_export, model_export
import numpy as np
import matplotlib.pyplot as plt
import tempfile
import io
import base64
from virall_dsl import create_equations as create


def define_equations(compartment, t, beta, gamma, mu, vacc, specified_compartments, population):
    differential_list = create.run_differentials(compartment, t, beta, gamma, mu, vacc, specified_compartments, population)

    if len(differential_list) == 1:
        return differential_list[0]
    if len(differential_list) == 2:
        return differential_list[0], differential_list[1]
    if len(differential_list) == 3:
        return differential_list[0], differential_list[1], differential_list[2]
    if len(differential_list) == 4:
        return differential_list[0], differential_list[1], differential_list[2], differential_list[3]
    if len(differential_list) == 5:
        return differential_list[0], differential_list[1], differential_list[2], differential_list[3], differential_list[4]
    if len(differential_list) == 6:
        return differential_list[0], differential_list[1], differential_list[2], differential_list[3], differential_list[4], differential_list[5]



def run_differentials(time, compartment_list, beta, gamma, mu, vacc, specified_compartments, population):
    differential_values = []
    t = np.linspace(0, time, 50)
    y0 = tuple(compartment_list)
    print(y0)
    ret = odeint(define_equations, y0, t, args=(beta,gamma, mu, vacc, specified_compartments, population))
    print(ret.T)
    print(specified_compartments)
    if len(compartment_list) == 1:
        S = ret.T
        differential_values.append(t)
        differential_values.append(S)
    if len(compartment_list) == 2:
        S, I = ret.T
        differential_values.append(t)
        differential_values.append(S)
        differential_values.append(I)
    if len(compartment_list) == 3:
        S, I, R = ret.T
        differential_values.append(t)
        differential_values.append(S)
        differential_values.append(I)
        differential_values.append(R)
    if len(compartment_list) == 4:
        S, I, R, other = ret.T
        differential_values.append(t)
        differential_values.append(S)
        differential_values.append(I)
        differential_values.append(R)
        differential_values.append(other)
    if len(compartment_list) == 5:
        S, I, R, other, other2 = ret.T
        differential_values.append(t)
        differential_values.append(S)
        differential_values.append(I)
        differential_values.append(R)
        differential_values.append(other)
        differential_values.append(other2)
    if len(compartment_list) == 6:
        S, I, R, other, other2, other3 = ret.T
        differential_values.append(t)
        differential_values.append(S)
        differential_values.append(I)
        differential_values.append(R)
        differential_values.append(other)
        differential_values.append(other2)
        differential_values.append(other3)
    
    return differential_values


def plotsir(differential_values, x_label, y_label, color_list, compartment_list):
    f, ax = plt.subplots(1, 1, figsize=(10, 4))
    print(color_list)
    for x in range(0, (len(compartment_list))):
        t = differential_values[0]
        color = color_list[x]
        compartment = compartment_list[x]
        ax.plot(t, differential_values[x+1], color, alpha=0.7, linewidth=2, label=compartment)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which="major", c="w", lw=2, ls="-")
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ("top", "right", "bottom", "left"):
        ax.spines[spine].set_visible(False)
    plt.show()

if __name__ == "__main__":
    specified_compartments = create.determine_compartments(True, True, True, True, False, False)

    # define_equations(2, 50, 0.1, 0.1, 0.1, 0.1, specified_compartments)
    compartment_list = [99, 1, 0, 0]
    beta = 0.9
    gamma = 0.2
    mu = 0.5
    vacc = 0.03

    x_label = "label"
    y_label = "label"
    population = 100
    color_list = ["yellow", "red", "blue", "green"]
    differential_values = run_differentials(50, compartment_list, beta, gamma, mu, vacc, specified_compartments, population)
    
    plotsir(differential_values, x_label, y_label, color_list, specified_compartments)
        
