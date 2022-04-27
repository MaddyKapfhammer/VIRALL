from scipy.integrate import odeint
from textx import metamodel_from_file
from textx.export import metamodel_export, model_export
import numpy as np
import matplotlib.pyplot as plt
import tempfile
import io
import base64


def define_equations(compartment, t, beta, gamma):
    S, I, R = compartment
    dS = -(S * I * beta)
    dI = (S * I * beta) - (I * gamma)
    dR = (I * R * gamma)

    return dS, dI, dR


def run_differentials(time, S0, I0, R0, beta, gamma):
    t = np.linspace(0, time, 50)
    y0 = S0, I0, R0
    ret = odeint(define_equations, y0, t, args=(beta,gamma))
    S, I, R = ret.T

    return t, S, I, R


def plotsir(t, S, I, R, x_label, y_label, susceptible_color, infected_color, recovered_color):
    f, ax = plt.subplots(1, 1, figsize=(10, 4))
    ax.plot(t, S, susceptible_color, alpha=0.7, linewidth=2, label="Susceptible")
    ax.plot(t, I, infected_color, alpha=0.7, linewidth=2, label="Infected")
    ax.plot(t, R, recovered_color, alpha=0.7, linewidth=2, label="Recovered")

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which="major", c="w", lw=2, ls="-")
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ("top", "right", "bottom", "left"):
        ax.spines[spine].set_visible(False)


class Virall(object):

    def __init__(self):
        self.susceptible = 1000
        self.infected = 1
        self.recovered = 0
        self.time_frame = 50
        self.beta = 0.1
        self.gamma = 0.2
        self.mu = 0.3

    def interpret(self, virall_model):
        for c in virall_model.commands:
            if c.__class__.__name__ == "CreateModel":
                print("Creating the Model!")
                # This will contain the code for generating equations once finished
            elif c.__class__.__name__ == "Calculate":
                print("Differentiating")

                compartment_list = c.compartment.x
                susceptible = compartment_list[0]
                infected = compartment_list[1]
                recovered = compartment_list[2]
                
                time_frame = c.time_frame.x

                rate_list = c.disease_rates.x
                beta = rate_list[0]
                gamma = rate_list[1]

                mu_value = rate_list[2]

                t, S, I, R = run_differentials(time_frame, susceptible, infected, recovered, beta, gamma)
             

            elif c.__class__.__name__ == "Plot":
                print("Plot")
                attributes = dir(c)
                # print(attributes)
                compartment_colors = c.compartment_color.compartment_color

                x_label = c.x_label

                y_label = c.y_label

                susceptible_color = compartment_colors[0]
                infected_color = compartment_colors[1]
                recovered_color = compartment_colors[2]

                return t, S, I, R, x_label, y_label, susceptible_color, infected_color, recovered_color
        
        return t, S, I, R, x_label, y_label, susceptible_color, infected_color, recovered_color


def main(text, debug=False):

    virall_meta = metamodel_from_file("virall_dsl/virall.tx", debug=False)
    # metamodel_export(virall_meta, "virall_dsl/virall_meta.dot")

    virall_model = virall_meta.model_from_str(text)
    # model_export(virall_model, "virall_dsl/program.dot")

    virall = Virall()

    t, S, I, R, x_label, y_label, susceptible_color, infected_color, recovered_color = virall.interpret(virall_model)

    return t, S, I, R, x_label, y_label, susceptible_color, infected_color, recovered_color
    

if __name__ == "__main__":
    text = "hello"
    file = text_to_file(text)
    main(file)
