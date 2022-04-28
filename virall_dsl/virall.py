from scipy.integrate import odeint
from textx import metamodel_from_file
from textx.export import metamodel_export, model_export
import numpy as np
import matplotlib.pyplot as plt
import tempfile
import io
import base64


def define_equations(compartment, t, beta, gamma, mu, vacc, vector, host, compartment_dict):
    # Compartment_dict --> dictionary with compartment as key and corresponding True or False as value  
    if compartment_dict["susceptible"] == True and compartment_dict["infected"] == True and compartment_dict["recovered"] == True and compartment_dict["exposed"] == False and compartment_dict["vector"] == False:
        S, I, R = compartment 
        dS = -(S * I * beta) - (S * vacc)
        dI = (S * I * beta) - (I * gamma) - (I * mu)
        dR = (I * gamma) - (R * mu)
        return dS, dI, dR
    if compartment_dict["susceptible"] == True and compartment_dict["infected"] == True and compartment_dict["recovered"] == True and compartment_dict["exposed"] == True and compartment_dict["vector"] == False:
        S, I, R, E = compartment
        dS = -(S * I * beta) + (S * E *beta) - (S * vacc)
        dI = (S * I * beta) - (I * gamma) + (I * E * beta) - (I * mu)
        dR = (I * gamma) - (R * mu)
        dE = (E * S * beta) + (E * I * beta) - (E * gamma) - (E * mu)
        return dS, dI, dR, dE
    if compartment_dict["susceptible"] == True and compartment_dict["infected"] == True and compartment_dict["recovered"] == True and compartment_dict["exposed"] == False and compartment_dict["vector"] == True:
        S, I, R, Vect = compartment
        dS = -(S * Vect * beta)
        dI = (S * Vect * beta)
        dR = (I * gamma) - (R * mu)
        dVect = (host*vector) - (host * mu)
        return dS, dI, dR, dVect


def run_differentials(time, compartments, compartment_dict, beta, gamma, mu, vacc, vector, host, used_compartments):
    differential_values = []
    t = np.linspace(0, time, 50)
    y0 = tuple(compartments)
    ret = odeint(define_equations, y0, t, args=(beta, gamma, mu, vacc, vector, host, compartment_dict))
    

    if "exposed" in used_compartments:
        S, I, R, E = ret.T
        differential_values.append(S)
        differential_values.append(I)
        differential_values.append(R)
        differential_values.append(E)
    elif "vector" in used_compartments:
        S, I, R, Vec = ret.T
        differential_values.append(S)
        differential_values.append(I)
        differential_values.append(R)
        differential_values.append(Vec)
    else:
        S, I, R = ret.T
        differential_values.append(S)
        differential_values.append(I)
        differential_values.append(R)

    return t, differential_values


def plotsir(differential_values, x_label, y_label, color_list, compartment_names):
    f, ax = plt.subplots(1, 1, figsize=(10, 4))
    for x in range(0, (len(compartment_names))):
        t = differential_values[0]
        color = color_list[x]
        compartment = compartment_names[x]
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


def convert_string_to_boolean(string):
    if string == "True":
        return True
    if string == "False":
        return False


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
                compartment_dict = {}
                compartment_names = c.compartment_definition.initial_compartment
                boolean_compartment_list = c.compartment_definition.boolean_value

                for x in range(0, len(compartment_names)):
                    key = compartment_names[x]
                    value = boolean_compartment_list[x]
                    compartment_dict[key] = convert_string_to_boolean(value)
        
            elif c.__class__.__name__ == "Calculate":
                used_compartments = c.compartment.initial_compartment
                compartment_values = c.compartment.x


                characteristic_list = c.disease_rates.disease_characteristic
                rate_list = c.disease_rates.x

                beta = 1
                gamma = 1
                vacc = 1
                mu = 1
                vector = 1
                host = 0

                for x in range(0, len(characteristic_list)):
                    if characteristic_list[x] == "transmission":
                        beta = rate_list[x]
                    if characteristic_list[x] == "recovery":
                        gamma = rate_list[x]
                    if characteristic_list[x] == "vaccination":
                        vacc = rate_list[x]
                    if characteristic_list[x] == "death":
                        mu = rate_list[x]
                    if characteristic_list[x] == "vectransmission":
                        vector = rate_list[x]
                    if characteristic_list[x] == "host":
                        host = rate_list[x]

                time_frame = c.time_frame.x
                population = c.population_size.x

                t, differential_values = run_differentials(time_frame, compartment_values, compartment_dict, beta, gamma, mu, vacc, vector, host, used_compartments)             
                print(len(differential_values))

            elif c.__class__.__name__ == "Plot":
                compartment_colors = c.compartment_color.compartment_color

                x_label = c.x_label

                y_label = c.y_label


        return t, differential_values, x_label, y_label, compartment_colors, compartment_names, used_compartments

def main(text, debug=False):

    virall_meta = metamodel_from_file("virall_dsl/virall.tx", debug=False)
    # metamodel_export(virall_meta, "virall_dsl/virall_meta.dot")

    virall_model = virall_meta.model_from_str(text)
    # model_export(virall_model, "virall_dsl/program.dot")

    virall = Virall()
        
    t, differential_values, x_label, y_label, compartment_colors, compartment_names, used_compartments = virall.interpret(virall_model)

    return t, differential_values, x_label, y_label, compartment_colors, compartment_names, used_compartments
    

if __name__ == "__main__":
    text = "hello"
    # t, diifferential_values, 