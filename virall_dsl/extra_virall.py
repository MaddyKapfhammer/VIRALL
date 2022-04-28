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
    ret = odeint(define_equations, y0, t, args=(beta,gamma, mu, vacc, specified_compartments, population), mxstep=1)
    # print(ret.T)
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
    
    return differential_values, population


def plotsir(differential_values, x_label, y_label, color_list, compartment_list):
    f, ax = plt.subplots(1, 1, figsize=(10, 4))
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
    # plt.show()


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
                compartment_list = c.compartment_definition.initial_compartment
                boolean_list = c.compartment_definition.boolean_value
                susceptible_value = False
                infected_value = False
                recovered_value = False
                exposed_value = False
                vaccinated_value = False
                dead_value = False
                
                for x in range(0, len(compartment_list)):
                    if compartment_list[x] == "susceptible":
                        susceptible_value = convert_string_to_boolean(boolean_list[x])
                    elif compartment_list[x] == "infected":
                        infected_value = convert_string_to_boolean(boolean_list[x])
                    elif compartment_list[x] == "recovered":
                        recovered_value = convert_string_to_boolean(boolean_list[x])
                    elif compartment_list[x] == "exposed":
                        exposed_value = convert_string_to_boolean(boolean_list[x])
                    elif compartment_list[x] == "vaccinated":
                        vaccinated_value = convert_string_to_boolean(boolean_list[x])
                    elif compartment_list[x] == "dead":
                        dead_value = convert_string_to_boolean(boolean_list[x])
                
                specified_compartments = create.determine_compartments(susceptible_value, infected_value, recovered_value, exposed_value, vaccinated_value, dead_value)

                # This will contain the code for generating equations once finished

            elif c.__class__.__name__ == "Calculate":
                final_value_list = []
                susceptible = 0
                infected = 0
                recovered = 0
                exposed = 0
                vaccinated = 0
                dead = 0
                compartment_list = c.compartment.initial_compartment
                value_list = c.compartment.x
                for x in range(0, len(compartment_list)):
                    if compartment_list[x] == "susceptible":
                        susceptible = value_list[x]
                        final_value_list.append(susceptible)
                    if compartment_list[x] == "infected":
                        infected = value_list[x]
                        final_value_list.append(infected)
                    if compartment_list[x] == "recovered":
                        recovered = value_list[x]
                        final_value_list.append(recovered)
                    if compartment_list[x] == "exposed":
                        exposed = value_list[x]
                        final_value_list.append(exposed)
                    if compartment_list[x] == "vaccinated":
                        vaccinated = value_list[x]
                        final_value_list.append(vaccinated)
                    if compartment_list[x] == "dead":
                        dead = value_list[x]
                        final_value_list.append(dead)

                characteristic_list = c.disease_rates.disease_characteristic
                rate_list = c.disease_rates.x

                beta = 1
                gamma = 1
                vacc = 1
                mu = 1
                arrival = 1

                for x in range(0, len(characteristic_list)):
                    if characteristic_list[x] == "transmission":
                        beta = rate_list[x]
                    if characteristic_list[x] == "recovery":
                        gamma = rate_list[x]
                    if characteristic_list[x] == "vaccination":
                        vacc = rate_list[x]
                    if characteristic_list[x] == "death":
                        mu = rate_list[x]
                    if characteristic_list[x] == "arrival":
                        arrival = rate_list[x]

                time_frame = c.time_frame.x
                population = c.population_size.x

                differential_values= run_differentials(time_frame, value_list, beta, gamma, mu, vacc, specified_compartments, population)                

            elif c.__class__.__name__ == "Plot":
                compartment_colors = c.compartment_color.compartment_color

                x_label = c.x_label

                y_label = c.y_label


        return differential_values, x_label, y_label, compartment_colors, compartment_list


def main(text, debug=False):

    virall_meta = metamodel_from_file("virall_dsl/virall.tx", debug=False)
    # metamodel_export(virall_meta, "virall_dsl/virall_meta.dot")

    virall_model = virall_meta.model_from_str(text)
    # model_export(virall_model, "virall_dsl/program.dot")

    virall = Virall()

    differential_values, x_label, y_label, compartment_colors, compartment_list = virall.interpret(virall_model)

    return differential_values, x_label, y_label, compartment_colors, compartment_list
    

if __name__ == "__main__":
    text = "hello"
    file = text_to_file(text)
    main(file)
