import pandas as pd


def define_relationships():
    relationships = pd.read_csv("virall_dsl/differential_relationships.csv")

    return relationships


def create_compartment_variables(input_string, compartment):
    if input_string == "S":
        S = compartment
        return S
    if input_string == "I":
        I = compartment
        return I
    if input_string == "R":
        R = compartment
        return R
    if input_string == "E":
        E = compartment
        return E
    if input_string == "D":
        D = compartment
        return D
    if input_string == "V":
        V = compartment
        return V


def create_rate_variables(input_string, mu, gamma, beta, vacc):
    if input_string == "mu":
        death_rate = mu
        return death_rate
    if input_string == "beta":
        transmission_rate = beta
        return transmission_rate
    if input_string == "gamma":
        recovery_rate = gamma
        return recovery_rate
    if input_string == "vaccination":
        vacc_rate = vacc
        return vacc_rate


def determine_compartments(S, I, R, E, V, D):
    specified_compartments = []
    if S == True:
        specified_compartments.append("S")
    if I == True:
        specified_compartments.append("I")
    if R == True:
        specified_compartments.append("R")
    if E == True:
        specified_compartments.append("E")
    if V == True:
        specified_compartments.append("V")
    if D == True:
        specified_compartments.append("D")

    # print(specified_compartments)
    return specified_compartments


def create_compartment_tables(specified_compartments, relationships):
    compartment_tables = {}
    for compartment in specified_compartments:
        compartment_dataframe = relationships.loc[
            (relationships["target_compartment"]) == compartment
        ]
        compartment_tables[compartment] = compartment_dataframe
    return compartment_tables


def create_equations(compartment_tables, compartment, gamma, mu, beta, vacc, population):
    equation_dictionary = {}
    for key, value in compartment_tables.items():
        equation_list = []
        for i, row in value.iterrows():
            target = row["target_compartment"]
            other = row["other_compartment"]
            take = row["take_away"]
            add = row["add_to"]
            rate = row["rate"]

            target_converted = create_compartment_variables(target, compartment)
            other_converted = create_compartment_variables(other, compartment)
            rate_converted = create_rate_variables(rate, mu, gamma, beta, vacc)

            if take == True and add == False:
                if other == "R" or other == "D" or other == "V":
                    first_equation = -(target_converted * rate_converted)
                else:
                    first_equation = -(target_converted * (other_converted/population) * rate_converted)

                equation_list.append(first_equation)
            
            if add == True and take == False:
                if other == "R" or other == "D" or other == "V":
                    first_equation = target_converted * rate_converted
                else:                     
                    first_equation = target_converted * (other_converted/population) * rate_converted

                equation_list.append(first_equation)
        equation_dictionary[key] = equation_list
    
    return equation_dictionary, population


def combine_equations(equation_dictionary, population):
    final_equation = 0
    dSdt = 0
    dIdt = 0
    dRdt = 0
    dEdt = 0
    dVdt = 0
    dDdt = 0
    # print(equation_dictionary)
    differentials_list = []
    for key, value in equation_dictionary.items():     
        # print(value)    
        for equation in value:
            for element in equation:
                final_equation = (final_equation + element)
            # print(final_equation)

        if key == "S":
            dSdt = final_equation/population
        if key == "I":
            dIdt = final_equation/population
        if key == "R":
            dRdt = final_equation/population
        if key == "E":
            dEdt = final_equation/population
        if key == "V":
            dVdt = final_equation/population
        if key == "D":
            dDdt = final_equation/population
        
    if dSdt !=0:
        differentials_list.append(dSdt)
    if dIdt != 0:
        differentials_list.append(dIdt)
    if dRdt != 0:
        differentials_list.append(dRdt)
    if dEdt != 0:
        differentials_list.append(dEdt)
    if dVdt != 0:
        differentials_list.append(dVdt)
    if dDdt != 0:
        differentials_list.append(dDdt)
    
    # print(len(differentials_list))
    return differentials_list


def run_differentials(compartment, t, beta, gamma, mu, vacc, specified_compartments, population):
    relationships = define_relationships()
    # print(specified_compartments)
    compartment_tables = create_compartment_tables(specified_compartments, relationships)

    # print(compartment_tables)
    equation_dictionary, population = create_equations(compartment_tables, compartment, gamma, mu, beta, vacc, population)
    # print(equation_dictionary)
    differentials_list = combine_equations(equation_dictionary, population)
    return differentials_list


if __name__ == "__main__":
    specified_compartments = determine_compartments(
        True, True, True, False, False, False
    )

    # print(specified_compartments)

    run_differentials(2, 50, 0.1, 0.1, 0.1, 0.1, specified_compartments)