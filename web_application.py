import streamlit as st
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from virall_dsl import virall

def read_markdown_file(file):
    return Path(file).read_text()


def original_setup():
    st.image("images/VIRALL.png")

    virall_dsl_button = st.sidebar.button("Use VIRALL's DSL")
    virall_basic_button = st.sidebar.button("VIRALL 2.0")
    documentation_button = st.sidebar.button("Documentation")
    about_button = st.sidebar.button("About VIRALL")
    transmission_button = st.sidebar.button("More on Infectious Disease")
    contribute_button = st.sidebar.button("Contributing")
    
    menu = st.selectbox("Check out VIRALL", ("Home", "Use VIRALL's DSL", "VIRALL 2.0", "Documentation", "About VIRALL", "More on Infectious Disease", "Contributing"))

    return virall_dsl_button, virall_basic_button, about_button, transmission_button, contribute_button, documentation_button, menu


def virall_dsl_page():
    st.title("Go ahead! Start writing your code!")
    st.subheader("Below is a sample model using the VIRALL DSL. Feel free to modify it, or even delete it and write your own code!")
    text_file = open("virall_dsl/SIR_model.virall")
    st.markdown("Press `ctrl+enter` to save any code modifications in the code window")
    sample_code = text_file.read()
    text = st.text_area("VIRALL Code Window", height=500, value = sample_code)
    st.subheader("Results")

    # temp_file = virall.text_to_file(text)

    # st.write(type(text))
    differential_values, x_label, y_label, compartment_colors, compartment_list = virall.main(text)


    f, ax = plt.subplots(1, 1, figsize=(10, 4))
    for x in range(0, (len(compartment_list))):
        t = differential_values[0]
        color = compartment_colors[x]
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
    
    st.pyplot(f)


def about_page():
    st.title("About the Developer")
    st.header("Meet Maddy Kapfhammer")
    st.image("images/profile.jpg")
    about = read_markdown_file("markdown_files/about.md")
    st.markdown(about, unsafe_allow_html=True)


def transmission_page():
    transmission = read_markdown_file("markdown_files/disease_transmission.md")
    st.markdown(transmission, unsafe_allow_html=True)

    st.header("Compartmental Modeling")
    st.image("images/SIRCompartmentDiagram.png")
    compartmental = read_markdown_file("markdown_files/compartmental_modeling.md")
    st.markdown(compartmental, unsafe_allow_html=True)

def contributing_page():
    contribute = read_markdown_file("markdown_files/contribute.md")
    st.markdown(contribute, unsafe_allow_html=True)


def determine_page(virall_dsl, virall_basic, about_button, transmission_button, contribute_button, documentation_button, menu):
    if virall_dsl == True or menu == "Use VIRALL's DSL":
        virall_dsl_page()
    if virall_basic == True or menu == "VIRALL 2.0":
        st.write("Coming soon!")
    if about_button == True or menu == "About VIRALL":
        about_page()
    if transmission_button == True or menu == "More on Infectious Disease":
        transmission_page()
    if contribute_button == True or menu == "Contributing":
        contributing_page()
    if documentation_button == True or menu == "Documentation":
        st.write("Documentation coming soon!")


if __name__ == "__main__":
    virall_dsl, virall_basic, about_button, transmission_button, contribute_button, documentation_button, menu = original_setup()
    determine_page(virall_dsl, virall_basic, about_button, transmission_button, contribute_button, documentation_button, menu)
