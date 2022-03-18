import streamlit as st

# Custom imports 
from multipage import MultiPage
from pages import presentation, Timeplots_2varieties, Timeplots_3varieties, virulence_complexity, Prevalence, AUDPC_n_var_ode

# Create an instance of the app 
app = MultiPage()

# Title of the main page
st.title("")

# Add all your applications (pages) here
app.add_page("Presentation", presentation.app)
app.add_page("Time plots n=2", Timeplots_2varieties.app)
app.add_page("Time plots n=3", Timeplots_3varieties.app)
app.add_page("Pathogen Fitness", virulence_complexity.app)
app.add_page("Prevalence of the disease", Prevalence.app)
app.add_page("AUDPC of the disease at the equilibrum", AUDPC_n_var_ode.app)
# app.add_page("AUDPC of the disease over time", AUDPC_transitoire_n_var_ode.app)
# The main app
app.run()
