import streamlit as st


def app():
    st.markdown("<h1 style='text-align: center;'> Host mixtures for plant disease control: benefits from pathogen selection and immune priming</h1>", unsafe_allow_html=True)
    st.markdown("<h3 style='text-align: center;'> Clin Pauline, Grognard Frédéric, Andrivon Didier, Mailleret Ludovic, Hamelin Frédéric</h3>", unsafe_allow_html=True)
    
    st.markdown("This application is linked to our study submitted for publication to the journal Evolutionary Applications.")
    st.markdown("This application allows to produce the different figures presented in the companion article for parameter values that can be chosen by the user.")
    st.markdown("Click on the left panel to choose the topic you are interested in:")
    
    st.write(r"""
    * **Pathogen fitness:** 
      * Maximization of the pathogen fitness as a function of the virulence complexity $k$.  
    * **Prevalence of the disease:**
      * Variation of the prevalence of the disease at the equilibrium as a function of the number of varieties $n$. 
    * **Plant-epidemic dynamics for $n=2$ varieties (Supplementary material):**
      * Infection dynamics over time of the full model and the reduce model for primed and infected hosts of the focal variety by each virulence complexity $k$.
    * **Plant-epidemic dynamics for $n=3$ varieties (Supplementary material):** 
      * Infection dynamics over time of the full model and the reduce model for primed and infected hosts of the focal variety by each virulence complexity $k$.
    * **The Area Under the Disease Progress Curve (AUDPC) of the disease over time (Supplementary material):**
      * This metric summarizes the epidemic size over the time as a function of the number of varieties $n$.
    """)
    
    st.write("\n")
    
  
    
    