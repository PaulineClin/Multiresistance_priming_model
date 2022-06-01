import streamlit as st


def app():
    st.markdown("<h1 style='text-align: center;'> Host mixtures for plant disease control: benefits from pathogen selection and immune priming</h1>", unsafe_allow_html=True)
    st.markdown("<h3 style='text-align: center;'> Pauline Clin, Frédéric Grognard, Didier Andrivon, Ludovic Mailleret, Frédéric Hamelin </h3>", unsafe_allow_html=True)
    
    st.markdown("This application is linked to our study submitted for publication to the journal Evolutionary Applications.")
    st.markdown("This application allows to produce the different figures presented in our article, for parameter values that can be chosen by the user.")
    st.markdown("Click on the left panel to choose the topic you are interested in:")
    
    st.write(r"""
    * **Pathogen fitness:** 
      * Maximization of the pathogen fitness as a function of the virulence complexity $k$.  
    * **Prevalence of the disease:**
      * Prevalence of the disease at equilibrium as a function of the number of varieties $n$. 
    * **Plant-epidemic dynamics for $n=2$ varieties (Supplementary material):**
      * Infection dynamics over time of the full and reduced models, for hosts of the focal variety that are primed or infected by pathogens of virulence complexity $k=1, ..., n$.
    * **Plant-epidemic dynamics for $n=3$ varieties (Supplementary material):** 
      * Infection dynamics over time of the full and reduced models, for hosts of the focal variety that are primed or infected by pathogens of virulence complexity $k=1, ..., n$.
    * **Area Under the Disease Progress Curve (AUDPC) over time (Supplementary material):**
      * This metric summarizes the epidemic size over time as a function of the number of varieties $n$.
    """)
    
    st.write("\n")
    
  
    
    