# -*- coding: utf-8 -*-
"""

Fitness fonction according to the virulence complexity of the pathogen. 

@author: clin
"""

# Importation des différents packages :
import streamlit as st
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

@st.cache(suppress_st_warning = True)
def app():

    st.markdown("# Pathogen fitness according the virulence complexity of the genotype pathogène")
    st.markdown("## Fitness function")
    # Figure de la fitness en fonction de i :
    
    R = st.slider('Transmission rate (R):', min_value=1, max_value=100, value = 20)
    c = st.slider('Virulence cost (c):', min_value=0.0, max_value=1.0, value = 0.3)
    i = np.arange(1, 21, 0.1)
    Fitness = R*(1-c)**i*i
    Fitness_max = np.where(Fitness == np.amax(Fitness))
    K_max = int(np.around(i[Fitness_max]))
    print(K_max)
    
    fig, ax = plt.subplots()
    plt.plot(i, Fitness, 'k')
    ax.set_ylim(0, np.max(Fitness+2))
    ax.set_xlim(1, np.max(i))
    ax.vlines(K_max, ymin=0, ymax=np.max(Fitness+2), colors = 'k', linestyles='dashed', label =' k = {} lines'.format(K_max))
    
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    
    ax.set_xlabel('Virulence complexity k')
    ax.set_ylabel('Pathogen Fitness')
    ax.legend(loc = 'upper right')
    
    # Show the pyplot figure in the app : 
    st.pyplot(fig)
