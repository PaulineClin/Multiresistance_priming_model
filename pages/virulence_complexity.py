# -*- coding: utf-8 -*-
"""

Fitness fonction according to the virulence complexity of the pathogen. 

@author: clin
"""

# Importation des différents packages :
import streamlit as st
import math
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def app():

    st.markdown("# Pathogen fitness according the virulence complexity of the pathogen genotype")
    st.markdown("## Fitness function")
    
    def FitnessMax(c):
        imax = -(1/(math.log(1-c)))
        i_inf = int(imax)
        i_sup = int(imax + 1)
        Fitness_i_inf = (1-c)**(i_inf)*R*i_inf
        Fitness_i_sup = (1-c)**(i_sup)*R*i_sup
        if Fitness_i_inf > Fitness_i_sup: 
            return(i_inf)
        else: 
            return(i_sup)
    
    # Figure de la fitness en fonction de i :
    
    R = st.slider('Transmission rate (R):', min_value=1, max_value=100, value = 20)
    c = st.slider('Virulence cost (c):', min_value=0.0, max_value=1.0, value = 0.3)
    i = np.arange(1, 21, 0.1)
    
    Fitness = R*(1-c)**i*i
    K_max = FitnessMax(c)
    
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
