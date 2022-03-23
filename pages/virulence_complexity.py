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

 
def app():    
    
    st.markdown("## Pathogen fitness according the virulence complexity of the genotype pathogen")
    st.markdown(r"""Pathogen fitness, $\phi_k=kR(1-c)^k$, as a function of virulence complexity $k$.""") 
                
    
    # Figure de la fitness en fonction de i :
    
    @st.cache(suppress_st_warning = True)  
    def fitness(R,c): 
        Fitness = R*(1-c)**i*i
        return Fitness
        
    @st.cache(suppress_st_warning = True)       
    def Fitness_max(R,c, Fitness):
        Fitness_max = np.where(Fitness == np.amax(Fitness))
        return int(Fitness_max[0][0])
    
    @st.cache(suppress_st_warning = True)      
    def K_max(R,c,Fitness_max):
        K_max = int(np.around(i[Fitness_max]))
        return K_max
           
    R = st.slider('Transmission rate (R):', min_value=1, max_value=100, value = 20)
    c = st.slider('Virulence cost (c):', min_value=0.0, max_value=1.0, value = 0.49)
    i = np.arange(0, 21, 0.1)

    path_fitness = fitness(R,c)
    path_fitness_max = Fitness_max(R,c,path_fitness)
    complex_max = K_max(R,c,path_fitness_max)
    
    #print(K_max)
    
    fig, ax = plt.subplots()
    plt.plot(i, path_fitness, 'k')
    ax.set_ylim(0, np.max(path_fitness+2))
    ax.set_xlim(0, np.max(i))
    ax.vlines(complex_max, ymin=0, ymax=np.max(path_fitness+2), colors = 'k', linestyles='dashed', label =' k = {} lines'.format(complex_max))

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    ax.set_xlabel('Virulence complexity k')
    ax.set_ylabel('Pathogen Fitness')
    ax.legend(loc = 'upper right')
    
    # Show the pyplot figure in the app : 
    st.pyplot(fig)
    
    st.caption(r"""Virulence complexity is the number of varieties a pathogen can infect, subject to a multiplicative cost $c$. 
    The maximum possible transmission rate is $R$. The pathogen fitness is maximized for an intermediate level of virulence complexity, $k^*$. 
    Note that the dotted line does not go through the maximum of the curve because $k$ can take only integer values.""")
