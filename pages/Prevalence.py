# -*- coding: utf-8 -*-
"""

Prevalence of the disease according to the number of host genotypes in the mixture. 

@author: clin
"""

# Importation des différents packages : 
import streamlit as st
import numpy as np 
import numpy.polynomial.polynomial as nppol
import matplotlib.pyplot as plt
import math

from matplotlib.ticker import MaxNLocator


plt.rc('axes', labelsize=16) 


def app():

    st.markdown("## Prevalence of the disease at the equilibrium: ")
    st.write(r""" Variation of the prevalence  $\mathcal{P}$ as a function of the number of varieties in the mixture $n$.""")
    

    # Calcul de la prévalence de la maladie : 
    # Paramètres 
    R = st.slider('Transmission rate (R):', min_value=1, max_value=100, value = 20)
    c = st.slider('Virulence cost (c):', min_value=0.0, max_value=1.0, value = 0.51)
    nu = 1
        
    # Calcul de imax : 
    # Selection de l'entier n pour le calcul de la prévalence à partir de imax: 
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
        
    i = FitnessMax(c)
    
    Prevalence_final= []
    
    for k in range(1,100):
        rho = k*0.01
        n_var = np.arange(1, 100)
        Prevalence = []
        for n in n_var:
            n = int(n)
            if n<i :
                Prevalence_i = 1-(n/((1-c)**n*R*n))
            else : 
                if n-1<i : 
                    A = 0
                else : 
                    A = math.factorial(n-1)/(math.factorial(i)*math.factorial((n-1)-i))
                    A = int(A)
                if n-1<i-1 : 
                    B = 0
                    C = 0
                else : 
                    B =  math.factorial(n-1)/(math.factorial(i-1)*math.factorial((n-1)-(i-1)))
                    B= int(B)
                    C = math.factorial(n)/(math.factorial(i)*math.factorial(n-i))
                    C = int(C)
                phi = i*R*(1-c)**i
                a2 = (B*phi*(1-rho)*(A+B))/rho
                a1 = (((nu-rho+1)*n-phi*(1-rho))*B+A*(n-phi*(1-rho)))/(n*rho)
                a0 = -(nu*(phi-n))/(n*phi*rho)
                roots = nppol.polyroots([a0,a1,a2])
                x = np.amax(roots)
                
                Prevalence_i = B*n*x # Float
                
            Prevalence.append(Prevalence_i) # List 
        Prevalence_final.append(Prevalence) # Array of float
    
    Prevalence_final = np.array(Prevalence_final)
    
    fig1, ax1 = plt.subplots()
    
    
    ax1.plot(n_var,Prevalence_final[98,:], color='r', label=r"$\rho=1$")
    ax1.plot(n_var,Prevalence_final[49,:], color='orange', label=r"$\rho = 0.5$")
    ax1.plot(n_var,Prevalence_final[0,:], color='k', label=r"$\rho = 0$")
    ax1.hlines(0.10, xmin=1, xmax=np.max(n), colors = 'grey', linestyles = 'dashed', label = '10% threshold') 
    ax1.vlines(i, ymin=0, ymax=1, colors = 'grey', linestyles = 'dashdot', label = ' n = {} lines'.format(i)) 
    ax1.set_ylim(0, 1)
    ax1.set_xlim(1, int(i*R*(1-c)**i)+1)
    ax1.set_xlabel('Number of varieties')
    ax1.set_ylabel('Disease prevalence')
    ax1.legend(loc = 'upper right')
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    
    # Show the pyplot figure in the app : 
    st.pyplot(fig1)
    
    st.caption(r""" Total equilibrium prevalence of the disease $\mathcal{P}$} as a function of the number of varieties in the mixture $n$ for $3$ values of $\rho$ and $\nu=1$.
               The $10\%$ prevalence threshold corresponds to a possible acceptable threshold in an agroecological context.  """)