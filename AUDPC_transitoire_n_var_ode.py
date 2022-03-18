#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calcul de l'AUDPC a partir du modèle général réduit en utilisant odeint en phase transitoire

@author: pauline
"""
import streamlit as st
import numpy as np # Numerical computing library
import matplotlib.pyplot as plt # Plotting library
import scipy.integrate #Integration library
from scipy.integrate import odeint
from matplotlib.ticker import MaxNLocator


def app(): 

    plt.rc('axes', labelsize=16) 
    
    st.markdown("# AUDPC according the number of varieties in the mixture")
    st.markdown("### There may be long calculation times for some values. ") 

    ######################################
    ### Etude l'AUDPC dans la phase transitoire
    ### pour l'ensemble des valeurs de p 
    ### (car on représente le diagramme de bifurcation)
    ######################################
    
    @st.cache(suppress_st_warning = True) 
    def inf_force_patho(k,x, R, c):
        inf_patho = k*R*(1-c)**k*x
        return(inf_patho)
        
        
    # Vector infection binomial coefficient:    
    @st.cache(suppress_st_warning = True) 
    def inf_combi(n,k):
        inf_comb = []
        for i in k:
            inf_comb.append(scipy.special.comb(n-1,i-1))  #math.comb(n-1,i-1)
        return(np.array(inf_comb))
    
    # Vector priming binomial coefficient:  
    @st.cache(suppress_st_warning = True) 
    def prim_combi(n,k):
        prim_comb = []
        for i in k:
            prim_comb.append(scipy.special.comb(n-1,i)) # math.comb(n-1,i)
        return(np.array(prim_comb))  
    
    # Time function:
    @st.cache(suppress_st_warning = True) 
    def duration(t_0, t_fin, pas_t):
        tspan = np.arange(t_0, t_fin, pas_t)
        return tspan
    
    # Parameters: 
    @st.cache(suppress_st_warning = True) 
    def parameters(R,c,rho,nu,n):
        params = np.array([n, R, c, rho, nu])
        return params
    
    # ode function:   
    @st.cache(suppress_st_warning = True) 
    def modele_1(EI,t,params):
        n, R, c, rho, nu = params
        
        etatdot = [k*R*(1-c)**k*EI[0:-1]*(((1/n)-sum(ICombi * EI[0:-1])))-EI[0:-1], 
                   sum(EI[0:-1])]
        
        etat_x = etatdot[0]
        etat_a = np.array([etatdot[1]])
        ETATDOT = np.concatenate((etat_x, etat_a))
        return ETATDOT
    
    # Définition du modèle : 
    @st.cache(suppress_st_warning = True) 
    def modele_2(EI,t,params):
        n, R, c, rho, nu = params
        
        etatdot = [k*R*(1-c)**k*EI[0:-2]*(((1/n)-EI[-2]-sum(ICombi * EI[0:-2])) +(1-rho)*EI[-2])-EI[0:-2],
                   ((1/n)-EI[-2]-sum(ICombi * EI[0:-2]))*sum(PCombi * k*R*(1-c)**k*EI[0:-2])-(1-rho)*EI[-2]*sum(ICombi * k*R*(1-c)**k*EI[0:-2])-nu*EI[-2], 
                   sum(EI[0:-2])]
        
        etat_x = etatdot[0]
        etat_m = np.array([etatdot[1]])
        etat_a = np.array([etatdot[2]])
        ETATDOT = np.concatenate((etat_x, etat_m, etat_a))
        return ETATDOT
    
    R = st.slider('Transmission rate (R):', min_value=1, max_value=100, value = 5)
    c = st.slider('Virulence cost (c):', min_value=0.0, max_value=1.0, value = 0.3)
    nu = st.slider('nu :', min_value=1.0, max_value=5.0, value = 1.0)
    rho = 1
    var = st.slider('Number of varieties (n):', min_value=1, max_value=100, value = 10)
    
    N = np.arange(1, var+1, 1)               # number of varieties in the mixture
    
    # Temps : 
    t_0 = 0.0
    t_fin = st.number_input('Temps max:', min_value=10, max_value = 100, value= 20)
    pas_t = 0.01
    
    # Définition du tspan (vecteur) via la fonction numpy :
    time = duration(t_0,t_fin,pas_t) # création d'un vecteur avec des valeurs uniformément espacées.
    
    # Maximum pathogen fitness 
    i = np.arange(1, 21, 0.1)
    Fitness = R*(1-c)**i*i
    Fitness_max = np.where(Fitness == np.amax(Fitness))
    K_max = int(np.around(i[Fitness_max]))
    
    
    Audpc = []
    
    for n in N:
        # Depends on n so is in the loop
        k = np.arange(1,n+1,1)  # vector of all possible virulence complexities
        
        # Combination for infection 
        ICombi = inf_combi(n,k)
        
        # Combination for priming
        PCombi = prim_combi(n, k)
    
        # Encapsulation des paramètres pour résoudre le système :
        params = np.array([n, R, c, rho, nu])
        
        if n <= K_max:
            # Initial conditions: 
            x_res = np.full((n), 0.01)  # for infected hote by each virulence complexities
            a = np.array([sum(x_res)])
            EI = np.concatenate((x_res, a))
    
    
            # Integration numerique 
            int_priming = odeint(modele_1, EI, time, args=(params,), hmax=pas_t)
            
            while np.any((int_priming[:,-1] < 0)):
                x_res = x_res*0.1  # for infected hote by each virulence complexities
                a = np.array([sum(x_res)])
                
                EI = np.concatenate((x_res, a))
            
                
                # Integration numerique 
                int_priming = odeint(modele_1, EI, time, args=(params,), hmax=pas_t)
            
        else:
            m_res = np.array([0.01])    # for priming hote   
            x_res = np.full((n), 0.01)  # for infected hote by each virulence complexities
            a = np.array([sum(x_res)])
            
            EI = np.concatenate((x_res, m_res, a))
        
        
            # Integration numerique 
            int_priming = odeint(modele_2, EI, time, args=(params,), hmax=pas_t)
            
            while np.any((int_priming[:,-2] < 0)):
                m_res = m_res*0.1    # for priming hote   
                x_res = x_res*0.1  # for infected hote by each virulence complexities
                a = np.array([sum(x_res)])
                
                EI = np.concatenate((x_res, m_res, a))
            
            
                # Integration numerique 
                int_priming = odeint(modele_2, EI, time, args=(params,), hmax=pas_t)
    
        Audpc.append(n*int_priming[:,-1])
    
    Finaux = np.vstack(Audpc)
    Finaux = np.transpose(np.array(Finaux))
    
     
    fig2, ax2 = plt.subplots() 
    
    ax2.plot(N, Finaux[100,:], label = 'T1', color='#fed976')
    ax2.plot(N, Finaux[499,:], label = 'T5', color ='#feb24c')
    ax2.plot(N, Finaux[999,:], label = 'T10', color='#fc4e2a')
    ax2.plot(N, Finaux[1999,:], label = 'T20', color ='black')
    
    ax2.set_xlabel('Proportion of resistant')
    ax2.set_ylabel('Audpc')
    
    # modification des bornes des axes
    ax2.set_ylim(0, np.max(Finaux))
    ax2.set_xlim(2, np.max(N))  
    
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True)) 
    ax2.yaxis.set_major_locator(MaxNLocator(integer=True)) 
    
    ax2.legend(loc='upper right', markerscale = 3, fontsize = 10)

    # Show the pyplot figure in the app : 
    st.pyplot(fig2)
