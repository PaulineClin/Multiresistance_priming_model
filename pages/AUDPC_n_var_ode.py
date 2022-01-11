#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calcul de l'AUDPC a partir du modèle général réduit en utilisant odeint. 

@author: pauline
"""

# Importation des différents packages :
import streamlit as st
import scipy
from scipy.integrate import odeint
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import math


def app(): 

    plt.rc('axes', labelsize=16) 
    
    st.markdown("# AUDPC according the number of varieties in the mixture")
    st.markdown("### There may be long calculation times for some values. ")
    
    ###############################################################################
    ## Model function for a focal variety:   
    # The infection force of a given pathogen genotype:
    def inf_force_patho(k,x):
        inf_patho = k*R*(1-c)**k*x
        return(inf_patho)
        
        
    # Vector infection binomial coefficient:
    def inf_combi(n,k):
        inf_comb = []
        for i in k:
            inf_comb.append(scipy.special.comb(n-1,i-1))  #math.comb(n-1,i-1)
        return(np.array(inf_comb))
    
    # Vector priming binomial coefficient:
    def prim_combi(n,k):
        prim_comb = []
        for i in k:
            prim_comb.append(scipy.special.comb(n-1,i)) # math.comb(n-1,i)
        return(np.array(prim_comb))  
    
    ###############################################################################
    ###############################################################################    
    # Fixed parameters: 
    R = st.slider('Transmission rate (R):', min_value=1, max_value=100, value = 20)
    c = st.slider('Virulence cost (c):', min_value=0.0, max_value=1.0, value = 0.49)
    nu = st.slider('nu :', min_value=1.0, max_value=5.0, value = 1.0)
    var = st.slider('Number of varieties (n):', min_value=1, max_value=100, value = 10)

    N = np.arange(1, var+1, 1)               # number of varieties in the mixture
    
    # Temps : 
    t_0 = 0.0
    t_fin = 100
    pas_t = 0.01
    
    # Définition du tspan (vecteur) via la fonction numpy :
    tspan = np.arange(t_0, t_fin, pas_t) # création d'un vecteur avec des valeurs uniformément espacées.
    
    # Maximum pathogen fitness 
    i = np.arange(1, 21, 0.1)
    Fitness = R*(1-c)**i*i
    Fitness_max = np.where(Fitness == np.amax(Fitness))
    K_max = int(np.around(i[Fitness_max]))
    
    
    Rho_AUDPC = []
    
    for x in range(1, 150, 50):
        rho = x*0.01
        AUDPC = []
        
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
                
                # Définition du modèle : 
                def modele_1(EI,t,params):
                    n, R, c, rho, nu = params
                    
                    etatdot = [k*R*(1-c)**k*EI[0:-1]*(((1/n)-sum(ICombi * EI[0:-1])))-EI[0:-1], 
                               sum(EI[0:-1])]
                    
                    etat_x = etatdot[0]
                    etat_a = np.array([etatdot[1]])
                    ETATDOT = np.concatenate((etat_x, etat_a))
                    return ETATDOT
        
        
                # Integration numerique 
                int_priming = odeint(modele_1, EI, tspan, args=(params,), hmax=pas_t)
                
                while np.any((int_priming[:,-1] < 0)):
                    x_res = x_res*0.1  # for infected hote by each virulence complexities
                    a = np.array([sum(x_res)])
                    
                    EI = np.concatenate((x_res, a))
                
                    # # Définition du modèle : 
                    def modele_1(EI,t,params):
                        n, R, c, rho, nu = params
                        
                        etatdot = [k*R*(1-c)**k*EI[0:-2]*(((1/n)-EI[-2]-sum(ICombi * EI[0:-2])) +(1-rho)*EI[-2])-EI[0:-2],
                                    ((1/n)-EI[-2]-sum(ICombi * EI[0:-2]))*sum(PCombi * k*R*(1-c)**k*EI[0:-2])-(1-rho)*EI[-2]*sum(ICombi * k*R*(1-c)**k*EI[0:-2])-nu*EI[-2], 
                                    sum(EI[0:-2])]
                        
                        etat_x = etatdot[0]
                        etat_m = np.array([etatdot[1]])
                        etat_a = np.array([etatdot[2]])
                        ETATDOT = np.concatenate((etat_x, etat_m, etat_a))
                        return ETATDOT
                
                
                    # Integration numerique 
                    int_priming = odeint(modele_1, EI, tspan, args=(params,), hmax=pas_t)
                
            else:
                m_res = np.array([0.01])    # for priming hote   
                x_res = np.full((n), 0.01)  # for infected hote by each virulence complexities
                a = np.array([sum(x_res)])
                
                EI = np.concatenate((x_res, m_res, a))
                
                
                # # Définition du modèle : 
                def modele_1(EI,t,params):
                    n, R, c, rho, nu = params
                    
                    etatdot = [k*R*(1-c)**k*EI[0:-2]*(((1/n)-EI[-2]-sum(ICombi * EI[0:-2])) +(1-rho)*EI[-2])-EI[0:-2],
                               ((1/n)-EI[-2]-sum(ICombi * EI[0:-2]))*sum(PCombi * k*R*(1-c)**k*EI[0:-2])-(1-rho)*EI[-2]*sum(ICombi * k*R*(1-c)**k*EI[0:-2])-nu*EI[-2], 
                               sum(EI[0:-2])]
                    
                    etat_x = etatdot[0]
                    etat_m = np.array([etatdot[1]])
                    etat_a = np.array([etatdot[2]])
                    ETATDOT = np.concatenate((etat_x, etat_m, etat_a))
                    return ETATDOT
            
            
                # Integration numerique 
                int_priming = odeint(modele_1, EI, tspan, args=(params,), hmax=pas_t)
                
                while np.any((int_priming[:,-2] < 0)):
                    m_res = m_res*0.1    # for priming hote   
                    x_res = x_res*0.1  # for infected hote by each virulence complexities
                    a = np.array([sum(x_res)])
                    
                    EI = np.concatenate((x_res, m_res, a))
                
                    # # Définition du modèle : 
                    def modele_1(EI,t,params):
                        n, R, c, rho, nu = params
                        
                        etatdot = [k*R*(1-c)**k*EI[0:-2]*(((1/n)-EI[-2]-sum(ICombi * EI[0:-2])) +(1-rho)*EI[-2])-EI[0:-2],
                                    ((1/n)-EI[-2]-sum(ICombi * EI[0:-2]))*sum(PCombi * k*R*(1-c)**k*EI[0:-2])-(1-rho)*EI[-2]*sum(ICombi * k*R*(1-c)**k*EI[0:-2])-nu*EI[-2], 
                                    sum(EI[0:-2])]
                        
                        etat_x = etatdot[0]
                        etat_m = np.array([etatdot[1]])
                        etat_a = np.array([etatdot[2]])
                        ETATDOT = np.concatenate((etat_x, etat_m, etat_a))
                        return ETATDOT
                
                
                    # Integration numerique 
                    int_priming = odeint(modele_1, EI, tspan, args=(params,), hmax=pas_t)
                
                 
            
            AUDPC.append(n*int_priming[-1,-1])
        Rho_AUDPC.append(AUDPC)
        
    Rho_AUDPC = np.transpose(np.array(Rho_AUDPC))
        
    fig1, ax1 = plt.subplots()
    
    # Tracer des simulations par rapport au temps : 
    ax1.plot(N, Rho_AUDPC[:,2], color='red', label=r'$\rho=1$')
    ax1.plot(N, Rho_AUDPC[:,1], color='orange', label=r'$\rho=0.5$')
    ax1.plot(N, Rho_AUDPC[:,0], color='black', label=r'$\rho=0$')
    ax1.hlines((np.max(Rho_AUDPC)*10)/100, xmin=1, xmax=np.max(N), colors = 'grey', linestyles = 'dashed', label = '10% threshold')
    ax1.legend(loc='upper right') # le label est exécuter via la fonction legend
    # ax1.vlines(2, ymin=0, ymax=np.max(Rho_AUDPC)+10, colors = 'grey', linestyles = 'dashdot', label = '$n=2$ line') 
    
    # labellisation des axes
    ax1.set_xlabel('Number of varieties')
    ax1.set_ylabel('AUDPC')
    
    # titre de la figure
    #fig1.suptitle(r'$R$ = 20, $c$ = 0.5, $\nu$ =1', fontsize = '14')
    #fig1.suptitle(r' $R = {}, c = {}, \nu = {} $'.format(R, c, nu), fontsize = '14') 
    
    # modification des bornes des axes
    # modification des bornes des axes:
    if np.any((Rho_AUDPC[-1,:] > 0.1)):
        x_max = np.max(N)
    else:
        x_max = np.where(Rho_AUDPC[:,2] < 0.1)
        x_max = x_max[0][0]
        
    ax1.set_xlim(np.min(N), x_max) 
    ax1.set_ylim(0, np.max(Rho_AUDPC)+10)
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True)) 
    
    # Show the pyplot figure in the app : 
    st.pyplot(fig1)
    


