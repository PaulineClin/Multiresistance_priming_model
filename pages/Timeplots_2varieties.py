import streamlit as st
from scipy.integrate import odeint
import numpy as np 
import matplotlib.pyplot as plt


def app():
    
    plt.rc('axes', labelsize=16) 
    
    st.markdown("## Plant-epidemic dynamics for a mixture with n=2 varieties")
    # st.markdown(r"""Infection dynamics over time for a model with $n=2$ varieties for the full and the reduced model.""") 
    
    @st.cache(suppress_st_warning = True)  
    def parameters(R,c,rho,nu,n):
        params = np.array([n, R, c, rho, nu])
        return params
        
    
    # Définition du modèle complet :
    @st.cache(suppress_st_warning = True)  
    def modele_priming_cpl(etats, t, params): 
        J1, J2, G1, G2, S1, S2 = etats          # Variables d'état
        n, R, c, rho, nu = params                  # Même ordre que dans le vecteur prédéfini
        
        etatdot = [ J1*(R*(1-c)*(1/2-S1-J1-G1+(1-rho)*S1)-1),
                    J2*(R*(1-c)*(1/2-S2-J2-G2+(1-rho)*S2)-1),
                    (1-c)**2*R*(G1+G2)*(1/2-S1-J1-G1)+(1-rho)*(1-c)**2*R*(G1+G2)*S1-G1,
                    (1-c)**2*R*(G1+G2)*(1/2-S2-J2-G2)+(1-rho)*(1-c)**2*R*(G1+G2)*S2-G2,
                    (1-c)*R*J1*(1/2-S1-J1-G1)-(1-rho)*S1*R*((1-c)*J1+(1-c)**2*(G1+G2))-nu*S1,   # on calcule la derivee de l'etat 
                    (1-c)*R*J2*(1/2-S2-J2-G2)-(1-rho)*S2*R*((1-c)*J2+(1-c)**2*(G1+G2))-nu*S2]
        return etatdot
    
    # Définition du modèle reduit : 
    @st.cache(suppress_st_warning = True) 
    def modele_priming_rdt(etats, t, params): 
        x1, x2, m = etats          # Variables d'état
        n, R, c, rho, nu = params                  # Même ordre que dans le vecteur prédéfini
        
        etatdot = [x1*(R*(1-c)*(1/n-m-x1-x2+(1-rho)*m)-1),
                    x2*(R*(1-c)**2*2*(1/n-m-x1-x2+(1-rho)*m)-1),
                    (1/n-m-x1-x2)*R*(1-c)*x1-(1-rho)*m*(R*(1-c)*x1+R*(1-c)**2*2*x2)-nu*m]
        return etatdot
    
    @st.cache(suppress_st_warning = True) 
    def duration(t_0,t_fin,pas_t):
        tspan = np.arange(t_0, t_fin, pas_t)
        return tspan
    
    # Paramètre du modèle :
    R = st.slider('Transmission rate (R):', min_value=1, max_value=100, value = 20)
    c = st.slider('Virulence cost (c):', min_value=0.0, max_value=1.0, value = 0.49)
    rho = st.slider(r'''Priming efficiency:''', min_value=0.0, max_value=1.0, value = 0.8)
    nu = 1
    n = 2
    
    # Encapsulation des paramètres du modèle dans une matrice  : 
    params = parameters(R,c,rho,nu,n)
    
    # Modèle complet :
    # Conditions initiales :
    EI = np.array([0.15, 0.01, 0.15, 0.01, 0.15, 0.01])
    
    # Temps : 
    t_0 = 0.0
    t_fin = st.number_input('Time max:', min_value=10, max_value = 500, value= 30)
    pas_t = 0.01
    
    # Définition du tspan (vecteur) via la fonction numpy :
    time = duration(t_0,t_fin,pas_t) # création d'un vecteur avec des valeurs uniformément espacées.
    
    # Integration numerique 
    int_priming_cpl = odeint(modele_priming_cpl, EI, time, args=(params,), hmax=pas_t)
    
    # ################################################################################################
    # Modèle réduit :
        
    # Conditions initiales :
    EI_bis = np.array([int_priming_cpl[2999,0], int_priming_cpl[2999,2], int_priming_cpl[2999,4]])
    
    
    # Integration numerique 
    int_priming_rdt = odeint(modele_priming_rdt, EI_bis, time, args=(params,), hmax=pas_t)
    
    # Représentation des trajectoires en fonction du temps : 
    # création d'une figure, et d'un (ou plusieurs) système d'axes
    fig, axes = plt.subplots(1,2, figsize=(15,5))
    
    # Tracer des simulations par rapport au temps pour le modèle complet : 
    axes[0].plot(time, int_priming_cpl[:,0], color='blue', label='$y_1$')
    axes[0].plot(time, int_priming_cpl[:,1], color='blue', linestyle='--', label='$y_2$')
    axes[0].plot(time, int_priming_cpl[:,2], color='red', label='$z_1$')
    axes[0].plot(time, int_priming_cpl[:,3], color='red', linestyle='--', label='$z_2$')
    axes[0].plot(time, int_priming_cpl[:,4], color='green', label='$m_1$')
    axes[0].plot(time, int_priming_cpl[:,5], color='green', linestyle='--', label='$m_2$')
    axes[0].legend(loc='upper right') # le label est exécuter via la fonction legend
    
    # labellisation des axes
    axes[0].set_xlabel('Time')
    axes[0].set_ylabel('Infected density')
    
    # modification des bornes des axes
    axes[0].set_ylim(0, np.max(int_priming_cpl+0.05))
    axes[0].set_xlim(0, 30) 
    
    # Tracer des simulations par rapport au temps modèle réduit : 
    axes[1].plot(time, int_priming_rdt[:,0], color='blue', label='$x_1$')
    axes[1].plot(time, int_priming_rdt[:,1], color='red', label='$x_2$')
    axes[1].plot(time, int_priming_rdt[:,2], color='green', label='$m$')
    axes[1].legend(loc='upper right') # le label est exécuter via la fonction legend
    
    # labellisation des axes
    axes[1].set_xlabel('Time')
    axes[1].set_ylabel('Infected density')
    
    # modification des bornes des axes
    axes[1].set_ylim(0, np.max(int_priming_cpl+0.05))
    axes[1].set_xlim(0, t_fin)

    # Show the pyplot figure in the app : 
    st.pyplot(fig)
    
    st.caption(r""" Infection dynamics over time for a model with $n =2$ varieties. Figure on the left shows a simulation from
             the full model. The full lines correspond to variety $1$, and the dashed lines correspond to variety $2$.
             The density of hosts of variety $i$ infected by the corresponding singly virulent pathogen genotype is
             $y_i$, $i =1,2$. The density of hosts of variety $i$ infected by the doubly virulent pathogen genotype is $z_i$,
             $i=1,2$. The density of hosts of variety $i$ that are primed is $m_i$, $i =1,2$. Figure on the right shows a simulation from the
             reduced model, i.e. the model keeping track of a single focal variety. The density of hosts of the
             focal variety infected by the corresponding singly virulent pathogen genotype is $x_1$. The density of
             hosts of the focal variety infected by the doubly virulent pathogen genotype is $x_2$. The density of hosts
             of the focal variety that are primed is $m$. """)
             

    
    