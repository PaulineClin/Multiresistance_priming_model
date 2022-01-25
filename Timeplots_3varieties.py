import streamlit as st
from scipy.integrate import odeint
import numpy as np 
import matplotlib.pyplot as plt


def app():
    st.markdown("# Epidemiological dynamics for a mixture with n=3 varieties")

    # Upload the dataset and save as csv
    st.markdown("## Time plots of the model.") 
    st.write("\n")
    st.latex(r''' \text{Let fixed } \nu=1 \, , \text{ which is a dimensionless parameter equal to }
                (\gamma+\alpha)/\alpha \\
                \text{ where } \gamma \text{ is the priming loss rate 
                and } \alpha \, , \text{ the harvest and replanting rate.} ''')
    
    
    plt.rc('axes', labelsize=16) 
    
    # Cas n = 3 variétés résistantes et priming efficace : 
    
    # Paramètre du modèle :
    R = st.slider('Transmission rate (R):', min_value=1, max_value=100, value = 20)
    c = st.slider('Virulence cost (c):', min_value=0.0, max_value=1.0, value = 0.49)
    rho = st.slider(r'Priming efficiency:', min_value=0.0, max_value=1.0, value=0.8)
    nu = 1
    n = 3
    
    # Encapsulation des paramètres du modèle dans une matrice  : 
    params = np.array([n, R, c, rho, nu])
    
    #############################################################
    # Modèle complet : 
        
        
    # Conditions initiales :
    EI = np.array([0.01, 0.03, 0.06,
                    0.01, 0.06, 
                    0.01, 0.06, 
                    0.01, 0.06,
                    0.01, 0.03, 0.06, 
                    0.01, 0.03, 0.06])
    
    # Temps : 
    t_0 = 0.0
    t_fin = st.number_input('Temps max:', min_value=10, max_value = 500, value= 200)
    pas_t = 0.01          
    
    # Définition du tspan (vecteur) via la fonction numpy :
    tspan = np.arange(t_0, t_fin, pas_t) # création d'un vecteur avec des valeurs uniformément espacées.
    
    # # Définition du modèle : 
    def modele_priming_cpl(etats, t, params): 
        W1, W2, W3, Y12, Y13, Y21, Y23, Y31, Y32, Z1, Z2, Z3, M1, M2, M3 = etats          # Variables d'état
        n, R, c, rho, nu = params                  # Même ordre que dans le vecteur prédéfini
        
        etatdot = [(1-c)*R*W1*(1/3-M1-W1-Y12-Y13-Z1)+(1-rho)*(1-c)*R*W1*M1-W1,
                    (1-c)*R*W2*(1/3-M2-W2-Y21-Y23-Z2)+(1-rho)*(1-c)*R*W2*M2-W2,
                    (1-c)*R*W3*(1/3-M3-W3-Y31-Y32-Z3)+(1-rho)*(1-c)*R*W3*M3-W3, 
                    (1-c)**2*R*(Y12+Y21)*(1/3-M1-W1-Y12-Y13-Z1)+(1-rho)*(1-c)**2*R*(Y12+Y21)*M1-Y12, 
                    (1-c)**2*R*(Y13+Y31)*(1/3-M1-W1-Y12-Y13-Z1)+(1-rho)*(1-c)**2*R*(Y13+Y31)*M1-Y13,
                    (1-c)**2*R*(Y12+Y21)*(1/3-M2-W2-Y21-Y23-Z2)+(1-rho)*(1-c)**2*R*(Y12+Y21)*M2-Y21,
                    (1-c)**2*R*(Y23+Y32)*(1/3-M2-W2-Y21-Y23-Z2)+(1-rho)*(1-c)**2*R*(Y23+Y32)*M2-Y23,
                    (1-c)**2*R*(Y13+Y31)*(1/3-M3-W3-Y31-Y32-Z3)+(1-rho)*(1-c)**2*R*(Y13+Y31)*M3-Y31,
                    (1-c)**2*R*(Y32+Y23)*(1/3-M3-W3-Y31-Y32-Z3)+(1-rho)*(1-c)**2*R*(Y32+Y23)*M3-Y32,
                    (1-c)**3*R*(Z1+Z2+Z3)*(1/3-M1-W1-Y12-Y13-Z1)+(1-rho)*(1-c)**3*R*Z1*M1-Z1,
                    (1-c)**3*R*(Z1+Z2+Z3)*(1/3-M2-W2-Y21-Y23-Z2)+(1-rho)*(1-c)**3*R*Z2*M2-Z2,
                    (1-c)**3*R*(Z1+Z2+Z3)*(1/3-M3-W3-Y31-Y32-Z3)+(1-rho)*(1-c)**3*R*Z3*M3-Z3,
                    (1-c)*R*(W2+W3)*(1/3-M1-W1-Y12-Y13-Z1)+(1-c)**2*R*(Y23+Y32)*(1/3-M1-W1-Y12-Y13-Z1)-(1-rho)*(1-c)*R*W1*M1-(1-rho)*(1-c)**2*R*(Y12+Y21)*M1-(1-rho)*(1-c)**2*R*(Y13+Y31)*M1-(1-rho)*(1-c)**3*R*(Z1+Z2+Z3)*M1-nu*M1,
                    (1-c)*R*(W1+W3)*(1/3-M2-W2-Y21-Y23-Z2)+(1-c)**2*R*(Y13+Y31)*(1/3-M2-W2-Y21-Y23-Z2)-(1-rho)*(1-c)*R*W2*M2-(1-rho)*(1-c)**2*R*(Y12+Y21)*M2-(1-rho)*(1-c)**2*R*(Y23+Y32)*M2-(1-rho)*(1-c)**3*R*(Z1+Z2+Z3)*M2-nu*M2,
                    (1-c)*R*(W1+W2)*(1/3-M3-W3-Y31-Y32-Z3)+(1-c)**2*R*(Y12+Y21)*(1/3-M3-W3-Y31-Y32-Z3)-(1-rho)*(1-c)*R*W3*M3-(1-rho)*(1-c)**2*R*(Y13+Y31)*M3-(1-rho)*(1-c)**2*R*(Y23+Y32)*M3-(1-rho)*(1-c)**3*R*(Z1+Z2+Z3)*M3-nu*M3]
        return etatdot
    
    # Integration numerique 
    int_priming_cpl = odeint(modele_priming_cpl, EI, tspan, args=(params,), hmax=pas_t) 
    
    # Modèle réduit :
    
    # Conditions initiales :
    EI = np.array([int_priming_cpl[-1,1], int_priming_cpl[-1,3], int_priming_cpl[-1,9], int_priming_cpl[-1,14]])
    # EI = np.array([0.02,0.02,0.02,0.02])

    # Temps : 
    t_0_bis = 0.0
    #Qt_fin_bis = 50
    pas_t_bis = 0.01
    
    # Définition du tspan (vecteur) via la fonction numpy :
    tspan_bis = np.arange(t_0_bis, t_fin, pas_t_bis) # création d'un vecteur avec des valeurs uniformément espacées.
    
    # Définition du modèle reduit : 
    def modele_priming_rdt(etats, t, params): 
        x1, x2, x3, m = etats          # Variables d'état
        n, R, c, rho, nu = params                  # Même ordre que dans le vecteur prédéfini
        
        etatdot = [x1*R*(1-c)*(1/n-m-x1-2*x2-x3)+(1-rho)*R*(1-c)*x1*m-x1,
                    2*x2*R*(1-c)**2*(1/n-m-x1-2*x2-x3)+(1-rho)*R*(1-c)**2*2*x2*m-x2,
                    3*x3*R*(1-c)**3*(1/n-m-x1-2*x2-x3)+(1-rho)*R*(1-c)**3*3*x3*m-x3,
                    2*R*(1/n-m-x1-2*x2-x3)*((1-c)*x1+(1-c)**2*x2)
                    -(1-rho)*R*m*((1-c)*x1+(1-c)**2*4*x2+(1-c)**3*3*x3)-nu*m]
        return etatdot
    
    
    # Integration numerique 
    int_priming_rdt = odeint(modele_priming_rdt, EI, tspan_bis, args=(params,), hmax=pas_t)
    
    # Représentation des trajectoires en fonction du temps : 
    # création d'une figure, et d'un (ou plusieurs) système d'axes
    fig, axes = plt.subplots(2,2, figsize=(15,10))
    
    # Tracer des simulations par rapport au temps pour le model complet en deux fois : 
    axes[0][0].plot(tspan, int_priming_cpl[:,0], color='blue', label='$y_1$')
    axes[0][0].plot(tspan, int_priming_cpl[:,1], color='blue', label='$y_2$', linestyle='--')
    axes[0][0].plot(tspan, int_priming_cpl[:,2], color='blue', label='$y_3$',linestyle='-.')
    axes[0][0].plot(tspan, int_priming_cpl[:,9], color = 'orange', label='$w_1$')
    axes[0][0].plot(tspan, int_priming_cpl[:,10], color = 'orange', label='$w_2$', linestyle='--')
    axes[0][0].plot(tspan, int_priming_cpl[:,11], color = 'orange', label='$w_3$', linestyle='-.')
    axes[0][0].plot(tspan, int_priming_cpl[:,12], color = 'green', label='$m_1$')
    axes[0][0].plot(tspan, int_priming_cpl[:,13], color = 'green', label='$m_2$', linestyle='--',)
    axes[0][0].plot(tspan, int_priming_cpl[:,14], color = 'green', label='$m_3$', linestyle='-.', )
    axes[0][0].legend(loc='upper right') # le label est exécuter via la fonction legend
    
    # labellisation des axes
    axes[0][0].set_xlabel('Time')
    axes[0][0].set_ylabel('Infected density')
        
    # modification des bornes des axes
    axes[0][0].set_ylim(0, np.max(int_priming_cpl+0.05))
    axes[0][0].set_xlim(0, np.max(tspan))   
    
    axes[1][0].plot(tspan, int_priming_cpl[:,3], color='darkred', label='$z_{12}$')
    axes[1][0].plot(tspan, int_priming_cpl[:,4], color='darkred', label='$z_{13}$', linestyle='--')
    axes[1][0].plot(tspan, int_priming_cpl[:,5], color='red', label = '$z_{21}$')
    axes[1][0].plot(tspan, int_priming_cpl[:,6], color= 'red', label='$z_{23}$', linestyle='--')
    axes[1][0].plot(tspan, int_priming_cpl[:,7], color='lightsalmon', label = '$z_{31}$')
    axes[1][0].plot(tspan, int_priming_cpl[:,8], color='lightsalmon', label='$z_{32}$', linestyle='--')
    axes[1][0].legend(loc='upper right') # le label est exécuter via la fonction legend
    
    # labellisation des axes
    axes[1][0].set_xlabel('Time')
    axes[1][0].set_ylabel('Infected density')
    
    # modification des bornes des axes
    axes[1][0].set_ylim(0,np.max(int_priming_cpl+0.05))
    axes[1][0].set_xlim(0, np.max(tspan))   
    
    # Tracer des simulations par rapport au temps pour le model réduit en deux fois :  
    axes[0][1].plot(tspan_bis, int_priming_rdt[:,0], color='blue', label='$x_1$')
    axes[0][1].plot(tspan_bis, int_priming_rdt[:,2], color='orange', label='$x_3$')
    axes[0][1].plot(tspan_bis, int_priming_rdt[:,3], color='green', label='$m$')
    axes[0][1].legend(loc='upper right') # le label est exécuter via la fonction legend
    
    # labellisation des axes
    axes[0][1].set_xlabel('Time')
    axes[0][1].set_ylabel('Infected density')

    # modification des bornes des axes
    axes[0][1].set_ylim(0, np.max(int_priming_cpl+0.05))
    axes[0][1].set_xlim(0, np.max(tspan_bis))   
    
    # Tracer des simulations par rapport au temps : 
    axes[1][1].plot(tspan_bis, int_priming_rdt[:,1], color='red', label='$x_2$')
    axes[1][1].legend(loc='upper right') # le label est exécuter via la fonction legend
    
    # labellisation des axes
    axes[1][1].set_xlabel('Time')
    axes[1][1].set_ylabel('Infected density')
    
    # modification des bornes des axes
    axes[1][1].set_ylim(0, np.max(int_priming_cpl+0.05))
    axes[1][1].set_xlim(0, np.max(tspan_bis)) 
    
    # Show the pyplot figure in the app : 
    st.pyplot(fig)