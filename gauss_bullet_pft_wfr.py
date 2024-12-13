import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("Gaussian Bullet with PFT and WFR in vacuum")

def initialize():
    var_string[0].set("0")    # z_{IF} / R_I
    var_string[1].set("0")    # \Pi_I
    var_string[2].set("5")    # \Omega_I
    var_string[3].set("0")     # C_I
    var_string[4].set("0")     # z / z_{IF}
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)

# Gaussian beam parameters
    
def w(z,RIinv): # in units of w_I, for z and R_I in units of z_{IF}
    return np.sqrt((1+z*RIinv)**2+z**2)

def A(z,RIinv): # in units of A_I, for z and R_I in units of z_{IF}
    return 1/w(z,RIinv)

def Rinv(z,RIinv): # in units of z_{IF}, for z and R_I in units of z_{IF}
    return (z+z*RIinv**2+RIinv)/(z**2+(z*RIinv+1)**2)

def PFR(z,RIinv,PFRI,WFRI): # for z and R_I in units of z_{IF}
    return (PFRI*(1+z*RIinv)+WFRI*z)/w(z,RIinv)**2

def WFR(z,RIinv,PFRI,WFRI): # for z and R_I in units of z_{IF}
    return (WFRI*(1+z*RIinv)-PFRI*z)/w(z,RIinv)**2

def TB(z,RIinv,PFRI,WFRI): # in units of T_I, for z and R_I in units of z_{IF}
    return 2/np.sqrt(4-PFRI**2+PFR(z,RIinv,PFRI,WFRI)**2*w(z,RIinv)**2)

def C(z,RIinv,PFRI,WFRI,CI): # for z and R_I in units of z_{IF}
    return CI+(PFRI**2*RIinv+(WFRI+PFRI*RIinv)**2*z-PFR(z,RIinv,PFRI,WFRI)**2*w(z,RIinv)**4*Rinv(z,RIinv))/4
    
def calculate():
    gui.change_cursor(root,"trek")
    try:  
        RIinv = float(var_string[0].get()) # in units of 1/z_{IF}
        PFRI = float(var_string[1].get())
        WFRI = float(var_string[2].get())
        CI = float(var_string[3].get())
        zmax = float(var_string[4].get()) # in units of z_{IF}
        
        if np.abs(PFRI) >= 2:
            gui.input_error("Pulse-front tilt too large for Gaussian bullet. Re-initializing with previous parameters...",reinitialize)
        else:

            fig.clf()
        
            a1 = fig.add_subplot(211)
        
            z = np.linspace(0,zmax,256,endpoint=False, retstep=False)
        
            a1.plot(z,A(z,RIinv),z,PFR(z,RIinv,PFRI,WFRI),z,TB(z,RIinv,PFRI,WFRI),z,C(z,RIinv,PFRI,WFRI,CI))
        
            factor = np.sqrt(4-PFR(zmax,RIinv,PFRI,WFRI)**2*TB(zmax,RIinv,PFRI,WFRI)**2*w(zmax,RIinv)**2)
            x = np.linspace(-4*w(zmax,RIinv)/factor,4*w(zmax,RIinv)/factor,256,endpoint=False, retstep=False)
            t = np.linspace(-4/np.sqrt(4-PFRI**2),4/np.sqrt(4-PFRI**2),256,endpoint=False, retstep=False)
            
            T,X = np.meshgrid(t,x, indexing='xy')
            u = A(zmax,RIinv)*np.exp(-(1/w(zmax,RIinv)**2-1j*Rinv(zmax,RIinv))*X**2+(PFR(zmax,RIinv,PFRI,WFRI)+1j*WFR(zmax,RIinv,PFRI,WFRI))*X*T-(1/TB(zmax,RIinv,PFRI,WFRI)**2+1j*C(zmax,RIinv,PFRI,WFRI,CI))*T**2) # without Gouy phase

            factor = np.sqrt((4*(1+CI**2)-2*PFRI*CI*WFRI-PFRI**2+WFRI**2)/(4-PFRI**2))
            omega = np.linspace(-4*factor,4*factor,256,endpoint=False, retstep=False)      
            
            OMEGA,X = np.meshgrid(omega,x, indexing='xy')
            ubar = (1+CI**2)**0.25*A(zmax,RIinv)*TB(zmax,RIinv,PFRI,WFRI)/(np.sqrt(1+1j*C(zmax,RIinv,PFRI,WFRI,CI)*TB(zmax,RIinv,PFRI,WFRI)**2))*np.exp(-(1/w(zmax,RIinv)**2-1j*Rinv(zmax,RIinv))*X**2+TB(zmax,RIinv,PFRI,WFRI)**2/(4*(1+1j*C(zmax,RIinv,PFRI,WFRI,CI)*TB(zmax,RIinv,PFRI,WFRI)**2))*((PFR(zmax,RIinv,PFRI,WFRI)+1j*WFR(zmax,RIinv,PFRI,WFRI))**2*X**2-2*(WFR(zmax,RIinv,PFRI,WFRI)-1j*PFR(zmax,RIinv,PFRI,WFRI))*X*OMEGA-OMEGA**2)) # without Gouy phase

            factor = np.sqrt((4*(1+RIinv**2)+2*PFRI*RIinv*WFRI-PFRI**2+WFRI**2)/(4-PFRI**2))
            kx = np.linspace(-4*factor,4*factor,256,endpoint=False, retstep=False)    

            T,KX = np.meshgrid(t,kx, indexing='xy')
            uhat = (1+RIinv**2)**0.5*A(zmax,RIinv)*w(zmax,RIinv)**2/(1-1j*Rinv(zmax,RIinv)*w(zmax,RIinv)**2)*np.exp(-(1/TB(zmax,RIinv,PFRI,WFRI)**2+1j*C(zmax,RIinv,PFRI,WFRI,CI))*T**2+w(zmax,RIinv)**2/(4*(1-1j*Rinv(zmax,RIinv)*w(zmax,RIinv)**2))*((PFR(zmax,RIinv,PFRI,WFRI)+1j*WFR(zmax,RIinv,PFRI,WFRI))**2*T**2+2*(WFR(zmax,RIinv,PFRI,WFRI)-1j*PFR(zmax,RIinv,PFRI,WFRI))*KX*T-KX**2)) # without Gouy phase
            
            a3 = fig.add_subplot(212)
#            im3 = a3.imshow(np.abs(u) ,extent=[t[0], t[-1], x[0], x[-1]], aspect='auto', origin='lower', vmin=0, cmap='jet')
#            a3.annotate(r'$|u(z)|/|u_0|_{\rm max}$', xy=(0.9*t[0],0.8*x[-1]),horizontalalignment='left', verticalalignment='bottom', color='w')
#            a3.set_xlabel(r'$\tau/T_{\rm I}$')
#            a3.set_ylabel(r'$x/w_{\rm I}$')
#            im3 = a3.imshow(np.abs(ubar) ,extent=[omega[0], omega[-1], x[0], x[-1]], aspect='auto', origin='lower', vmin=0, cmap='jet')
#            a3.annotate(r'$|\bar{u}(z)|/|\bar{u}_0|_{\rm max}$', xy=(0.9*omega[0],0.8*x[-1]),horizontalalignment='left', verticalalignment='bottom', color='w')
#            a3.set_xlabel(r'$\bar{\omega}T_{\rm I}$')
#            a3.set_ylabel(r'$x/w_{\rm I}$')
            im3 = a3.imshow(np.abs(uhat) ,extent=[t[0], t[-1], kx[0], kx[-1]], aspect='auto', origin='lower', vmin=0, cmap='jet')
            a3.annotate(r'$|\hat{u}(z)|/|\hat{u}_0|_{\rm max}$', xy=(0.9*t[0],0.8*kx[-1]),horizontalalignment='left', verticalalignment='bottom', color='w')
            a3.set_xlabel(r'$\tau/T_{\rm I}$')
            a3.set_ylabel(r'$k_x w_{\rm I}$')
            plt.colorbar(im3,location='top',shrink=0.75)


            
            plt.tight_layout()
                
            plt.savefig('gauss_bullet_pfr_wfr.pdf',bbox_inches='tight',dpi=300, transparent=True)

            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()       
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

fig = plt.figure(1,[4,7])

canvas = gui.create_canvas(root,fig)
canvas.draw() # for faster feedback to user on startup
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(5)
var_save = gui.create_stringvar_vector(5)

initialize()

row = 1
row = gui.create_formula_with_latex(mainframe,r'$u_0 \propto \exp\!\left[-\left(1-\mathrm{i} \frac{ z_{\rm IF} }{R_{\rm I}}\right)\frac{x^2+y^2}{w_{\rm I}^{2}} \right.$',r'$\left. + \left(\Pi_{\rm I}+\mathrm{i}\Omega_{\rm I}\right)\frac{x\tau}{w_{\rm I}T_{\rm I}}  - \left(1+\mathrm{i}C_{\rm I}\right)\frac{\tau^2}{T_{\rm I}^2} \right]$',row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r'Inverse phase curvature $z_{\rm IF} / R_{\rm I}=$',var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r'Pulse-front tilt $\Pi_{\rm I} =$',var_string[1],row)
row = gui.create_entry_with_latex(mainframe,r'Wave-front rotation $\Omega_{\rm I} =$',var_string[2],row)
row = gui.create_entry_with_latex(mainframe,r'Chirp $C_{\rm I} =$',var_string[3],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r'Propagation distance $z/z_{\rm IF}=$',var_string[4],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)