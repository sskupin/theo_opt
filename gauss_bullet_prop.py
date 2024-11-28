import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("Gaussian Bullet Propagation in in cubic approximation")

def initialize():
    var_double[0].set(100)    # w_0 k_0
    var_double[1].set(100)    # T_p \omega_0
    var_double[2].set(1)      # v_g k_0 / \omega_0 = 1 / (k_0^{(1)} \omega_0 / k_0) 
    var_double[3].set(0)      # k_0^{(2)} \omega_0^2 / k_0 
    var_double[4].set(0)      # k_0^{{3)} \omega_0^3 / k_0 
    var_double[5].set(40000)  # z k_0  
    calculate()
    
def reinitialize():
    calculate() # because sliders may have changed

# Gaussian beam parameters
    
def w(z): # in units of w_0, for z in units of z_0
    return np.sqrt(1+z**2)

def A(z): # in units of A_0, for z in units of z_0
    return 1/np.sqrt(1+z**2)

def Rinv(z): # in units of 1/z_0, for z in units of z_0
    return z/(1+z**2)

def phi(z): # for z in units of z_0
    return -np.arctan(z)
    
def calculate():
    gui.change_cursor(root,"trek")
    try:  
        k0 = var_double[0].get() # in units of 1/w_0
        omega0 = var_double[1].get() # in units of 1/T_p
        vg = var_double[2].get()*omega0/k0 # group velocity in units of w_0 / T_p
        k02 = var_double[3].get()*k0**2/omega0**2/2 # GVD coefficient in units of T_p^2 / z_{0F} = 2 T_p^2 / (k_0 w_0^2)
        k03 = var_double[4].get()*k0**2/omega0**3/2 # TOD coefficient in units of T_p^3 / z_{0F} = 2 T_p^3 / (k_0 w_0^2)
        z = var_double[5].get()*2/k0**2 # in units of z_{0F} = k_0 w_0^2 / 2
        zLF = z/2 # in units of L_{F}
        zLF_string.set(zLF)
        zLSTC = z/(k0*vg) # in units of L_{STC}
        zLSTC_string.set(zLSTC)
        zLD = z*np.abs(k02) # in units of L_{D}
        zLD_string.set(zLD)
        zLDT = z*np.abs(k03) # in units of L_{DT}
        zLDT_string.set(zLDT)
            
        f.clf()
        
        rperp = np.linspace(-2*w(z),2*w(z),256,endpoint=False, retstep=False)
        Nt = 4096
        factor = w(np.amax([zLSTC,zLD,zLDT])*2)
        t,delta_t = np.linspace(-4*factor,4*factor,Nt,endpoint=False, retstep=True)
        omega = 2*np.pi*np.fft.fftfreq(Nt,delta_t) # in units of 1/T_p, shifted
        
        T,Rperp = np.meshgrid(t,rperp, indexing='xy')
        u0 = np.exp(-Rperp**2-T**2)
        
        OMEGA,Rperp = np.meshgrid(omega,rperp, indexing='xy')
        ZOMEGA = z*(1 - OMEGA/(k0*vg))
        U = A(ZOMEGA)*np.exp(-Rperp**2/w(ZOMEGA)**2 + 1j*Rperp**2*Rinv(ZOMEGA) + 1j*phi(ZOMEGA) - OMEGA**2/4 + 1j*k02*OMEGA**2*z/2 + 1j*k03*OMEGA**3*z/6)/(2*np.sqrt(np.pi))*omega[1]
        u = np.fft.fftshift(np.fft.fft(U,axis=-1),axes=-1) # because its temporal backtransform

        a1 = f.add_subplot(211)
        im1 = a1.imshow(np.abs(u0) ,extent=[t[0], t[-1], rperp[0], rperp[-1]], aspect='auto', origin='lower', vmin=0, cmap='jet')
        a1.annotate(r'$|u_0|/|u_0|_{\rm max}$', xy=(0.9*t[0],0.8*rperp[-1]),horizontalalignment='left', verticalalignment='bottom', color='w')
        a1.set_xlabel(r'$\tau/T_{\rm p}$')
        a1.set_ylabel(r'$r/w_0$')
        plt.colorbar(im1,location='top',shrink=0.75)
                        
        a3 = f.add_subplot(212)
        im3 = a3.imshow(np.abs(u) ,extent=[t[0], t[-1], rperp[0], rperp[-1]], aspect='auto', origin='lower', vmin=0, cmap='jet')
        a3.annotate(r'$|u(z)|/|u_0|_{\rm max}$', xy=(0.9*t[0],0.8*rperp[-1]),horizontalalignment='left', verticalalignment='bottom', color='w')
        a3.set_xlabel(r'$\tau/T_{\rm p}$')
        a3.set_ylabel(r'$r/w_0$')
        plt.colorbar(im3,location='top',shrink=0.75)
            
        plt.tight_layout()
                
#        plt.savefig('gaussian_bullet_prop.pdf',bbox_inches='tight',dpi=300, transparent=True)

        canvas.draw()       
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[4,7])

canvas = gui.create_canvas(root,f)
canvas.draw() # for faster feedback to user on startup
mainframe = gui.create_mainframe(root)


var_double = gui.create_doublevar_vector(6)
zLF_string = Tk.StringVar()
zLSTC_string = Tk.StringVar()
zLD_string = Tk.StringVar()
zLDT_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_formula_with_latex(mainframe,r'$u_0 \propto $',r'$\exp\!\left(-\frac{r_\perp^2}{w_0^{2}} - \frac{\tau^2}{T_{\rm p}^2} \right)$',row)
row = gui.create_slider_with_latex(mainframe,r'Bullet radius $w_0 k_0=$',var_double[0],20,200,row)
row = gui.create_slider_with_latex(mainframe,r'Bullet length $T_{\rm p}\omega_0=$',var_double[1],20,200,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r'Group velocity $v_g k_0/\omega_0=$',var_double[2],0.9,1.1,row)
row = gui.create_slider_with_latex(mainframe,r'GVD coefficient $k_0^{(2)}\omega_0^2/k_0=$',var_double[3],-1,1,row)
row = gui.create_slider_with_latex(mainframe,r'TOD coefficient $k_0^{(3)}\omega_0^3/k_0=$',var_double[4],-2,2,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r'Propagation distance $z k_0=$',var_double[5],0,50000,row)
row = gui.create_label_with_latex(mainframe,r'$z/L_{\rm F}=$',zLF_string,row)
row = gui.create_label_with_latex(mainframe,r'$z/L_{\rm STC}=$',zLSTC_string,row)
row = gui.create_label_with_latex(mainframe,r'$z/L_{\rm D}=$',zLD_string,row)
row = gui.create_label_with_latex(mainframe,r'$z/L_{\rm DT}=$',zLDT_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)