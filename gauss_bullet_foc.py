import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("Focusing of Gaussian Bullet in vacuum")

def initialize():
    var_double[0].set(2000)   # w_I k_0
    var_double[1].set(100)    # T_I \omega_0
    var_double[2].set(100)    # f / w_I
    var_double[3].set(1)      # z / f
    calculate()
    
def reinitialize():
    calculate() # because sliders may have changed

# Gaussian beam parameters
    
def w(zI,f,z): # in units of w_I, for z, f, z_I in units of w_I
    return np.sqrt((1-z/f)+(z/zI)**2)

def A(zI,f,z): # in units of A_I, for z, f, z_I in units of w_I
    return 1/np.sqrt((1-z/f)+(z/zI)**2)

def Rinv(zI,f,z): # in units of 1/w_I, for z, f, z_I in units of w_I
    return (z+(zI/f)**2*(z-f))/(z**2+(zI/f)**2*(z-f)**2)

def phi(zI,f,z): # for z, f, z_I in units of w_I
    return np.angle(1-z/f-1j*z/zI)
    
def calculate():
    gui.change_cursor(root,"trek")
    try:  
        k0 = var_double[0].get() # in units of 1/w_I
        omega0 = var_double[1].get() # in units of 1/T_I
        f = var_double[2].get() # in units of w_I
        z = var_double[3].get()*f # in units of w_I
        zI = k0/2 # in units of w_I
        LSTC = 2*omega0*f**2/k0 # in units of w_I

        fig.clf()
        
        rperp = np.linspace(-2*w(zI,f,z),2*w(zI,f,z),256,endpoint=False, retstep=False)
        Nt = 4096
        factor = np.sqrt(1+f**2/LSTC**2)
        t,delta_t = np.linspace(-4*factor,4*factor,Nt,endpoint=False, retstep=True)
        omegabar = 2*np.pi*np.fft.fftfreq(Nt,delta_t) # in units of 1/T_I, shifted
        
        OMEGAbar,Rperp = np.meshgrid(omegabar,rperp, indexing='xy')
        OMEGA = omega0+OMEGAbar
        zIOMEGA = zI*OMEGA/omega0
        U = A(zIOMEGA,f,z)*np.exp(-Rperp**2/w(zIOMEGA,f,z)**2 + 1j*zIOMEGA*Rperp**2*Rinv(zIOMEGA,f,z) + 1j*phi(zIOMEGA,f,z) - OMEGAbar**2/4)/(2*np.sqrt(np.pi))*omegabar[1]
        u = np.fft.fftshift(np.fft.fft(U,axis=-1),axes=-1) # because its temporal backtransform
           
        a3 = fig.add_subplot(111)
        im3 = a3.imshow(np.abs(u) ,extent=[t[0], t[-1], rperp[0], rperp[-1]], aspect='auto', origin='lower', vmin=0, cmap='jet')
        a3.annotate(r'$|u(z)|/|u_0|_{\rm max}$', xy=(0.9*t[0],0.8*rperp[-1]),horizontalalignment='left', verticalalignment='bottom', color='w')
        a3.set_xlabel(r'$\tau/T_{\rm I}$')
        a3.set_ylabel(r'$r_\perp/w_{\rm I}$')
        plt.colorbar(im3,location='top',shrink=0.75)
            
        plt.tight_layout()
                
#        plt.savefig('gauss_bullet_foc.pdf',bbox_inches='tight',dpi=300, transparent=True)

        canvas.draw()       
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

fig = plt.figure(1,[4,3.5])

canvas = gui.create_canvas(root,fig)
canvas.draw() # for faster feedback to user on startup
mainframe = gui.create_mainframe(root)


var_double = gui.create_doublevar_vector(4)

initialize()

row = 1
row = gui.create_formula_with_latex(mainframe,r'$\bar{u}_0 \propto $',r'$\exp\!\left[-\frac{r_\perp^2}{w_{\rm I}^{2}}\left(1+\mathrm{i} \frac{\omega}{c}\frac{ w_{\rm I}^{2} }{2f}\right) - \frac{T_{\rm I}^2(\omega-\omega_0)^2}{4} \right]$',row)
row = gui.create_slider_with_latex(mainframe,r'Bullet radius $w_{\rm I} k_0=$',var_double[0],500,5000,row)
row = gui.create_slider_with_latex(mainframe,r'Bullet length $T_{\rm I}\omega_0=$',var_double[1],10,200,row)
row = gui.create_slider_with_latex(mainframe,r'Focal length $f / w_{\rm I}=$',var_double[2],50,500,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r'Propagation distance $z/f=$',var_double[3],0,1,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)