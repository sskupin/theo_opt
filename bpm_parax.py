import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import bpm_stuff as bpm

gui.set_rcParams()
root = Tk.Tk()
root.title("Paraxial scalar BPM in homogeneous media")

def initialize():
    global var_save
    var_string[0].set("128")   # Nx
    var_string[1].set("25")    # Lx in w0
    var_string[2].set("100")   # Nz
    var_string[3].set("2")     # Lz in LF 
    var_string[4].set("0")     # Nabs
    var_string[5].set("1")     # alpha
    gui.copy_stringvar_vector(var_string,var_save)
    calculate()
    
def reinitialize():
    global var_string
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()
    
def calculate():
    global var_save
    try:
        Nx = int(float(var_string[0].get()))
        var_string[0].set(Nx)
        Lx = float(var_string[1].get())
        Nz = int(float(var_string[2].get())) + 1 # Number of points = Number of steps + 1
        var_string[2].set(Nz-1)
        Lz = float(var_string[3].get())
        Nabs = int(float(var_string[4].get()))
        var_string[4].set(Nabs)
        alpha = float(var_string[5].get())    
        
        if Nx <= 1 or Lx <= 0 or Nz <= 1 or Lz <= 0: 
            gui.input_error("Box lengths or numbers of points invalid. Re-initializing ...",reinitialize)
        elif Nabs < 0 or 20*Nabs > Nx: 
            gui.input_error("Number of absober points must be nonzero and must not exceed Nx/20. Re-initializing ...",reinitialize) 
        elif alpha < 1: 
            gui.input_error("Order of Super-Gaussian should be larger or equal one. Re-initializing ...",reinitialize)  
        else:

            x, delta_x = np.linspace(-Lx/2,Lx/2,Nx,endpoint=False, retstep=True)
            v0 =  np.exp(-(x**2)**alpha)
            z, delta_z = np.linspace(0,Lz,Nz,endpoint=True, retstep=True)
    
            v = bpm.propagation_parax(Nx,v0,delta_x,Nz,delta_z*2*np.pi,Nabs) # multiplication of delta_z by 2pi to take into account scaling to LF
    
            f.clf()
    
            a1 = plt.subplot2grid((1, 4), (0, 1), colspan=2)
            a1.imshow(np.abs(np.transpose(v)) ,extent=[z[0], z[-1], x[0], x[-1]] , aspect='auto', origin='lower', cmap='jet')
            a1.set_xlabel(r'$z/L_{\rm F}$')
            a1.set_ylabel(r'$x/w_0$')
            
            a2 = plt.subplot2grid((1, 4), (0, 0))       
            kx = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(Nx,delta_x))
            V0 = np.fft.fftshift(np.fft.fft(np.fft.fftshift(v0)))
            lns1 = a2.plot(np.abs(V0)/np.max(np.abs(V0)),kx,'r--', label=r'$|V_0(k_x)|$')
            a2.set_xlabel(r'$|v_0|$, $|V_0|/\mathrm{max}|V_0|$')
            a2.set_ylabel(r'$k_x w_0$')
            a2.set_ylim([kx[0],kx[-1]])
            a2.invert_yaxis()
            a2bis = a2.twinx()
            lns2 = a2bis.plot(np.abs(v0)/np.max(np.abs(v0)),x, label=r'$|v_0(x)|$',color='b')
            a2bis.set_ylim([x[0],x[-1]])
            a2bis.invert_yaxis()
            a2bis.set(yticklabels=[]) 
            lns = lns2+lns1
            labs = [l.get_label() for l in lns]
            a2.legend(lns, labs, loc=1)
    
            a3 = plt.subplot2grid((1, 4), (0, 3), sharey=a1)
            lns1 = a3.plot(np.abs(v[-1,:])/np.max(np.abs(v0)),x,'b', label=r'$|v(x,z=L_z)|$')
            a3.set_xlabel(r'$|v|, |V|/\mathrm{max}|V_0|$')
            a3.set_ylabel(r'$x/w_0$')
            a3.invert_yaxis()   
            a3bis = a3.twinx()
            V = np.fft.fftshift(np.fft.fft(np.fft.fftshift(v[-1,:])))
            lns2 = a3bis.plot(np.abs(V)/np.max(np.abs(V0)),kx,'r--', label=r'$|V(k_x,z=L_z))|$')
            a3bis.set_ylim([kx[0],kx[-1]])
            a3bis.set_ylabel(r'$k_x w_0$')
            a3bis.invert_yaxis()
            lns = lns1+lns2
            labs = [l.get_label() for l in lns]
            a3.legend(lns, labs, loc=1)
            
            plt.tight_layout()
                
#            plt.savefig('bpm_parax.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()       
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", initialize)

f = plt.figure(1,[10,4])

canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(6)
var_save = gui.create_stringvar_vector(6)

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"Transverse number of points $N_x=$",var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r"Transverse box length $L_x/w_0=$",var_string[1],row)
row = gui.create_entry_with_latex(mainframe,r"Longitudinal number of steps $N_z=$",var_string[2],row)
row = gui.create_entry_with_latex(mainframe,r"Longitudinal box length $L_z/L_{\rm F}=$",var_string[3],row)
row = gui.create_entry_with_latex(mainframe,r"Number of absorber points $N_{\rm abs}=$",var_string[4],row)
row = gui.create_formula_with_latex(mainframe,r'$\partial_z v = $',r'$\mathrm{i}\frac{1}{2k}\partial^2_x v $',row)
row = gui.create_formula_with_latex(mainframe,r'$v_0=$',r'$\exp\!\left[ -\left(\frac{x}{w_0}\right)^{2\alpha} \right]$',row)
row = gui.create_entry_with_latex(mainframe,r"Degree of super-Gaussian $\alpha=$",var_string[5],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)