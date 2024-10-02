import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import bpm_stuff as bpm

gui.set_rcParams()
root = Tk.Tk()
root.title("Modulational Instability (NLS)")

def initialize():
    LZ_double.set(2)
    N_double.set(1.15)
    var_string[0].set("+1")     # D
    var_string[1].set("+1")     # Gamma    
    calculate()
    
def calculate():
    gui.change_cursor(root,"trek")
    Neta = 1024
    Leta = 40 
    NZ = 20001
    LZ = LZ_double.get()
    Nabs = 0
    N = N_double.get()
    D = float(var_string[0].get())    
    Gamma = float(var_string[1].get())
        
    eta, delta_eta = np.linspace(-Leta/2,Leta/2,Neta,endpoint=False, retstep=True)
    rng = np.random.default_rng()
    u0 = np.ones(Neta) + 0.001 * rng.standard_normal(Neta)
    Z, delta_Z = np.linspace(0,LZ,NZ,endpoint=True, retstep=True)
    
    u = bpm.propagation_nls(Neta,u0,delta_eta,NZ,delta_Z*2*np.pi**2,Nabs,N,D,Gamma,0,0) # multiplication of delta_z by 2pi^2 to take into account scaling of LZ
    
    f.clf()
    
    a1 = plt.subplot2grid((5, 4), (1, 1), rowspan=4, colspan=2)
    im = a1.imshow(np.abs(np.transpose(u))**2 ,extent=[Z[0], Z[-1], eta[0], eta[-1]] , aspect='auto', origin='lower', vmin=0, cmap='jet')
    a1.set_xlabel(r'$Z/\pi$')
    a1.set_ylabel(r'$\eta$')
    ac = plt.subplot2grid((5, 4), (0, 1), colspan=2)
    f.colorbar(im, cax=ac, orientation='horizontal')
    
    a2 = plt.subplot2grid((5, 4), (1, 0), rowspan=4)       
    keta = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(Neta,delta_eta))
    U0 = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(u0))) # ifft because we assume time dependent problem
    lns1 = a2.plot(np.abs(U0)/np.max(np.abs(U0)),keta,'r', label=r'$|U_0(\mu)|$')
    a2.set_xlabel(r'$|u_0|/N$, $|U_0|/\mathrm{max}|U_0|$')
    a2.set_ylabel(r'$\mu$')
    a2.set_ylim([-20,20])
    a2.invert_yaxis()
    a2bis = a2.twinx()
    lns2 = a2bis.plot(np.abs(u0)/np.max(np.abs(u0)),eta, label=r'$|u_0(\eta)|$',color='b')
    a2bis.set_ylim([eta[0],eta[-1]])
    a2bis.invert_yaxis()
    a2bis.set(yticklabels=[]) 
    lns = lns2+lns1
    labs = [l.get_label() for l in lns]
    a2.legend(lns, labs, bbox_to_anchor=(0, 1.05, 1, 0), loc="lower right")
    
    a3 = plt.subplot2grid((5, 4), (1, 3), rowspan=4, sharey=a1)
    lns1 = a3.plot(np.abs(u[-1,:])/np.max(np.abs(u0)),eta,'b', label=r'$|u(\eta,Z_L)|$')
    a3.set_xlabel(r'$|u|/N, |U|/\mathrm{max}|U_0|$')
    a3.set_ylabel(r'$\eta$')
    a3.invert_yaxis()   
    a3bis = a3.twinx()
    U = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(u[-1,:]))) # ifft because we assume time dependent problem
    lns2 = a3bis.plot(np.abs(U)/np.max(np.abs(U0)),keta,'r', label=r'$|U(\mu,Z_L))|$')
    a3bis.set_ylim([-20,20])
    a3bis.set_ylabel(r'$\mu$')
    a3bis.invert_yaxis()
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    a3.legend(lns, labs, bbox_to_anchor=(0, 1.05, 1, 0), loc="lower right")
        
    plt.tight_layout()

    ac.set_position([0.35,0.825,0.3,0.025])  
    ac.xaxis.tick_top()
    ac.set_xlabel(r'$|u|^2/N^2$')
    ac.xaxis.set_label_position('top') 
     
#    plt.savefig('mi.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
    canvas.draw()
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[10,5])

canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

LZ_double = Tk.DoubleVar()
N_double = Tk.DoubleVar()
var_string = gui.create_stringvar_vector(2)

initialize()

row = 1
row = gui.create_slider_with_latex(mainframe,r'Propagation length $Z_L/\pi=$',LZ_double,0.5,3,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_formula_with_latex(mainframe,r'$\partial_Z u - \mathrm{i}\frac{D}{2}\partial^2_\eta u = $',r'$\mathrm{i} \Gamma |u|^2 u$',row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r'Amplitude (Soliton order) $N=$',N_double,0.5,2,row)
row = gui.create_radiobutton(mainframe,['Sign of D:','+1','-1'],var_string[0],2,row)
row = gui.create_radiobutton(mainframe,[u'Sign of \u0393:','+1','-1'],var_string[1],2,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)