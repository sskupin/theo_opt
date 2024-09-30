import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import bpm_stuff as bpm

gui.set_rcParams()
root = Tk.Tk()
root.title("Pulse Envelope Propagation")

def initialize():
    k2_double.set(1)
    k3_double.set(0.1)
    kloss_double.set(0.1)
    pulse_string.set("Gaussian")   # profile  
    calculate()
    
def calculate():
    gui.change_cursor(root,"trek")
    Neta = 4096*4
    Leta = 1000 
    NZ = 201
    LZ = 1
    Nabs = 128
    k2 = k2_double.get()
    k3 = k3_double.get()
    kloss = kloss_double.get()
    profile = pulse_string.get()
        
    eta, delta_eta = np.linspace(-Leta/2,Leta/2,Neta,endpoint=False, retstep=True)
    if profile == 'sech':
        u0 =  1/np.cosh(eta)
    elif profile == 'Gaussian':
        u0 = np.exp(-eta**2)
    else:
        u0 = np.exp(-eta**8)
    Z, delta_Z = np.linspace(0,LZ,NZ,endpoint=True, retstep=True)
    
    u = bpm.propagation_nls(Neta,u0,delta_eta,NZ,delta_Z,Nabs,1,-2*np.pi*k2,0,np.pi*k3/3,4*np.pi*kloss) # factors adjusted do meet scaled nls, see bpm_stuff.py
    
    f.clf()
    
    a1 = plt.subplot2grid((5, 4), (1, 1), rowspan=4, colspan=2)
    im = a1.imshow(np.abs(np.transpose(u))**2 ,extent=[Z[0], Z[-1], eta[0], eta[-1]] , aspect='auto', origin='lower', vmin=0, cmap='jet')
    a1.set_xlabel(r'$z/L$')
    a1.set_ylabel(r'$\tau/T_{\rm p}$')
    ac = plt.subplot2grid((5, 4), (0, 1), colspan=2)
    f.colorbar(im, cax=ac, orientation='horizontal')
            
    a2 = plt.subplot2grid((5, 4), (1, 0), rowspan=4)       
    keta = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(Neta,delta_eta))
    U0 = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(u0))) # ifft because we assume time dependent problem
    lns1 = a2.plot(np.abs(U0)/np.max(np.abs(U0)),keta,'r', label=r'$|V_0(\overline\omega)|$')
    a2.set_xlabel(r'$|v_0|/\mathrm{max}|v_0|$, $|V_0|/\mathrm{max}|V_0|$')
    a2.set_ylabel(r'$\overline\omega T_{\rm p}$')
    a2.set_ylim([-20,20])
    a2.invert_yaxis()
    a2bis = a2.twinx()
    lns2 = a2bis.plot(np.abs(u0)/np.max(np.abs(u0)),eta, label=r'$|v_0(\tau)|$',color='b')
    a2bis.set_ylim([-10,10])
    a2bis.invert_yaxis()
    a2bis.set(yticklabels=[]) 
    lns = lns2+lns1
    labs = [l.get_label() for l in lns]
    a2.legend(lns, labs, bbox_to_anchor=(0, 1.05, 1, 0), loc="lower right")
    
    a3 = plt.subplot2grid((5, 4), (1, 3), rowspan=4, sharey=a1)
    lns1 = a3.plot(np.abs(u[-1,:])/np.max(np.abs(u0)),eta,'b', label=r'$|v(L,\tau)|$')
    a3.set_xlabel(r'$|v|/\mathrm{max}|v_0|, |V|/\mathrm{max}|V_0|$')
    a3.set_ylabel(r'$\tau/T_{\rm p}$')
    a3.set_ylim([-10,10])
    a3.invert_yaxis()   
    a3bis = a3.twinx()
    U = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(u[-1,:]))) # ifft because we assume time dependent problem
    lns2 = a3bis.plot(np.abs(U)/np.max(np.abs(U0)),keta,'r', label=r'$|V(L,\overline\omega)|$')
    a3bis.set_ylim([-20,20])
    a3bis.set_ylabel(r'$\overline\omega T_{\rm p}$')
    a3bis.invert_yaxis()
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    a3.legend(lns, labs, bbox_to_anchor=(0, 1.05, 1, 0), loc="lower right")
        
    plt.tight_layout()

    ac.set_position([0.35,0.825,0.3,0.025])  
    ac.xaxis.tick_top()
    ac.set_xlabel(r'$|v|^2$')
    ac.xaxis.set_label_position('top') 
     
#    plt.savefig('scalar_pulse_prop.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
    canvas.draw()
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[10,5])

canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

k2_double = Tk.DoubleVar()
k3_double = Tk.DoubleVar()
kloss_double = Tk.DoubleVar()
pulse_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_formula_with_latex(mainframe,r'$\mathrm{i}\partial_z v - \frac{k_0^{(2)}}{2}\partial^2_\tau v  = $',r'$\mathrm{i} \frac{k_0^{(3)}}{6}\partial^3_\tau v - \mathrm{i}k_0'' v$',row)
row = gui.create_radiobutton_single_column(mainframe,[u'Input beam profile:','sech','Gaussian','super-Gaussian'],pulse_string,3,row)
row = gui.create_slider_with_latex(mainframe,r'Normalized GVD $k_0^{(2)}L/T_{\rm p}^2=$',k2_double,-2,2,row,increment=0.05)
row = gui.create_slider_with_latex(mainframe,r'Normalized TOD $k_0^{(3)}L/T_{\rm p}^3=$',k3_double,-2,2,row,increment=0.05)
row = gui.create_slider_with_latex(mainframe,r'Normalized linear losses $k_0''L=$',kloss_double,0,1,row,increment=0.05)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)