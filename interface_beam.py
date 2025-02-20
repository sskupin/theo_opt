import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import strat_stuff as strat

gui.set_rcParams()
title = "Reflection of a 1D beam at an Interface"
root = Tk.Tk()
root.title(title)

def tilted_Gaussian_beam(epsilon_s,beam_width,phi,Nx,Lx): # compute envelope E_i in the z=0 plane for tilted incident beam
    Nx_tilted = 512
    Lx_tilted = 128*beam_width
    x_tilted,delta_x_tilted = np.linspace(-Lx_tilted/2, Lx_tilted/2, num=Nx_tilted, endpoint=False, retstep=True) # x_tilted in units of lambda
    kx_tilted = np.fft.fftshift(2*np.pi*np.fft.fftfreq(Nx_tilted,delta_x_tilted)) # kx_tilted in units of 1/lambda (not shifted)
    ks=2*np.pi*np.sqrt(epsilon_s)
    FTE_tilted = np.fft.fftshift(np.fft.fft(np.fft.fftshift(np.exp(-x_tilted**2/beam_width**2) + 0j))) # spectrum in titled system (not shifted)
    FTE_tilted = np.where(kx_tilted**2<ks**2,FTE_tilted,0+0j) # remove any evanescent amplitudes
    x,delta_x = np.linspace(-Lx/2, Lx/2, num=Nx, endpoint=False, retstep=True) # x in units of lambda
    E_i = np.dot(np.exp(1j*np.sin(phi)*np.outer(x,np.sqrt(ks**2-kx_tilted**2+0j))+1j*np.cos(phi)*np.outer(x,kx_tilted)),FTE_tilted) # field in the z=0 plane
    E_i = E_i*np.exp(-1j*ks*np.sin(phi)*x) # beam envelope 
    E_i = np.where(x**2<x_tilted[0]**2/np.cos(phi),E_i,0+0j) # remove next period
    E_i = E_i/np.amax(np.abs(E_i))
    return x,delta_x,E_i

def reflection(epsilon_s,epsilon_c,kx): # computing coefficients of reflection for all frequencies
    ksz = np.sqrt(epsilon_s*(2*np.pi)**2+0j-kx**2)
    kcz = np.sqrt(epsilon_c*(2*np.pi)**2-kx**2)    
    RTE,RTM,TTE,TTM,tauTE,tauTM = strat.RTAU(ksz,kcz,epsilon_s,epsilon_c,np.identity(2),np.identity(2))
    RTE = np.where(epsilon_s*(2*np.pi)**2-kx**2>0,RTE,0)
    RTM = np.where(epsilon_s*(2*np.pi)**2-kx**2>0,RTM,0)
    return RTE,RTM
    
def plot_subplot(ax,t,curves,labels,colors):
    for index in range(len(labels)):
        ax.plot(t,curves[index],colors[index],label=labels[index])
    ax.set_xlabel(r'$x$ [$\lambda$]')
    ax.set_ylabel(','.join(labels)+' [norm. u.]')
    ax.legend()

def initialize():
    var_string[0].set("2.73") # epsilon_s
    var_string[1].set("2.11") # epsilon_c_real
    var_string[2].set("0") # epsilon_c_imag
    beam_width_double.set(5) # beam width in the plane orthogonal to mean wave vector at x=z=0 in units of lambda
    phi_double.set(0.375) # angle of incidence in units of pi
    gui.copy_stringvar_vector(var_string,var_save)  
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()

def calculate():
    gui.change_cursor(root,"trek")
    try:
        epsilon_s = float(var_string[0].get())
        epsilon_c_real = float(var_string[1].get())
        epsilon_c_imag = float(var_string[2].get())
        beam_width = beam_width_double.get()
        phi = phi_double.get()*np.pi
        
        if epsilon_s <= 1:
            gui.input_error("Substrate epsilon must not be smaller than one. Re-initializing with previous parameters...",reinitialize)
        elif epsilon_c_real == 0 and epsilon_c_imag == 0: 
            gui.input_error("Cladding epsilon must not be zero. Re-initializing with previous parameters...",reinitialize)
        elif beam_width < 0 or np.abs(phi)>=np.pi/2:
            gui.input_error("bla",reinitialize)
        else:
            f.clf()
            epsilon_c = epsilon_c_real + 1j*epsilon_c_imag
            beam_width_z=beam_width/np.cos(phi) # approximate beam width in the z=0 plane
            Nx=2096
            Lx=64*beam_width_z
            x,delta_x,E_i = tilted_Gaussian_beam(epsilon_s,beam_width,phi,Nx,Lx)
            delta = 2*np.pi*np.fft.fftfreq(Nx,delta_x) # k_x-k_{x0} in units of 1/lambda (shifted)
            kx=2*np.pi*np.sqrt(epsilon_s)*np.sin(phi)+delta # k_x in units of 1/lambda (shifted)
            RTE,RTM = reflection(epsilon_s,epsilon_c,kx)
            E_r = np.fft.ifft(RTE*np.fft.fft(E_i))  
            E_rx = np.fft.ifft(RTM*np.fft.fft(E_i)) 
     
            a1 = f.add_subplot(221)
            plot_subplot(a1,x,[np.abs(E_i),np.abs(E_r)],[r'$\left|\tilde E_{\rm i}\right|$',r'$\left|\tilde E_{\rm r}\right|$'],['b','r'])
            a1.set_xlim([-4*beam_width_z, 4*beam_width_z])
            a1.set_title(r'TE polarization, $\varphi_{\rm i0}=$'+str(round(phi/np.pi,3))+r'$\pi$')

            a2 = f.add_subplot(222)
            plot_subplot(a2,x,[np.angle(E_i),np.angle(E_r)],[r'$\arg \tilde E_{\rm i}$',r'$\arg \tilde E_{\rm r}$'],['b','r'])
            a2.set_xlim([-4*beam_width_z, 4*beam_width_z])
            a2.set_title(r'TE polarization, $\varphi_{\rm i0}=$'+str(round(phi/np.pi,3))+r'$\pi$')
            
            a3 = f.add_subplot(223)
            plot_subplot(a3,x,[np.abs(E_i),np.abs(E_rx)],[r'$\left|\tilde E_{{\rm i}x}\right|$',r'$\left|\tilde E_{{\rm r}x}\right|$'],['b','r'])
            a3.set_xlim([-4*beam_width_z, 4*beam_width_z])
            a3.set_title(r'TM polarization, $\varphi_{\rm i0}=$'+str(round(phi/np.pi,3))+r'$\pi$')
            
            a4 = f.add_subplot(224)
            plot_subplot(a4,x,[np.angle(E_i),np.angle(E_rx)],[r'$\arg \tilde E_{{\rm i}x}$',r'$\arg \tilde E_{{\rm r}x}$'],['b','r'])
            a4.set_xlim([-4*beam_width_z, 4*beam_width_z])
            a4.set_title(r'TM polarization, $\varphi_{\rm i0}=$'+str(round(phi/np.pi,3))+r'$\pi$')
            
            plt.tight_layout()  
            
#            plt.savefig('interface_1Dbeam.pdf',bbox_inches='tight',dpi=300, transparent=True)

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[8,6])

canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(3)
var_save = gui.create_stringvar_vector(3)
phi_double = Tk.DoubleVar()
beam_width_double = Tk.DoubleVar()

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"substrate $\varepsilon_{\rm s} =$",var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r"cladding $\varepsilon_{\rm c}' =$",var_string[1],row)
row = gui.create_entry_with_latex(mainframe,r"cladding $\varepsilon_{\rm c}'' =$",var_string[2],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r"angle of incidence $\varphi_{\rm i}$ [$\pi$] =",phi_double,-0.4,0.4,row,increment=0.025)
row = gui.create_slider_with_latex(mainframe,r"beam width [$\lambda$] =",beam_width_double,2,10,row,increment=0.1)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)