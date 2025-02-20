import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import media as media
import strat_stuff as strat

gui.set_rcParams()
title = "Reflection of a Pulse at an Interface"
root = Tk.Tk()
root.title(title)

def reflection(epsilon_s,epsilon_c,phi): # computing coefficients of reflection for all frequencies
    kx,ksz,kcz = strat.KSC(epsilon_s,epsilon_c,phi)
    RTE,RTM,TTE,TTM,tauTE,tauTM = strat.RTAU(ksz,kcz,epsilon_s,epsilon_c,np.identity(2),np.identity(2))
    return RTE,RTM
    
def plot_subplot(ax,t,curves,labels,colors):
    for index in range(len(labels)):
        ax.plot(t,curves[index],colors[index],label=labels[index])
    ax.set_xlabel(r'$t$ [fs]')
    ax.set_ylabel(','.join(labels))
    ax.legend()

def initialize():
    duration_double.set(1.5) # pulse duration in optical cycles
    N_int.set(10)
    cladding_string.set("Ag")
    substrate_string.set("Vacuum")
    phi_double.set(0.25) # angle of incidence in units of pi
    lambda_double.set(800) # center wavelength in nm
        
    calculate()

def calculate():
    gui.change_cursor(root,"trek")
    lambda0 = lambda_double.get()
    duration = duration_double.get()
    phi = phi_double.get()*np.pi
    N = N_int.get()
    cladding = cladding_string.get()
    substrate = substrate_string.get()
        
    f.clf()
    omega0=media.omega2lambda(lambda0) # center frequency in 1/fs
    Tp = duration*2*np.pi/omega0 # pulse duration in fs
    Nt=4096
    Lt=64*Tp
    t,delta_t = np.linspace(-Lt/2, Lt/2, num=Nt, endpoint=False, retstep=True) # t in units of fs
    E_i = np.exp(-t**2/Tp**2) + 0j
    delta = 2*np.pi*np.fft.fftfreq(Nt,delta_t) # \omega-\omega_0 in units of 1/fs (shifted)
    omega=omega0+delta # (shifted)
    lambdav=media.omega2lambda(omega) # (shifted)
    lambdav=np.where(lambdav > 250, lambdav, 250)
    lambdav=np.where(lambdav < 12000, lambdav, 12000)
    epsilon_s=media.epsilon(substrate,lambdav)
    epsilon_c=media.epsilon(cladding,lambdav)
    RTE,RTM = reflection(epsilon_s,epsilon_c,phi)
    E_r = np.fft.fft(RTE**N*np.fft.ifft(E_i)) # inverse because temporal FFT   
    E_rx = np.fft.fft(RTM**N*np.fft.ifft(E_i)) # inverse because temporal FFT   
     
    a1 = f.add_subplot(221)
    plot_subplot(a1,t,[np.abs(E_i),np.abs(E_r)],[r'$\left|\tilde E_{\rm i}\right|$',r'$\left|\tilde E_{\rm r}\right|$'],['b','r'])
    a1.plot(t,np.real(E_i*np.exp(-1j*omega0*t)),'b:')
    a1.plot(t,np.real(E_r*np.exp(-1j*omega0*t)),'r:')
    a1.set_xlim([-4*Tp, 4*Tp])
    a1.set_title('TE polarization, '+str(N)+r' reflections, $\varphi_{\rm i}=$'+str(round(phi/np.pi,3))+r'$\pi$')

    a2 = f.add_subplot(222)
    plot_subplot(a2,t,[np.angle(E_i),np.angle(E_r)],[r'$\arg \tilde E_{\rm i}$',r'$\arg \tilde E_{\rm r}$'],['b','r'])
    a2.set_xlim([-4*Tp, 4*Tp])
    a2.set_title('TE polarization, '+str(N)+r' reflections, $\varphi_{\rm i}=$'+str(round(phi/np.pi,3))+r'$\pi$')
            
    a3 = f.add_subplot(223)
    plot_subplot(a3,t,[np.abs(E_i),np.abs(E_rx)],[r'$\left|\tilde E_{{\rm i}x}\right|$',r'$\left|\tilde E_{{\rm r}x}\right|$'],['b','r'])
    a3.plot(t,np.real(E_i*np.exp(-1j*omega0*t)),'b:')
    a3.plot(t,np.real(E_rx*np.exp(-1j*omega0*t)),'r:')
    a3.set_xlim([-4*Tp, 4*Tp])
    a3.set_title('TM polarization, '+str(N)+r' reflections, $\varphi_{\rm i}=$'+str(round(phi/np.pi,3))+r'$\pi$')
            
    a4 = f.add_subplot(224)
    plot_subplot(a4,t,[np.angle(E_i),np.angle(E_rx)],[r'$\arg \tilde E_{{\rm i}x}$',r'$\arg \tilde E_{{\rm r}x}$'],['b','r'])
    a4.set_xlim([-4*Tp, 4*Tp])
    a4.set_title('TM polarization, '+str(N)+r' reflections, $\varphi_{\rm i}=$'+str(round(phi/np.pi,3))+r'$\pi$')
            
    plt.tight_layout()  
            
#    plt.savefig('interface_pulse.pdf',bbox_inches='tight',dpi=300, transparent=True)

    canvas.draw()
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[8,6])

canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

cladding_string = Tk.StringVar()
substrate_string = Tk.StringVar()
phi_double = Tk.DoubleVar()
lambda_double = Tk.DoubleVar()
duration_double = Tk.DoubleVar()
N_int = Tk.IntVar()

initialize()

row = 1
row = gui.create_label(mainframe,r'substrate medium:',substrate_string,row)
row = gui.create_radiobutton(mainframe,['cladding medium:','Al','Ag'],cladding_string,2,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r"angle of incidence $\varphi_{\rm i}$ [$\pi$] =",phi_double,-0.475,0.475,row,increment=0.025)
row = gui.create_slider_with_latex(mainframe,r"center wavelength $\lambda$ [nm] =",lambda_double,450,1500,row,increment=25)
row = gui.create_slider_with_latex(mainframe,r"pulse duration $T_{\rm p}$ [$\lambda/c$] =",duration_double,1,5,row,increment=0.1)
row = gui.create_intslider_with_latex(mainframe,r"number of reflections =",N_int,1,10,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)