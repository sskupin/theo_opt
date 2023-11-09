import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
title = "Reflection and Transmission at Interface with Bloch medium"
root = Tk.Tk()
root.title(title)

def mTE(kfz,z):
    if np.abs(kfz) == 0:
        m = np.array([[1,z],[0,1]])
    else:
        m = np.array([[np.cos(kfz*z),np.sin(kfz*z)/kfz],[-np.sin(kfz*z)*kfz,np.cos(kfz*z)]])
    return m

def mTM(kfz,epsilon_f,z):
    if np.abs(kfz) == 0:
        m = np.array([[1,z*epsilon_f],[0,1]])
    else:
        m = np.array([[np.cos(kfz*z),np.sin(kfz*z)*epsilon_f/kfz],[-np.sin(kfz*z)*kfz/epsilon_f,np.cos(kfz*z)]])
    return m    

def reflection_transmission(epsilon_s,epsilon_f1,epsilon_f2,d1,d2,phi): # computing coefficients of reflection and transmission
    kx = np.sqrt(epsilon_s)*np.sin(phi)
    ksz = np.sqrt(epsilon_s-kx**2)*2*np.pi
    kf1z = np.sqrt(epsilon_f1+0j-kx**2)*2*np.pi
    kf2z = np.sqrt(epsilon_f2+0j-kx**2)*2*np.pi
    m1TE = mTE(kf1z,d1)
    m2TE = mTE(kf2z,d2)
    m1TM = mTM(kf1z,epsilon_f1,d1)
    m2TM = mTM(kf2z,epsilon_f2,d2)
    MTE = np.matmul(m2TE,m1TE)
    MTM = np.matmul(m2TM,m1TM)
    kzTE = np.arccos((MTE[1,1]+MTE[0,0])/2)/(d1+d2)
    kzTM = np.arccos((MTM[1,1]+MTM[0,0])/2)/(d1+d2)
    kappaTE = -1j*(np.exp(1j*kzTE*(d1+d2))-MTE[0,0])/MTE[0,1]
    if np.real(kappaTE) < 0:
        kappaTE = -1j*(np.exp(-1j*kzTE*(d1+d2))-MTE[0,0])/MTE[0,1]
    kappaTM = -1j*(np.exp(1j*kzTM*(d1+d2))-MTM[0,0])/MTM[0,1]
    if np.real(kappaTM) < 0:
        kappaTM = -1j*(np.exp(-1j*kzTM*(d1+d2))-MTM[0,0])/MTM[0,1]
    RTE = (ksz-kappaTE)/(ksz+kappaTE)
    RTM = (kappaTM-ksz/epsilon_s)/(ksz/epsilon_s+kappaTM) # for electric field (negative of magnetic coeff.)
        
    return RTE,RTM

def plot_subplot(ax,phi,curve_1,curve_2,label_1,label_2):
    ax.plot(phi,curve_1,'b',label=label_1)
    ax.plot(phi,curve_2,'r',label=label_2)
    ax.set_xticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2])
    ax.set_xticklabels([r'$0$', r'$\pi/8$', r'$\pi/4$', r'$3\pi/8$', r'$\pi/2$'])
    ax.set_xlabel(r'$\varphi_{\rm i}$')
    ax.set_xlim([0,np.pi/2])
    ax.set_ylabel(label_1+','+label_2)
    ax.legend()
    
def initialize():
    epsilon_s_string.set("1")
    d1_string.set("1")
    d2_string.set("1")
    epsilon_f1_string.set("2.25")
    epsilon_f2_string.set("2.2")
        
    calculate()

def calculate():
    try:
        epsilon_s = float(epsilon_s_string.get())
        d1 = float(d1_string.get())
        d2 = float(d2_string.get())
        D_string.set(d1+d2)
        epsilon_f1 = float(epsilon_f1_string.get())
        epsilon_f2 = float(epsilon_f2_string.get())
        
        if epsilon_s <= 0 or d1 <= 0 or d2 <= 0 or epsilon_f1 <= 0 or epsilon_f2 <= 0 or epsilon_f1 == epsilon_f2: 
            gui.input_error(initialize)
        else:
            f.clf()
            phi = np.linspace(0, np.pi/2, num=401, endpoint=False) # angle of incidence
            vreflection_transmission = np.vectorize(reflection_transmission)
            RTE,RTM = vreflection_transmission(epsilon_s,epsilon_f1,epsilon_f2,d1,d2,phi)
            a1 = f.add_subplot(131)
            plot_subplot(a1,phi,np.abs(RTE)**2,1-np.abs(RTE)**2,r'$\rho_{\rm TE}$',r'$\tau_{\rm TE}$')
            a1.set_ylim([-0.025,1.025])
            a2 = f.add_subplot(132)
            plot_subplot(a2,phi,np.abs(RTM)**2,1-np.abs(RTM)**2,r'$\rho_{\rm TM}$',r'$\tau_{\rm TM}$')
            a2.set_ylim([-0.025,1.025])
            a3 = f.add_subplot(133)
            plot_subplot(a3,phi,np.angle(RTE),np.angle(RTM),r'$\theta_{\rm TE}$',r'$\theta_{\rm TM}$')  
            a3.set_yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
            a3.set_yticklabels([r'$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])
            a3.set_ylim([-1.05*np.pi,1.05*np.pi])
            
            plt.tight_layout()
            
#            plt.savefig('Bloch_interface.pdf',bbox_inches='tight',dpi=300, transparent=True)

            canvas.draw()
    except ValueError: gui.input_error(initialize)

f = plt.figure(1,[7,2])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

epsilon_s_string = Tk.StringVar()
d1_string = Tk.StringVar()
d2_string = Tk.StringVar()
D_string = Tk.StringVar()
epsilon_f1_string = Tk.StringVar()
epsilon_f2_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_entry(mainframe,u"substrate: \u03B5 =",epsilon_s_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_title(mainframe,"unit cell parameters",row)
row = gui.create_entry(mainframe,u"film 1 thickness: d/\u03BB =",d1_string,row)
row = gui.create_entry(mainframe,u"film 1: \u03B5 =",epsilon_f1_string,row)
row = gui.create_entry(mainframe,u"film 2 thickness: d/\u03BB =",d2_string,row)
row = gui.create_entry(mainframe,u"film 2: \u03B5 =",epsilon_f2_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)