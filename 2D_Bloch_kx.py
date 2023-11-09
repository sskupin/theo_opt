import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
title = "2D DR of Bloch modes for fixed frequency"
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

def DR(d1,epsilon_f1,d2,epsilon_f2,kx,polarization): # computing dispersion relation for Bloch modes
    kf1z = np.sqrt(epsilon_f1+0j-kx**2)*2*np.pi
    kf2z = np.sqrt(epsilon_f2+0j-kx**2)*2*np.pi
    if polarization == 'TE':
        m1 = mTE(kf1z,d1)
        m2 = mTE(kf2z,d2)
    else:
        m1 = mTM(kf1z,epsilon_f1,d1)
        m2 = mTM(kf2z,epsilon_f2,d2)
    M = np.matmul(m2,m1)
    Kz = np.arccos((M[1,1]+M[0,0])/2)/(2*np.pi)
    Kx = kx*(d1+d2)
    
    return Kz,Kx
    
def initialize():
    d1_string.set("1")
    d2_string.set("1")
    epsilon_f1_string.set("2.25")
    epsilon_f2_string.set("2.2")
    polarization_string.set("TE")
    cone_string.set("n")
   
    calculate()

def calculate():
    try:
        d1 = float(d1_string.get())
        d2 = float(d2_string.get())
        D_string.set(d1+d2)
        epsilon_f1 = float(epsilon_f1_string.get())
        epsilon_f2 = float(epsilon_f2_string.get())
        polarization = polarization_string.get()
        cone = cone_string.get()     
        
        if d1 <= 0 or d2 < 0 or epsilon_f1 <= 0 or epsilon_f2 <= 0 or np.maximum(np.sqrt(epsilon_f1),np.sqrt(epsilon_f2))*(d1+d2) > 3:
            gui.input_error(initialize)
        else:
            f.clf()
            kx = np.linspace(0, np.maximum(np.sqrt(epsilon_f1),np.sqrt(epsilon_f2)), num=10001, endpoint=True) # transverse wavevector in 2\pi/\lambda
            vDR = np.vectorize(DR)
            Kz,Kx = vDR(d1,epsilon_f1,d2,epsilon_f2,kx,polarization)
            a1 = f.add_subplot(111, aspect='equal')
            a1.plot(np.real(Kz),Kx,'b',np.real(-Kz),Kx,'b')
            a1.plot(np.real(Kz-1),Kx,'b',np.real(-Kz+1),Kx,'b')
            a1.set_xlim([-1, 1])
            a1.set_ylim([Kx[0], Kx[-1]])
            a1.fill_betweenx(Kx, -1, 1, where=np.imag(Kz)!=0., facecolor='red', alpha=1, zorder=100, interpolate=True)
            a1.set_xlabel(r'$K_z^{\prime}$')
            a1.set_ylabel(r'$K_x$')
            if cone=='y':
                a1.fill_betweenx(Kx, -1, 1, where=Kx>=d1+d2, facecolor='gray', alpha=.5, zorder=100, interpolate=True)
                    
            plt.tight_layout()  
            
#            plt.savefig('2D_Bloch_kx.pdf',bbox_inches='tight',dpi=300, transparent=True)

            canvas.draw()
    except ValueError: gui.input_error(initialize)

f = plt.figure(1,[6,4])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

d1_string = Tk.StringVar()
d2_string = Tk.StringVar()
D_string = Tk.StringVar()
epsilon_f1_string = Tk.StringVar()
epsilon_f2_string = Tk.StringVar()
polarization_string = Tk.StringVar()
cone_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_title(mainframe,"unit cell parameters",row)
row = gui.create_entry(mainframe,u"film 1 thickness: d/\u03BB =",d1_string,row)
row = gui.create_entry(mainframe,u"film 1: \u03B5 =",epsilon_f1_string,row)
row = gui.create_entry(mainframe,u"film 2 thickness: d/\u03BB =",d2_string,row)
row = gui.create_entry(mainframe,u"film 2: \u03B5 =",epsilon_f2_string,row)
row = gui.create_label(mainframe,u"\u03A9 = \u039b/\u03BB =",D_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_radiobutton(mainframe,['polarization:','TE','TM'],polarization_string,2,row)
row = gui.create_checkbutton(mainframe,"vacuum light cone",'n','y',cone_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)