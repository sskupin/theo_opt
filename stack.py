# merge with FP_mirror

import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import strat_stuff as strat

gui.set_rcParams()
title = "Reflection and Transmission at a Stack"
root = Tk.Tk()
root.title(title)

def mTE(kfz,z):
    return np.array([[np.cos(kfz*2*np.pi*z),np.sin(kfz*2*np.pi*z)/kfz],[-np.sin(kfz*2*np.pi*z)*kfz,np.cos(kfz*2*np.pi*z)]])

def mTM(kfz,epsilon_f,z):
    return np.array([[np.cos(kfz*2*np.pi*z),np.sin(kfz*2*np.pi*z)*epsilon_f/kfz],[-np.sin(kfz*2*np.pi*z)*kfz/epsilon_f,np.cos(kfz*2*np.pi*z)]])

def reflection_transmission(epsilon_s,d1,epsilon_f1,d2,epsilon_f2,N,epsilon_c,phi): # computing coefficients of reflection and transmission
    kx,ksz,kcz = strat.KSC(epsilon_s,epsilon_c,phi)
    MTE,MTM = strat.MP(kx,d1,epsilon_f1,d2,epsilon_f2,N)
    RTE,RTM,TTE,TTM,tauTE,tauTM = strat.RTAU(ksz,kcz,epsilon_s,epsilon_c,MTE,MTM)   
    return RTE,RTM,tauTE,tauTM

def plot_subplot(ax,phi,curves,labels,colors,phi_min, phi_max):
    for index in range(len(labels)):
        ax.plot(phi,curves[index],colors[index],label=labels[index])
    if np.floor(8*phi_max/np.pi)-np.ceil(8*phi_min/np.pi) >= 1:
        ax.set_xticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2])
        ax.set_xticklabels([r'$0$', r'$\pi/8$', r'$\pi/4$', r'$3\pi/8$', r'$\pi/2$'])
    ax.set_xlabel(r'$\varphi_{\rm i}$')
    ax.set_xlim([phi_min, phi_max])
    ax.set_ylabel(','.join(labels))
    ax.legend()

def plot_amplitude(ax,M,F,G,z):
    Fplot = M[0,0]*F+M[0,1]*G
    Gplot = M[1,0]*F+M[1,1]*G
    FGabs = np.abs(Fplot)+np.abs(Gplot)
    if np.any(np.less(FGabs,1.e-6)):
        Fplot[np.where(FGabs<1.e-6)[0][0]:]=0
        Gplot[np.where(FGabs<1.e-6)[0][0]:]=0
    F=Fplot[-1]
    G=Gplot[-1]
    ax.plot(z,np.abs(Fplot),'b')
    
    return F,G
        
def initialize():
    d1_string.set("1")
    d2_string.set("1")
    epsilon_s_string.set("1")
    epsilon_f1_real_string.set("2.25")
    epsilon_f1_imag_string.set("0")
    epsilon_f2_real_string.set("2.2")
    epsilon_f2_imag_string.set("0")
    N_string.set("10")
    epsilon_c_real_string.set("1")
    epsilon_c_imag_string.set("0")
    polarization_string.set("TE")
    phi_0_string.set("0")
    phi_min_string.set("0")
    phi_max_string.set(".5")
        
    calculate()

def calculate():
    try:
        epsilon_s = float(epsilon_s_string.get())
        d1 = float(d1_string.get())
        epsilon_f1_real = float(epsilon_f1_real_string.get())
        epsilon_f1_imag = float(epsilon_f1_imag_string.get())
        d2 = float(d2_string.get())
        epsilon_f2_real = float(epsilon_f2_real_string.get())
        epsilon_f2_imag = float(epsilon_f2_imag_string.get())        
        N = int(N_string.get())
        epsilon_c_real = float(epsilon_c_real_string.get())
        epsilon_c_imag = float(epsilon_c_imag_string.get())
        polarization = polarization_string.get()
        phi_0 = float(phi_0_string.get())*np.pi
        phi_min = float(phi_min_string.get())*np.pi
        phi_max = float(phi_max_string.get())*np.pi
        
        if epsilon_c_imag < 0 or epsilon_c_real == 0 or epsilon_f1_real == 0 or epsilon_f2_real == 0 or epsilon_s <= 0\
                              or d1 <= 0 or d2 < 0 or N <= 0 or N > 50 or phi_max > np.pi/2 or phi_min < 0 or phi_min >= phi_max:
            gui.input_error(initialize)
        else:
            f.clf()
            phi = np.linspace(phi_min, phi_max, num=10001, endpoint=False) # angle of incidence
            if phi_0 < phi_min or phi_0 > phi_max:
                phi_0 = phi_min
                phi_0_string.set(phi_0/np.pi)
            epsilon_c = epsilon_c_real + 1j*epsilon_c_imag
            epsilon_f1 = epsilon_f1_real + 1j*epsilon_f1_imag
            epsilon_f2 = epsilon_f2_real + 1j*epsilon_f2_imag
            vreflection_transmission = np.vectorize(reflection_transmission)
            RTE,RTM,tauTE,tauTM = vreflection_transmission(epsilon_s,d1,epsilon_f1,d2,epsilon_f2,N,epsilon_c,phi)
            a1 = f.add_subplot(221)
            a2 = f.add_subplot(222)
            if epsilon_f1_imag == 0 and epsilon_f2_imag == 0:
                plot_subplot(a1,phi,[np.abs(RTE)**2,tauTE],[r'$\rho_{\rm TE}$',r'$\tau_{\rm TE}$'],['b','r'],phi_min, phi_max)
                plot_subplot(a2,phi,[np.abs(RTM)**2,tauTM],[r'$\rho_{\rm TM}$',r'$\tau_{\rm TM}$'],['b','r'],phi_min, phi_max)
            else:
                plot_subplot(a1,phi,[np.abs(RTE)**2,tauTE,1-np.abs(RTE)**2-tauTE],[r'$\rho_{\rm TE}$',r'$\tau_{\rm TE}$',r'$a_{\rm TE}$'],['b','r','g'],phi_min, phi_max)
                plot_subplot(a2,phi,[np.abs(RTM)**2,tauTM,1-np.abs(RTM)**2-tauTM],[r'$\rho_{\rm TM}$',r'$\tau_{\rm TM}$',r'$a_{\rm TM}$'],['b','r','g'],phi_min, phi_max)
                if polarization == 'TE':
                    a1.plot(phi_0,np.interp(phi_0,phi,1-np.abs(RTE)**2-tauTE),'go')
                else:
                    a2.plot(phi_0,np.interp(phi_0,phi,1-np.abs(RTM)**2-tauTM),'go')
            if polarization == 'TE':
                a1.plot(phi_0,np.interp(phi_0,phi,np.abs(RTE)**2),'bo')
                a1.plot(phi_0,np.interp(phi_0,phi,tauTE),'ro')
            else:
                a2.plot(phi_0,np.interp(phi_0,phi,np.abs(RTM)**2),'bo')
                a2.plot(phi_0,np.interp(phi_0,phi,tauTM),'ro')
            left, right = a1.get_ylim()
            a1.set_ylim([min(left,-0.025),max(right,1.025)])
            left, right = a2.get_ylim()
            a2.set_ylim([min(left,-0.025),max(right,1.025)])
            
            a3 = f.add_subplot(212)
            zs = np.linspace(-np.maximum(0.75*N*(d1+d2)/2,1), 0, num=501, endpoint=True)
            zf1 = np.linspace(0, d1, num=101, endpoint=True)
            zf2 = np.linspace(0, d2, num=101, endpoint=True) 
            zc = np.linspace(0, np.maximum(0.75*N*(d1+d2)/2,1), num=501, endpoint=True)
            for index in range(N):
                a3.axvspan(index*(d1+d2), (index+d1/(d1+d2))*(d1+d2), color='0.75')
                a3.axvspan((index+d1/(d1+d2))*(d1+d2), (index+1)*(d1+d2), color='0.875')
            a3.set_xlabel(r'$z/\lambda$')
            a3.set_xlim([zs[0],N*(d1+d2)+zc[-1]])
            RTE,RTM,tauTE,tauTM = reflection_transmission(epsilon_s,d1,epsilon_f1,d2,epsilon_f2,N,epsilon_c,phi_0)
            if polarization == 'TE':
                rho = np.abs(RTE)**2
                tau = tauTE
                theta = np.angle(RTE)
                F = 1+RTE
                G = 1j*np.sqrt(epsilon_s)*np.cos(phi_0)*(1-RTE)
                MTE = mTE(np.sqrt(epsilon_s)*np.cos(phi_0),zs)
                a3.plot(zs,np.abs(MTE[0,0]*F+MTE[0,1]*G),'b')
                for index in range(N):
                    MTE = mTE(np.sqrt(epsilon_f1-epsilon_s*np.sin(phi_0)**2),zf1)
                    F,G = plot_amplitude(a3,MTE,F,G,index*(d1+d2)+zf1)
                    MTE = mTE(np.sqrt(epsilon_f2-epsilon_s*np.sin(phi_0)**2),zf2)
                    F,G = plot_amplitude(a3,MTE,F,G,(index+d1/(d1+d2))*(d1+d2)+zf2)
                a3.plot(N*(d1+d2)+zc,np.abs(F)*np.exp(-np.imag(np.sqrt(epsilon_c-epsilon_s*np.sin(phi_0)**2))*2*np.pi*zc),'b')
                a3.set_ylabel(r'$| E_\mathrm{TE} | / | E_\mathrm{TE}^\mathrm{i} |$ for $\varphi_\mathrm{i}=$ '+str(round(phi_0/np.pi,4))+r'$\pi$')
            else:
                rho = np.abs(RTM)**2
                tau = tauTM
                theta = np.angle(RTM)
                F = 1-RTM
                G = 1j*np.sqrt(epsilon_s)*np.cos(phi_0)*(1+RTM)/epsilon_s
                MTM = mTM(np.sqrt(epsilon_s)*np.cos(phi_0),epsilon_s,zs)
                a3.plot(zs,np.abs(MTM[0,0]*F+MTM[0,1]*G),'b')
                for index in range(N):
                    MTM = mTM(np.sqrt(epsilon_f1-epsilon_s*np.sin(phi_0)**2),epsilon_f1,zf1)
                    F,G = plot_amplitude(a3,MTM,F,G,index*(d1+d2)+zf1)
                    MTM = mTM(np.sqrt(epsilon_f2-epsilon_s*np.sin(phi_0)**2),epsilon_f2,zf2)
                    F,G = plot_amplitude(a3,MTM,F,G,(index+d1/(d1+d2))*(d1+d2)+zf2)
                a3.plot(N*(d1+d2)+zc,np.abs(F)*np.exp(-np.imag(np.sqrt(epsilon_c-epsilon_s*np.sin(phi_0)**2))*2*np.pi*zc),'b')                
                a3.set_ylabel(r'$| H_\mathrm{TM} | / | H_\mathrm{TM}^\mathrm{i} |$ for $\varphi_\mathrm{i}=$ '+str(round(phi_0/np.pi,4))+r'$\pi$')
            ylimits = a3.get_ylim()
            a3.set_ylim([ylimits[0],ylimits[1]+0.1*(ylimits[1]-ylimits[0])])
            a3.text(0.5, 0.925, r'stack', verticalalignment='center', horizontalalignment='center', transform=a3.transAxes)
            a3.text(0.1, 0.925, r'substrate', verticalalignment='center', horizontalalignment='center', transform=a3.transAxes)
            a3.text(0.9, 0.925, r'cladding', verticalalignment='center', horizontalalignment='center', transform=a3.transAxes)
            a3.set_title(r'$\rho=$ '+str(round(rho,3))+r', $\tau=$ '+str(round(tau,3))+r', $\theta=$ '+str(round(theta/np.pi,3))+r'$\pi$')
            
            plt.tight_layout()  
            
#            plt.savefig('stack.pdf',bbox_inches='tight',dpi=300, transparent=True)

            canvas.draw()
    except ValueError: gui.input_error(initialize)

f = plt.figure(1,[8,5])
canvas = gui.create_canvas(root,f)
canvas.draw() # for faster feedback to user on startup
mainframe = gui.create_mainframe(root)

epsilon_s_string = Tk.StringVar()
d1_string = Tk.StringVar()
epsilon_f1_real_string = Tk.StringVar()
epsilon_f1_imag_string = Tk.StringVar()
d2_string = Tk.StringVar()
epsilon_f2_real_string = Tk.StringVar()
epsilon_f2_imag_string = Tk.StringVar()
N_string = Tk.StringVar()
epsilon_c_real_string = Tk.StringVar()
epsilon_c_imag_string = Tk.StringVar()
polarization_string = Tk.StringVar()
phi_0_string = Tk.StringVar()
phi_min_string = Tk.StringVar()
phi_max_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_entry(mainframe,u"substrate: \u03B5 =",epsilon_s_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"film 1 thickness: d/\u03BB =",d1_string,row)
row = gui.create_double_entry(mainframe,u"film 1: \u03B5' =",epsilon_f1_real_string,u"\u03B5'' =",epsilon_f1_imag_string,row)
row = gui.create_entry(mainframe,u"film 2 thickness: d/\u03BB =",d2_string,row)
row = gui.create_double_entry(mainframe,u"film 2: \u03B5' =",epsilon_f2_real_string,u"\u03B5'' =",epsilon_f2_imag_string,row)
row = gui.create_entry(mainframe,u"Number of periods =",N_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_entry(mainframe,u"cladding: \u03B5' =",epsilon_c_real_string,u"\u03B5'' =",epsilon_c_imag_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'angles of incidence:',row)
row = gui.create_double_entry_with_latex(mainframe,r'$\varphi_\mathrm{i}/\pi>$',phi_min_string,r'$\varphi_\mathrm{i}/\pi<$',phi_max_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_title(mainframe,"field parameters",row)
row = gui.create_radiobutton(mainframe,['polarization:','TE','TM'],polarization_string,2,row)
row = gui.create_entry_with_latex(mainframe,r'angle of incidence: $\varphi_\mathrm{i}/\pi=$',phi_0_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)