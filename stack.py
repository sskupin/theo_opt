import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import strat_stuff as strat

gui.set_rcParams()
title = "Reflection and Transmission at a Stack"
root = Tk.Tk()
root.title(title)

def reflection_transmission(epsilon_s,d1,epsilon_f1,d2,epsilon_f2,N,epsilon_c,phi): # computing coefficients of reflection and transmission
    kx,ksz,kcz = strat.KSC(epsilon_s,epsilon_c,phi)
    MTE,MTM = strat.MP(kx,d1,epsilon_f1,d2,epsilon_f2,N)
    RTE,RTM,TTE,TTM,tauTE,tauTM = strat.RTAU(ksz,kcz,epsilon_s,epsilon_c,MTE,MTM)   
    return RTE,RTM,tauTE,tauTM

def plot_subplot(ax,phi,curves,labels,colors,phi_min, phi_max):
    if np.floor(8*phi_max/np.pi)-np.ceil(8*phi_min/np.pi) >= 1:
        for index in range(len(labels)):
            ax.plot(phi,curves[index],colors[index],label=labels[index])
        ax.set_xticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2])
        ax.set_xticklabels([r'$0$', r'$\pi/8$', r'$\pi/4$', r'$3\pi/8$', r'$\pi/2$'])
        ax.set_xlabel(r'$\varphi_{\rm i}$')
        ax.set_xlim([phi_min, phi_max])
    else:
        for index in range(len(labels)):
            ax.plot(phi/np.pi,curves[index],colors[index],label=labels[index])
        ax.set_xlabel(r'$\varphi_{\rm i}/\pi$')
        ax.set_xlim([phi_min/np.pi, phi_max/np.pi])       
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

def Nz(d,epsilon,epsilon_s,phi_0):
    return np.rint(np.maximum(101,d*np.real(50*np.sqrt(epsilon-epsilon_s*np.sin(phi_0)**2)))).astype(np.int64)

def initialize():
    var_string[0].set("1") # epsilon_s
    var_string[1].set("1.129") # d1 in units of lambda
    var_string[2].set("2.12") # Re epsilon_f1
    var_string[3].set("0") # Im epsilon_f1
    var_string[4].set("0") # d2 in units of lambda
    var_string[5].set("2.2") # Re epsilon_f2
    var_string[6].set("0") # Im epsilon_f2
    var_string[7].set("1") # Number of periods
    var_string[8].set("15.1")
    var_string[9].set("0")
    var_string[10].set("0")
    var_string[11].set(".5")
    var_string[12].set("TE")
    var_string[13].set("0.395")
    gui.copy_stringvar_vector(var_string,var_save)
    calculate()

def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()  
        
def calculate():
    gui.change_cursor(root,"trek")
    try:
        epsilon_s = float(var_string[0].get())
        d1 = float(var_string[1].get())
        epsilon_f1_real = float(var_string[2].get())
        epsilon_f1_imag = float(var_string[3].get())
        d2 = float(var_string[4].get())
        epsilon_f2_real = float(var_string[5].get())
        epsilon_f2_imag = float(var_string[6].get())        
        N = int(var_string[7].get())
        epsilon_c_real = float(var_string[8].get())
        epsilon_c_imag = float(var_string[9].get())
        phi_min = float(var_string[10].get())*np.pi
        phi_max = float(var_string[11].get())*np.pi
        polarization = var_string[12].get()
        phi_0 = float(var_string[13].get())*np.pi
        
        if epsilon_c_imag < 0 or epsilon_c_real == 0 or epsilon_f1_real == 0 or epsilon_f2_real == 0 or epsilon_s <= 0\
                              or d1 < 0 or d2 < 0 or N < 0 or N > 50 or phi_max > np.pi/2 or phi_min < 0 or phi_min >= phi_max:
            gui.input_error("Values out of range. Re-initializing ...", reinitialize)
        else:
            f.clf()
            phi = np.linspace(phi_min, phi_max, num=10001, endpoint=False) # angle of incidence
            if phi_0 < phi_min or phi_0 > phi_max:
                phi_0 = phi_min
                var_string[13].set(phi_0/np.pi)
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
            zs = np.linspace(-np.maximum(0.75*N*(d1+d2)/2,1), 0, num=Nz(np.maximum(0.75*N*(d1+d2)/2,1),epsilon_s,epsilon_s,phi_0), endpoint=True)
            zf1 = np.linspace(0, d1, num=Nz(d1,epsilon_f1,epsilon_s,phi_0), endpoint=True)
            zf2 = np.linspace(0, d2, num=Nz(d2,epsilon_f2,epsilon_s,phi_0), endpoint=True) 
            zc = np.linspace(0, np.maximum(0.75*N*(d1+d2)/2,1), num=Nz(np.maximum(0.75*N*(d1+d2)/2,1),epsilon_c,epsilon_s,phi_0), endpoint=True)
            if d1 != 0 or d2 != 0:
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
                G = 1j*2*np.pi*np.sqrt(epsilon_s)*np.cos(phi_0)*(1-RTE)
                MTE = strat.mTE(2*np.pi*np.sqrt(epsilon_s)*np.cos(phi_0),zs)
                a3.plot(zs,np.abs(MTE[0,0]*F+MTE[0,1]*G),'b')
                if d1 != 0 or d2 != 0:
                    for index in range(N):
                        MTE = strat.mTE(2*np.pi*np.sqrt(epsilon_f1-epsilon_s*np.sin(phi_0)**2),zf1)
                        F,G = plot_amplitude(a3,MTE,F,G,index*(d1+d2)+zf1)
                        MTE = strat.mTE(2*np.pi*np.sqrt(epsilon_f2-epsilon_s*np.sin(phi_0)**2),zf2)
                        F,G = plot_amplitude(a3,MTE,F,G,(index+d1/(d1+d2))*(d1+d2)+zf2)
                a3.plot(N*(d1+d2)+zc,np.abs(F)*np.exp(-np.imag(np.sqrt(epsilon_c-epsilon_s*np.sin(phi_0)**2))*2*np.pi*zc),'b')
                a3.set_ylabel(r'$| E_\mathrm{TE} | / | E_\mathrm{TE}^\mathrm{i} |$ for $\varphi_\mathrm{i}=$ '+str(round(phi_0/np.pi,4))+r'$\pi$')
            else:
                rho = np.abs(RTM)**2
                tau = tauTM
                theta = np.angle(RTM)
                F = 1-RTM
                G = 1j*2*np.pi*np.sqrt(epsilon_s)*np.cos(phi_0)*(1+RTM)/epsilon_s
                MTM = strat.mTM(2*np.pi*np.sqrt(epsilon_s)*np.cos(phi_0),epsilon_s,zs)
                a3.plot(zs,np.abs(MTM[0,0]*F+MTM[0,1]*G),'b')
                if d1 != 0 or d2 != 0:
                    for index in range(N):
                        MTM = strat.mTM(2*np.pi*np.sqrt(epsilon_f1-epsilon_s*np.sin(phi_0)**2),epsilon_f1,zf1)
                        F,G = plot_amplitude(a3,MTM,F,G,index*(d1+d2)+zf1)
                        MTM = strat.mTM(2*np.pi*np.sqrt(epsilon_f2-epsilon_s*np.sin(phi_0)**2),epsilon_f2,zf2)
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
            
            plt.savefig('stack.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[8,5])
canvas = gui.create_canvas(root,f)
canvas.draw() # for faster feedback to user on startup
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(14)
var_save = gui.create_stringvar_vector(14)

initialize()

row = 1
row = gui.create_entry(mainframe,u"substrate: \u03B5 =",var_string[0],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"film 1 thickness: d/\u03BB =",var_string[1],row)
row = gui.create_double_entry(mainframe,u"film 1: \u03B5' =",var_string[2],u"\u03B5'' =",var_string[3],row)
row = gui.create_entry(mainframe,u"film 2 thickness: d/\u03BB =",var_string[4],row)
row = gui.create_double_entry(mainframe,u"film 2: \u03B5' =",var_string[5],u"\u03B5'' =",var_string[6],row)
row = gui.create_entry(mainframe,u"Number of periods =",var_string[7],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_entry(mainframe,u"cladding: \u03B5' =",var_string[8],u"\u03B5'' =",var_string[9],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'angles of incidence:',row)
row = gui.create_double_entry_with_latex(mainframe,r'$\varphi_\mathrm{i}/\pi>$',var_string[10],r'$\varphi_\mathrm{i}/\pi<$',var_string[11],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_title(mainframe,"field parameters",row)
row = gui.create_radiobutton(mainframe,['polarization:','TE','TM'],var_string[12],2,row)
row = gui.create_entry_with_latex(mainframe,r'angle of incidence: $\varphi_\mathrm{i}/\pi=$',var_string[13],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)