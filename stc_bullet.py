import numpy as np
import tkinter as Tk
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
import gui_stuff as gui

gui.set_rcParams()
title = "STC Bullet"
root = Tk.Tk()
root.title(title)

def delta(RbarIinv,Z): # focusing strength
    return 1/np.sqrt((1+Z*RbarIinv)**2+Z**2)
vdelta = np.vectorize(delta)

def Rbarinv(RbarIinv,Z): # normalized inverse radius spatial phase
    if RbarIinv == 0:
        Rbarinv = Z/(Z**2+1)
    else:
        Rbarinv = (Z+RbarIinv**2*(Z+1/RbarIinv))/(Z**2+RbarIinv**2*(Z+1/RbarIinv)**2)
    return Rbarinv

def b(bI,OmegabarI,RbarIinv,Z): # normalized rotation parameter
    return bI-(OmegabarI-bI*RbarIinv)*Z
vb = np.vectorize(b)

def Omegabar(bI,OmegabarI,RbarIinv,Z): # normalized linear chirp parameter
    return OmegabarI-bI*RbarIinv+b(bI,OmegabarI,RbarIinv,Z)*Rbarinv(RbarIinv,Z)
vOmegabar = np.vectorize(Omegabar)

def deltat(bI,OmegabarI,RbarIinv,Z): # temporal focusing strength
    return np.sqrt((1+delta(RbarIinv,Z)**2*b(bI,OmegabarI,RbarIinv,Z)**2)/(1+bI**2))
vdeltat = np.vectorize(deltat)

def Xu(bI,OmegabarI,RbarIinv,Z,T): # normalized spatial coordinate for temporal ellipse (upper)
    return np.sqrt(1-T**2)/delta(RbarIinv,Z)+b(bI,OmegabarI,RbarIinv,Z)*T
vXu = np.vectorize(Xu)

def Xl(bI,OmegabarI,RbarIinv,Z,T): # normalized spatial coordinate for temporal ellipse (lower)
    return -np.sqrt(1-T**2)/delta(RbarIinv,Z)+b(bI,OmegabarI,RbarIinv,Z)*T
vXl = np.vectorize(Xl)

def pq(bI,OmegabarI,RbarIinv,freqnorm,Z): # compute parameters for frequency space ellipse   
    b = bI-(OmegabarI-bI*RbarIinv)*Z
    Chirp = (OmegabarI-bI*RbarIinv)**2*Z+bI**2*RbarIinv-b**2*Rbarinv(RbarIinv,Z)
    Omegabar = OmegabarI-bI*RbarIinv+b*Rbarinv(RbarIinv,Z)
    delta = 1/np.sqrt((1+Z*RbarIinv)**2+Z**2)
    denom = 4*Omegabar**2*b**2*delta**2+4*b**2*delta**4+8*Chirp*Omegabar*b*delta**2+4*Chirp**2*delta**2+4*Omegabar**2+4*delta**2
    p = (-4*Omegabar*freqnorm*b**2*delta**2-4*Chirp*freqnorm*b*delta**2-4*Omegabar*freqnorm)/denom
    q = (-4*b**4*delta**4+freqnorm**2*b**2*delta**2-8*b**2*delta**2-4*Chirp**2+freqnorm**2-4)/denom
    return p,q

def XFu(bI,OmegabarI,RbarIinv,freqnorm,Z): # normalized spatial coordinate for frequency space ellipse (upper)
    p,q = pq(bI,OmegabarI,RbarIinv,freqnorm,Z)
    return -p/2+np.sqrt(np.abs(p**2/4-q))
vXFu = np.vectorize(XFu)

def XFl(bI,OmegabarI,RbarIinv,freqnorm,Z): # normalized spatial coordinate for frequency space ellipse (lower)
    p,q = pq(bI,OmegabarI,RbarIinv,freqnorm,Z)
    return -p/2-np.sqrt(np.abs(p**2/4-q))
vXFl = np.vectorize(XFl)

def freqnormmax(OmegabarI,bI):
    return 2*np.sqrt(OmegabarI**2*bI**2+OmegabarI**2+bI**2+1)   

def initialize():
    wI_string.set("0.2") # width at input in mm
    lambda0_string.set("0.8") # wavelength in microns
    aI_string.set("0.") # rotation parameter at input in mm/fs
    OmegaI_string.set("0.025") # spatial chirp parameter at input in 1/ps
    RIinv_string.set("-0.25") # inverse radius spatial phase at input in 1/mm
    Tp_string.set("250") # pulse duration at input without rotation in fs
    normalization_string.set("physical")
    
    calculate()

def getparameters():
    wI = float(wI_string.get())*1.e-3
    Tp = int(Tp_string.get())*1.e-15 
    normalization = normalization_string.get()
    if normalization == 'physical':
        lambda0 = float(lambda0_string.get())*1.e-6
        aI = float(aI_string.get())*1.e12
        OmegaI = float(OmegaI_string.get())*1.e12
        RIinv = float(RIinv_string.get())*1.e3
        k0wI = 2*np.pi/lambda0*wI
        k0wI_string.set(np.around(k0wI,decimals=8))
        bI = aI*Tp/wI
        bI_string.set(np.around(bI,decimals=8))
        OmegabarI = OmegaI*Tp*wI*np.pi/lambda0
        OmegabarI_string.set(np.around(OmegabarI,decimals=8))
        RbarIinv = RIinv*wI**2*np.pi/lambda0
        RbarIinv_string.set(np.around(RbarIinv,decimals=8))
    elif normalization == 'normalized':
        k0wI = float(k0wI_string.get())
        bI = float(bI_string.get())
        OmegabarI = float(OmegabarI_string.get())
        RbarIinv = float(RbarIinv_string.get())
        lambda0 = 2*np.pi/k0wI*wI
        lambda0_string.set(np.around(lambda0*1.e6,decimals=8))
        aI = bI*wI/Tp
        aI_string.set(np.around(aI*1.e-12,decimals=8))
        OmegaI = OmegabarI*lambda0/(Tp*wI*np.pi)
        OmegaI_string.set(np.around(OmegaI*1.e-12,decimals=8))
        RIinv = RbarIinv*lambda0/(wI**2*np.pi)
        RIinv_string.set(np.around(RIinv*1.e-3,decimals=8))
    zIF_string.set(np.around(np.pi*wI**2/lambda0*1.e3,decimals=8))
    TI_string.set(np.around(Tp/np.sqrt(1+bI**2)*1.e15,decimals=8))
    if RbarIinv == 0:
        ZD = 0
    else:
        ZD = -1/RbarIinv/(1/RbarIinv**2+1) # normalized distance to focus
    w0 = wI/delta(RbarIinv,ZD) # width at focus
    w0_string.set(np.around(w0*1.e3,decimals=8))         
    z0F_string.set(np.around(np.pi*w0**2/lambda0*1.e3,decimals=8))    
    Z = np.linspace(0, 2*np.maximum(np.abs(ZD),2*w0**2/wI**2), num=500, endpoint=True)
    T = np.linspace(-1, 1, num=500, endpoint=True)
    freqnorm = np.linspace(-freqnormmax(OmegabarI,bI), freqnormmax(OmegabarI,bI), num=500, endpoint=True)
    return wI,lambda0,Tp,bI,OmegabarI,RbarIinv,Z,T,freqnorm

def plotellipse(color,label,ax5,ax6,bI,OmegabarI,RbarIinv,freqnorm,Z,T):
    el1, = ax5.plot(vXu(bI,OmegabarI,RbarIinv,Z,T),T,color,label=label)
    el2, = ax5.plot(vXl(bI,OmegabarI,RbarIinv,Z,T),T,color) 
    el3, = ax6.plot(vXFu(bI,OmegabarI,RbarIinv,freqnorm,Z),freqnorm,color,label=label)
    el4, = ax6.plot(vXFl(bI,OmegabarI,RbarIinv,freqnorm,Z),freqnorm,color)
    return el1,el2,el3,el4

def Zaxes(ax,Z):
    axbis = ax.twiny()        
    axbis.set_xlim(ax.get_xlim())
    axbis.set_xticks([Z[0], Z[-1]/2, Z[-1]])
    axbis.set_xticklabels(['0',r'$Z_{\rm f}$',r'$2Z_{\rm f}$'])
    ax.set_xlabel(r'$z/z_{\rm IF}$')

def calculate():
    try:
        wI,lambda0,Tp,bI,OmegabarI,RbarIinv,Z,T,freqnorm = getparameters()        
        if wI <= 0 or lambda0 <= 0 or Tp <= 0:
            gui.input_error(initialize)
        else:      
            f.clf()
            ax1 = f.add_subplot(321)
            ax2 = f.add_subplot(322)
            ax3 = f.add_subplot(323)
            ax4 = f.add_subplot(324)
            ax5 = f.add_subplot(325)
            ax6 = f.add_subplot(326)
            cu1, = ax1.plot(Z,vdelta(RbarIinv,Z),'b')
            Zaxes(ax1,Z)
            ax1.set_ylabel(r'$\delta=w_{\rm I}/w$')
            cu2, = ax2.plot(Z,vdeltat(bI,OmegabarI,RbarIinv,Z),'b')
            Zaxes(ax2,Z)
            ax2.set_ylabel(r'$T_{\rm I}/T_{\rm B}$') 
            cu3, = ax3.plot(Z,vb(bI,OmegabarI,RbarIinv,Z),'b')
            Zaxes(ax3,Z)
            ax3.set_ylabel(r'$b=a T_{\rm p}/w_{\rm I}$')
            cu4, = ax4.plot(Z,vOmegabar(bI,OmegabarI,RbarIinv,Z),'b')
            Zaxes(ax4,Z)
            ax4.set_ylabel(r'$\overline{\Omega} = \Omega T_{\rm p} z_{\rm IF}/w_{\rm I}$')
            plotellipse('w','',ax5,ax6,bI,OmegabarI,RbarIinv,freqnorm,Z[-1]/2,T)
            plotellipse('w','',ax5,ax6,bI,OmegabarI,RbarIinv,freqnorm,Z[-1],T)
            el1,el2,el3,el4, = plotellipse('b',r'$z/z_{\rm IF}=2Z_{\rm f}$',ax5,ax6,bI,OmegabarI,RbarIinv,freqnorm,Z[0],T)
            ax5.set_ylabel(r'$t/T_{\rm p}$')
            ax5.set_xlabel(r'$x/w_{\rm I}$')
            ax6.set_ylabel(r'$\overline{\omega}T_{\rm p}$')
            ax6.set_xlabel(r'$x/w_{\rm I}$')  
            plt.tight_layout()
            canvas.draw()
            for index in range(int(Z.size/20)-1,Z.size,int(Z.size/20)):
                cu1.set_data(Z[0:index],vdelta(RbarIinv,Z[0:index]))
                cu2.set_data(Z[0:index],vdeltat(bI,OmegabarI,RbarIinv,Z[0:index]))
                cu3.set_data(Z[0:index],vb(bI,OmegabarI,RbarIinv,Z[0:index]))
                cu4.set_data(Z[0:index],vOmegabar(bI,OmegabarI,RbarIinv,Z[0:index]))
                el1.set_xdata(vXu(bI,OmegabarI,RbarIinv,Z[index],T))
                el2.set_xdata(vXl(bI,OmegabarI,RbarIinv,Z[index],T))
                el3.set_xdata(vXFu(bI,OmegabarI,RbarIinv,freqnorm,Z[index]))
                el4.set_xdata(vXFl(bI,OmegabarI,RbarIinv,freqnorm,Z[index]))                
                canvas.draw()
                time.sleep(.1)
            canvas.draw()
            plotellipse('r',r'$z/z_{\rm IF}=Z_{\rm f}$',ax5,ax6,bI,OmegabarI,RbarIinv,freqnorm,Z[-1]/2,T)
            plotellipse('k',r'$z/z_{\rm IF}=0$',ax5,ax6,bI,OmegabarI,RbarIinv,freqnorm,Z[0],T)
            ax5.legend()
            ax6.legend()
            plt.tight_layout()
            
            canvas.draw()
    except ValueError: gui.input_error(initialize)

f = plt.figure(1,[8,10])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

wI_string = Tk.StringVar()
lambda0_string = Tk.StringVar()
aI_string = Tk.StringVar()
OmegaI_string = Tk.StringVar()
RIinv_string = Tk.StringVar()
Tp_string = Tk.StringVar()
w0_string = Tk.StringVar()
z0F_string = Tk.StringVar()
normalization_string = Tk.StringVar()
k0wI_string = Tk.StringVar()
bI_string = Tk.StringVar()
OmegabarI_string = Tk.StringVar()
RbarIinv_string = Tk.StringVar()
zIF_string = Tk.StringVar()
TI_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"$w_{\rm I}$ [mm] $=$",wI_string,row)
row = gui.create_entry_with_latex(mainframe,r"$T_{\rm P}$ [fs] $=$",Tp_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_checkbutton(mainframe,"physical parameters",'normalized','physical',normalization_string,row)
row = gui.create_entry_with_latex(mainframe,r"$\lambda_0$ [$\mu$m] $=$",lambda0_string,row)
row = gui.create_entry_with_latex(mainframe,r"$a_{\rm I}$ [mm/fs] $=$",aI_string,row)
row = gui.create_entry_with_latex(mainframe,r"$\Omega_{\rm I}$ [1/ps] $=$",OmegaI_string,row)
row = gui.create_entry_with_latex(mainframe,r"$1/R_{\rm I}$ [1/mm] $=$",RIinv_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_checkbutton(mainframe,"normalized parameters",'physical','normalized',normalization_string,row)
row = gui.create_entry_with_latex(mainframe,r"$k_0 w_{\rm I} =$",k0wI_string,row)
row = gui.create_entry_with_latex(mainframe,r"$b_{\rm I} = a_{\rm I} T_{\rm P} / w_{\rm I}=$",bI_string,row)
row = gui.create_entry_with_latex(mainframe, r"$\overline{\Omega}_{\rm I} = \Omega_{\rm I} T_{\rm P} z_{\rm IF}/w_{\rm I}=$",OmegabarI_string,row)
row = gui.create_entry_with_latex(mainframe, r"$1/\overline{R}_{\rm I} = z_{\rm IF} / R_{\rm I}$",RbarIinv_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_label_with_latex(mainframe, r"$T_{\rm I}$ [fs] $=$",TI_string,row)
row = gui.create_label_with_latex(mainframe,r"$z_{\rm IF}$ [mm] $=$",zIF_string,row)
row = gui.create_label_with_latex(mainframe,r"$w_{0}$ [mm] $=$",w0_string,row)
row = gui.create_label_with_latex(mainframe,r"$z_{\rm 0F}$ [mm] $=$",z0F_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)