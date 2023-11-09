import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
title = "Reflection and Transmission of a Stack (Bloch)"
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

def MP(kx,d1,epsilon_f1,d2,epsilon_f2,N):
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
    kappafTE = -1j*(np.exp(1j*kzTE*(d1+d2))-MTE[0,0])/MTE[0,1]
    if np.real(kappafTE) < 0:
        kzTE = -kzTE
        kappafTE = -1j*(np.exp(1j*kzTE*(d1+d2))-MTE[0,0])/MTE[0,1]
    kappafTM = -1j*(np.exp(1j*kzTM*(d1+d2))-MTM[0,0])/MTM[0,1]
    if np.real(kappafTM) < 0:
        kzTM = -kzTM
        kappafTM = -1j*(np.exp(1j*kzTM*(d1+d2))-MTM[0,0])/MTM[0,1]
    kappabTE = 1j*(np.exp(1j*kzTE*(d1+d2))-MTE[1,1])/MTE[0,1]
    kappabTM = 1j*(np.exp(1j*kzTM*(d1+d2))-MTM[1,1])/MTM[0,1]
    if kappabTE == kappafTE:
        MTEP = np.array([[0,1],[0,0]])
    else:
        MTEP = 1/(kappabTE-kappafTE)*np.array([[kappabTE*np.exp(1j*kzTE*N*(d1+d2))-kappafTE*np.exp(-1j*kzTE*N*(d1+d2)),-2*np.sin(kzTE*N*(d1+d2))],
                                              [-2*kappafTE*kappabTE*np.sin(kzTE*N*(d1+d2)),-kappafTE*np.exp(1j*kzTE*N*(d1+d2))+kappabTE*np.exp(-1j*kzTE*N*(d1+d2))]])
    if kappabTM == kappafTM:
        MTMP = np.array([[0,1],[0,0]])
    else:
        MTMP = 1/(kappabTM-kappafTM)*np.array([[kappabTM*np.exp(1j*kzTM*N*(d1+d2))-kappafTM*np.exp(-1j*kzTM*N*(d1+d2)),-2*np.sin(kzTM*N*(d1+d2))],
                                              [-2*kappafTM*kappabTM*np.sin(kzTM*N*(d1+d2)),-kappafTM*np.exp(1j*kzTM*N*(d1+d2))+kappabTM*np.exp(-1j*kzTM*N*(d1+d2))]])        
    return MTEP,MTMP

def reflection_transmission(epsilon_s,da,epsilon_fa,db,epsilon_fb,N,epsilon_c,phi): # computing coefficients of reflection and transmission
    kx = np.sqrt(epsilon_s)*np.sin(phi)
    ksz = np.sqrt(epsilon_s-kx**2)*2*np.pi
    MTE = np.array([[1,0],[0,1]])
    MTM = np.array([[1,0],[0,1]])
    for index in range(len(N)):
        MTEP,MTMP = MP(kx,da[index],epsilon_fa[index],db[index],epsilon_fb[index],N[index])
        MTE = np.matmul(MTEP,MTE)
        MTM = np.matmul(MTMP,MTM)
    kcz = np.sqrt(epsilon_c+0j-kx**2)*2*np.pi
    NTE = ksz*MTE[1,1]+kcz*MTE[0,0]+1j*MTE[1,0]-1j*ksz*kcz*MTE[0,1]
    RTE = (ksz*MTE[1,1]-kcz*MTE[0,0]-1j*MTE[1,0]-1j*ksz*kcz*MTE[0,1])/NTE
    NTM = ksz*MTM[1,1]/epsilon_s+kcz*MTM[0,0]/epsilon_c+1j*MTM[1,0]-1j*ksz*kcz*MTM[0,1]/(epsilon_s*epsilon_c)
    RTM = -(ksz*MTM[1,1]/epsilon_s-kcz*MTM[0,0]/epsilon_c-1j*MTM[1,0]-1j*ksz*kcz*MTM[0,1]/(epsilon_s*epsilon_c))/NTM # for electric field (negative of magnetic coeff.)
    tauTE = np.real(kcz)/ksz*np.abs(2*ksz/NTE)**2
    tauTM = np.real(kcz/epsilon_c)*epsilon_s/ksz*(np.abs(2*ksz/epsilon_s/NTM))**2
    
    return RTE,RTM,tauTE,tauTM

def plot_subplot(ax,phi,curve_1,curve_2,label_1,label_2,phi_min, phi_max):
    ax.plot(phi,curve_1,'b',label=label_1)
    ax.plot(phi,curve_2,'r',label=label_2)
    if np.floor(8*phi_max/np.pi)-np.ceil(8*phi_min/np.pi) >= 1:
        ax.set_xticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2])
        ax.set_xticklabels([r'$0$', r'$\pi/8$', r'$\pi/4$', r'$3\pi/8$', r'$\pi/2$'])
    ax.set_xlabel(r'$\varphi_{\rm i}$')
    ax.set_xlim([phi_min, phi_max])
    ax.set_ylabel(label_1+','+label_2)
    ax.legend()
        
def initialize():
    da_string[0].set("1")
    db_string[0].set("1")
    da_string[1].set("1")
    db_string[1].set("1")
    da_string[2].set("1")
    db_string[2].set("1")
    epsilon_s_string.set("1")
    epsilon_fa_string[0].set("2.25")
    epsilon_fb_string[0].set("2.2")
    epsilon_fa_string[1].set("2.35")
    epsilon_fb_string[1].set("2.3")
    epsilon_fa_string[2].set("2.45")
    epsilon_fb_string[2].set("2.4")
    N_string[0].set("10")
    N_string[1].set("10")
    N_string[2].set("10")
    epsilon_c_string.set("1")
    phi_min_string.set("0")
    phi_max_string.set(".5")
       
    calculate()

def calculate():
    try:
        epsilon_s = float(epsilon_s_string.get())
        da = np.zeros(len(N_string))
        db = np.zeros(len(N_string))
        epsilon_fa = np.zeros(len(N_string))
        epsilon_fb = np.zeros(len(N_string))
        N = np.zeros(len(N_string),dtype=int)
        for index in range(len(N_string)):
            da[index] = float(da_string[index].get())
            epsilon_fa[index] = float(epsilon_fa_string[index].get())
            db[index] = float(db_string[index].get())
            epsilon_fb[index] = float(epsilon_fb_string[index].get())        
            N[index] = int(N_string[index].get())
        epsilon_c = float(epsilon_c_string.get())
        phi_min = float(phi_min_string.get())*np.pi
        phi_max = float(phi_max_string.get())*np.pi
        
        if epsilon_c <= 0 or (epsilon_fa <= 0).any() or (epsilon_fb <= 0).any() or epsilon_s <= 0\
            or (da <= 0).any() or (db < 0).any() or (N < 0).any() or N[0] == 0 or (N > 50).any() or phi_max > np.pi/2 or phi_min < 0 or phi_min >= phi_max:
            gui.input_error(initialize)
        else:
            f.clf()
            phi = np.linspace(phi_min, phi_max, num=10001, endpoint=False) # angle of incidence
            vreflection_transmission = np.vectorize(reflection_transmission, excluded=[1,2,3,4,5])
            RTE,RTM,tauTE,tauTM = vreflection_transmission(epsilon_s,da,epsilon_fa,db,epsilon_fb,N,epsilon_c,phi)
            a1 = f.add_subplot(131)
            plot_subplot(a1,phi,np.abs(RTE)**2,tauTE,r'$\rho_{\rm TE}$',r'$\tau_{\rm TE}$',phi_min, phi_max)
            a1.set_ylim([-0.025,1.025])
            a2 = f.add_subplot(132)
            plot_subplot(a2,phi,np.abs(RTM)**2,tauTM,r'$\rho_{\rm TM}$',r'$\tau_{\rm TM}$',phi_min, phi_max)
            a2.set_ylim([-0.025,1.025])
            a3 = f.add_subplot(133)
            plot_subplot(a3,phi,np.angle(RTE),np.angle(RTM),r'$\theta_{\rm TE}$',r'$\theta_{\rm TM}$',phi_min, phi_max)  
            a3.set_yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
            a3.set_yticklabels([r'$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])
            a3.set_ylim([-1.05*np.pi,1.05*np.pi])
                     
            plt.tight_layout()  
            
#            plt.savefig('Bloch_stack.pdf',bbox_inches='tight',dpi=300, transparent=True)

            canvas.draw()
    except ValueError: gui.input_error(initialize)

f = plt.figure(1,[7,2])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

epsilon_s_string = Tk.StringVar()
da_string = gui.create_stringvar_vector(3)
epsilon_fa_string = gui.create_stringvar_vector(3)
db_string = gui.create_stringvar_vector(3)
epsilon_fb_string = gui.create_stringvar_vector(3)
N_string = gui.create_stringvar_vector(3)
epsilon_c_string = Tk.StringVar()
phi_min_string = Tk.StringVar()
phi_max_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_entry(mainframe,u"substrate: \u03B5 =",epsilon_s_string,row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"medium 1: N =",N_string[0],row)
row = gui.create_double_entry(mainframe,u"d\u2081/\u03BB =",da_string[0],u"\u03B5\u2081 =",epsilon_fa_string[0],row)
row = gui.create_double_entry(mainframe,u"d\u2082/\u03BB =",db_string[0],u"\u03B5\u2082 =",epsilon_fb_string[0],row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"medium 2: N =",N_string[1],row)
row = gui.create_double_entry(mainframe,u"d\u2081/\u03BB =",da_string[1],u"\u03B5\u2081 =",epsilon_fa_string[1],row)
row = gui.create_double_entry(mainframe,u"d\u2082/\u03BB =",db_string[1],u"\u03B5\u2082 =",epsilon_fb_string[1],row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"medium 3: N =",N_string[2],row)
row = gui.create_double_entry(mainframe,u"d\u2081/\u03BB =",da_string[2],u"\u03B5\u2081 =",epsilon_fa_string[2],row)
row = gui.create_double_entry(mainframe,u"d\u2082/\u03BB =",db_string[2],u"\u03B5\u2082 =",epsilon_fb_string[2],row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"cladding: \u03B5 =",epsilon_c_string,row)
row = gui.create_double_entry_with_latex(mainframe,r'$\varphi_\mathrm{i}/\pi>$',phi_min_string,r'$\varphi_\mathrm{i}/\pi<$',phi_max_string,row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)