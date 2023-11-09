import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
title = "Transmission of a Fabry-Perot Resonator (monochromatic)"
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
    kf1z = np.sqrt(epsilon_f1-kx**2)*2*np.pi
    kf2z = np.sqrt(epsilon_f2-kx**2)*2*np.pi
    m1TE = mTE(kf1z,d1)
    m2TE = mTE(kf2z,d2)
    m1TM = mTM(kf1z,epsilon_f1,d1)
    m2TM = mTM(kf2z,epsilon_f2,d2)
    MTEP = np.linalg.matrix_power(np.matmul(m2TE,m1TE),N)
    MTMP = np.linalg.matrix_power(np.matmul(m2TM,m1TM),N)
    return MTEP,MTMP

def reflection_transmission(epsilon_s,da,epsilon_fa,db,epsilon_fb,N,d,epsilon_f,epsilon_c,phi): # computing coefficients of reflection and transmission
    kx = np.sqrt(epsilon_s)*np.sin(phi)
    ksz = np.sqrt(epsilon_s-kx**2)*2*np.pi
    kfz = np.sqrt(epsilon_f-kx**2)*2*np.pi
    kcz = np.sqrt(epsilon_c+0j-kx**2)*2*np.pi
    MTEP,MTMP = MP(kx,da[0],epsilon_fa[0],db[0],epsilon_fb[0],N[0])
    MTE = np.matmul(mTE(kfz,d),MTEP)
    MTM = np.matmul(mTM(kfz,epsilon_f,d),MTMP)
    MTEP,MTMP = MP(kx,da[1],epsilon_fa[1],db[1],epsilon_fb[1],N[1])
    MTE = np.matmul(MTEP,MTE)
    MTM = np.matmul(MTMP,MTM)
    
    NTE = ksz*MTE[1,1]+kcz*MTE[0,0]+1j*MTE[1,0]-1j*ksz*kcz*MTE[0,1]
    RTE = (ksz*MTE[1,1]-kcz*MTE[0,0]-1j*MTE[1,0]-1j*ksz*kcz*MTE[0,1])/NTE
    NTM = ksz*MTM[1,1]/epsilon_s+kcz*MTM[0,0]/epsilon_c+1j*MTM[1,0]-1j*ksz*kcz*MTM[0,1]/(epsilon_s*epsilon_c)
    RTM = -(ksz*MTM[1,1]/epsilon_s-kcz*MTM[0,0]/epsilon_c-1j*MTM[1,0]-1j*ksz*kcz*MTM[0,1]/(epsilon_s*epsilon_c))/NTM # for electric field (negative of magnetic coeff.)
    tauTE = np.real(kcz)/ksz*np.abs(2*ksz/NTE)**2
    tauTM = np.real(kcz/epsilon_c)*epsilon_s/ksz*(np.abs(2*ksz/epsilon_s/NTM))**2
    
    return RTE,RTM,tauTE,tauTM

def initialize():
    da_string[0].set("0.147")
    db_string[0].set("0.13")
    da_string[1].set("0.13")
    db_string[1].set("0.147")
    d_string.set("3")
    epsilon_s_string.set("1")
    epsilon_fa_real_string[0].set("2.13")
    epsilon_fa_imag_string[0].set("0")
    epsilon_fb_real_string[0].set("6.68")
    epsilon_fb_imag_string[0].set("0")
    epsilon_fa_real_string[1].set("6.68")
    epsilon_fa_imag_string[1].set("0")
    epsilon_fb_real_string[1].set("2.13")
    epsilon_fb_imag_string[1].set("0")
    epsilon_f_real_string.set("1")
    epsilon_f_imag_string.set("0")
    N_string[0].set("7")
    N_string[1].set("7")
    epsilon_c_string.set("1")
    sinphi_min_string.set("0.7")
    sinphi_max_string.set("0.77")
        
    calculate()

def calculate():
    try:
        epsilon_s = float(epsilon_s_string.get())
        da = np.zeros(len(N_string))
        db = np.zeros(len(N_string))
        epsilon_fa_real = np.zeros(len(N_string))
        epsilon_fa_imag = np.zeros(len(N_string))
        epsilon_fb_real = np.zeros(len(N_string))
        epsilon_fb_imag = np.zeros(len(N_string))
        N = np.zeros(len(N_string),dtype=int)
        for index in range(len(N_string)):
            da[index] = float(da_string[index].get())
            epsilon_fa_real[index] = float(epsilon_fa_real_string[index].get())
            epsilon_fa_imag[index] = float(epsilon_fa_imag_string[index].get())
            db[index] = float(db_string[index].get())
            epsilon_fb_real[index] = float(epsilon_fb_real_string[index].get())     
            epsilon_fb_imag[index] = float(epsilon_fb_imag_string[index].get())     
            N[index] = int(N_string[index].get())
        d = float(d_string.get())
        epsilon_f_real = float(epsilon_f_real_string.get())
        epsilon_f_imag = float(epsilon_f_imag_string.get())
        epsilon_c = float(epsilon_c_string.get())
        sinphi_min = float(sinphi_min_string.get())
        sinphi_max = float(sinphi_max_string.get())

        if epsilon_c <= 0 or (epsilon_fa_imag < 0).any() or (epsilon_fb_imag < 0).any() or epsilon_s <= 0 or sinphi_min<0 or sinphi_max>1\
           or sinphi_min>=sinphi_max or (da <= 0).any() or (db < 0).any() or (N < 0).any() or N[0] == 0 or (d<0) or (epsilon_f_real <= 0):
            gui.input_error(initialize)
        else:
            f.clf()
            phi = np.linspace(np.arcsin(sinphi_min), np.arcsin(sinphi_max), num=10001, endpoint=False) # angle of incidence
            epsilon_fa = epsilon_fa_real + 1j*epsilon_fa_imag
            epsilon_fb = epsilon_fb_real + 1j*epsilon_fb_imag
            epsilon_f = epsilon_f_real + 1j*epsilon_f_imag
            vreflection_transmission = np.vectorize(reflection_transmission, excluded=[1,2,3,4,5])
            RTE,RTM,tauTE,tauTM = vreflection_transmission(epsilon_s,da,epsilon_fa,db,epsilon_fb,N,d,epsilon_f,epsilon_c,phi)
            a1 = f.add_subplot(211)
            a1.plot(np.sqrt(epsilon_s)*np.sin(phi),tauTE,'r')
            a1.set_xlabel(r'$\sqrt{\varepsilon_{\rm s}}\sin\varphi_{\rm i}=k_x \lambda/(2\pi)$')
            a1.set_xlim([np.sqrt(epsilon_s)*sinphi_min,np.sqrt(epsilon_s)*sinphi_max])
            a1.set_ylim([-0.025,1.025])
            a1.set_ylabel(r'$\tau_{\rm TE}$')
            a2 = f.add_subplot(212)
            a2.plot(np.sqrt(epsilon_s)*np.sin(phi),tauTM,'r')
            a2.set_xlabel(r'$\sqrt{\varepsilon_{\rm s}}\sin\varphi_{\rm i}=k_x \lambda/(2\pi)$')
            a2.set_xlim([np.sqrt(epsilon_s)*sinphi_min,np.sqrt(epsilon_s)*sinphi_max])
            a2.set_ylim([-0.05,1.05])
            a2.set_ylabel(r'$\tau_{\rm TM}$')            
                     
            plt.tight_layout()  
            
#            plt.savefig('FP_kx.pdf',bbox_inches='tight',dpi=300, transparent=True)

            canvas.draw()
    except ValueError: gui.input_error(initialize)

f = plt.figure(1,[5,5])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

epsilon_s_string = Tk.StringVar()
da_string = gui.create_stringvar_vector(2)
epsilon_fa_real_string = gui.create_stringvar_vector(2)
epsilon_fa_imag_string = gui.create_stringvar_vector(2)
db_string = gui.create_stringvar_vector(2)
epsilon_fb_real_string = gui.create_stringvar_vector(2)
epsilon_fb_imag_string = gui.create_stringvar_vector(2)
N_string = gui.create_stringvar_vector(2)
d_string = Tk.StringVar()
epsilon_f_real_string = Tk.StringVar()
epsilon_f_imag_string = Tk.StringVar()
epsilon_c_string = Tk.StringVar()
sinphi_min_string = Tk.StringVar()
sinphi_max_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_entry(mainframe,u"substrate: \u03B5 =",epsilon_s_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"mirror 1: N =",N_string[0],row)
row = gui.create_triple_entry(mainframe,u"d\u2081/\u03BB =",da_string[0],u"\u03B5\u2081' =",epsilon_fa_real_string[0],u"\u03B5\u2081'' =",epsilon_fa_imag_string[0],row)
row = gui.create_triple_entry(mainframe,u"d\u2082/\u03BB =",db_string[0],u"\u03B5\u2082' =",epsilon_fb_real_string[0],u"\u03B5\u2082'' =",epsilon_fb_imag_string[0],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_triple_entry(mainframe,u"cavity: d/\u03BB =",d_string,u"\u03B5' =",epsilon_f_real_string,u"\u03B5'' =",epsilon_f_imag_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"mirror 2: N =",N_string[1],row)
row = gui.create_triple_entry(mainframe,u"d\u2081/\u03BB =",da_string[1],u"\u03B5\u2081' =",epsilon_fa_real_string[1],u"\u03B5\u2081'' =",epsilon_fa_imag_string[1],row)
row = gui.create_triple_entry(mainframe,u"d\u2082/\u03BB =",db_string[1],u"\u03B5\u2082' =",epsilon_fb_real_string[1],u"\u03B5\u2082'' =",epsilon_fb_imag_string[1],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"cladding: \u03B5 =",epsilon_c_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_entry_with_latex(mainframe,r'$\sin\varphi_\mathrm{i}>$',sinphi_min_string,r'$\sin\varphi_\mathrm{i}<$',sinphi_max_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)