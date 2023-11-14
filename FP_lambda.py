import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import tkinter as Tk
import gui_stuff as gui
import media as media

gui.set_rcParams()
title = "Transmission of a Fabry-Perot Resonator (normal incidence)"
root = Tk.Tk()
root.title(title)

def m(kfz,z):
    return np.array([[np.cos(kfz*z),np.sin(kfz*z)/kfz],[-np.sin(kfz*z)*kfz,np.cos(kfz*z)]])

def MP(d1,epsilon_f1,d2,epsilon_f2,N,lambdav):
    kf1z = np.sqrt(epsilon_f1)*2*np.pi/lambdav
    kf2z = np.sqrt(epsilon_f2)*2*np.pi/lambdav
    m1 = m(kf1z,d1)
    m2 = m(kf2z,d2)

    return np.linalg.matrix_power(np.matmul(m2,m1),N)

def reflection_transmission(epsilon_s,da,epsilon_fa0,epsilon_fa1,db,epsilon_fb0,epsilon_fb1,N,d,epsilon_f,epsilon_c,lambdav): # computing coefficients of reflection and transmission
    ksz = np.sqrt(epsilon_s)*2*np.pi/lambdav
    kfz = np.sqrt(epsilon_f)*2*np.pi/lambdav
    kcz = np.sqrt(epsilon_c)*2*np.pi/lambdav
    MP1 = MP(da[0],epsilon_fa0,db[0],epsilon_fb0,N[0],lambdav)
    M = np.matmul(m(kfz,d),MP1)
    MP2 = MP(da[1],epsilon_fa1,db[1],epsilon_fb1,N[1],lambdav)
    M = np.matmul(MP2,M)
    DENOM = ksz*M[1,1]+kcz*M[0,0]+1j*M[1,0]-1j*ksz*kcz*M[0,1]
    R = (ksz*M[1,1]-kcz*M[0,0]-1j*M[1,0]-1j*ksz*kcz*M[0,1])/DENOM
    tau = np.real(kcz)/np.real(ksz)*np.abs(2*ksz/DENOM)**2
    
    return R,tau

def epsilon(medium,lambdav):
    if medium == "Vacuum":
        epsilon_medium = np.ones_like(lambdav)
    elif medium == "AlAs":
        epsilon_medium = media.AlAs(lambdav)
    elif medium == "BaSF":
        epsilon_medium = media.BaSF(lambdav)
    elif medium == "AlGaAs (31.5% Al)":
        epsilon_medium = media.AlGaAs31(lambdav)
    elif medium == "TiO2":
        epsilon_medium = media.TiO2(lambdav)
    elif medium == "fused silica":
        epsilon_medium = media.silica(lambdav)
    elif medium == "Ag":
        epsilon_medium = media.Ag(lambdav)
    else:
        print("Oops! Medium not known")
        
    return epsilon_medium

def plot_subplot(ax,lambdav,curves,labels,colors):
    for index in range(len(labels)):
        ax.plot(lambdav,curves[index],colors[index],label=labels[index])
    ax.set_xlabel(r'$\lambda$ [nm]')
    ax.set_ylabel(','.join(labels))
    ax.legend()

def initialize():
    da_string[0].set("10")
    db_string[0].set("20")
    da_string[1].set("20")
    db_string[1].set("10")
    d_string.set("3000")
    epsilon_s_string.set("Vacuum")
    epsilon_fa_string[0].set("TiO2")
    epsilon_fb_string[0].set("fused silica")
    epsilon_fa_string[1].set("fused silica")
    epsilon_fb_string[1].set("TiO2")
    epsilon_f_string.set("Vacuum")
    N_string[0].set("50")
    N_string[1].set("50")
    epsilon_c_string.set("Vacuum")
    lambda_min_string.set("400")
    lambda_max_string.set("2400")
        
    calculate()
    
def calculate():
    try:
        lambda_min = float(lambda_min_string.get())
        lambda_max = float(lambda_max_string.get())
        cladding = epsilon_c_string.get()
        substrate = epsilon_s_string.get()
        da = np.zeros(len(N_string))
        db = np.zeros(len(N_string))
        filma = np.empty(len(N_string), dtype=object)
        filmb = np.empty(len(N_string), dtype=object)
        N = np.zeros(len(N_string),dtype=int)
        for index in range(len(N_string)):
            da[index] = float(da_string[index].get())
            filma[index] = epsilon_fa_string[index].get()
            db[index] = float(db_string[index].get())
            filmb[index] = epsilon_fb_string[index].get()        
            N[index] = int(N_string[index].get())
        d = float(d_string.get())
        film = epsilon_f_string.get()

        if (da <= 0).any() or (db < 0).any() or (N < 0).any() or N[0] == 0 or (d<0) or lambda_min < 400 or lambda_max > 2400 or lambda_min >= lambda_max:
            gui.input_error(initialize)
        else:
            f.clf()
            lambdav = np.linspace(lambda_min, lambda_max, num=10001, endpoint=False) # vacuum wavelength in nm
            epsilon_s = epsilon(substrate,lambdav)
            epsilon_c = epsilon(cladding,lambdav)
            epsilon_fa = np.zeros([len(N),len(lambdav)])+0j
            epsilon_fb = np.zeros([len(N),len(lambdav)])+0j
            for index in range(len(N)):
                epsilon_fa[index,:] = epsilon(filma[index],lambdav)
                epsilon_fb[index,:] = epsilon(filmb[index],lambdav)
            epsilon_f = epsilon(film,lambdav)
            vreflection_transmission = np.vectorize(reflection_transmission, excluded=[1,4,7])
            R,tau = vreflection_transmission(epsilon_s,da,epsilon_fa[0,:],epsilon_fa[1,:],db,epsilon_fb[0,:],epsilon_fb[1,:],N,d,epsilon_f,epsilon_c,lambdav)
            a1 = f.add_subplot(111)
            a1.plot(lambdav,tau,'r')
            a1.set_xlabel(r'$\lambda$ [nm]')
            a1.set_ylabel(r'$\tau$')
            a1.set_xlim([lambda_min, lambda_max])
            a1.set_ylim([-0.025,1.025])
            plt.tight_layout()  
            
#            plt.savefig('FP_lambda.pdf',bbox_inches='tight',dpi=300, transparent=True)

            canvas.draw()
    except ValueError: gui.input_error(initialize)
    
f = plt.figure(1,[5,2.5])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

epsilon_s_string = Tk.StringVar()
da_string = gui.create_stringvar_vector(2)
epsilon_fa_string = gui.create_stringvar_vector(2)
db_string = gui.create_stringvar_vector(2)
epsilon_fb_string = gui.create_stringvar_vector(2)
N_string = gui.create_stringvar_vector(2)
d_string = Tk.StringVar()
epsilon_f_string = Tk.StringVar()
epsilon_c_string = Tk.StringVar()
lambda_min_string = Tk.StringVar()
lambda_max_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_radiobutton(mainframe,['substrate medium:','Vacuum','fused silica'],epsilon_s_string,2,row)
row = gui.create_triple_entry(mainframe,u"mirror 1: N =",N_string[0],u"d\u2081 [nm] =",da_string[0],u"d\u2082 [nm] =",db_string[0],row)
row = gui.create_radiobutton(mainframe,['film 1 medium:','Ag','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica','BaSF'],epsilon_fa_string[0],6,row)
row = gui.create_radiobutton(mainframe,['film 2 medium:','Ag','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica','BaSF'],epsilon_fb_string[0],6,row)
row = gui.create_radiobutton_with_entry(mainframe,u"cavity: d [nm] =",d_string,['Vacuum','fused silica'],epsilon_f_string,2,row)
row = gui.create_triple_entry(mainframe,u"mirror 2: N =",N_string[1],u"d\u2081 [nm] =",da_string[1],u"d\u2082 [nm] =",db_string[1],row)
row = gui.create_radiobutton(mainframe,['film 1 medium:','Ag','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica','BaSF'],epsilon_fa_string[1],6,row)
row = gui.create_radiobutton(mainframe,['film 2 medium:','Ag','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica','BaSF'],epsilon_fb_string[1],6,row)
row = gui.create_radiobutton(mainframe,['cladding medium:','Vacuum','fused silica'],epsilon_c_string,2,row)
row = gui.create_double_entry(mainframe,u"\u03bb [nm] >",lambda_min_string,u"\u03bb [nm] <",lambda_max_string,row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)