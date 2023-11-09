import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import media as media

gui.set_rcParams()
title = "Fabry-Perot Resonator with 2f Imaging (incident Gaussian)"
root = Tk.Tk()
root.title(title)

def m(kfz,z):
    return np.array([[np.cos(kfz*z),np.sin(kfz*z)/kfz],[-np.sin(kfz*z)*kfz,np.cos(kfz*z)]])

def MP(d1,epsilon_f1,d2,epsilon_f2,N,lambdav,kx):
    kf1z = np.sqrt(epsilon_f1-kx**2)*2*np.pi/lambdav
    kf2z = np.sqrt(epsilon_f2-kx**2)*2*np.pi/lambdav
    m1 = m(kf1z,d1)
    m2 = m(kf2z,d2)

    return np.linalg.matrix_power(np.matmul(m2,m1),N)

def reflection_transmission(epsilon_s,da,epsilon_fa,db,epsilon_fb,N,d,epsilon_f,epsilon_c,lambdav,kx): # computing coefficients of reflection and transmission
    ksz = np.sqrt(epsilon_s-kx**2)*2*np.pi/lambdav
    kfz = np.sqrt(epsilon_f-kx**2)*2*np.pi/lambdav
    kcz = np.sqrt(epsilon_c-kx**2)*2*np.pi/lambdav
    MP1 = MP(da[0],epsilon_fa[0],db[0],epsilon_fb[0],N[0],lambdav,kx)
    M = np.matmul(m(kfz,d),MP1)
    MP2 = MP(da[1],epsilon_fa[1],db[1],epsilon_fb[1],N[1],lambdav,kx)
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
    w0_string.set("15") # in microns
    da_string[0].set("30")
    db_string[0].set("0")
    da_string[1].set("30")
    db_string[1].set("0")
    d_string.set("10") # in mm
    epsilon_s_string.set("Vacuum")
    epsilon_fa_string[0].set("Ag")
    epsilon_fb_string[0].set("fused silica")
    epsilon_fa_string[1].set("Ag")
    epsilon_fb_string[1].set("TiO2")
    epsilon_f_string.set("Vacuum")
    N_string[0].set("1")
    N_string[1].set("1")
    epsilon_c_string.set("Vacuum")
    lambdav_string.set("600")
    f_string.set("100") # in mm
        
    calculate()
    
def calculate():
    try:
        lambdav = float(lambdav_string.get())
        w0 = float(w0_string.get())
        flense = float(f_string.get())
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
        d = float(d_string.get())*1000000
        film = epsilon_f_string.get()

        if (da <= 0).any() or (db < 0).any() or (N < 0).any() or N[0] == 0 or (d<0) or lambdav < 400 or lambdav > 2400 or w0*1000 < 10*lambdav:
            gui.input_error(initialize)
        else:
            f.clf()
            k0 = 2*np.pi/lambdav
            kx = np.linspace(0, 2.5/(w0*1000*k0), num=10001, endpoint=False) # kx in k0
            epsilon_s = epsilon(substrate,lambdav)
            epsilon_c = epsilon(cladding,lambdav)
            epsilon_fa = np.zeros(len(N))+0j
            epsilon_fb = np.zeros(len(N))+0j
            for index in range(len(N)):
                epsilon_fa[index] = epsilon(filma[index],lambdav)
                epsilon_fb[index] = epsilon(filmb[index],lambdav)
            epsilon_f = epsilon(film,lambdav)
            vreflection_transmission = np.vectorize(reflection_transmission, excluded=[1,2,3,4,5])
            R,tau = vreflection_transmission(epsilon_s,da,epsilon_fa,db,epsilon_fb,N,d,epsilon_f,epsilon_c,lambdav,kx)
            azimuths = np.radians(np.linspace(0, 360, 200))
            zeniths = kx*flense
            r, theta = np.meshgrid(zeniths, azimuths)
            values = np.zeros((azimuths.size, zeniths.size))
            for index in range(azimuths.size):
                values[index,:] = tau*np.exp(-2*kx**2*k0**2*w0**2*1000**2/4)
            a1 = f.subplots(subplot_kw=dict(projection='polar'))
            a1.contourf(theta, r, values, levels=100, cmap="gray")
            a1.set_rlabel_position(0)
            a1.set_thetagrids([])
            a1.grid(False)
            label_position=a1.get_rlabel_position()
            a1.text(np.radians(label_position-5),a1.get_rmax()/2,'r [mm]',
                    rotation=label_position,ha='center',va='center',color='w')
            a1.tick_params(axis='y', colors='w')
            
            plt.tight_layout()  
            
#            plt.savefig('FP_image.pdf',bbox_inches='tight',dpi=300, transparent=True)

            canvas.draw()
    except ValueError: gui.input_error(initialize)
    
f = plt.figure(1,[6,6])
canvas = gui.create_canvas(root,f)
canvas.draw()
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
w0_string = Tk.StringVar()
lambdav_string = Tk.StringVar()
f_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_radiobutton(mainframe,['substrate medium:','Vacuum','fused silica'],epsilon_s_string,2,row)
row = gui.create_triple_entry(mainframe,u"mirror 1: N =",N_string[0],u"d\u2081 [nm] =",da_string[0],u"d\u2082 [nm] =",db_string[0],row)
row = gui.create_radiobutton(mainframe,['film 1 medium:','Ag','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica','BaSF'],epsilon_fa_string[0],6,row)
row = gui.create_radiobutton(mainframe,['film 2 medium:','Ag','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica','BaSF'],epsilon_fb_string[0],6,row)
row = gui.create_radiobutton_with_entry(mainframe,u"cavity: d [mm] =",d_string,['Vacuum','fused silica'],epsilon_f_string,2,row)
row = gui.create_triple_entry(mainframe,u"mirror 2: N =",N_string[1],u"d\u2081 [nm] =",da_string[1],u"d\u2082 [nm] =",db_string[1],row)
row = gui.create_radiobutton(mainframe,['film 1 medium:','Ag','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica','BaSF'],epsilon_fa_string[1],6,row)
row = gui.create_radiobutton(mainframe,['film 2 medium:','Ag','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica','BaSF'],epsilon_fb_string[1],6,row)
row = gui.create_radiobutton(mainframe,['cladding medium:','Vacuum','fused silica'],epsilon_c_string,2,row)
row = gui.create_triple_entry(mainframe,u"w\u2080 [\u03bcm] =",w0_string,u"\u03bb [nm] =",lambdav_string,u"f [mm] =",f_string,row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)