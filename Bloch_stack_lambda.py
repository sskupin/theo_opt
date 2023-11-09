import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
import tkinter as Tk
import gui_stuff as gui
import media as media

gui.set_rcParams()
title = "Reflection and Transmission of a Stack vs. Wavelength (Bloch)"
root = Tk.Tk()
root.title(title)

def m(kfz,z):
    if np.abs(kfz) == 0:
        m = np.array([[1,z],[0,1]])
    else:
        m = np.array([[np.cos(kfz*z),np.sin(kfz*z)/kfz],[-np.sin(kfz*z)*kfz,np.cos(kfz*z)]])
    return m

def MP(lambdav,d1,epsilon_f1,d2,epsilon_f2,N):
    kf1z = np.sqrt(epsilon_f1+0j)*2*np.pi/lambdav
    kf2z = np.sqrt(epsilon_f2+0j)*2*np.pi/lambdav
    m1 = m(kf1z,d1)
    m2 = m(kf2z,d2)
    M = np.matmul(m2,m1)
    kz = np.arccos((M[1,1]+M[0,0])/2)/(d1+d2)
    kappaf = -1j*(np.exp(1j*kz*(d1+d2))-M[0,0])/M[0,1]
    if np.real(kappaf) < 0:
        kz = -kz
        kappaf = -1j*(np.exp(1j*kz*(d1+d2))-M[0,0])/M[0,1]
    kappab = 1j*(np.exp(1j*kz*(d1+d2))-M[1,1])/M[0,1]
    if kappab == kappaf:
        MP = np.array([[0,1],[0,0]])
    else:
        MP = 1/(kappab-kappaf)*np.array([[kappab*np.exp(1j*kz*N*(d1+d2))-kappaf*np.exp(-1j*kz*N*(d1+d2)),-2*np.sin(kz*N*(d1+d2))],
                                              [-2*kappaf*kappab*np.sin(kz*N*(d1+d2)),-kappaf*np.exp(1j*kz*N*(d1+d2))+kappab*np.exp(-1j*kz*N*(d1+d2))]])
    return MP

def reflection_transmission(epsilon_s,da,epsilon_fa,db,epsilon_fb,N,epsilon_c,lambdav): # computing coefficients of reflection and transmission
    ksz = np.sqrt(epsilon_s)*2*np.pi/lambdav
    kcz = np.sqrt(epsilon_c)*2*np.pi/lambdav
    M = np.array([[1,0],[0,1]])
    for index in range(len(N)):
        M_P = MP(lambdav,da[index],epsilon_fa,db[index],epsilon_fb,N[index])
        if N[index] > 0:
            M = np.matmul(M_P,M)
    ND = ksz*M[1,1]+kcz*M[0,0]+1j*M[1,0]-1j*ksz*kcz*M[0,1]
    R = (ksz*M[1,1]-kcz*M[0,0]-1j*M[1,0]-1j*ksz*kcz*M[0,1])/ND
    tau = np.real(kcz)/np.real(ksz)*np.abs(2*ksz/ND)**2
    
    return R,tau

def plot_subplot(ax,lambdav,curve_1,curve_2,label_1,label_2,lambda_min, lambda_max):
    ax.plot(lambdav,curve_1,'b',label=label_1)
    ax.plot(lambdav,curve_2,'r',label=label_2)
    ax.set_xlabel(r'$\lambda$ [nm]')
    ax.set_xlim([lambda_min, lambda_max])
    ax.set_ylabel(label_1+','+label_2)
    ax.legend()
    
def epsilon(medium,lambdav):
    if medium == "Vacuum":
        epsilon_medium = (1+0j)*np.ones_like(lambdav)
    elif medium == "TiO2":
        epsilon_medium = media.TiO2(lambdav)
    elif medium == "fused silica":
        epsilon_medium = media.silica(lambdav)
    else:
        print("Oops! Medium not known")
        
    return epsilon_medium
        
def initialize():
    da_string[0].set("1")
    db_string[0].set("1")
    da_string[1].set("2")
    db_string[1].set("2")
    da_string[2].set("3")
    db_string[2].set("3")
    da_string[3].set("4")
    db_string[3].set("4")
    da_string[4].set("5")
    db_string[4].set("5")
    da_string[5].set("6")
    db_string[5].set("6")
    da_string[6].set("7")
    db_string[6].set("7")
#    da_string[7].set("8")
#    db_string[7].set("8")
#    da_string[8].set("9")
#    db_string[8].set("9")
#    da_string[9].set("10")
#    db_string[9].set("10")
    substrate_string.set("Vacuum")
    cladding_string.set("fused silica")
    filma_string.set("fused silica")
    filmb_string.set("TiO2")
    N_string[0].set("10")
    N_string[1].set("10")
    N_string[2].set("10")
    N_string[3].set("10")
    N_string[4].set("10")
    N_string[5].set("10")
    N_string[6].set("10")
#    N_string[7].set("10")
#    N_string[8].set("10")
#    N_string[9].set("10")
    lambda_min_string.set("400")
    lambda_max_string.set("2400")
       
    calculate()

def calculate():
    try:
        lambda_min = float(lambda_min_string.get())
        lambda_max = float(lambda_max_string.get())
        substrate = substrate_string.get()
        da = np.zeros(len(N_string))
        db = np.zeros(len(N_string))
        filma = filma_string.get()
        filmb = filmb_string.get()
        cladding = cladding_string.get()
        N = np.zeros(len(N_string),dtype=int)
        for index in range(len(N_string)):
            da[index] = float(da_string[index].get())
            db[index] = float(db_string[index].get())
            N[index] = int(N_string[index].get())
        
        if (da <= 0).any() or (db < 0).any() or (N < 0).any() or N[0] == 0 or (N > 50).any() or lambda_min < 400 or lambda_max > 2400 or lambda_min >= lambda_max:
            gui.input_error()
        else:
            f.clf()
            omega_max = 2*np.pi*spc.c/lambda_min*1.e-6 # maximum frequency in fs^-1
            omega_min = 2*np.pi*spc.c/lambda_max*1.e-6 # minimum frequency in fs^-1
            omegav, delta_omegav = np.linspace(omega_min, omega_max, num=10001, retstep=True, endpoint=False) # frequency in fs^-1
            lambdav = 2*np.pi*spc.c/omegav*1.e-6 # vacuum wavelength in nm
#            lambdav, delta_lambdav = np.linspace(lambda_min, lambda_max, num=10001, retstep=True, endpoint=False) # vacuum wavelength in nm
            epsilon_s = epsilon(substrate,lambdav)
            epsilon_c = epsilon(cladding,lambdav)
            epsilon_fa = epsilon(filma,lambdav)
            epsilon_fb = epsilon(filmb,lambdav)

            vreflection_transmission = np.vectorize(reflection_transmission, excluded=[1,3,5])
            R,tau = vreflection_transmission(epsilon_s,da,epsilon_fa,db,epsilon_fb,N,epsilon_c,lambdav)
            a1 = f.add_subplot(211)
            plot_subplot(a1,lambdav,np.abs(R)**2,tau,r'$\rho$',r'$\tau$',lambda_min, lambda_max)
            a1.set_ylim([-0.025,1.025])
            a2 = f.add_subplot(212)
            theta = np.unwrap(np.angle(R))
#            thetad = (theta[2:]-theta[:-2])/(2*delta_lambdav)
#            thetad2 = (theta[2:]+theta[:-2]-2*theta[1:-1])/delta_lambdav**2
#            a2.plot(lambdav[1:-1],.5*(thetad2*lambdav[1:-1]**4/(4*np.pi**2*spc.c**2)+thetad*lambdav[1:-1]**3/(2*np.pi**2*spc.c**2))*1e12,'b')
#            a2.plot(lambdav[1:-1],.5*(thetad2*lambdav[1:-1]**4/(4*np.pi**2*spc.c**2))*1e12,lambdav[1:-1],0.5*(thetad*lambdav[1:-1]**3/(2*np.pi**2*spc.c**2))*1e12)
            a2.plot(lambdav[1:-1],(theta[2:]+theta[:-2]-2*theta[1:-1])/delta_omegav**2,'b')
#            a2.plot(lambdav,theta,'b')
            a2.set_xlabel(r'$\lambda$ [nm]')
            a2.set_xlim([lambda_min, lambda_max])
            a2.set_ylabel(r'$\frac{\partial^2\theta}{\partial \omega^2}$ [fs$^2$]')
                     
            plt.tight_layout()  
            
#            plt.savefig('Bloch_stack_lambda.pdf',bbox_inches='tight',dpi=300, transparent=True)

            canvas.draw()
    except ValueError: gui.input_error()

f = plt.figure(1,[5,5])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

substrate_string = Tk.StringVar()
da_string = gui.create_stringvar_vector(7)
filma_string = Tk.StringVar()
db_string = gui.create_stringvar_vector(7)
filmb_string = Tk.StringVar()
N_string = gui.create_stringvar_vector(7)
cladding_string = Tk.StringVar()
lambda_min_string = Tk.StringVar()
lambda_max_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_radiobutton(mainframe,['substrate medium:','Vacuum','fused silica'],substrate_string,2,row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"medium 1: N =",N_string[0],row)
row = gui.create_double_entry(mainframe,u"fused silica: d [nm] =",da_string[0],u"TiO2: d [nm] =",db_string[0],row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"medium 2: N =",N_string[1],row)
row = gui.create_double_entry(mainframe,u"fused silica: d [nm] =",da_string[1],u"TiO2: d [nm] =",db_string[1],row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"medium 3: N =",N_string[2],row)
row = gui.create_double_entry(mainframe,u"fused silica: d [nm] =",da_string[2],u"TiO2: d [nm] =",db_string[2],row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"medium 4: N =",N_string[3],row)
row = gui.create_double_entry(mainframe,u"fused silica: d [nm] =",da_string[3],u"TiO2: d [nm] =",db_string[3],row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"medium 5: N =",N_string[4],row)
row = gui.create_double_entry(mainframe,u"fused silica: d [nm] =",da_string[4],u"TiO2: d [nm] =",db_string[4],row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"medium 6: N =",N_string[5],row)
row = gui.create_double_entry(mainframe,u"fused silica: d [nm] =",da_string[5],u"TiO2: d [nm] =",db_string[5],row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"medium 7: N =",N_string[6],row)
row = gui.create_double_entry(mainframe,u"fused silica: d [nm] =",da_string[6],u"TiO2: d [nm] =",db_string[6],row)
#row = gui.create_spacer(mainframe,row)
#row = gui.create_entry(mainframe,u"medium 8: N =",N_string[7],row)
#row = gui.create_double_entry(mainframe,u"fused silica: d [nm] =",da_string[7],u"TiO2: d [nm] =",db_string[7],row)
#row = gui.create_spacer(mainframe,row)
#row = gui.create_entry(mainframe,u"medium 9: N =",N_string[8],row)
#row = gui.create_double_entry(mainframe,u"fused silica: d [nm] =",da_string[8],u"TiO2: d [nm] =",db_string[8],row)
#row = gui.create_spacer(mainframe,row)
#row = gui.create_entry(mainframe,u"medium 10: N =",N_string[9],row)
#row = gui.create_double_entry(mainframe,u"fused silica: d [nm] =",da_string[9],u"TiO2: d [nm] =",db_string[9],row)
row = gui.create_radiobutton(mainframe,['cladding medium:','Vacuum','fused silica'],cladding_string,2,row)
row = gui.create_double_entry(mainframe,u"\u03bb [nm] >",lambda_min_string,u"\u03bb [nm] <",lambda_max_string,row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)