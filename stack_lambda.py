import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import media as media

gui.set_rcParams()
title = "Reflection and Transmission of a Stack vs. Wavelength"
root = Tk.Tk()
root.title(title)

def m(kfz,z):
    return np.array([[np.cos(kfz*z),np.sin(kfz*z)/kfz],[-np.sin(kfz*z)*kfz,np.cos(kfz*z)]])

def reflection_transmission(epsilon_s,d1,epsilon_f1,d2,epsilon_f2,N,epsilon_c,lambdav): # computing coefficients of reflection and transmission
    ksz = np.sqrt(epsilon_s)*2*np.pi/lambdav
    kf1z = np.sqrt(epsilon_f1)*2*np.pi/lambdav
    kf2z = np.sqrt(epsilon_f2)*2*np.pi/lambdav
    kcz = np.sqrt(epsilon_c)*2*np.pi/lambdav
    m1 = m(kf1z,d1)
    m2 = m(kf2z,d2)
    M = np.linalg.matrix_power(np.matmul(m2,m1),N)
    DENOM = ksz*M[1,1]+kcz*M[0,0]+1j*M[1,0]-1j*ksz*kcz*M[0,1]
    R = (ksz*M[1,1]-kcz*M[0,0]-1j*M[1,0]-1j*ksz*kcz*M[0,1])/DENOM
    tau = np.real(kcz)/np.real(ksz)*np.abs(2*ksz/DENOM)**2
    
    return R,tau

def plot_subplot_twinx(ax,lambdav,curves,labels,colors):
    axbis = ax.twinx()
    axbis.yaxis.get_major_formatter().set_powerlimits((0, 3))
    lns1 = ax.plot(lambdav,curves[0],colors[0],label=labels[0])
    ax.set_xlabel(r'$\lambda$ [nm]')
    ax.set_ylabel(labels[0])
    lns2 = axbis.plot(lambdav,curves[1],colors[1],label=labels[1])
    axbis.set_ylabel(labels[1])
    ax.legend(lns1+lns2,labels)
    return axbis
    
def plot_subplot(ax,lambdav,curves,labels,colors):
    for index in range(len(labels)):
        ax.plot(lambdav,curves[index],colors[index],label=labels[index])
    ax.set_xlabel(r'$\lambda$ [nm]')
    ax.set_ylabel(','.join(labels))
    ax.legend()
    
def epsilon(medium,lambdav):
    if medium == "Vacuum":
        epsilon_medium = (1+0j)*np.ones_like(lambdav)
    elif medium == "Si":
        epsilon_medium = media.Si(lambdav)
    elif medium == "AlAs":
        epsilon_medium = media.AlAs(lambdav)
    elif medium == "GaAs":
        epsilon_medium = media.GaAs(lambdav)
    elif medium == "AlGaAs (70% Al)":
        epsilon_medium = media.AlGaAs70(lambdav)
    elif medium == "AlGaAs (31.5% Al)":
        epsilon_medium = media.AlGaAs31(lambdav)
    elif medium == "TiO2":
        epsilon_medium = media.TiO2(lambdav)
    elif medium == "Ag":
        epsilon_medium = media.Ag(lambdav)
    elif medium == "fused silica":
        epsilon_medium = media.silica(lambdav)
    elif medium == "BaSF":
        epsilon_medium = media.BaSF(lambdav)
    else:
        print("Oops! Medium not known")
        
    return epsilon_medium

def initialize():
    lambda_min_string.set("400")
    lambda_max_string.set("800")
    d1_string.set("543")
    d2_string.set("0")
    film1_string.set("fused silica")
    film2_string.set("fused silica")
    N_string.set("1")
    cladding_string.set("Si")
    substrate_string.set("Vacuum")
        
    calculate()
    
def reinitialize():
    lambda_min_string.set("400")
    lambda_max_string.set("2400")
    d1_string.set("75")
    d2_string.set("90")
    film1_string.set("GaAs")
    film2_string.set("AlAs")
    N_string.set("10")
    cladding_string.set("GaAs")
    substrate_string.set("Vacuum")
    
    calculate()
    
def calculate():
    try:
        lambda_min = float(lambda_min_string.get())
        lambda_max = float(lambda_max_string.get())
        substrate = substrate_string.get()
        d1 = float(d1_string.get())
        d2 = float(d2_string.get())
        N = int(N_string.get())
        film1 = film1_string.get()
        film2 = film2_string.get()
        cladding = cladding_string.get()
        
        if d1 < 0 or d2 < 0 or N < 0 or lambda_min < 400 or lambda_max > 2400 or lambda_min >= lambda_max:
            gui.input_error("bla",reinitialize)
        else:
            f.clf()
            lambdav = np.linspace(lambda_min, lambda_max, num=10001, endpoint=False) # vacuum wavelength in nm
            epsilon_s = epsilon(substrate,lambdav)
            epsilon_c = epsilon(cladding,lambdav)
            epsilon_f1 = epsilon(film1,lambdav)
            epsilon_f2 = epsilon(film2,lambdav)
            epsilonp_min = min(np.amin(np.real(epsilon_f1)),np.amin(np.real(epsilon_f2)),np.amin(np.real(epsilon_c)))
            epsilonp_min = epsilonp_min - .1*np.abs(epsilonp_min)
            epsilonp_max = max(np.amax(np.real(epsilon_f1)),np.amax(np.real(epsilon_f2)),np.amax(np.real(epsilon_c)))
            epsilonp_max = epsilonp_max + .1*np.abs(epsilonp_max)
            epsilonpp_max = 1.1*max(np.amax(np.imag(epsilon_f1)),np.amax(np.imag(epsilon_f2)),np.amax(np.imag(epsilon_c)))
            vreflection_transmission = np.vectorize(reflection_transmission)
            R,tau = vreflection_transmission(epsilon_s,d1,epsilon_f1,d2,epsilon_f2,N,epsilon_c,lambdav)
            a1 = f.add_subplot(221)
#            a1 = f.add_subplot(231)
            plot_subplot(a1,lambdav,[np.abs(R)**2,tau,1-np.abs(R)**2-tau],[r'$\rho$',r'$\tau$',r'$a$'],['b','r','g'])
            a1.set_xlim([lambda_min, lambda_max])
            a1.set_ylim([-0.025,1.025])
            a2 = f.add_subplot(222)
#            a2 = f.add_subplot(232)
            a2bis = plot_subplot_twinx(a2,lambdav,[np.real(epsilon_c),np.imag(epsilon_c)],[r'$\varepsilon_{\rm c}^{\prime}$',r'$\varepsilon_{\rm c}^{\prime\prime}$'],['b','r'])
            a2.set_xlim([lambda_min, lambda_max])
            a2.set_ylim([epsilonp_min,epsilonp_max])
            a2bis.set_ylim([0,max(1.e-3,epsilonpp_max)])
            a2.text(0.5, 0.925, r'cladding', verticalalignment='center', horizontalalignment='center', transform=a2.transAxes)
#            a2.text(0.5, 0.925, r'Si cladding', verticalalignment='center', horizontalalignment='center', transform=a2.transAxes)
            a3 = f.add_subplot(223)
            a3bis = plot_subplot_twinx(a3,lambdav,[np.real(epsilon_f1),np.imag(epsilon_f1)],[r'$\varepsilon_{\rm f1}^{\prime}$',r'$\varepsilon_{\rm f1}^{\prime\prime}$'],['b','r'])
            a3.set_xlim([lambda_min, lambda_max])
            a3.set_ylim([epsilonp_min,epsilonp_max])
            a3bis.set_ylim([0,max(1.e-3,epsilonpp_max)])
            a3.text(0.5, 0.925, r'film 1', verticalalignment='center', horizontalalignment='center', transform=a3.transAxes)
            a4 = f.add_subplot(224)
            a4bis = plot_subplot_twinx(a4,lambdav,[np.real(epsilon_f2),np.imag(epsilon_f2)],[r'$\varepsilon_{\rm f2}^{\prime}$',r'$\varepsilon_{\rm f2}^{\prime\prime}$'],['b','r'])
            a4.set_xlim([lambda_min, lambda_max])
            a4.set_ylim([epsilonp_min,epsilonp_max])
            a4bis.set_ylim([0,max(1.e-3,epsilonpp_max)])
            a4.text(0.5, 0.925, r'film 2', verticalalignment='center', horizontalalignment='center', transform=a4.transAxes)

#            a3 = f.add_subplot(233)
#            plot_subplot_twinx(a3,lambdav,[np.real(epsilon_f1),np.imag(epsilon_f1)],[r'$\varepsilon_{\rm f}^{\prime}$',r'$\varepsilon_{\rm f}^{\prime\prime}$'],['b','r'])
#            a3.set_xlim([lambda_min, lambda_max])
#            a3.text(0.5, 0.925, r'BaSF glass film', verticalalignment='center', horizontalalignment='center', transform=a3.transAxes)

            plt.tight_layout()  
            
#            plt.savefig('stack_lambda.pdf',bbox_inches='tight',dpi=300, transparent=True)

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[8,6])
#f = plt.figure(1,[19,10])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

lambda_min_string = Tk.StringVar()
lambda_max_string = Tk.StringVar()
d1_string = Tk.StringVar()
film1_string = Tk.StringVar()
d2_string = Tk.StringVar()
film2_string = Tk.StringVar()
N_string = Tk.StringVar()
cladding_string = Tk.StringVar()
substrate_string = Tk.StringVar()

initialize()

row = 1
#row = gui.create_entry(mainframe,u"min. wavelength: \u03bb [nm] =",lambda_min_string,row)
#row = gui.create_entry(mainframe,u"max. wavelength: \u03bb [nm] =",lambda_max_string,row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_radiobutton(mainframe,['substrate medium:','Vacuum','fused silica'],substrate_string,2,row)
#row = gui.create_spacer(mainframe,row)
#row = gui.create_title(mainframe,"stack parameters",row)
row = gui.create_entry(mainframe,u"film 1 thickness: d [nm] =",d1_string,row)
row = gui.create_radiobutton(mainframe,['film 1 medium:','Vacuum','Si','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','GaAs','TiO2','Ag','fused silica','BaSF'],film1_string,10,row)
row = gui.create_entry(mainframe,u"film 2 thickness: d [nm] =",d2_string,row)
row = gui.create_radiobutton(mainframe,['film 2 medium:','Vacuum','Si','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','GaAs','TiO2','Ag','fused silica','BaSF'],film2_string,10,row)
row = gui.create_entry(mainframe,u"Number of periods =",N_string,row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_radiobutton(mainframe,['cladding medium:','Vacuum','Si','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','GaAs','TiO2','Ag','fused silica','BaSF'],cladding_string,10,row)
#row = gui.create_spacer(mainframe,row)
row = gui.create_double_entry(mainframe,u"\u03bb [nm] >",lambda_min_string,u"\u03bb [nm] <",lambda_max_string,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)