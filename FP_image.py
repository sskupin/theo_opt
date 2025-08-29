import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import strat_stuff as strat
import media as media

gui.set_rcParams()
title = "Far-Field Image of the Fabry-Perot Output (Paraxial)"
root = Tk.Tk()
root.title(title)

def MP(d1,epsilon_f1,d2,epsilon_f2,N,lambdav,kx):
    kf1z = np.sqrt(epsilon_f1-kx**2)*2*np.pi/lambdav
    kf2z = np.sqrt(epsilon_f2-kx**2)*2*np.pi/lambdav
    m1 = strat.mTE(kf1z,d1)
    m2 = strat.mTE(kf2z,d2)

    return np.linalg.matrix_power(np.matmul(m2,m1),N)

def reflection_transmission(epsilon_s,da,epsilon_fa,db,epsilon_fb,N,d,epsilon_f,epsilon_c,lambdav,kx): # computing coefficients of reflection and transmission
    ksz = np.sqrt(epsilon_s-kx**2)*2*np.pi/lambdav
    kfz = np.sqrt(epsilon_f-kx**2)*2*np.pi/lambdav
    kcz = np.sqrt(epsilon_c-kx**2)*2*np.pi/lambdav
    MP1 = MP(da[0],epsilon_fa[0],db[0],epsilon_fb[0],N[0],lambdav,kx)
    M = np.matmul(strat.mTE(kfz,d),MP1)
    MP2 = MP(da[1],epsilon_fa[1],db[1],epsilon_fb[1],N[1],lambdav,kx)
    M = np.matmul(MP2,M)
    DENOM = ksz*M[1,1]+kcz*M[0,0]+1j*M[1,0]-1j*ksz*kcz*M[0,1]
    R = (ksz*M[1,1]-kcz*M[0,0]-1j*M[1,0]-1j*ksz*kcz*M[0,1])/DENOM
    tau = np.real(kcz)/np.real(ksz)*np.abs(2*ksz/DENOM)**2
    return R,tau

def initialize():
    var_string[0].set("Vacuum") # substrate medium
    var_string[1].set("93") # mirror 1: thickness layer 1 in nm
    var_string[2].set("fused silica") # mirror 1: layer 1 medium
    var_string[3].set("82") # mirror 1: thickness layer 2 in nm
    var_string[4].set("TiO2") # mirror 1: layer 2 medium
    var_string[5].set("5") # mirror 1: Number of periods
    var_string[6].set("10") # cavity thickness D in mm
    var_string[7].set("fused silica") # cavity medium
    var_string[8].set("82") # mirror 2: thickness layer 1 in nm
    var_string[9].set("TiO2") # mirror 2: layer 1 medium
    var_string[10].set("93") # mirror 2: thickness layer 2 in nm
    var_string[11].set("fused silica") # mirror 2: layer 2 medium
    var_string[12].set("5") # mirror 2: Number of periods
    var_string[13].set("Vacuum") # cladding medium
    var_string[14].set("10") # in microns
    var_string[15].set("633") # in nm
    var_string[16].set("100") # in mm
    gui.copy_stringvar_vector(var_string,var_save) 
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()  

def calculate():
    gui.change_cursor(root,"trek")
    try:
        lambdav = float(var_string[15].get())
        w0 = float(var_string[14].get())
        flense = float(var_string[16].get())
        cladding = var_string[13].get()
        substrate = var_string[0].get()
        da = np.zeros(2)
        db = np.zeros(2)
        filma = np.empty(2, dtype=object)
        filmb = np.empty(2, dtype=object)
        N = np.zeros(2,dtype=int)
        da[0] = float(var_string[1].get())
        da[1] = float(var_string[8].get())
        filma[0] = var_string[2].get()
        filma[1] = var_string[9].get()
        db[0] = float(var_string[3].get())
        db[1] = float(var_string[10].get())
        filmb[0] = var_string[4].get()
        filmb[1] = var_string[11].get()
        N[0] = int(var_string[5].get())
        N[1] = int(var_string[12].get())
        D = float(var_string[6].get())*1000000
        film = var_string[7].get()

        if (da <= 0).any() or (db < 0).any() or (D<0) or flense <= 0 or w0 <= 0:
            gui.input_error("Values out of range. Re-initializing ...", reinitialize)
        elif lambdav < 400 or lambdav > 2400:
            gui.input_error("Wavelength range between 400 and 2400 nm. Re-initializing ...", reinitialize)
        elif (N < 0).any() or (N > 50).any():
            gui.input_error("Number of periods must be between 0 and 50. Re-initializing ...", reinitialize)
        else:
            f.clf()
            k0 = 2*np.pi/lambdav
            kx = np.linspace(0, 2.5/(w0*1000*k0), num=10001, endpoint=False) # kx in k0
            epsilon_s = media.epsilon(substrate,lambdav)
            epsilon_c = media.epsilon(cladding,lambdav)
            epsilon_fa = np.zeros(len(N))+0j
            epsilon_fb = np.zeros(len(N))+0j
            for index in range(len(N)):
                epsilon_fa[index] = media.epsilon(filma[index],lambdav)
                epsilon_fb[index] = media.epsilon(filmb[index],lambdav)
            epsilon_f = media.epsilon(film,lambdav)
            vreflection_transmission = np.vectorize(reflection_transmission, excluded=[1,2,3,4,5])
            R,tau = vreflection_transmission(epsilon_s,da,epsilon_fa,db,epsilon_fb,N,D,epsilon_f,epsilon_c,lambdav,kx)
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

            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")
    
f = plt.figure(1,[6,6])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(17)
var_save = gui.create_stringvar_vector(17)

initialize()

row = 1
row = gui.create_radiobutton(mainframe,['substrate medium:','Vacuum','fused silica'],var_string[0],2,row)
row = gui.create_triple_entry(mainframe,u"mirror 1: N =",var_string[5],u"d\u2081 [nm] =",var_string[1],u"d\u2082 [nm] =",var_string[3],row)
row = gui.create_radiobutton(mainframe,['film 1 medium:','Ag','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica','BaSF'],var_string[2],6,row)
row = gui.create_radiobutton(mainframe,['film 2 medium:','Ag','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica','BaSF'],var_string[4],6,row)
row = gui.create_radiobutton_with_entry(mainframe,u"cavity: D [mm] =",var_string[6],['Vacuum','fused silica'],var_string[7],2,row)
row = gui.create_triple_entry(mainframe,u"mirror 2: N =",var_string[12],u"d\u2081 [nm] =",var_string[8],u"d\u2082 [nm] =",var_string[10],row)
row = gui.create_radiobutton(mainframe,['film 1 medium:','Ag','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica','BaSF'],var_string[9],6,row)
row = gui.create_radiobutton(mainframe,['film 2 medium:','Ag','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica','BaSF'],var_string[11],6,row)
row = gui.create_radiobutton(mainframe,['cladding medium:','Vacuum','fused silica'],var_string[13],2,row)
row = gui.create_triple_entry(mainframe,u"w\u2080 [\u03bcm] =",var_string[14],u"\u03bb [nm] =",var_string[15],u"f [mm] =",var_string[16],row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)