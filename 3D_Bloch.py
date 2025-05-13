import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import strat_stuff as strat
import media as media
from matplotlib.colors import LogNorm

gui.set_rcParams()
title = "3D DR of Bloch modes"
root = Tk.Tk()
root.title(title)

def DR(d1,epsilon_f1,d2,epsilon_f2,polarization,Kx,Omega): # computing dispersion relation for Bloch modes
    kf1z = np.sqrt(epsilon_f1+0j-(Kx/Omega)**2)*2*np.pi
    kf2z = np.sqrt(epsilon_f2+0j-(Kx/Omega)**2)*2*np.pi
    if polarization == 'TE':
        m1 = strat.mTE(kf1z,d1*Omega/(d1+d2))
        m2 = strat.mTE(kf2z,d2*Omega/(d1+d2))
    else:
        m1 = strat.mTM(kf1z,epsilon_f1,d1*Omega/(d1+d2))
        m2 = strat.mTM(kf2z,epsilon_f2,d2*Omega/(d1+d2))
    M = np.matmul(m2,m1)
    Kz = np.arccos((M[1,1]+M[0,0])/2)/(2*np.pi)
    
    return Kz

def initialize():
    var_string[0].set("140") # thickness layer 1 in nm
    var_string[1].set("fused silica") # layer 1 medium
    var_string[2].set("180") # thickness layer 2 in nm
    var_string[3].set("TiO2") # layer 2 medium
    var_string[5].set("400") # minimum vacuum wavelength in nm
    var_string[6].set("2400") # maximum vacuum wavelength in nm
    var_string[7].set("TE")
    var_string[8].set("band")
    var_string[9].set("n")
    gui.copy_stringvar_vector(var_string,var_save)
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate() 

def calculate():
    gui.change_cursor(root,"trek")
    try:
        lambda_min = float(var_string[5].get())
        lambda_max = float(var_string[6].get())
        d1 = float(var_string[0].get())
        d2 = float(var_string[2].get())
        var_string[4].set(d1+d2)
        film1 = var_string[1].get()
        film2 = var_string[3].get()
        polarization = var_string[7].get()
        band = var_string[8].get()
        cone = var_string[9].get()        
        
        if d1 <= 0 or d2 < 0 or (d1+d2)/lambda_min- (d1+d2)/lambda_max > 1 or lambda_min < 400 or lambda_max > 2400 or lambda_min >= lambda_max:
            gui.input_error("Values out of range. Re-initializing ...", reinitialize)
        else:
            if (film1 == "AlAs" or film1 == "AlGaAs (70% Al)" or film1 == "AlGaAs (31.5% Al)" or film2 == "AlAs" or film2 == "AlGaAs (70% Al)" or film2 == "AlGaAs (31.5% Al)" ) and lambda_min < 800:
                gui.input_error("Oops! Losses not negigible in this wavelength range ... adjusting ...")
                lambda_min = 800
                lambda_max = 2400
                var_string[5].set(lambda_min)
                var_string[6].set(lambda_max)
            f.clf()
            Omega = np.linspace((d1+d2)/lambda_max, (d1+d2)/lambda_min, num=501, endpoint=True) # normalized frequency
            lambdav = (d1+d2)/Omega # vacuum wavelength in nm
            epsilon_f1 = np.real(media.epsilon(film1,lambdav))
            epsilon_f2 = np.real(media.epsilon(film2,lambdav))
            Kx = np.linspace(0, np.amax(np.maximum(np.sqrt(epsilon_f1),np.sqrt(epsilon_f2))*Omega), num=501, endpoint=True) # normalized transverse wavevector
            Kz = np.zeros([Omega.size,Kx.size]) + 0j
            vDR = np.vectorize(DR)
            for index in range(Omega.size):
                Kz[index,:] = vDR(d1,epsilon_f1[index],d2,epsilon_f2[index],polarization,Kx,Omega[index])

            a1 = plt.subplot2grid((1, 1), (0, 0))            
            if band == 'gap':
                for index in range(Omega.size):
                    Kz[index,np.where(Kx>np.maximum(np.sqrt(epsilon_f1[index]),np.sqrt(epsilon_f2[index]))*Omega[index])] = np.nan
                im = a1.imshow(np.abs(np.imag(Kz)), origin='lower', extent=[Kx[0], Kx[-1], Omega[0], Omega[-1]], cmap='jet', aspect='auto', norm=LogNorm(vmin=.01*np.amax(np.abs(np.imag(Kz))), vmax=np.amax(np.abs(np.imag(Kz)))))
                if polarization == 'TM':
                    a1.plot(np.sqrt(epsilon_f1*epsilon_f2/(epsilon_f1+epsilon_f2))*Omega,Omega,'k--')
                plt.colorbar(im, pad=0.15).set_label(r'$K_z^{\prime\prime}$')
            else:
                im = a1.imshow(np.real(Kz), origin='lower', extent=[Kx[0], Kx[-1], Omega[0], Omega[-1]], aspect='auto', cmap='jet')
                if polarization == 'TM':
                    a1.plot(np.sqrt(epsilon_f1*epsilon_f2/(epsilon_f1+epsilon_f2))*Omega,Omega,'w--')
                plt.colorbar(im, pad=0.15).set_label(r'$K_z^{\prime}$')
            if cone=='y':
                if Omega[-1]<Kx[-1]:
                    a1.fill_between(Kx,np.ones_like(Kx)*Omega[0], Kx, facecolor='gray', alpha=.5, interpolate=True)
                    a1.set_ylim([Omega[0],Omega[-1]])
                else:
                    a1.fill_between(Omega, np.ones_like(Omega)*Omega[0], Omega, facecolor='gray', alpha=.5, interpolate=True)
                    a1.set_xlim([Kx[0],Kx[-1]])
            a1.fill_between(np.maximum(np.sqrt(epsilon_f1),np.sqrt(epsilon_f2))*Omega, np.ones_like(Omega)*Omega[0], Omega, facecolor='gray', alpha=1, zorder=100, interpolate=True)
            a1.set_xlabel(r'$K_x$')
            a1.set_ylabel(r'$\Omega$')
            a1bis = a1.twinx()
            a1bis.set_ylim(a1.get_ylim())
            a1bis.set_ylabel(r'$\lambda$ [nm]')
            a1bis.set_yticks([Omega[-1], 2*Omega[-1]/3+Omega[0]/3, Omega[-1]/3+2*Omega[0]/3, Omega[0]])
            a1bis.set_yticklabels([str(round(lambdav[-1],1)),str(round((d1+d2)/(2*Omega[-1]/3+Omega[0]/3),1)),str(round((d1+d2)/(Omega[-1]/3+2*Omega[0]/3),1)),str(round(lambdav[0],1))])                                                

            plt.tight_layout()  
            
#            plt.savefig('3D_Bloch.pdf',bbox_inches='tight',dpi=300, transparent=True)

            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[6,4])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(10)
var_save = gui.create_stringvar_vector(10)

initialize()

row = 1

row = gui.create_title(mainframe,"unit cell parameters",row)
row = gui.create_entry(mainframe,u"film 1 thickness: d [nm] =",var_string[0],row)
row = gui.create_radiobutton_with_latex(mainframe,[r'film 1 medium:',r'Vacuum',r'Al$_{0.7}$Ga$_{0.3}$GAs',r'AlAs',r'Al$_{0.315}$Ga$_{0.685}$As',r'TiO$_2$ ($\varepsilon_{\rm or}$)',r'fused silica'],['film 1 medium:','Vacuum','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica'],var_string[1],6,row)
row = gui.create_entry(mainframe,u"film 2 thickness: d [nm] =",var_string[2],row)
row = gui.create_radiobutton_with_latex(mainframe,[r'film 2 medium:',r'Vacuum',r'Al$_{0.7}$Ga$_{0.3}$GAs',r'AlAs',r'Al$_{0.315}$Ga$_{0.685}$As',r'TiO$_2$ ($\varepsilon_{\rm or}$)',r'fused silica'],['film 2 medium:','Vacuum','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica'],var_string[3],6,row)
row = gui.create_label(mainframe,u"\u039b [nm] =",var_string[4],row)
row = gui.create_double_entry(mainframe,u"\u03bb [nm] >",var_string[5],u"\u03bb [nm] <",var_string[6],row)
row = gui.create_radiobutton(mainframe,['polarization:','TE','TM'],var_string[7],2,row)
row = gui.create_radiobutton(mainframe,['show:','band','gap'],var_string[8],2,row)
row = gui.create_checkbutton(mainframe,"mark area outside vacuum light cone",'n','y',var_string[9],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)