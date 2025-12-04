import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import sys
sys.path.append('./aux')
import gui_stuff as gui
import strat_stuff as strat
import media as media

gui.set_rcParams()
title = "Dispersion Relation of Bloch Modes - Normal Incidence"
root = Tk.Tk()
root.title(title)

def initialize():
    var_string[0].set("217.2") # thickness layer 1 in nm
    var_string[1].set("fused silica") # layer 1 medium
    var_string[2].set("122.5") # thickness layer 2 in nm
    var_string[3].set("TiO2") # layer 2 medium
    var_string[5].set("400") # minimum vacuum wavelength in nm
    var_string[6].set("2400") # maximum vacuum wavelength in nm
    var_string[7].set("y")
    gui.copy_stringvar_vector(var_string,var_save)        
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate() 

def show_manual():
    gui.show_manual("man/bloch_dr_lambda.png",title)

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
        foldback = var_string[7].get()
        
        if d1 <= 0 or d2 < 0 or (d1+d2)/lambda_min- (d1+d2)/lambda_max > 5 or lambda_min < 400 or lambda_max > 2400 or lambda_min >= lambda_max:
            gui.input_error("Values out of range. Re-initializing ...", reinitialize)
        else:
            if (film1 == "AlAs" or film1 == "AlGaAs (70% Al)" or film1 == "AlGaAs (31.5% Al)" or film2 == "AlAs" or film2 == "AlGaAs (70% Al)" or film2 == "AlGaAs (31.5% Al)" ) and lambda_min < 800:
                gui.input_error("Oops! Losses not negigible in this wavelength range ... adjusting ...")
                lambda_min = 800
                lambda_max = 2400
                var_string[5].set(lambda_min)
                var_string[6].set(lambda_max)
            f.clf()
            lambdav = np.linspace(lambda_min, lambda_max, num=10001, endpoint=True) # vacuum wavelength in nm
            epsilon_f1 = np.real(media.epsilon(film1,lambdav))
            epsilon_f2 = np.real(media.epsilon(film2,lambdav))
            Omega = (d1+d2)/lambdav
            vDR = np.vectorize(strat.DR_Bloch)
            Kz,dummy1,dummy2 = vDR(d1,epsilon_f1+0j,d2,epsilon_f2+0j,'TE',0,Omega)
            
            if foldback == 'n':
                a1 = plt.subplot2grid((1, 1), (0, 0))
                a1.plot(np.real(Kz),Omega,'b',np.real(-Kz),Omega,'b')
                a1.plot(np.real(Kz-1),Omega,'b',np.real(-Kz+1),Omega,'b')
                a1.set_xlim([-1, 1])
                a1.set_ylim([Omega[-1], Omega[0]])
                a1.fill_betweenx(Omega, -1, 1, where=np.imag(Kz)!=0., facecolor='red', alpha=1, zorder=100, interpolate=True)
                a1.set_xlabel(r'$K_z^{\prime}$')
                a1.set_ylabel(r'$\Omega$')
                a1bis = a1.twinx()
                a1bis.set_ylim(a1.get_ylim())
                a1bis.set_ylabel(r'$\lambda$ [nm]')
                a1bis.set_yticks([Omega[-1], 2*Omega[-1]/3+Omega[0]/3, Omega[-1]/3+2*Omega[0]/3, Omega[0]])
                a1bis.set_yticklabels([str(round(lambdav[-1],1)),str(round((d1+d2)/(2*Omega[-1]/3+Omega[0]/3),1)),str(round((d1+d2)/(Omega[-1]/3+2*Omega[0]/3),1)),str(round(lambdav[0],1))])
            else:               
                a1 = plt.subplot2grid((2, 3), (0, 0), rowspan=2)
                a1.plot(np.real(Kz),Omega,'b')
                a1.set_xlim([0, .5])
                a1.xaxis.set_ticks([0, .25, .5])
                a1.set_ylim([Omega[-1], Omega[0]])
                a1.fill_betweenx(Omega, 0, .5, where=np.imag(Kz)!=0., facecolor='red', alpha=1, zorder=100, interpolate=True)
                a1.set_xlabel(r'$K_z^{\prime}$')
                a1.set_ylabel(r'$\Omega$')
                
                a2 = plt.subplot2grid((2, 3), (0, 1), colspan=2)
                a2.plot(Omega,np.real(Kz),'b')
                a2.yaxis.set_ticks([0, .25, .5])
                a2.set_xlim([Omega[-1], Omega[0]])
                a2.set_ylim([0, 0.5])
                a2.set_ylabel(r'$K_z^{\prime}$')
                a2.set_xlabel(r'$\Omega$')
                a2bis = a2.twiny()
                a2bis.set_xlim(a2.get_xlim())
                a2.fill_between(Omega, 0, .5, where=np.imag(Kz)!=0., facecolor='red', alpha=1, zorder=0, interpolate=True)
                a2bis.set_xlabel(r'$\lambda$ [nm]', labelpad=10)
                a2bis.set_xticks([Omega[-1], 2*Omega[-1]/3+Omega[0]/3, Omega[-1]/3+2*Omega[0]/3, Omega[0]])
                a2bis.set_xticklabels([str(round(lambdav[-1],1)),str(round((d1+d2)/(2*Omega[-1]/3+Omega[0]/3),1)),str(round((d1+d2)/(Omega[-1]/3+2*Omega[0]/3),1)),str(round(lambdav[0],1))])
                a2bis.tick_params(direction='out', pad=0)
                
                a3 = plt.subplot2grid((2, 3), (1, 1), colspan=2)
                a3.plot(Omega,np.abs(np.imag(Kz)),'b')
                a3.set_ylim([0,np.amax(np.abs(np.imag(Kz)))*1.05])
                a3.set_xlim([Omega[-1], Omega[0]])
                a3.fill_between(Omega, 0, np.amax(np.abs(np.imag(Kz)))*1.05, where=np.imag(Kz)!=0., facecolor='red', alpha=1, zorder=0, interpolate=True)
                a3.set_ylabel(r'$K_z^{\prime\prime}$')
                a3.set_xlabel(r'$\Omega$')

            plt.tight_layout()  
            
#            plt.savefig('2D_Bloch_lambda.pdf',bbox_inches='tight',dpi=300, transparent=True)

            gui.copy_stringvar_vector(var_string,var_save) 

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[6,4])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(8)
var_save = gui.create_stringvar_vector(8)

initialize()

row = 1
row = gui.create_title(mainframe,"unit cell parameters",row)
row = gui.create_entry(mainframe,u"film 1 thickness: d [nm] =",var_string[0],row)
row = gui.create_radiobutton_with_latex(mainframe,[r'film 1 medium:',r'Vacuum',r'Al$_{0.7}$Ga$_{0.3}$GAs',r'AlAs',r'Al$_{0.315}$Ga$_{0.685}$As',r'TiO$_2$ ($\varepsilon_{\rm or}$)',r'fused silica'],['film 1 medium:','Vacuum','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica'],var_string[1],6,row)
row = gui.create_entry(mainframe,u"film 2 thickness: d [nm] =",var_string[2],row)
row = gui.create_radiobutton_with_latex(mainframe,[r'film 2 medium:',r'Vacuum',r'Al$_{0.7}$Ga$_{0.3}$GAs',r'AlAs',r'Al$_{0.315}$Ga$_{0.685}$As',r'TiO$_2$ ($\varepsilon_{\rm or}$)',r'fused silica'],['film 2 medium:','Vacuum','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica'],var_string[3],6,row)
row = gui.create_label(mainframe,u"\u039b [nm] =",var_string[4],row)
row = gui.create_double_entry(mainframe,u"\u03bb [nm] >",var_string[5],u"\u03bb [nm] <",var_string[6],row)
row = gui.create_checkbutton(mainframe,"fold back",'n','y',var_string[7],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)