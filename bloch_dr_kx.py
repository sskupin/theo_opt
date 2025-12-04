import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import sys
sys.path.append('./aux')
import strat_stuff as strat
import gui_stuff as gui

gui.set_rcParams()
title = "Dispersion Relation of Bloch Modes - Fixed Frequency"
root = Tk.Tk()
root.title(title)
    
def initialize():
    var_string[0].set("0.31") # d1 in units of lambda
    var_string[1].set("2.15") # Re epsilon_f1
    var_string[2].set("0.4") # d2 in units of lambda
    var_string[3].set("7.91") # Re epsilon_f2
    var_string[5].set("TE")
    var_string[6].set("n")
    var_string[7].set("y")
    var_string[8].set("n")
    gui.copy_stringvar_vector(var_string,var_save) 
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate() 

def show_manual():
    gui.show_manual("man/bloch_dr_kx.png",title)

def calculate():
    gui.change_cursor(root,"trek")
    try:
        d1 = float(var_string[0].get())
        d2 = float(var_string[2].get())
        var_string[4].set(d1+d2)
        epsilon_f1 = float(var_string[1].get())
        epsilon_f2 = float(var_string[3].get())
        polarization = var_string[5].get()
        cone = var_string[6].get()     
        foldback = var_string[7].get()
        diffcoeff = var_string[8].get()
        
        if d1 <= 0 or d2 < 0 or epsilon_f1 <= 0 or epsilon_f2 <= 0 or np.maximum(np.sqrt(epsilon_f1),np.sqrt(epsilon_f2))*(d1+d2) > 3:
            gui.input_error("Values out of range. Re-initializing ...", reinitialize)
        else:
            f.clf()
            Kx = np.linspace(0, np.maximum(np.sqrt(epsilon_f1),np.sqrt(epsilon_f2)), num=10001, endpoint=True)*(d1+d2) # normalized transverse wavevector
            vDR = np.vectorize(strat.DR_Bloch)
            Kz,dummy1,dummy2 = vDR(d1,epsilon_f1+0j,d2,epsilon_f2+0j,polarization,Kx,d1+d2)
            
            if foldback == 'n':
                a1 = f.add_subplot(111, aspect='equal')
                a1.plot(np.real(Kz),Kx,'b',np.real(-Kz),Kx,'b')
                a1.plot(np.real(Kz-1),Kx,'b',np.real(-Kz+1),Kx,'b')
                a1.set_xlim([-1, 1])
                a1.set_ylim([Kx[0], Kx[-1]])
                a1.fill_betweenx(Kx, -1, 1, where=np.imag(Kz)!=0., facecolor='red', alpha=1, zorder=100, interpolate=True)
                a1.set_xlabel(r'$K_z^{\prime}$')
                a1.set_ylabel(r'$K_x$')
                if cone=='y':
                    a1.fill_betweenx(Kx, -1, 1, where=Kx>=d1+d2, facecolor='gray', alpha=.5, zorder=100, interpolate=True)
            else:
                a1 = plt.subplot2grid((2, 3), (0, 0), rowspan=2)
                a1.plot(np.real(Kz),Kx,'b')
                a1.set_xlim([0, .5])
                a1.xaxis.set_ticks([0, .25, .5])
                a1.set_ylim([Kx[0], Kx[-1]])
                a1.fill_betweenx(Kx, 0, 0.5, where=np.imag(Kz)!=0., facecolor='red', alpha=1, zorder=100, interpolate=True)
                a1.set_xlabel(r'$K_z^{\prime}$')
                a1.set_ylabel(r'$K_x$')                
                if cone=='y':
                    a1.fill_betweenx(Kx, 0, 0.5, where=Kx>=d1+d2, facecolor='gray', alpha=.5, zorder=100, interpolate=True)

                a2 = plt.subplot2grid((2, 3), (0, 1), colspan=2)
                a2.plot(Kx,np.real(Kz),'b')
                a2.yaxis.set_ticks([0, .25, .5])
                a2.set_xlim([Kx[0], Kx[-1]])
                a2.set_ylim([0, 0.5])
                a2.set_ylabel(r'$K_z^{\prime}$')
                a2.set_xlabel(r'$K_x$')
                a2.fill_between(Kx, 0, .5, where=np.imag(Kz)!=0., facecolor='red', alpha=1, zorder=0, interpolate=True)
                if cone=='y':
                    a2.fill_between(Kx, 0, 0.5, where=Kx>=d1+d2, facecolor='gray', alpha=.5, zorder=100, interpolate=True)

                a3 = plt.subplot2grid((2, 3), (1, 1), colspan=2)
                if diffcoeff=='y':
                    DIFF = np.gradient(np.gradient(np.real(Kz),Kx,edge_order=2),Kx,edge_order=2)
                    a3.plot(Kx,DIFF,'b')
                    a3.set_ylim([-2,2])
                    a3.set_ylabel(r'$\partial^2K_z^{\prime}/\partial K_x^2$')
                else:
                    a3.plot(Kx,np.abs(np.imag(Kz)),'b')
                    a3.set_ylim([0,np.amax(np.abs(np.imag(Kz)))*1.05])
                    a3.set_ylabel(r'$K_z^{\prime\prime}$')
                a3.set_xlim([Kx[0], Kx[-1]])
                a3.fill_between(Kx, a3.get_ylim()[0], a3.get_ylim()[1], where=np.imag(Kz)!=0., facecolor='red', alpha=1, zorder=0, interpolate=True)
                a3.set_xlabel(r'$K_x$')
                if cone=='y':
                    a3.fill_between(Kx, a3.get_ylim()[0], a3.get_ylim()[1], where=Kx>=d1+d2, facecolor='gray', alpha=.5, zorder=100, interpolate=True)
                    
            plt.tight_layout()  
            
#            plt.savefig('2D_Bloch_kx.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            gui.copy_stringvar_vector(var_string,var_save) 

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[6,4])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(9)
var_save = gui.create_stringvar_vector(9)

initialize()

row = 1
row = gui.create_title(mainframe,"unit cell parameters",row)
row = gui.create_entry(mainframe,u"film 1 thickness: d/\u03BB =",var_string[0],row)
row = gui.create_entry(mainframe,u"film 1: \u03B5 =",var_string[1],row)
row = gui.create_entry(mainframe,u"film 2 thickness: d/\u03BB =",var_string[2],row)
row = gui.create_entry(mainframe,u"film 2: \u03B5 =",var_string[3],row)
row = gui.create_label(mainframe,u"\u03A9 = \u039b/\u03BB =",var_string[4],row)
row = gui.create_radiobutton(mainframe,['polarization:','TE','TM'],var_string[5],2,row)
row = gui.create_checkbutton(mainframe,"vacuum light cone",'n','y',var_string[6],row)
row = gui.create_checkbutton(mainframe,"fold back",'n','y',var_string[7],row)
row = gui.create_checkbutton(mainframe,"show diff. coeff.",'n','y',var_string[8],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)