import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import media as media

gui.set_rcParams()
title = "2D DR of Bloch modes for normal incidence"
root = Tk.Tk()
root.title(title)

def m(kfz,z):
    return np.array([[np.cos(kfz*z),np.sin(kfz*z)/kfz],[-np.sin(kfz*z)*kfz,np.cos(kfz*z)]])

def DR(d1,epsilon_f1,d2,epsilon_f2,lambdav): # computing dispersion relation for Bloch modes
    kf1z = np.sqrt(epsilon_f1+0j)*2*np.pi/lambdav
    kf2z = np.sqrt(epsilon_f2+0j)*2*np.pi/lambdav
    m1 = m(kf1z,d1)
    m2 = m(kf2z,d2)
    M = np.matmul(m2,m1)
    Kz = np.arccos((M[1,1]+M[0,0])/2)/(2*np.pi)
    Omega = (d1+d2)/lambdav
    
    return Kz,Omega
    
def epsilon(medium,lambdav):
    if medium == "Vacuum":
        epsilon_medium = (1+0j)*np.ones_like(lambdav)
    elif medium == "AlAs":
        epsilon_medium = media.AlAs(lambdav)
    elif medium == "AlGaAs (70% Al)":
        epsilon_medium = media.AlGaAs70(lambdav)
    elif medium == "AlGaAs (31.5% Al)":
        epsilon_medium = media.AlGaAs31(lambdav)
    elif medium == "TiO2":
        epsilon_medium = media.TiO2(lambdav)
    elif medium == "fused silica":
        epsilon_medium = media.silica(lambdav)
    else:
        print("Oops! Medium not known")
        
    return epsilon_medium

def initialize():
    lambda_min_string.set("400")
    lambda_max_string.set("2400")
    d1_string.set("75")
    d2_string.set("90")
    film1_string.set("fused silica")
    film2_string.set("TiO2")
    foldback_string.set("n")
        
    calculate()

def calculate():
    try:
        lambda_min = float(lambda_min_string.get())
        lambda_max = float(lambda_max_string.get())
        d1 = float(d1_string.get())
        d2 = float(d2_string.get())
        D_string.set(d1+d2)
        film1 = film1_string.get()
        film2 = film2_string.get()
        foldback = foldback_string.get()
        
        if d1 <= 0 or d2 < 0 or (d1+d2)/lambda_min- (d1+d2)/lambda_max > 5 or lambda_min < 400 or lambda_max > 2400 or lambda_min >= lambda_max:
            gui.input_error(initialize)
        else:
            if (film1 == "AlAs" or film1 == "AlGaAs (70% Al)" or film1 == "AlGaAs (31.5% Al)" or film2 == "AlAs" or film2 == "AlGaAs (70% Al)" or film2 == "AlGaAs (31.5% Al)" ) and lambda_min < 800:
                print("Oops! Losses not negigible in this wavelength range ... adjusting ...")
                lambda_min = 800
                lambda_max = 2400
                lambda_min_string.set(lambda_min)
                lambda_max_string.set(lambda_max)
            f.clf()
            lambdav = np.linspace(lambda_min, lambda_max, num=10001, endpoint=True) # vacuum wavelength in nm
            epsilon_f1 = np.real(epsilon(film1,lambdav))
            epsilon_f2 = np.real(epsilon(film2,lambdav))
            vDR = np.vectorize(DR)
            Kz,Omega = vDR(d1,epsilon_f1,d2,epsilon_f2,lambdav)
            
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
                a2.set_ylabel(r'$K_z^{\prime}$')
                a2.set_xlabel(r'$\Omega$')
                a2bis = a2.twiny()
                a2bis.set_xlim(a2.get_xlim())
                a2bis.set_xlabel(r'$\lambda$ [nm]', labelpad=10)
                a2bis.set_xticks([Omega[-1], 2*Omega[-1]/3+Omega[0]/3, Omega[-1]/3+2*Omega[0]/3, Omega[0]])
                a2bis.set_xticklabels([str(round(lambdav[-1],1)),str(round((d1+d2)/(2*Omega[-1]/3+Omega[0]/3),1)),str(round((d1+d2)/(Omega[-1]/3+2*Omega[0]/3),1)),str(round(lambdav[0],1))])
                a2bis.tick_params(direction='out', pad=0)
                
                a3 = plt.subplot2grid((2, 3), (1, 1), colspan=2)
                a3.plot(Omega,np.abs(np.imag(Kz)),'b')
                a3.set_ylim([0,np.amax(np.abs(np.imag(Kz)))*1.05])
                a3.set_xlim([Omega[-1], Omega[0]])
                a3.set_ylabel(r'$K_z^{\prime\prime}$')
                a3.set_xlabel(r'$\Omega$')

            plt.tight_layout()  
            
#            plt.savefig('2D_Bloch_lambda.pdf',bbox_inches='tight',dpi=300, transparent=True)

            canvas.draw()
    except ValueError: gui.input_error(initialize)

f = plt.figure(1,[6,4])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

lambda_min_string = Tk.StringVar()
lambda_max_string = Tk.StringVar()
d1_string = Tk.StringVar()
film1_string = Tk.StringVar()
d2_string = Tk.StringVar()
film2_string = Tk.StringVar()
D_string = Tk.StringVar()
foldback_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_title(mainframe,"unit cell parameters",row)
row = gui.create_entry(mainframe,u"film 1 thickness: d [nm] =",d1_string,row)
row = gui.create_radiobutton(mainframe,['film 1 medium:','Vacuum','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica'],film1_string,6,row)
row = gui.create_entry(mainframe,u"film 2 thickness: d [nm] =",d2_string,row)
row = gui.create_radiobutton(mainframe,['film 2 medium:','Vacuum','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica'],film2_string,6,row)
row = gui.create_label(mainframe,u"\u039b [nm] =",D_string,row)
row = gui.create_double_entry(mainframe,u"\u03bb [nm] >",lambda_min_string,u"\u03bb [nm] <",lambda_max_string,row)
row = gui.create_checkbutton(mainframe,"fold back",'n','y',foldback_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)