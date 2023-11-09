import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import media as media
from matplotlib.colors import LogNorm

gui.set_rcParams()
title = "3D DR of Bloch modes"
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

def DR(d1,epsilon_f1,d2,epsilon_f2,polarization,Kx,Omega): # computing dispersion relation for Bloch modes
    kf1z = np.sqrt(epsilon_f1+0j-(Kx/Omega)**2)*2*np.pi
    kf2z = np.sqrt(epsilon_f2+0j-(Kx/Omega)**2)*2*np.pi
    if polarization == 'TE':
        m1 = mTE(kf1z,d1*Omega/(d1+d2))
        m2 = mTE(kf2z,d2*Omega/(d1+d2))
    else:
        m1 = mTM(kf1z,epsilon_f1,d1*Omega/(d1+d2))
        m2 = mTM(kf2z,epsilon_f2,d2*Omega/(d1+d2))
    M = np.matmul(m2,m1)
    Kz = np.arccos((M[1,1]+M[0,0])/2)/(2*np.pi)
    
    return Kz
    
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
    polarization_string.set("TE")
    band_string.set("band")
    cone_string.set("n")
        
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
        polarization = polarization_string.get()
        band = band_string.get()
        cone = cone_string.get()        
        
        if d1 <= 0 or d2 < 0 or (d1+d2)/lambda_min- (d1+d2)/lambda_max > 1 or lambda_min < 400 or lambda_max > 2400 or lambda_min >= lambda_max:
            gui.input_error(initialize)
        else:
            if (film1 == "AlAs" or film1 == "AlGaAs (70% Al)" or film1 == "AlGaAs (31.5% Al)" or film2 == "AlAs" or film2 == "AlGaAs (70% Al)" or film2 == "AlGaAs (31.5% Al)" ) and lambda_min < 800:
                print("Oops! Losses not negigible in this wavelength range ... adjusting ...")
                lambda_min = 800
                lambda_max = 2400
                lambda_min_string.set(lambda_min)
                lambda_max_string.set(lambda_max)
            f.clf()
            Omega = np.linspace((d1+d2)/lambda_max, (d1+d2)/lambda_min, num=501, endpoint=True) # normalized frequency
            lambdav = (d1+d2)/Omega # vacuum wavelength in nm
            epsilon_f1 = np.real(epsilon(film1,lambdav))
            epsilon_f2 = np.real(epsilon(film2,lambdav))
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
polarization_string = Tk.StringVar()
band_string = Tk.StringVar()
cone_string = Tk.StringVar()

initialize()

row = 1

row = gui.create_title(mainframe,"unit cell parameters",row)
row = gui.create_entry(mainframe,u"film 1 thickness: d [nm] =",d1_string,row)
row = gui.create_radiobutton(mainframe,['film 1 medium:','Vacuum','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica'],film1_string,6,row)
row = gui.create_entry(mainframe,u"film 2 thickness: d [nm] =",d2_string,row)
row = gui.create_radiobutton(mainframe,['film 2 medium:','Vacuum','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica'],film2_string,6,row)
row = gui.create_label(mainframe,u"\u039b [nm] =",D_string,row)
row = gui.create_double_entry(mainframe,u"\u03bb [nm] >",lambda_min_string,u"\u03bb [nm] <",lambda_max_string,row)
row = gui.create_radiobutton(mainframe,['polarization:','TE','TM'],polarization_string,2,row)
row = gui.create_radiobutton(mainframe,['show:','band','gap'],band_string,2,row)
row = gui.create_checkbutton(mainframe,"mark area outside vacuum light cone",'n','y',cone_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)