import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import sys
sys.path.append('./aux')
import gui_stuff as gui
import strat_stuff as strat

gui.set_rcParams()
title = "Reflection and Transmission at Graded-Index Films"
root = Tk.Tk()
root.title(title)

def epsilon_graded(profile,epsilon_l,epsilon_r,dL,z):
    if dL==0:
        epsgrad = epsilon_l
    elif profile == 'linear':
        epsgrad = epsilon_l + (epsilon_r-epsilon_l)*z/dL
    else:
        epsgrad = epsilon_l + (epsilon_r-epsilon_l)*np.exp(-6*(z-dL)**2/dL**2)
    return epsgrad

def reflection_transmission(epsilon_s,d,epsilon_f,epsilon_c,phi): # computing coefficients of reflection and transmission
    kx,ksz,kcz = strat.KSC(epsilon_s,epsilon_c,phi)
    kfz = np.sqrt(epsilon_f*(2*np.pi)**2+1j*0.-kx**2)
    MTE = strat.mTE(kfz[0],d[0])
    MTM = strat.mTM(kfz[0],epsilon_f[0],d[0])
    for index in range(1,len(d)):
        MTE = np.matmul(strat.mTE(kfz[index],d[index]),MTE)
        MTM = np.matmul(strat.mTM(kfz[index],epsilon_f[index],d[index]),MTM)
    RTE,RTM,TTE,TTM,tauTE,tauTM = strat.RTAU(ksz,kcz,epsilon_s,epsilon_c,MTE,MTM)
    return RTE,RTM,tauTE,tauTM
        
def initialize():
    var_string[0].set("1") # epsilon_s
    var_string[1].set("0") # overlay thickness in units of lambda
    var_string[2].set("1") # epsilon_f (overlay)
    var_string[3].set("4.1") # Re epsilon_c
    var_string[4].set("0") # Im epsilon_c
    var_string[5].set("2.3") # GRIN film thickness in units of lambda
    var_string[6].set("1") # epsilon left
    var_string[7].set("4.1") # espilon right
    var_string[8].set("Gaussian") # profile
    var_string[9].set("40") # number of layers
    var_string[10].set("TE")
    var_string[11].set("0.375")
    gui.copy_stringvar_vector(var_string,var_save)        
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()  

def show_manual():
    gui.show_manual("man/rt_graded_index.png",title)

def calculate():
    gui.change_cursor(root,"trek")
    try:
        epsilon_s = float(var_string[0].get())
        df = float(var_string[1].get())
        epsilon_f = float(var_string[2].get())
        epsilon_c_real = float(var_string[3].get())
        epsilon_c_imag = float(var_string[4].get())
        dL = float(var_string[5].get())
        epsilon_l = float(var_string[6].get())
        epsilon_r = float(var_string[7].get())
        profile = var_string[8].get()
        NL = int(var_string[9].get())
        polarization = var_string[10].get()
        phi_0 = float(var_string[11].get())*np.pi
        
        if epsilon_s <= 0 or (epsilon_c_imag == 0 and epsilon_c_real == 0) or dL < 0 or df < 0 or phi_0 < 0 or phi_0 >= np.pi/2:
            gui.input_error("Values out of range. Re-initializing ...", reinitialize)
        elif NL <= 0 or NL > 100:
            gui.input_error("Number of layers must be between 1 and 1000. Re-initializing ...", reinitialize)
        else:
            f.clf()
            phi = np.linspace(0, np.pi/2, num=1001, endpoint=False) # angle of incidence
            epsilon_c = epsilon_c_real + 1j*epsilon_c_imag
            vepsilon_graded = np.vectorize(epsilon_graded)
            if df > 0:
                d = np.ones(NL+1)*dL/NL
                d[0] = df
                zf = np.cumsum(d)
                epsilon_fL = np.ones(NL+1)
                epsilon_fL[0] = epsilon_f
                epsilon_fL[1:] = vepsilon_graded(profile,epsilon_l,epsilon_r,dL,zf[1:]-d[1:]/2-df)
            else:
                d = np.ones(NL)*dL/NL
                zf = np.cumsum(d)
                epsilon_fL = vepsilon_graded(profile,epsilon_l,epsilon_r,dL,zf-d/2)
            vreflection_transmission = np.vectorize(reflection_transmission,excluded=[1,2])
            RTE,RTM,tauTE,tauTM = vreflection_transmission(epsilon_s,d,epsilon_fL,epsilon_c,phi)
            a1 = f.add_subplot(221)
            a2 = f.add_subplot(222)
            strat.plot_curves_vs_angle(a1,phi,[np.abs(RTE)**2,tauTE],[r'$\rho_{\rm TE}$',r'$\tau_{\rm TE}$'],['b','r'],phi[0],phi[-1])
            strat.plot_curves_vs_angle(a2,phi,[np.abs(RTM)**2,tauTM],[r'$\rho_{\rm TM}$',r'$\tau_{\rm TM}$'],['b','r'],phi[0],phi[-1])
            if polarization == 'TE':
                a1.plot(phi_0,np.interp(phi_0,phi,np.abs(RTE)**2),'bo')
                a1.plot(phi_0,np.interp(phi_0,phi,tauTE),'ro')
            else:
                a2.plot(phi_0,np.interp(phi_0,phi,np.abs(RTM)**2),'bo')
                a2.plot(phi_0,np.interp(phi_0,phi,tauTM),'ro')
            left, right = a1.get_ylim()
            a1.set_ylim([min(left,-0.025),max(right,1.025)])
            left, right = a2.get_ylim()
            a2.set_ylim([min(left,-0.025),max(right,1.025)])
            
            a3 = f.add_subplot(212)
            zs = np.linspace(-np.maximum(0.75*zf[-1]/2,1), 0, num=501, endpoint=True)
            zc = np.linspace(0, np.maximum(0.75*zf[-1]/2,1), num=501, endpoint=True)
            for index in range(len(d)):
                if np.ceil((index)/2) == index/2:
                    color = '0.875'
                else:
                    color='0.75'
                a3.axvspan(zf[index]-d[index], zf[index], color=color)
            a3.set_xlabel(r'$z/\lambda$')
            a3.set_xlim([zs[0],zf[-1]+zc[-1]])
            RTE,RTM,tauTE,tauTM = reflection_transmission(epsilon_s,d,epsilon_fL,epsilon_c,phi_0)
            if polarization == 'TE':
                F = 1+RTE
                G = 1j*2*np.pi*np.sqrt(epsilon_s)*np.cos(phi_0)*(1-RTE)
                MTE = strat.mTE(2*np.pi*np.sqrt(epsilon_s)*np.cos(phi_0),zs)
                a3.plot(zs,np.abs(MTE[0,0]*F+MTE[0,1]*G),'b')
                if zf[-1]>0:
                    for index in range(len(d)):
                        zfi = np.linspace(0, d[index], num=max(1,int(np.ceil(d[index]/zf[-1]*1001))), endpoint=True)
                        MTE = strat.mTE(2*np.pi*np.sqrt(epsilon_fL[index]+1j*0-epsilon_s*np.sin(phi_0)**2),zfi)
                        F,G = strat.plot_amplitude(a3,MTE,F,G,zf[index]-d[index]+zfi)
                lns1 = a3.plot(zf[-1]+zc,np.abs(F)*np.exp(-np.imag(np.sqrt(epsilon_c-epsilon_s*np.sin(phi_0)**2))*2*np.pi*zc),'b',label=r'$|E|$')
                a3.set_ylabel(r'$| E_\mathrm{TE} | / | E_\mathrm{TE}^\mathrm{i} |$ for $\varphi_\mathrm{i}=$ '+str(round(phi_0/np.pi,4))+r'$\pi$')
            else:
                F = 1-RTM
                G = 1j*2*np.pi*np.sqrt(epsilon_s)*np.cos(phi_0)*(1+RTM)/epsilon_s
                MTM = strat.mTM(2*np.pi*np.sqrt(epsilon_s)*np.cos(phi_0),epsilon_s,zs)
                a3.plot(zs,np.abs(MTM[0,0]*F+MTM[0,1]*G),'b')
                if zf[-1]>0:
                    for index in range(len(d)):
                        zfi = np.linspace(0, d[index], num=max(1,int(np.ceil(d[index]/zf[-1]*1001))), endpoint=True)
                        MTM = strat.mTM(2*np.pi*np.sqrt(epsilon_fL[index]+1j*0-epsilon_s*np.sin(phi_0)**2),epsilon_fL[index],zfi)
                        F,G = strat.plot_amplitude(a3,MTM,F,G,zf[index]-d[index]+zfi)
                lns1 = a3.plot(zf[-1]+zc,np.abs(F)*np.exp(-np.imag(np.sqrt(epsilon_c-epsilon_s*np.sin(phi_0)**2))*2*np.pi*zc),'b',label=r'$|H|$')                
                a3.set_ylabel(r'$| H_\mathrm{TM} | / | H_\mathrm{TM}^\mathrm{i} |$ for $\varphi_\mathrm{i}=$ '+str(round(phi_0/np.pi,4))+r'$\pi$')
            ylimits = a3.get_ylim()
            a3.set_ylim([ylimits[0],ylimits[1]+0.1*(ylimits[1]-ylimits[0])])
            a3bis = a3.twinx()
            zfi = np.linspace(df, zf[-1], num=1001, endpoint=True)
            a3bis.plot([zs[0],zs[-1]],[epsilon_s,epsilon_s],'r',[0,df],[epsilon_f,epsilon_f],'r',[zf[-1]+zc[0],zf[-1]+zc[-1]],[np.real(epsilon_c),np.real(epsilon_c)],'r',label=r'$\epsilon$')
            lns2 = a3bis.plot(zfi,vepsilon_graded(profile,epsilon_l,epsilon_r,dL,zfi-df),'r',label=r'$\epsilon$')
            a3bis.set_ylabel(r'$\epsilon$') 
            if polarization == 'TE':
                a3.legend(lns1+lns2,[r'$|E|$',r'$\epsilon$'])
            else:
                a3.legend(lns1+lns2,[r'$|H|$',r'$\epsilon$'])
            
            plt.tight_layout()  
            
#            plt.savefig('graded_index.pdf',bbox_inches='tight',dpi=300, transparent=True)

            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[8,5])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(12)
var_save = gui.create_stringvar_vector(12)

initialize()

row = 1
row = gui.create_title(mainframe,"environment parameters",row)
row = gui.create_entry(mainframe,u"substrate: \u03B5 =",var_string[0],row)
row = gui.create_entry(mainframe,u"overlay thickness: d/\u03BB =",var_string[1],row)
row = gui.create_entry(mainframe,u"overlay: \u03B5 =",var_string[2],row)
row = gui.create_entry(mainframe,u"cladding: \u03B5' =",var_string[3],row)
row = gui.create_entry(mainframe,u"cladding: \u03B5'' =",var_string[4],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_title(mainframe,"graded index parameters",row)
row = gui.create_entry(mainframe,u"film thickness: d/\u03BB =",var_string[5],row)
row = gui.create_entry(mainframe,u"lhs: \u03B5 =",var_string[6],row)
row = gui.create_entry(mainframe,u"rhs: \u03B5 =",var_string[7],row)
row = gui.create_radiobutton(mainframe,['profile:','linear','Gaussian'],var_string[8],2,row)# add quintic from ARC
row = gui.create_entry(mainframe,u"number of layers =",var_string[9],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_title(mainframe,"field parameters",row)
row = gui.create_radiobutton(mainframe,['polarization:','TE','TM'],var_string[10],2,row)
row = gui.create_entry_with_latex(mainframe,r'angle of incidence: $\varphi_\mathrm{i}/\pi=$',var_string[11],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)