import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
title = "Reflection and Transmission of a Graded Index Structure"
root = Tk.Tk()
root.title(title)

def epsilon_graded(profile,epsilon_l,epsilon_r,dL,z):
    if profile == 'linear':
        epsgrad = epsilon_l + (epsilon_r-epsilon_l)*z/dL
    else:
        epsgrad = epsilon_l + (epsilon_r-epsilon_l)*np.exp(-6*(z-dL)**2/dL**2)
    return epsgrad

def mTE(kfz,z):
    return np.array([[np.cos(kfz*2*np.pi*z),np.sin(kfz*2*np.pi*z)/kfz],[-np.sin(kfz*2*np.pi*z)*kfz,np.cos(kfz*2*np.pi*z)]])

def mTM(kfz,epsilon_f,z):
    return np.array([[np.cos(kfz*2*np.pi*z),np.sin(kfz*2*np.pi*z)*epsilon_f/kfz],[-np.sin(kfz*2*np.pi*z)*kfz/epsilon_f,np.cos(kfz*2*np.pi*z)]])

def reflection_transmission(epsilon_s,d,epsilon_f,epsilon_c,phi): # computing coefficients of reflection and transmission
    ksx = np.sqrt(np.real(epsilon_s))*np.sin(phi)
    ksz = np.sqrt(epsilon_s-ksx**2)
    kfz = np.sqrt(epsilon_f+1j*0.-ksx**2)
    kcz = np.sqrt(epsilon_c-ksx**2)
    MTE = mTE(kfz[0],d[0])
    MTM = mTM(kfz[0],epsilon_f[0],d[0])
    for index in range(1,len(d)):
        MTE = np.matmul(mTE(kfz[index],d[index]),MTE)
        MTM = np.matmul(mTM(kfz[index],epsilon_f[index],d[index]),MTM)
    NTE = ksz*MTE[1,1]+kcz*MTE[0,0]+1j*MTE[1,0]-1j*ksz*kcz*MTE[0,1]
    RTE = (ksz*MTE[1,1]-kcz*MTE[0,0]-1j*MTE[1,0]-1j*ksz*kcz*MTE[0,1])/NTE
    NTM = ksz*MTM[1,1]/epsilon_s+kcz*MTM[0,0]/epsilon_c+1j*MTM[1,0]-1j*ksz*kcz*MTM[0,1]/(epsilon_s*epsilon_c)
    RTM = -(ksz*MTM[1,1]/epsilon_s-kcz*MTM[0,0]/epsilon_c-1j*MTM[1,0]-1j*ksz*kcz*MTM[0,1]/(epsilon_s*epsilon_c))/NTM # for electric field (negative of magnetic coeff.)
    tauTE = np.real(kcz)/ksz*np.abs(2*ksz/NTE)**2
    tauTM = np.real(kcz/epsilon_c)*epsilon_s/ksz*(np.abs(2*ksz/epsilon_s/NTM))**2
    
    return RTE,RTM,tauTE,tauTM

def plot_subplot(ax,phi,curves,labels,colors):
    for index in range(len(labels)):
        ax.plot(phi,curves[index],colors[index],label=labels[index])
    ax.set_xticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2])
    ax.set_xticklabels([r'$0$', r'$\pi/8$', r'$\pi/4$', r'$3\pi/8$', r'$\pi/2$'])
    ax.set_xlabel(r'$\varphi_{\rm i}$')
    ax.set_xlim([0,np.pi/2])
    ax.set_ylabel(','.join(labels))
    ax.legend()

def plot_amplitude(ax,M,F,G,z):
    Fplot = M[0,0]*F+M[0,1]*G
    Gplot = M[1,0]*F+M[1,1]*G
    FGabs = np.abs(Fplot)+np.abs(Gplot)
    if np.any(np.less(FGabs,1.e-6)):
        Fplot[np.where(FGabs<1.e-6)[0][0]:]=0
        Gplot[np.where(FGabs<1.e-6)[0][0]:]=0
    F=Fplot[-1]
    G=Gplot[-1]
    ax.plot(z,np.abs(Fplot),'b')
    
    return F,G
        
def initialize():
    dL_string.set("2.3")
    df_string.set("0")
    epsilon_s_string.set("1")
    epsilon_f_string.set("1")
    epsilon_c_real_string.set("4.1")
    epsilon_c_imag_string.set("0")
    epsilon_l_string.set("1")
    epsilon_r_string.set("4.1")
    NL_string.set("40")
    profile_string.set("linear")
    polarization_string.set("TE")
    phi_0_string.set("0")
        
    calculate()

def calculate():
    try:
        epsilon_s = float(epsilon_s_string.get())
        epsilon_f = float(epsilon_f_string.get())
        epsilon_l = float(epsilon_l_string.get())
        epsilon_r = float(epsilon_r_string.get())
        epsilon_c_real = float(epsilon_c_real_string.get())
        epsilon_c_imag = float(epsilon_c_imag_string.get())
        df = float(df_string.get())
        dL = float(dL_string.get())
        NL = int(NL_string.get())
        profile = profile_string.get()
        polarization = polarization_string.get()
        phi_0 = float(phi_0_string.get())*np.pi
        
        if epsilon_s <= 0 or epsilon_c_imag < 0 or epsilon_c_real == 0 or dL <= 0 or df < 0 or NL <= 0 or NL >= 1000 or phi_0 < 0 or phi_0 >= np.pi/2:
            gui.input_error(initialize)
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
                epsilon_fL[1:] = vepsilon_graded(profile,epsilon_l,epsilon_r,dL,zf[1:]-d[1:]/2)
            else:
                d = np.ones(NL)*dL/NL
                zf = np.cumsum(d)
                epsilon_fL = vepsilon_graded(profile,epsilon_l,epsilon_r,dL,zf-d/2-df)
            vreflection_transmission = np.vectorize(reflection_transmission,excluded=[1,2])
            RTE,RTM,tauTE,tauTM = vreflection_transmission(epsilon_s,d,epsilon_fL,epsilon_c,phi)
            a1 = f.add_subplot(221)
            a2 = f.add_subplot(222)
            plot_subplot(a1,phi,[np.abs(RTE)**2,tauTE],[r'$\rho_{\rm TE}$',r'$\tau_{\rm TE}$'],['b','r'])
            plot_subplot(a2,phi,[np.abs(RTM)**2,tauTM],[r'$\rho_{\rm TM}$',r'$\tau_{\rm TM}$'],['b','r'])
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
                G = 1j*np.sqrt(epsilon_s)*np.cos(phi_0)*(1-RTE)
                MTE = mTE(np.sqrt(epsilon_s)*np.cos(phi_0),zs)
                a3.plot(zs,np.abs(MTE[0,0]*F+MTE[0,1]*G),'b')
                for index in range(len(d)):
                    zfi = np.linspace(0, d[index], num=max(1,int(np.ceil(d[index]/zf[-1]*1001))), endpoint=True)
                    MTE = mTE(np.sqrt(epsilon_fL[index]+1j*0-epsilon_s*np.sin(phi_0)**2),zfi)
                    F,G = plot_amplitude(a3,MTE,F,G,zf[index]-d[index]+zfi)
                lns1 = a3.plot(zf[-1]+zc,np.abs(F)*np.exp(-np.imag(np.sqrt(epsilon_c-epsilon_s*np.sin(phi_0)**2))*2*np.pi*zc),'b',label=r'$|E|$')
                a3.set_ylabel(r'$| E_\mathrm{TE} | / | E_\mathrm{TE}^\mathrm{i} |$ for $\varphi_\mathrm{i}=$ '+str(round(phi_0/np.pi,4))+r'$\pi$')
            else:
                F = 1-RTM
                G = 1j*np.sqrt(epsilon_s)*np.cos(phi_0)*(1+RTM)/epsilon_s
                MTM = mTM(np.sqrt(epsilon_s)*np.cos(phi_0),epsilon_s,zs)
                a3.plot(zs,np.abs(MTM[0,0]*F+MTM[0,1]*G),'b')
                for index in range(len(d)):
                    zfi = np.linspace(0, d[index], num=max(1,int(np.ceil(d[index]/zf[-1]*1001))), endpoint=True)
                    MTM = mTM(np.sqrt(epsilon_fL[index]+1j*0-epsilon_s*np.sin(phi_0)**2),epsilon_fL[index],zfi)
                    F,G = plot_amplitude(a3,MTM,F,G,zf[index]-d[index]+zfi)
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

            canvas.draw()
    except ValueError: gui.input_error(initialize)

f = plt.figure(1,[8,5])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

epsilon_s_string = Tk.StringVar()
epsilon_f_string = Tk.StringVar()
epsilon_l_string = Tk.StringVar()
epsilon_r_string = Tk.StringVar()
epsilon_c_real_string = Tk.StringVar()
epsilon_c_imag_string = Tk.StringVar()
df_string = Tk.StringVar()
dL_string = Tk.StringVar()
NL_string = Tk.StringVar()
profile_string = Tk.StringVar()
polarization_string = Tk.StringVar()
phi_0_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_title(mainframe,"stack parameters",row)
row = gui.create_entry(mainframe,u"substrate: \u03B5 =",epsilon_s_string,row)
row = gui.create_entry(mainframe,u"film thickness: d/\u03BB =",df_string,row)
row = gui.create_entry(mainframe,u"film: \u03B5 =",epsilon_f_string,row)
row = gui.create_entry(mainframe,u"cladding: \u03B5' =",epsilon_c_real_string,row)
row = gui.create_entry(mainframe,u"cladding: \u03B5'' =",epsilon_c_imag_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_title(mainframe,"graded index parameters",row)
row = gui.create_entry(mainframe,u"domain thickness: L/\u03BB =",dL_string,row)
row = gui.create_entry(mainframe,u"lhs: \u03B5 =",epsilon_l_string,row)
row = gui.create_entry(mainframe,u"rhs: \u03B5 =",epsilon_r_string,row)
row = gui.create_radiobutton(mainframe,['profile:','linear','Gaussian'],profile_string,2,row)# add quintic from ARC
row = gui.create_entry(mainframe,u"Number of layers =",NL_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_title(mainframe,"field parameters",row)
row = gui.create_radiobutton(mainframe,['polarization:','TE','TM'],polarization_string,2,row)
row = gui.create_entry_with_latex(mainframe,r'angle of incidence: $\varphi_\mathrm{i}/\pi=$',phi_0_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)