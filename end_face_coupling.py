import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import bpm_stuff as bpm
import film_stuff as film

gui.set_rcParams()
root = Tk.Tk()
root.title("End-face coupling")

def initialize():
    w0_double.set(15)
    x0_double.set(5)
    d_string.set(25)
    z0_double.set(-500)
    epsilon_f_string.set(2.250)
    epsilon_s_string.set(2.249)
    epsilon_c_string.set(2.249)
    
    calculate()
    
def calculate():
    gui.change_cursor(root,"trek")
    w0 = w0_double.get()
    x0 = x0_double.get()
    d = float(d_string.get())
    epsilon_f = float(epsilon_f_string.get())
    epsilon_s = float(epsilon_s_string.get())
    epsilon_c = float(epsilon_c_string.get())

    N = 4096*2
    x, delta_x = np.linspace(-1000,1000,N,endpoint=False, retstep=True)

    n_eff = np.linspace(np.sqrt(epsilon_s), np.sqrt(epsilon_f), num=500, endpoint=False)
    mu_max = film.number_of_modes(epsilon_f,epsilon_s,epsilon_s,d,'TE')
    kz = np.zeros(mu_max)
    modes = np.zeros([mu_max,N])
    for index in range(mu_max):
        kz[index] = film.n_eff_d(epsilon_f,epsilon_s,epsilon_c,n_eff,d,index,'TE')[0]
        modes[index,:],dummy1,dummy2 = film.mode_profile(epsilon_f,epsilon_s,epsilon_c,kz[index],d,x/d,'TE')
        modes[index,:] = modes[index,:]/np.sqrt(np.sum(np.abs(modes[index,:])**2))
        kz[index] = 2*np.pi*kz[index] 
        
    uf = np.exp(-(x-x0)**2/w0**2)
    uf = uf/np.sqrt(np.sum(np.abs(uf)**2))
    u = bpm.propagation_wg(N,uf,np.ones([2,N])*epsilon_s,delta_x,2,-z0_double.get(),16)
    u0 = u[-1,:]
    z0 = -2.e3
    u = bpm.propagation_wg(N,u0,np.ones([2,N])*epsilon_s,delta_x,2,z0,16)
    u1 = u[-1,:]
    Nz = 5000
    L = 8.e3
    z, delta_z = np.linspace(z0,L,Nz,endpoint=True, retstep=True)  
    epsilon = np.ones([Nz,N])*epsilon_s 
    epsilon[:,np.where(x>-d/2)] = epsilon_f  
    epsilon[:,np.where(x>=d/2)] = epsilon_c
    epsilon[np.where(z<=0),:] = epsilon_s
    u = bpm.propagation_wg(N,u1,epsilon,delta_x,Nz,delta_z,16)
    
    f.clf()
    
    a1 = f.add_subplot(211)
    a1.imshow(np.abs(np.transpose(u[:,15*256:17*256]))**2 ,extent=[z0, L, x[17*256], x[15*256]] , aspect=L/(3*(x[17*256]-x[15*256])), origin='upper', cmap='jet')
    a1.plot([0,L],[-d/2,-d/2],'w:')
    a1.plot([0,L],[d/2,d/2],'w:')  
    a1.plot([0,0],[x[15*256],x[17*256]],'w:')  
    a1.set_xlabel(r'$z/\lambda$')
    a1.set_ylabel(r'$x/\lambda$')
        
    a3 = f.add_subplot(223)
    a3.plot(x,np.abs(u0),'b', label=r'input beam')
    a3.plot(x,np.abs(np.sum(u0*modes[0,:])*modes[0,:]),'r', label='mode $\mu=0$')
    a3.set_xlim([x[15*256], x[17*256]])
    a3.set_xlabel(r'$x/\lambda$')
    a3.set_ylabel(r'$|E_y|$ [arb.u.]') 
    a3.set_title(r'conversion '+str(round(100*np.abs(np.sum(u0*modes[0,:]))**2,3))+r' \%')
    a3.legend()
    ylim = a3.get_ylim()
    a3.plot([-d/2,-d/2],ylim,'k:')
    a3.plot([d/2,d/2],ylim,'k:')
    a3.set_ylim(ylim)

    a4 = f.add_subplot(224)
    a4.plot(x,np.abs(u0),'b', label=r'input beam')
    a4.plot(x,np.abs(np.sum(u0*modes[1,:])*modes[1,:]),'r', label='mode $\mu=1$')
    a4.set_xlim([x[15*256], x[17*256]])
    a4.set_xlabel(r'$x/\lambda$')
    a4.set_ylabel(r'$|E_y|$ [arb.u.]') 
    a4.set_title(r'conversion '+str(round(100*np.abs(np.sum(u0*modes[1,:]))**2,3))+r' \%')
    a4.legend()
    ylim = a4.get_ylim()
    a4.plot([-d/2,-d/2],ylim,'k:')
    a4.plot([d/2,d/2],ylim,'k:')
    a4.set_ylim(ylim)

    plt.tight_layout()
                
#    plt.savefig('end_face_coupling.pdf',bbox_inches='tight',dpi=300, transparent=True)

    canvas.draw()       
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[7,5])

canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

w0_double = Tk.DoubleVar()
x0_double = Tk.DoubleVar()
d_string = Tk.StringVar()
z0_double = Tk.DoubleVar()
epsilon_f_string = Tk.StringVar()
epsilon_s_string = Tk.StringVar()
epsilon_c_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_slider_with_latex(mainframe,r'beam width $[\lambda]~=$',w0_double,10,25,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r'transversal focus position $[\lambda]~=$',x0_double,-25,25,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r'longitudinal focus position $[\lambda]~=$',z0_double,-1000,1000,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_label_with_latex(mainframe,r'$\varepsilon^{\rm left}~=$',epsilon_s_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_label_with_latex(mainframe,r'film: $\varepsilon^{\rm right}_{\rm f}~=$',epsilon_f_string,row)
row = gui.create_label_with_latex(mainframe,r'substrate: $\varepsilon^{\rm right}_{\rm s}~=$',epsilon_s_string,row)
row = gui.create_label_with_latex(mainframe,r'cladding: $\varepsilon^{\rm right}_{\rm c}~=$',epsilon_c_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_label_with_latex(mainframe,r'film thickness: $d/\lambda~=$',d_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)