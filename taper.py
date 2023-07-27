import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import bpm_stuff as bpm
import film_stuff as film

gui.set_rcParams()
root = Tk.Tk()
root.title("Taper")

def initialize():
    TL_double.set(5000)
    
    calculate()
    
def calculate():
    gui.change_cursor(root,"trek")
    TL = TL_double.get()

    N = 4096
    x, delta_x = np.linspace(-500,500,N,endpoint=False, retstep=True)
    
    epsilon_f = 2.25
    epsilon_s = 2.249
    d1 = 10
    d2 = 40
    n_eff = np.linspace(np.sqrt(epsilon_s), np.sqrt(epsilon_f), num=500, endpoint=False)
    kz1 = film.n_eff_d(epsilon_f,epsilon_s,epsilon_s,n_eff,d1,0,'TE')
    mode1,dummy1,dummy2 = film.mode_profile(epsilon_f,epsilon_s,epsilon_s,kz1,d1,x/d1,'TE')
    mode1 = mode1/np.sqrt(np.sum(np.abs(mode1)**2))
    kz20 = film.n_eff_d(epsilon_f,epsilon_s,epsilon_s,n_eff,d2,0,'TE')
    mode20,dummy1,dummy2 = film.mode_profile(epsilon_f,epsilon_s,epsilon_s,kz20,d2,x/d2,'TE')
    mode20 = mode20/np.sqrt(np.sum(np.abs(mode20)**2))
    kz22 = film.n_eff_d(epsilon_f,epsilon_s,epsilon_s,n_eff,d2,2,'TE')
    mode22,dummy1,dummy2 = film.mode_profile(epsilon_f,epsilon_s,epsilon_s,kz22,d2,x/d2,'TE')
    mode22 = mode22/np.sqrt(np.sum(np.abs(mode22)**2))

    u0 = mode1
    u0 = u0/np.sqrt(np.sum(np.abs(u0)**2))
    Nz = 10000
    L = 2.e4
    z, delta_z = np.linspace(0,L,Nz,endpoint=True, retstep=True)
    epsilon = np.ones([Nz,N])*epsilon_s
    
    epsilon[:,np.where((x>-d1/2) & (x<d1/2))] = epsilon_f
    
    epsilon[np.argwhere(z>=(L+TL)/2),np.where((x>-d2/2) & (x<d2/2))] = epsilon_f
        
    for index in np.argwhere((z>(L-TL)/2) & (z<(L+TL)/2)):
        d = (z[index]-(L-TL)/2)/TL*d2 - (z[index]-(L+TL)/2)/TL*d1
        epsilon[index,np.where((x>-d/2) & (x<d/2))] = epsilon_f
    
    u = bpm.propagation_wg(N,u0,epsilon,delta_x,Nz,delta_z,16)
    
    f.clf()
    
    a1 = plt.subplot2grid((1, 4), (0, 1), colspan=2)
    a1.imshow(np.abs(np.transpose(u[:,7*256:9*256]))**2 ,extent=[0, L, x[7*256], x[9*256]] , aspect='auto', origin='upper', cmap='jet') #L/(1.5*(x[3*1024]-x[7*256])), origin='upper', cmap='jet')
    a1.plot([0,(L-TL)/2],[-d1/2,-d1/2],'w:')
    a1.plot([0,(L-TL)/2],[d1/2,d1/2],'w:') 
    a1.plot([(L+TL)/2,L],[-d2/2,-d2/2],'w:')
    a1.plot([(L+TL)/2,L],[d2/2,d2/2],'w:') 
    a1.plot([(L-TL)/2,(L+TL)/2],[-d1/2,-d2/2],'w:') 
    a1.plot([(L-TL)/2,(L+TL)/2],[d1/2,d2/2],'w:') 
    a1.set_xlabel(r'$z/\lambda$')
    a1.set_ylabel(r'$x/\lambda$')
    a1bis = a1.twiny()
    a1bis.set_xlim(a1.get_xlim())
    a1bis.annotate(text='', xy=((L-TL)/2,x[7*266]), xytext=((L+TL)/2,x[7*266]), arrowprops=dict(arrowstyle='<->', color='w'))
    if TL>0:
        a1bis.annotate(r'T', xy=(L/2,x[7*263]),horizontalalignment='center', verticalalignment='bottom', color='w')
    a1bis.set_xticks([(L-TL)/2,(L+TL)/2])
    a1bis.set_xticklabels([r'', r''])
    a1bis.tick_params(direction='out', pad=0)
        
    a2 = plt.subplot2grid((1, 4), (0, 0), sharey=a1)
    lns1 = a2.plot(epsilon[0,:],x,'r', label=r'$\varepsilon(x,z=0)$')
    a2.set_ylim([x[7*256],x[9*256]])
    a2.set_xlabel(r'$\varepsilon$')
    a2.set_xticks([epsilon_s,epsilon_f])
    a2.set_xticklabels([round(epsilon_s,3),round(epsilon_f,3)])
    a2.invert_yaxis()
    a2bis = a2.twiny()
    lns2 = a2bis.plot(np.abs(u0),x, label=r'input beam',color='b')
    lns3 = a2bis.plot(np.abs(mode1),x,'g--', label='mode $\mu=0$')
    a2bis.set_xlabel(r'$|E_y|$ [arb.u.]')
    lns = lns1+lns2+lns3
    labs = [l.get_label() for l in lns]
    a2.legend(lns, labs, loc=0)
    
    a3 = plt.subplot2grid((1, 4), (0, 3), sharey=a1)
    lns1 = a3.plot(epsilon[-1,:],x,'r', label=r'$\varepsilon(x,z=L)$')
    a3.set_ylim([x[7*256],x[9*256]])
    a3.set_xlabel(r'$\varepsilon$')
    a3.set_xticks([epsilon_s,epsilon_f])
    a3.set_xticklabels([round(epsilon_s,3),round(epsilon_f,3)])
    a3.invert_yaxis()
    a3bis = a3.twiny()
    lns2 = a3bis.plot(np.abs(u[-1,:]),x, label=r'output beam',color='b')
    lns3 = a3bis.plot(np.abs(np.sum(u[-1,:]*mode20))*mode20,x,'g--', label='mode $\mu=0$')
    lns4 = a3bis.plot(np.abs(np.sum(u[-1,:]*mode22))*np.abs(mode22),x,'y--', label='mode $\mu=2$')
    a3bis.set_xlabel(r'$|E_y|$ [arb.u.]')
    lns = lns1+lns2+lns3+lns4
    labs = [l.get_label() for l in lns]
    a3.legend(lns, labs, loc=0)
    
    a1.set_title(r'coupling to mode $\mu=0$ : '+str(round(100*np.abs(np.sum(u[-1,:]*mode20))**2,3))+r' \%; mode $\mu=2$ : '+str(round(100*np.abs(np.sum(u[-1,:]*mode22))**2,3))+r' \%')

    plt.tight_layout()
                
#    plt.savefig('taper.pdf',bbox_inches='tight',dpi=300, transparent=True)

    canvas.draw() 
    gui.change_cursor(root,"arrow")      

f = plt.figure(1,[10,5])

canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

TL_double = Tk.DoubleVar()

initialize()

row = 1
row = gui.create_slider_with_latex(mainframe,r"Taper length $T/\lambda =$",TL_double,0,1.e4,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)