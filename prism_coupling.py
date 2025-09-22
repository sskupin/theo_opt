import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import bpm_stuff as bpm
import film_stuff as film

gui.set_rcParams()
title = "Prism Coupling"
root = Tk.Tk()
root.title(title)

def initialize():
    LP_double.set(10000)
    tilt_double.set(0.1)
    gap_double.set(5)
    calculate()

def show_manual():
    gui.show_manual("taylor_series.png",title) 
    
def calculate():
    gui.change_cursor(root,"trek")
    LP = LP_double.get()
    tilt = tilt_double.get()
    gap = gap_double.get()

    N = 4096
    x, delta_x = np.linspace(-500,500,N,endpoint=False, retstep=True)
    
    epsilon_f = 2.25
    epsilon_s = 2.249
    d = 10
    n_eff = np.linspace(np.sqrt(epsilon_s), np.sqrt(epsilon_f), num=500, endpoint=False)
    kz = film.n_eff_d(epsilon_f,epsilon_s,epsilon_s,n_eff,d,0,'TE')
    mode,dummy1,dummy2 = film.mode_profile(epsilon_f,epsilon_s,epsilon_s,kz,d,x/d,'TE')
    mode = mode/np.sqrt(np.sum(np.abs(mode)**2))

    u0 = np.exp(-(x+100)**2/50**2+1j*x*tilt)
    norm = np.sum(np.abs(u0)**2)
    Nz = 10000
    L = 2.e4
    z, delta_z = np.linspace(0,L,Nz,endpoint=True, retstep=True)
    epsilon = np.ones([Nz,N])*epsilon_s
    
    epsilon[np.argwhere(z<=LP),np.where(x<=-d/2-gap)] = epsilon_f 
    
    epsilon[:,np.where((x>-d/2) & (x<d/2))] = epsilon_f
    
    u = bpm.propagation_wg(N,u0,epsilon,delta_x,Nz,delta_z,16)
    
    f.clf()
    
    a1 = plt.subplot2grid((1, 4), (0, 1), colspan=2)
    im1 = a1.imshow(np.abs(np.transpose(u[:,4*256:9*256]))**2 ,extent=[0, L, x[4*256], x[9*256]] , aspect='auto', origin='lower', cmap='jet') #L/(1.5*(x[3*1024]-x[7*256])), origin='upper', cmap='jet')
    a1.annotate(r'$|E_y|^2/|E_{0y}|^2_{\rm max}$', xy=(0.05*L,40),horizontalalignment='left', verticalalignment='bottom', color='w')
    a1.plot([0,L],[-d/2,-d/2],'w:')
    a1.plot([0,L],[d/2,d/2],'w:') 
    a1.plot([0,LP],[-d/2-gap,-d/2-gap],'w:') 
    a1.plot([LP,LP],[-d/2-gap,x[4*256]],'w:') 
    a1.set_xlabel(r'$z/\lambda$')
    a1.set_ylabel(r'$x/\lambda$')
    ca1 = a1.inset_axes([0.3*L, 50, 0.65*L, 10], transform=a1.transData)
    cb = plt.colorbar(im1, cax=ca1, orientation='horizontal')
    cb.ax.xaxis.set_tick_params(color='w', labelcolor='w')
    cb.outline.set_edgecolor('w')
    a1bis = a1.twiny()
    a1bis.set_xlim(a1.get_xlim())
    a1bis.set_xticks([0,L])
    a1bis.annotate(text='', xy=(0,x[1047]), xytext=(LP,x[1047]), arrowprops=dict(arrowstyle='<->', color='w'))
    if LP>0:
        a1bis.annotate(r'$P$', xy=(LP/2,x[1070]),horizontalalignment='center', verticalalignment='bottom', color='w')
    a1bis.set_xticklabels([r'0', r'L'])
    a1bis.tick_params(direction='out', pad=0)
        
    a2 = plt.subplot2grid((1, 4), (0, 0), sharey=a1)
    lns1 = a2.plot(epsilon[0,:],x,'r', label=r'$\varepsilon(x,z=0)$')
    a2.set_ylim([x[4*256],x[9*256]])
    a2.set_xlabel(r'$\varepsilon$')
    a2.set_xticks([epsilon_s,epsilon_f])
    a2.set_xticklabels([round(epsilon_s,3),round(epsilon_f,3)])
#    a2.invert_yaxis()
    a2bis = a2.twiny()
    lns2 = a2bis.plot(np.abs(u0),x, label=r'input beam',color='b')
    a2bis.set_xlabel(r'$|E_y|/|E_{0y}|_{\rm max}$')
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    a2.legend(lns, labs, loc=0)
    
    a3 = plt.subplot2grid((1, 4), (0, 3), sharey=a1)
    lns1 = a3.plot(epsilon[-1,:],x,'r', label=r'$\varepsilon(x,z=L)$')
    a3.set_ylim([x[4*256],x[9*256]])
    a3.set_xlabel(r'$\varepsilon$')
    a3.set_xticks([epsilon_s,epsilon_f])
    a3.set_xticklabels([round(epsilon_s,3),round(epsilon_f,3)])
#    a3.invert_yaxis()
    a3bis = a3.twiny()
    lns2 = a3bis.plot(np.abs(u[-1,:]),x, label=r'output beam',color='b')
    lns3 = a3bis.plot(np.abs(np.sum(u[-1,:]*mode))*mode,x, label='mode $\mu=0$',color='g')
    a3bis.set_xlabel(r'$|E_y|/|E_{0y}|_{\rm max}$')
    lns = lns1+lns2+lns3
    labs = [l.get_label() for l in lns]
    a3.legend(lns, labs, loc=0)
    
    a1.set_title(r'coupling efficiency '+str(round(100*np.abs(np.sum(u[-1,:]*mode))**2/norm,1))+r' \%')

    plt.tight_layout()
                
#    plt.savefig('prism_coupling.pdf',bbox_inches='tight',dpi=300, transparent=True)

    canvas.draw()       
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[10,5])

canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

LP_double = Tk.DoubleVar()
tilt_double = Tk.DoubleVar()
gap_double = Tk.DoubleVar()

initialize()

row = 1
row = gui.create_slider_with_latex(mainframe,r'Prism length $P/\lambda=$',LP_double,0,2.e4,row,increment=500)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r'Beam tilt $k_{x0}\lambda=$',tilt_double,0,0.2,row,increment=.005)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r'Gap width $d_{\mathrm{g}}/\lambda=$',gap_double,0,20,row,increment=.5)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)