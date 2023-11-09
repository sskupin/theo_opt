import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("4f image of a finite periodic structure with Fourier mask (paraxial)")

def initialize():
    kx0_double.set(0)
    b_double.set(6)
    N_int.set(4)
    Rin_double.set(0.1)
    Rout_double.set(0.3)
    log_RA_double.set(np.log(0.6))
    calculate()
    
def calculate():
    gui.change_cursor(root,"trek")
    a = 1
    b = b_double.get()*a
    kx0 = kx0_double.get()*2*np.pi/b
    N = N_int.get()
    Rin = Rin_double.get()
    Rout = Rout_double.get()
    RA = np.exp(log_RA_double.get())
        
    Nx = 2**22
    AO = (N-1)*b+2*a
    Lx = a*Nx/10.01
    x, delta_x = np.linspace(-Lx/2,Lx/2,Nx,endpoint=False, retstep=True)
    v = np.where(np.abs(x+a-AO/2)<a,np.exp(1j*kx0*x),0)
    for index in range(1,N):
        v = v + np.where(np.abs(x+a-AO/2+index*b)<a,np.exp(1j*kx0*x),0)  
    f.clf()
    a1 = f.add_subplot(311)
    a1.plot(x/a,np.abs(v)**2,'k')
    a1.set_xlim([-np.max([AO/a,10]),np.max([AO/a,10])])
    a1.set_xlabel(r'$x/a$')
    a1.set_ylabel(r'$|u(z=0)|^2$ [norm. u.]')     
    
    V = np.fft.fft(v)
    k = 2*np.pi
    a2 = f.add_subplot(312)
    kx = 2*np.pi*np.fft.fftfreq(Nx,delta_x)
    V_mask = np.where(np.logical_or(np.logical_and(np.abs(kx/k*a)>=Rin,np.abs(kx/k*a)<=Rout),np.abs(kx/k*a)>=RA),0,V)
    a2.plot(np.fft.fftshift(kx)/k*a,(delta_x*np.abs(np.fft.fftshift(V)))**2/a**2,'b:')
    a2.plot(np.fft.fftshift(kx)/k*a,(delta_x*np.abs(np.fft.fftshift(V_mask)))**2/a**2,'b')
    a2.set_xlim([-1,1])
    a2.set_xlabel(r'$xa/(\lambda f)$')
    a2.set_ylabel(r'$|u_+(z=2f)|^2\lambda f/a^2$ [norm. u.]')
    a2.axvspan(RA, np.fft.fftshift(kx)[-1]/k, color='0.75')
    a2.axvspan(np.fft.fftshift(kx)[0]/k, -RA, color='0.75')
    if Rout>Rin:
        a2.axvspan(Rin, Rout, color='0.75')
        a2.axvspan(-Rout, -Rin, color='0.75')
    ylim = a2.get_ylim()
    if RA>0:
        a2.annotate(r'$\tilde{R}_{\rm A}$', xy=(RA-0.01,ylim[1]*0.9),horizontalalignment='left', verticalalignment='center')
        a2.plot([RA,RA],[ylim[1]*0.975,ylim[1]],'k')
    if Rout>Rin+0.05 and Rout<RA-0.11:
        a2.annotate(r'$\tilde{R}_{\rm out}$', xy=(Rout-0.01,ylim[1]*0.9),horizontalalignment='left', verticalalignment='center')
    if Rout>Rin and Rout<RA:
        a2.plot([Rout,Rout],[ylim[1]*0.975,ylim[1]],'k')
    if Rin<Rout and Rin<RA-0.05:
        a2.annotate(r'$\tilde{R}_{\rm in}$', xy=(Rin-0.01,ylim[1]*0.9),horizontalalignment='left', verticalalignment='center')
    if Rin<Rout and Rin<RA:
        a2.plot([Rin,Rin],[ylim[1]*0.975,ylim[1]],'k')
    a2.set_ylim(ylim)
    
    v = np.fft.ifft(V_mask)
    a3 = f.add_subplot(313)
    a3.plot(-x/a,np.abs(v)**2,'b')
    a3.set_xlim([-np.max([AO/a,10]),np.max([AO/a,10])])
    a3.set_xlabel(r'$x/a$')
    a3.set_ylabel(r'$|u(z=4f)|^2$ [norm. u.]')   
    
    plt.tight_layout()
    
#    plt.savefig('fps_4f_mask.pdf',bbox_inches='tight',dpi=300, transparent=True)
    
    canvas.draw()       
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[6,6.75])

canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

kx0_double = Tk.DoubleVar()
b_double = Tk.DoubleVar()
N_int = Tk.IntVar()
Rin_double = Tk.DoubleVar()
Rout_double = Tk.DoubleVar()
log_RA_double = Tk.DoubleVar()

initialize()

row = 1
row = gui.create_slider_with_latex(mainframe,r"period length $b/a =$",b_double,4,10,row)
row = gui.create_intslider_with_latex(mainframe,r"Number of periods $N =$",N_int,1,7,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r"angle of incidence $b \varphi /\lambda =$",kx0_double,-0.6,0.6,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r"inner radius of Fourier mask $\tilde{R}_{\rm in} = R_{\rm in}a/(\lambda f) =$",Rin_double,0,0.9,row)
row = gui.create_slider_with_latex(mainframe,r"outer radius of Fourier mask $\tilde{R}_{\rm out} = R_{\rm out}a/(\lambda f) =$",Rout_double,0,0.9,row)
row = gui.create_logslider_with_latex(mainframe,r"radius of Fourier aperture $\tilde{R}_{\rm A} = R_{\rm A}a/(\lambda f) =$",log_RA_double,0.01,10,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)