import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
title = "Diffraction at a Finite Periodic Structure"
root = Tk.Tk()
root.title(title)

def initialize():
    kx0_double.set(0)
    a_double.set(2)
    b_double.set(6)
    N_int.set(4)
    log_N_F_double.set(np.log(.1))
    propagation_string.set("exact")
    calculate()
    
def show_manual():
    gui.show_manual("taylor_series.png",title)
    
def calculate():
    gui.change_cursor(root,"trek")
    a = a_double.get()
    b = b_double.get()*a
    kx0 = kx0_double.get()*2*np.pi/b
    N = N_int.get()
    A = ((N-1)*b+2*a)/2
    DA_string.set(round(2*A,2))
    N_F = np.exp(log_N_F_double.get())
    z_B = A**2/N_F
    z_B_string.set(round(z_B,2))
    propagation = propagation_string.get()
    
    Nx = 2**22
    Lx = np.max([a*Nx/200.0/N_F,1000])
    x, delta_x = np.linspace(-Lx/2,Lx/2,Nx,endpoint=False, retstep=True)
    v = np.where(np.abs(x+a-A)<a,np.exp(1j*kx0*x),0)
    for index in range(1,N):
        v = v + np.where(np.abs(x+a-A+index*b)<a,np.exp(1j*kx0*x),0)    
    f.clf()  
    a1 = f.add_subplot(311)    
    a1.plot(x,np.abs(v)**2,'k')
    a1.set_xlim([-np.max([2*A,10*a]),np.max([2*A,10*a])])
    a1.set_ylim([-0.05,1.15])
    a1.annotate(r'', xy=(-A,1.075), xytext=(-A-0.5*np.max([2*A,10*a]),1.075), arrowprops=dict(arrowstyle='->'))
    a1.annotate(r'', xy=(A,1.075), xytext=(A+0.5*np.max([2*A,10*a]),1.075), arrowprops=dict(arrowstyle='->'))
    a1.annotate(r'2A', xy=(0,1.06),horizontalalignment='center', verticalalignment='center')
    a1.set_xlabel(r'$x/\lambda$')
    a1.set_ylabel(r'$|u(z=0)|^2$ [norm. u.]')

    V = np.fft.fft(v)
    kx = 2*np.pi*np.fft.fftfreq(Nx,delta_x)
    k = 2*np.pi  
    a2 = f.add_subplot(312) 
    a2.plot(np.fft.fftshift(kx)/k,(np.abs(np.fft.fftshift(V))/np.max(np.abs(V)))**2,'b')
    a2.set_xlim([-1/a,1/a])
    a2.set_xlabel(r'$k_x/k$')
    a2.set_ylabel(r'$|U(z=0)|^2$ [norm. u.]')

    if propagation == 'paraxial':
        V = V*np.exp(-1j*kx**2/(2*k)*z_B)
    else:
        V = V*np.exp(1j*np.sqrt(k**2-kx**2+0j)*z_B)
    v = np.fft.ifft(V)
    a3 = f.add_subplot(313) 
    a3.plot(x/z_B,np.abs(v)**2,'b')
    a3.set_xlim([-np.max([1/a,a1.get_xlim()[1]/z_B]),np.max([1/a,a1.get_xlim()[1]/z_B])])
    a3.set_xlabel(r'$x/z_{\rm B}$')
    a3.set_ylabel(r'$|u(z=z_{\rm B})|^2$ [norm. u.]')
    
    plt.tight_layout()
    
#    plt.savefig('FPS.pdf',bbox_inches='tight',dpi=300, transparent=True)
    
    canvas.draw()       
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[6,6.75])

canvas = gui.create_canvas(root,f)
canvas.draw() # for faster feedback to user on startup
mainframe = gui.create_mainframe(root)

kx0_double = Tk.DoubleVar()
a_double = Tk.DoubleVar()
b_double = Tk.DoubleVar()
N_int = Tk.IntVar()
DA_string = Tk.StringVar()
log_N_F_double = Tk.DoubleVar()
z_B_string = Tk.StringVar()
propagation_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_slider_with_latex(mainframe,r"half slit width $a/\lambda =$",a_double,.25,4,row,increment=.25)
row = gui.create_slider_with_latex(mainframe,r"period length $b/a =$",b_double,4,10,row,increment=.25)
row = gui.create_intslider_with_latex(mainframe,r"Number of periods $N =$",N_int,1,7,row)
row = gui.create_label_with_latex(mainframe,r"Aperture size $2A/\lambda =$",DA_string,row)
row = gui.create_logslider_with_latex(mainframe,r"Fresnel Number $N_{\rm F} = $",log_N_F_double,0.1,1000,row)
row = gui.create_label_with_latex(mainframe,r"Distance to screen $z_{\rm B} / \lambda =$",z_B_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r"angle of incidence $b\sin(\varphi)/\lambda =$",kx0_double,-1,1,row,increment=.1)
row = gui.create_spacer(mainframe,row)
row = gui.create_radiobutton(mainframe,['vacuum propagation:','paraxial','exact'],propagation_string,2,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)