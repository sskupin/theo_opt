import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import tkinter as Tk
import gui_stuff as gui

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rc('text', usetex=True)
mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')
mpl.rcParams.update({'font.size': 10})

root = Tk.Tk()
root.title("A finite periodic structure as a spectrometer")

def initialize():
    kx0_double.set(.5)
    log_delta_double.set(np.log(0.01))
    a_double.set(.5)
    N_int.set(70)
    propagation_string.set("exact")
    log_window_double.set(np.log(.2))
    
    calculate()
    
def G(x,N):
    if np.sin(x) == 0:
        G = N
    else:
        G = np.sin(N*x)/np.sin(x)
    return G
    
def calculate():
    delta = np.exp(log_delta_double.get())
    a = a_double.get()
    b = 4*a
    b_string.set(round(b,2))
    kx0 = kx0_double.get()/b*2*np.pi
    N = N_int.get()
    A = (N-1)*b+2*a
    A_string.set(round(A,2))
    N_F = 0.1
    N_F_string.set(round(N_F,2))
    z_B = A**2/N_F
    z_B_string.set(round(z_B,2))
    propagation = propagation_string.get()
    window = np.exp(log_window_double.get())
    
    Nx = 2**16
    Lx = 2*150*a
    x, delta_x = np.linspace(-Lx/2,Lx/2,Nx,retstep=True) # x normalized to lambda
    v = np.where(np.abs(x+a-A/2)<a,1,0)
    for index in range(1,N):
        v = v + np.where(np.abs(x+a-A/2+index*b)<a,1,0)    
    a1.cla()        
    a1.plot(x,np.abs(v)**2,'k', linewidth=.5)
    a1.set_xlim([-150*a,150*a])
    a1.set_ylim([-0.05,1.15])
    a1.annotate(r'', xy=(-A/2,1.075), xytext=(np.max([-A/2-0.5*150*a,-150*a]),1.075), arrowprops=dict(arrowstyle='->'))
    a1.annotate(r'', xy=(A/2,1.075), xytext=(np.min([A/2+0.5*150*a,150*a]),1.075), arrowprops=dict(arrowstyle='->'))
    a1.annotate(r'A', xy=(0,1.06),horizontalalignment='center', verticalalignment='center')
    a1.set_xlabel(r'$x/\lambda$')
    a1.set_ylabel(r'$|u(z=0)|^2$ [norm. u.]')

    vG = np.vectorize(G)

    k1 = 2*np.pi/(1-delta)
    k2 = 2*np.pi/(1+delta)
    a2.cla()    
    xzB = np.linspace(-1/a,1/a,Nx)*window # x normalized to zB, or kx normalized to k
    xbarzB = xzB-kx0/(2*np.pi) # shifted normalized x or kx due to oblique incidence
    a2.plot(xzB,(np.sinc(k1*a*xbarzB/np.pi)*vG(k1*b*xbarzB/2,N)/N)**2,'b',label=r'$\lambda = \lambda_0-\Delta\lambda$')
    a2.plot(xzB,(np.sinc(k2*a*xbarzB/np.pi)*vG(k2*b*xbarzB/2,N)/N)**2,'r',label=r'$\lambda = \lambda_0+\Delta\lambda$')
    a2.set_xlim([-window/a,window/a])
    a2.set_xlabel(r'$k_x/k_0$')
    a2.set_ylabel(r'$|U(z=0)|^2$ [norm. u.]')
    a2.legend(loc='upper right')
    
    if propagation == 'paraxial':
        I1FH = (2*a)**2/z_B/(1-delta)*(np.sinc(k1*a*xbarzB/np.pi)*vG(k1*b*xbarzB/2,N))**2
        I2FH = (2*a)**2/z_B/(1+delta)*(np.sinc(k2*a*xbarzB/np.pi)*vG(k2*b*xbarzB/2,N))**2
    else:
        rBzB = np.sqrt(xzB**2+1) # rB normalized to zB
        rB = rBzB*z_B # rB normalized to lambda
        xbarzB = xzB-kx0*rBzB/(2*np.pi)
        I1FH = (2*a)**2*z_B**2/rB**3/(1-delta)*(np.sinc(k1*a*xbarzB/rBzB/np.pi)*vG(k1*b*xbarzB/rBzB/2,N))**2
        I2FH = (2*a)**2*z_B**2/rB**3/(1+delta)*(np.sinc(k2*a*xbarzB/rBzB/np.pi)*vG(k2*b*xbarzB/rBzB/2,N))**2
    a3.cla()
    a3.plot(xzB,I1FH,'b')
    a3.plot(xzB,I2FH,'r')
    a3.set_xlim([-window/a,window/a])
    a3.set_xlabel(r'$x/z_{\rm B}$')
    a3.set_ylabel(r'$|u(z=z_{\rm B})|^2$ [norm. u.]')
    
    plt.tight_layout()
    
#    plt.savefig('FPS_spectro.pdf',bbox_inches='tight',dpi=300, transparent=True)
    
    canvas.draw()       

f = plt.figure(1,[6,6.75])
a1 = f.add_subplot(311)
a2 = f.add_subplot(312)
a3 = f.add_subplot(313)

canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

kx0_double = Tk.DoubleVar()
log_delta_double = Tk.DoubleVar()
a_double = Tk.DoubleVar()
b_string = Tk.StringVar()
N_int = Tk.IntVar()
A_string = Tk.StringVar()
N_F_string = Tk.StringVar()
z_B_string = Tk.StringVar()
propagation_string = Tk.StringVar()
log_window_double = Tk.DoubleVar()

initialize()

row = 1
row = gui.create_logslider_with_latex(mainframe,r"$\Delta\lambda/\lambda_0 =$",log_delta_double,0.001,0.1,row)
row = gui.create_slider_with_latex(mainframe,r"half slit width $a/\lambda_0 =$",a_double,.1,4,row)
row = gui.create_label_with_latex(mainframe,r"period length $b/\lambda_0 =$",b_string,row)
row = gui.create_intslider_with_latex(mainframe,r"Number of periods $N =$",N_int,1,1000,row)
row = gui.create_label_with_latex(mainframe,r"Aperture $A/\lambda_0 =$",A_string,row)
row = gui.create_label_with_latex(mainframe,r"Fresnel Number $N_{\rm F} = A^2/(\lambda_0 z_{\rm B}) <$",N_F_string,row)
row = gui.create_label_with_latex(mainframe,r"Distance to screen $z_{\rm B} / \lambda_0 >$",z_B_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r"angle of incidence $b\sin(\varphi)/\lambda_0 =$",kx0_double,-1,1,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_logslider_with_latex(mainframe,r"window size $|k_x|a/(2 \pi) = |x|a/(z_{\rm B}\lambda_0) <$",log_window_double,0.001,1,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_radiobutton(mainframe,['vacuum propagation:','paraxial','exact'],propagation_string,2,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)