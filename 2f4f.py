import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("2f/4f system")
#root.geometry("1280x800")

def plot_spect(ax,xlim,ylabel):
    global kx,k,V,foufac
    ax.plot(np.fft.fftshift(kx)/k,(np.abs(np.fft.fftshift(V))*foufac)**2,'b')
    ax.set_xlim(xlim)
    ax.set_xlabel(r'$k_x/k$')
    ax.set_ylabel(ylabel)
    
def plot_conf(ax,xlim,ylabel):
    global x,f,v
    ax.plot(x/f,np.abs(v)**2,'b')
    ax.set_xlim(xlim)
    ax.set_ylim([-np.max((np.abs(v))**2)*0.05,np.max((np.abs(v))**2)*1.15])
    ax.set_xlabel(r'$x/f$')
    ax.set_ylabel(ylabel)
    
def lens(ax,xlim,ylabel):
    global propagation,k,x,v,AL
    if propagation == 'paraxial':
        v = np.where(np.abs(x)<AL,v,0)*np.exp(-1j*k*x**2/(2*f))
    else:
        v = np.where(np.abs(x)<AL,v,0)*np.exp(-1j*k*np.sqrt(f**2+x**2))
    ax.plot(x/f,np.abs(v)**2,'b')
    ax.set_xlim(xlim)
    ax.set_ylim([-np.max((np.abs(v))**2)*0.05,np.max((np.abs(v))**2)*1.15])
    if AL/(f)<xlim[1]:
        if AL/(f)>0.04*xlim[1]:
           ax.annotate(r'2A', xy=(0,np.max((np.abs(v))**2)*1.075),horizontalalignment='center', verticalalignment='center')
        if AL/(f)>0.1*xlim[1]:
           ax.annotate(r'', xy=(-AL/(f),np.max((np.abs(v))**2)*1.075), xytext=(-0.05*xlim[1],np.max((np.abs(v))**2)*1.075), arrowprops=dict(arrowstyle='->'))
           ax.annotate(r'', xy=(AL/(f),np.max((np.abs(v))**2)*1.075), xytext=(0.05*xlim[1],np.max((np.abs(v))**2)*1.075), arrowprops=dict(arrowstyle='->'))
    ax.plot([-AL/f,-AL/f],ax.get_ylim(),'k:',[AL/f,AL/f],ax.get_ylim(),'k:')
    ax.set_xlabel(r'$x/f$')
    ax.set_ylabel(ylabel)
    
def prop():
    global propagation,kx,k,f,V,v
    if propagation == 'paraxial':
        V = V*np.exp(-1j*kx**2/(2*k)*f)
    else:
        V = V*np.exp(1j*np.sqrt(k**2-kx**2+0j)*f) 
    v = np.fft.ifft(V) 
    
def initialize():
    kx0_double.set(0)
    b_double.set(2)
    w_double.set(6)
    f_double.set(1000)
    log_A_double.set(np.log(.1)) 
    propagation_string.set("exact")
    setup_string.set("2f")
    
    calculate()
    
def calculate():
    global propagation,kx,k,f,V,foufac,v,x,AL
    
    gui.change_cursor(root,"trek")
    b = b_double.get()
    w = w_double.get()*b
    kx0 = kx0_double.get()*2*np.pi/b
    f = f_double.get()*b
    AL = np.exp(log_A_double.get())*f
    propagation = propagation_string.get()
    setup = setup_string.get()
    
    Nx = 2**22
    Lx = b*Nx/50.01
    x, delta_x = np.linspace(-Lx/2,Lx/2,Nx,endpoint=False, retstep=True)
    v = np.exp(-x**2/w**2)*(1+np.cos(2*np.pi*x/b))*np.exp(1j*kx0*x)/2   
    fig.clf() 
    a1 = fig.add_subplot(421) 
    a1.plot(x,np.abs(v)**2,'k')
    a1.set_xlim([-3*w,3*w])
    a1.set_xlabel(r'$x/\lambda$')
    a1.set_ylabel(r'$|u(z=0)|^2$ [norm. u.]')
    
    V = np.fft.fft(v)
    kx = 2*np.pi*np.fft.fftfreq(Nx,delta_x)
    k = 2*np.pi     
    a2 = fig.add_subplot(422)    
    foufac = 1/np.max(np.abs(V))
    plot_spect(a2,[-2/b,2/b],r'$|U(z=0)|^2$ [norm. u.]')
    
    prop()
    xlim = [-np.max([3*w/f,2/b]),np.max([3*w/f,2/b])]
    a3 = fig.add_subplot(423)  
    a4 = fig.add_subplot(424) 
    a5 = fig.add_subplot(425) 
    if setup == '2f':
        plot_conf(a3,xlim,r'$|u_-(z=f)|^2$ [norm. u.]')
        plot_spect(a4,[-2/b,2/b],r'$|U_-(z=f)|^2$ [norm. u.]')
        lens(a5,xlim,r'$|u_+(z=f)|^2$ [norm. u.]')
    else:
        lens(a3,xlim,r'$|u_+(z=f)|^2$ [norm. u.]')
    
    V = np.fft.fft(v)     
    a6 = fig.add_subplot(426) 
    if setup == '2f':
        plot_spect(a6,[-3*w/f,3*w/f],r'$|U_+(z=f)|^2$ [norm. u.]')
        
    prop()    
    a7 = fig.add_subplot(427)
    a8 = fig.add_subplot(428) 
    if setup == '2f':
        plot_conf(a7,xlim,r'$|u(z=2f)|^2$ [norm. u.]')
        plot_spect(a8,[-3*w/f,3*w/f],r'$|U(z=2f)|^2$ [norm. u.]')
    else:    
        plot_conf(a4,xlim,r'$|u(z=2f)|^2$ [norm. u.]')
        plot_spect(a5,[-3*w/f,3*w/f],r'$|U(z=2f)|^2$ [norm. u.]')
        prop()
        lens(a6,xlim,r'$|u_+(z=3f)|^2$ [norm. u.]')
        V = np.fft.fft(v) 
        prop()
        a7.plot(x,np.abs(v)**2,'b')
        a7.set_xlim([-3*w,3*w])
        a7.set_xlabel(r'$x/\lambda$')
        a7.set_ylabel(r'$|u(z=4f)|^2$ [norm. u.]')
        plot_spect(a8,[-2/b,2/b],r'$|U(z=4f)|^2$ [norm. u.]')
    
    plt.tight_layout()
    
#    plt.savefig('2f4f.pdf',bbox_inches='tight',dpi=300, transparent=True)
    
    canvas.draw()       
    gui.change_cursor(root,"arrow")

fig = plt.figure(1,[9,9])

canvas = gui.create_canvas(root,fig)
mainframe = gui.create_mainframe(root)

kx0_double = Tk.DoubleVar()
b_double = Tk.DoubleVar()
w_double = Tk.DoubleVar()
f_double = Tk.DoubleVar()
log_A_double = Tk.DoubleVar()
propagation_string = Tk.StringVar()
setup_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_slider_with_latex(mainframe,r"period length $b/\lambda =$",b_double,.8,8,row)
row = gui.create_slider_with_latex(mainframe,r"Gaussian width $w/b =$",w_double,1,20,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r"angle of incidence $b\sin(\varphi)/\lambda =$",kx0_double,-0.6,0.6,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r"focal length of lens $f/b =$",f_double,10,2000,row)
row = gui.create_logslider_with_latex(mainframe,r"half aperture size of lens $A/f =$",log_A_double,0.01,100,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_radiobutton(mainframe,['vacuum propagation:','paraxial','exact'],propagation_string,2,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_radiobutton(mainframe,['setup:','2f','4f'],setup_string,2,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)