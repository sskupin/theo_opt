import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import sys
sys.path.append('./aux')
import gui_stuff as gui
import ff_stuff as ff

gui.set_rcParams()
title = "Image Formation by Thin Lenses - 4f System"
root = Tk.Tk()
root.title(title)

def initialize():
    kx0_double.set(0)
    b_double.set(2)
    w_double.set(6)
    f_double.set(1000)
    log_A_double.set(np.log(.1)) 
    propagation_string.set("exact")
    calculate()
    
def show_manual():
    gui.show_manual("man/image_4f.png",title)
    
def calculate():    
    gui.change_cursor(root,"trek")
    b = b_double.get()
    w = w_double.get()*b
    kx0 = kx0_double.get()*2*np.pi/b
    f = f_double.get()*b
    AL = np.exp(log_A_double.get())*f
    propagation = propagation_string.get()
    
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
    ff.plot_spect(a2,[-2/b,2/b],r'$|U(z=0)|^2$ [norm. u.]',kx,k,V,foufac)
    
    v,V = ff.prop(propagation,kx,k,f,V,v)
    xlim = [-np.max([3*w/f,2/b]),np.max([3*w/f,2/b])]
    a3 = fig.add_subplot(423)  
    a4 = fig.add_subplot(424) 
    a5 = fig.add_subplot(425) 
    v = ff.lens(a3,xlim,r'$|u_+(z=f)|^2$ [norm. u.]',propagation,k,x,v,AL,f)
    
    V = np.fft.fft(v)     
    a6 = fig.add_subplot(426) 
        
    v,V = ff.prop(propagation,kx,k,f,V,v)
    a7 = fig.add_subplot(427)
    a8 = fig.add_subplot(428) 
    ff.plot_conf(a4,xlim,r'$|u(z=2f)|^2$ [norm. u.]',x,f,v)
    ff.plot_spect(a5,[-3*w/f,3*w/f],r'$|U(z=2f)|^2$ [norm. u.]',kx,k,V,foufac)
    v,V = ff.prop(propagation,kx,k,f,V,v)
    v = ff.lens(a6,xlim,r'$|u_+(z=3f)|^2$ [norm. u.]',propagation,k,x,v,AL,f)
    V = np.fft.fft(v) 
    v,V = ff.prop(propagation,kx,k,f,V,v)
    a7.plot(x,np.abs(v)**2,'b')
    a7.set_xlim([-3*w,3*w])
    a7.set_xlabel(r'$x/\lambda$')
    a7.set_ylabel(r'$|u(z=4f)|^2$ [norm. u.]')
    ff.plot_spect(a8,[-2/b,2/b],r'$|U(z=4f)|^2$ [norm. u.]',kx,k,V,foufac)
    
    plt.tight_layout()
    
#    plt.savefig('4f.pdf',bbox_inches='tight',dpi=300, transparent=True)
    
    canvas.draw()       
    gui.change_cursor(root,"arrow")

fig = plt.figure(1,[9,9])

canvas = gui.create_canvas(root,fig)
canvas.draw() 
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
row = gui.create_slider_with_latex(mainframe,r"period length $b/\lambda =$",b_double,.8,8,row,increment=.25)
row = gui.create_slider_with_latex(mainframe,r"Gaussian width $w/b =$",w_double,1,20,row,increment=.25)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r"angle of incidence $b\sin(\varphi)/\lambda =$",kx0_double,-0.6,0.6,row,increment=.05)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r"focal length of lens $f/b =$",f_double,10,2000,row,increment=10)
row = gui.create_logslider_with_latex(mainframe,r"half aperture size of lens $A/f =$",log_A_double,0.01,100,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_radiobutton(mainframe,['vacuum propagation:','paraxial','exact'],propagation_string,2,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)