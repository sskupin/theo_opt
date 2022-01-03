import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("Taylor Series")

def taylor_sin(x,x0,N):
    taylor = np.sin(x0)*np.ones_like(x) # compute zero order
    fac = 1
    for index in range(1,N+1):
        fac = fac*index # compute factorial
        taylor = taylor + np.sin(x0+index*np.pi/2)*(x-x0)**index/fac # add term of order "index"
    
    return taylor

def initialize():
    global x0_save,N_save,lx_save
    x0_string.set("1")
    N_string.set("10")
    lx_string.set("10")
    
    x0_save = x0_string.get()
    N_save = N_string.get()
    lx_save = lx_string.get()
    
    calculate()
    
def reinitialize():
    global x0_save,N_save,lx_save
    x0_string.set(x0_save)
    N_string.set(N_save)
    lx_string.set(lx_save)

def calculate():
    global x0_save,N_save,lx_save
    try:
        x0 = float(x0_string.get())
        N = int(float(N_string.get()))
        N_string.set(N)
        lx = float(lx_string.get())
        
        if N < 0 or lx <= 0: 
            gui.input_error("N and L should be positive. Re-initializing ...", reinitialize)
        elif N > 100:
            gui.input_error("N is too large. Re-initializing ...", reinitialize)
        else:
            x = np.linspace(x0-lx/2,x0+lx/2,1000)           

            f.clf()
            a = f.add_subplot(111)
            a.plot(x, np.sin(x), 'b', x, taylor_sin(x,x0,N) ,'r')
            plt.xlabel(r'$x$')
            plt.ylabel(r'$f(x), \mathcal{T}^{(N)}_{f}(x;x_0)$')
            a.legend((r'$\sin(x)$',r'$\mathcal{T}^{(N)}_{\sin}(x;x_0)$'))
            plt.tight_layout()
            
#            plt.savefig('taylor_series.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            x0_save = x0_string.get()
            N_save = N_string.get()
            lx_save = lx_string.get()

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)

f = plt.figure(1,[7,3.5])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

x0_string = Tk.StringVar()
N_string = Tk.StringVar()
lx_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_title(mainframe,"Taylor polynome",row)
row = gui.create_entry_with_latex(mainframe,r"center $x_0 =$",x0_string,row)
row = gui.create_entry_with_latex(mainframe,r"order $N =$",N_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r" window size $L =$",lx_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)