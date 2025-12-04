import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import sys
sys.path.append('./aux')
import gui_stuff as gui

gui.set_rcParams()
title = "Taylor Series"
root = Tk.Tk()
root.title(title)

def taylor_sin(x,x0,N):
    taylor = np.sin(x0)*np.ones_like(x) # compute zero order
    fac = 1
    for index in range(1,N+1):
        fac = fac*index # compute factorial
        taylor = taylor + np.sin(x0+index*np.pi/2)*(x-x0)**index/fac # add term of order "index"
    
    return taylor

def initialize():
    var_string[0].set("1")
    var_string[1].set("10")
    var_string[2].set("10")
    
    gui.copy_stringvar_vector(var_string,var_save)
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()
    
def show_manual():
    gui.show_manual("man/taylor_series.png",title) 

def calculate():
    try:
        x0 = float(var_string[0].get())
        N = int(float(var_string[1].get()))
        var_string[1].set(N)
        lx = float(var_string[2].get())
        
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
            
            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)

f = plt.figure(1,[7,3.5])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(3)
var_save = gui.create_stringvar_vector(3)

initialize()

row = 1
row = gui.create_title(mainframe,"Taylor polynome",row)
row = gui.create_entry_with_latex(mainframe,r"center $x_0 =$",var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r"order $N =$",var_string[1],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r" window size $L =$",var_string[2],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)