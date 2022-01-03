import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("Fourier Series")

def fourier_sawtooth(x,N):
    fourier = np.zeros_like(x) 
    for index in range(1,N+1):
         fourier = fourier + 2*(-1)**(index+1)/(index*np.pi)*np.sin(index*x) # add term of order "index"
    
    return fourier

def sawtooth(x):
    sawtooth = np.zeros_like(x) 
    sawtooth[np.where((x>-np.pi) & (x<np.pi))[0]] = x[np.where((x>-np.pi) & (x<np.pi))[0]]/np.pi
    sawtooth[np.where((x>-3*np.pi) & (x<-np.pi))[0]] = x[np.where((x>-3*np.pi) & (x<-np.pi))[0]]/np.pi+2
    sawtooth[np.where((x>np.pi) & (x<3*np.pi))[0]] = x[np.where((x>np.pi) & (x<3*np.pi))[0]]/np.pi-2
    sawtooth[np.where((x>-5*np.pi) & (x<-3*np.pi))[0]] = x[np.where((x>-5*np.pi) & (x<-3*np.pi))[0]]/np.pi+4
    sawtooth[np.where((x>3*np.pi) & (x<5*np.pi))[0]] = x[np.where((x>3*np.pi) & (x<5*np.pi))[0]]/np.pi-4
    
    return sawtooth

def initialize():
    N_string.set("5")
    
    calculate()

def calculate():
    try:
        N = int(float(N_string.get()))
        N_string.set(N)
        
        if N < 0: 
            gui.input_error("N should be positive. Re-initializing ...", initialize)
        elif N > 1000:
            gui.input_error("N is too large. Re-initializing ...", initialize) 
        else:
            x = np.linspace(-10,10,100000)           

            f.clf()
            a = f.add_subplot(111)
            a.plot(x, sawtooth(x), 'b', x, fourier_sawtooth(x,N) ,'r')
            plt.xlabel(r'$x$')
            plt.ylabel(r'$f(x), \mathcal{F}^{(N)}_{f}(x)$')
            a.legend((r'sawtooth$(x)$',r'$\mathcal{F}^{(N)}_{\mathrm{sawtooth}}(x)$'))
            plt.tight_layout()
            
#            plt.savefig('fourier_series.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", initialize)

f = plt.figure(1,[7,3.5])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

x0_string = Tk.StringVar()
N_string = Tk.StringVar()
lx_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_title(mainframe,"Fourier sum",row)
row = gui.create_entry_with_latex(mainframe,r"order $N = $",N_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)