import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("Fast Fourier Transform of Gaussian")
    
def initialize():
    N_string.set("64")
    L_string.set("16")
    
    calculate()

def calculate():
    try:
        N = int(float(N_string.get()))
        N_string.set(N)
        L = float(L_string.get())
        
        if N <= 0 or L <= 0: 
            gui.input_error("N and L should be positive. Re-initializing ...",initialize)
        elif np.floor(N/2) < N/2:
            gui.input_error("N should be even. Re-initializing ...",initialize)
        else:
            f.clf()
            x,delta_x = np.linspace(-L/2, L/2, num=N, endpoint=False, retstep=True) # x in units of w_0
            xa = np.linspace(-L/2, L/2, num=int(16*L), endpoint=False)
            ka = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(int(128/delta_x),delta_x))
            Gauss = np.exp(-xa**2) 
            Gauss_diskret = np.exp(-x**2) + 0j
            FTGauss = np.exp(-ka**2/4)/(2*np.sqrt(np.pi))
            k = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(N,delta_x))
            FFTGauss_diskret = delta_x/(2*np.pi)*np.fft.fftshift(np.fft.fft(np.fft.fftshift(Gauss_diskret)))
            a1 = f.add_subplot(131)
            a1.plot(xa,Gauss,'b', label=r'$f$')
            a1.plot(x,np.real(Gauss_diskret),'r*', label=r'$f_n$')
            plt.legend()
            a1.set_xlabel(r'$x$')
            a1.set_ylabel(r'$f$, $f_n$')
            a2 = f.add_subplot(132)
            a2.plot(ka,FTGauss,'b',label=r'$\hat{f}$')
            a2.plot(k,np.real(FFTGauss_diskret),'r*', label=r'$F_m$')
            a2.set_xlabel(r'$k$')
            a2.set_ylabel(r'$\hat{f}$, $F_m$')
            plt.legend()
            a3 = f.add_subplot(133)
            a3.semilogy(ka,FTGauss,'b',label=r'$\hat{f}$')
            a3.semilogy(k,np.real(FFTGauss_diskret),'r*', label=r'$F_m$')
            a3.set_xlabel(r'$k$')
            a3.set_ylabel(r'$\hat{f}$, $F_m$')
            plt.legend()
              
            plt.tight_layout()
            Lk_string.set(2*np.pi*N/L)
            
#            plt.savefig('fft_gauss.pdf',bbox_inches='tight',dpi=300, transparent=True)

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", initialize)

f = plt.figure(1,[8,2])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

N_string = Tk.StringVar()
L_string = Tk.StringVar()
Lk_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"sample points $N =$",N_string,row)
row = gui.create_entry_with_latex(mainframe,r"window size $L_x =$",L_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_label_with_latex(mainframe,r'$L_k=2 \pi N / L_x =$',Lk_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)