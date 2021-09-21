import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.integrate as spi
import tkinter as Tk
import gui_stuff as gui

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rc('text', usetex=True)
mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')
mpl.rcParams.update({'font.size': 10})

root = Tk.Tk()
root.title("SHG with pulse")

def plot_2D(ax, Z, T, AMP, title): # plot 2D amplitude on non-equidistant ZxT grid
    im = mpl.image.NonUniformImage(ax, extent=(Z[0], Z[-1], T[0], T[-1]))
    im.set_data(Z, T, AMP)
    ax.images.append(im)
    ax.set_xlim([Z[0], Z[-1]])
    ax.set_ylim([T[0], T[-1]])
    ax.set_xlabel(r'$Z$')
    ax.set_ylabel(r'$T$')
    ax.set_title(title)
    plt.colorbar(mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=np.amin(AMP), vmax=np.amax(AMP))))
    
def plot_1D(ax, T, labelT, AMP, labelAMP, PHASE, labelPHASE, Tmin, Tmax): # plot 1D amplitude and phase 
    ax.plot(T, AMP, 'b')
    ax.set_xlim([Tmin, Tmax])
    ax.set_xlabel(labelT)
    ax.set_ylabel(labelAMP, color='b')
    ax.tick_params(axis='y', labelcolor='b')
    axbis = ax.twinx()
    axbis.plot(T, PHASE, 'r:')
    axbis.set_ylabel(labelPHASE, color='r')
    axbis.tick_params(axis='y', labelcolor='r')  

def initialize():
    LDeltak_double.set(1)
    LTpDeltakp_double.set(-1)
    A1_double.set(1)
    
    calculate()
    
def calculate():
    LDeltak = LDeltak_double.get()*np.pi
    LTpDeltakp = LTpDeltakp_double.get()
    A1 = A1_double.get()
    
    N = 1024 # number of grid points in T (even)
    LT = 16 # size of temporal window (delay) in units of T_p
    A = np.zeros(2*N) + 0j 
    T,delta_T = np.linspace(-LT/2, LT/2, num=N, endpoint=False, retstep=True) # delay in units of T_p
    A[0:N] = A1 * np.exp(-T**2) # set A_1 Gaussian pulse amplitude at z = 0
    A[N:2*N] = 0 # * np.exp(-T**2) # set A_2 amplitude at z = 0

    def compute_rhs(Z, A): # computes rhs of ode system, A[0:N] = A_1, A[N:2*N] = A_2
        rhs1 = 1j * A[N:2*N] * np.conj(A[0:N])
        rhs2 = 1j * (A[0:N]**2 - LDeltak * A[N:2*N]) - LTpDeltakp * np.gradient(A[N:2*N],delta_T)
        return np.concatenate((rhs1, rhs2))

    sol = spi.solve_ivp(compute_rhs, [0, 1], A, max_step = 1.e-3)
    
    f.clf() 

    # plot amplitudes and phases
    a1 = f.add_subplot(231)
    plot_2D(a1, sol.t, T, np.abs(sol.y[0:N,:]), r'$|A_1|$')
    a2 = f.add_subplot(232)
    plot_1D(a2, T, r'$T$', np.abs(sol.y[0:N,-1]), r'$|A_1(Z=1)|$', np.angle(sol.y[0:N,-1]), r'arg$A_1(Z=1)$', -6, 6)
    a2.set_title(r'$L \Delta k=$ '+str(round(LDeltak/np.pi,4))+r'$\pi \qquad L \Delta k^\prime / T_{\rm p} =$'+str(round(LTpDeltakp,4)))
    a3 = f.add_subplot(234)
    plot_2D(a3, sol.t, T, np.abs(sol.y[N:2*N,:]), r'$|A_2|$')
    a4 = f.add_subplot(235)
    plot_1D(a4, T, r'$T$', np.abs(sol.y[N:2*N,-1]), r'$|A_2(Z=1)|$', np.angle(sol.y[N:2*N,-1]), r'arg$\,A_2(Z=1)$', -6, 6)
    
    # plot spectral amplitudes and phases
    Omega = 2 * np.pi * np.fft.fftshift(np.fft.fftfreq(N, d=delta_T)) # centered frequency in units of 1/T_p
    FTA_1 = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(sol.y[0:N,-1])))
    a5 = f.add_subplot(233)
    plot_1D(a5, Omega, r'$\Omega+\Omega_0$', np.abs(FTA_1), r'$|FT[A_1(Z=1)]|$', np.angle(FTA_1), r'arg$FT[A_1(Z=1)]$', -20, 20)
    FTA_2 = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(sol.y[N:2*N,-1])))
    a6 = f.add_subplot(236)
    plot_1D(a6, Omega, r'$\Omega+2\Omega_0$', np.abs(FTA_2), r'$|FT[A_2(Z=1)]|$', np.angle(FTA_2), r'arg$FT[A_2(Z=1)]$', -20, 20)

    plt.tight_layout()  

              
#    plt.savefig('shg_pulse.pdf',bbox_inches='tight',dpi=300, transparent=True)

    canvas.draw()       

f = plt.figure(1,[12,6])

canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

LDeltak_double = Tk.DoubleVar()
LTpDeltakp_double = Tk.DoubleVar()
A1_double = Tk.DoubleVar()

initialize()

row = 1
row = gui.create_slider_with_latex(mainframe,r'phase mismatch $L \Delta k / \pi =$',LDeltak_double,-10,10,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r'group velocity mismatch $L \Delta k^\prime / T_{\rm p} =$',LTpDeltakp_double,-5,5,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r'fundamental amplitude $|A_1|=$',A1_double,0.1,10,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)