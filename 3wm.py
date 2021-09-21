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
root.title("3-wave mixing")

def initialize():
    global LDeltak_save,A1abs_save,A1phase_save,A2abs_save,A2phase_save,A3abs_save,A3phase_save
    LDeltak_string.set("5")
    A1abs_string.set("1")
    A1phase_string.set("0")
    A2abs_string.set("10")
    A2phase_string.set("0")
    A3abs_string.set("0")
    A3phase_string.set("0")
    A1show_string.set("showA1")
    A2show_string.set("noshow")
    A3show_string.set("showA3")
    MRshow_string.set("noshow")
    UDPA2show_string.set("showUDPA2")
    UDPA3show_string.set("noshow")
    
    LDeltak_save = LDeltak_string.get()
    A1abs_save = A1abs_string.get()
    A1phase_save = A1phase_string.get()
    A2abs_save = A2abs_string.get()
    A2phase_save = A2phase_string.get()
    A3abs_save = A3abs_string.get()
    A3phase_save = A3phase_string.get()
    
    calculate()
    
def reinitialize():
    global LDeltak_save,A1abs_save,A1phase_save,A2abs_save,A2phase_save,A3abs_save,A3phase_save
    LDeltak_string.set(LDeltak_save)
    A1abs_string.set(A1abs_save)
    A1phase_string.set(A1phase_save)
    A2abs_string.set(A2abs_save)
    A2phase_string.set(A2phase_save)
    A3abs_string.set(A3abs_save)
    A3phase_string.set(A3phase_save)

def calculate():
    global LDeltak_save,A1abs_save,A1phase_save,A2abs_save,A2phase_save,A3abs_save,A3phase_save
    try:
        LDeltak = float(LDeltak_string.get())*np.pi
        A = np.zeros(3) + 0j 
        A[0] = float(A1abs_string.get())*np.exp(1j*float(A1phase_string.get()))
        A[1] = float(A2abs_string.get())*np.exp(1j*float(A2phase_string.get()))
        A[2] = float(A3abs_string.get())*np.exp(1j*float(A3phase_string.get()))
        
        def compute_rhs(Z,A): # computes rhs of ode system, A[j-1] = A_j
            rhs1 = 1j * A[2] * np.conj(A[1])
            rhs2 = 1j * A[2] * np.conj(A[0])
            rhs3 = 1j * (A[0] * A[1] - LDeltak * A[2])
            return np.array([rhs1, rhs2, rhs3])
        
        sol = spi.solve_ivp(compute_rhs, [0, 1], A, max_step = 1.e-3)
        
        if UDPA2show_string.get() == 'showUDPA2' and UDPA3show_string.get() == 'showUDPA3':
            gui.input_warning("You are about to plot the undepleted pump approximation for both sum and difference frequency generation. Let us pick the larger pump only...")
            if np.abs(A[1]) >= np.abs(A[2]):
                UDPA3show_string.set('noshow')
        f.clf() 
        
        a1 = f.add_subplot(111)
        if A1show_string.get() == 'showA1':
            a1.plot(sol.t, np.abs(sol.y[0,:])**2, 'r', label=r'$|A_1|^2$')
        if A2show_string.get() == 'showA2':
            a1.plot(sol.t, np.abs(sol.y[1,:])**2, 'g', label=r'$|A_2|^2$')
        if A3show_string.get() == 'showA3':            
            a1.plot(sol.t, np.abs(sol.y[2,:])**2, 'b', label=r'$|A_3|^2$')
        if MRshow_string.get() == 'showMR':             
            a1.plot(sol.t, np.abs(sol.y[1,:])**2+np.abs(sol.y[2,:])**2, 'k:', label=r'$|A_2|^2+|A_3|^2$')
            a1.plot(sol.t, np.abs(sol.y[0,:])**2+np.abs(sol.y[2,:])**2, 'k:', label=r'$|A_1|^2+|A_3|^2$')
            a1.plot(sol.t, np.abs(sol.y[0,:])**2-np.abs(sol.y[1,:])**2, 'k:', label=r'$|A_1|^2-|A_2|^2$')
        a1.set_xlabel(r'$Z$')
        a1.set_ylabel(r'Normalized Intensity')
        a1.set_title(r'$L \Delta k=$ '+str(round(LDeltak/np.pi,4))+r'$\pi$')
        if UDPA2show_string.get() == 'showUDPA2':
            Delta = LDeltak/2
            g = np.sqrt(LDeltak**2/4+A[1]**2)
            a1.plot(sol.t, np.abs(A[0]*np.cos(g*sol.t)+1j*(A[0]*Delta+A[1]*A[2])*np.sin(g*sol.t)/g)**2, color='orange', linestyle='--', label=r'$|A_1^{\rm undepleted}|^2$')
            a1.plot(sol.t, np.abs(A[2]*np.cos(g*sol.t)+1j*(A[1]*A[0]-A[2]*Delta)*np.sin(g*sol.t)/g)**2, color='cyan', linestyle='--',label=r'$|A_3^{\rm undepleted}|^2$')
        if UDPA3show_string.get() == 'showUDPA3': 
            Delta = LDeltak/2
            g = np.sqrt(-LDeltak**2/4+A[2]**2)
            a1.plot(sol.t, np.abs(A[0]*np.cosh(g*sol.t)+1j*(A[0]*Delta+np.conj(A[1])*A[2])*np.sinh(g*sol.t)/g)**2, color='orange', linestyle='--', label=r'$|A_1^{\rm undepleted}|^2$')
            a1.plot(sol.t, np.abs(A[1]*np.cosh(g*sol.t)+1j*(A[1]*Delta+np.conj(A[0])*A[2])*np.sinh(g*sol.t)/g)**2, color='lime', linestyle='--', label=r'$|A_2^{\rm undepleted}|^2$')
        a1.legend()
        
        LDeltak_save = LDeltak_string.get()
        A1abs_save = A1abs_string.get()
        A1phase_save = A1phase_string.get()
        A2abs_save = A2abs_string.get()
        A2phase_save = A2phase_string.get()
        A3abs_save = A3abs_string.get()
        A3phase_save = A3phase_string.get()        

        canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)

f = plt.figure(1,[7,4])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

LDeltak_string = Tk.StringVar()
A1abs_string = Tk.StringVar()
A1phase_string = Tk.StringVar()
A1show_string = Tk.StringVar()
A2abs_string = Tk.StringVar()
A2phase_string = Tk.StringVar()
A2show_string = Tk.StringVar()
A3abs_string = Tk.StringVar()
A3phase_string = Tk.StringVar()
A3show_string = Tk.StringVar()
MRshow_string = Tk.StringVar()
UDPA2show_string = Tk.StringVar()
UDPA3show_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"phase mismatch $L \Delta k / \pi =$",LDeltak_string,row)
row = gui.create_double_entry_with_latex(mainframe,r"$|A_1(Z=0)|=$",A1abs_string,r"arg$[A_1(Z=0)]=$",A1phase_string,row)
row = gui.create_double_entry_with_latex(mainframe,r"$|A_2(Z=0)|=$",A2abs_string,r"arg$[A_2(Z=0)]=$",A2phase_string,row)
row = gui.create_double_entry_with_latex(mainframe,r"$|A_3(Z=0)|=$",A3abs_string,r"arg$[A_3(Z=0)]=$",A3phase_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'$|A_1|^2$','noshow','showA1',A1show_string,r'Manley-Rowe','noshow','showMR',MRshow_string,row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'$|A_2|^2$','noshow','showA2',A2show_string,r'undepleted $|A_2|^2$','noshow','showUDPA2',UDPA2show_string,row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'$|A_3|^2$','noshow','showA3',A3show_string,r'undepleted $|A_3|^2$','noshow','showUDPA3',UDPA3show_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)