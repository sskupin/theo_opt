import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import sys
sys.path.append('./stuff')
import gui_stuff as gui
import bpm_stuff as bpm
import aniso_stuff as ani

gui.set_rcParams()
title = "Beam Propagation in Homogeneous and Anisotropic Media"
root = Tk.Tk()
root.title(title)

#K. Kato and E. Takaoka. Sellmeier and thermo-optic dispersion formulas for KTP, Appl. Opt. 41, 5040-5044 (2002).
# 587.6 nm

def initialize():
    var_string[0].set("3.1258")   # epsilon1
    var_string[1].set("3.1606")   # epsilon2
    var_string[2].set("3.5108")   # epsilon3
    var_string[3].set("angles")   # propagation along optical axis
    var_string[4].set("0.3")  # theta0/pi
    var_string[5].set("0.4")  # phi0/pi    
    var_double[0].set(10)    # w_0/\lambda
    var_double[1].set(2)   # z\delta_{max}/w_0 
    var_double[2].set(0.25)     # polarization angle
    var_double[3].set(0)     # polarization ellipticity
    calculate()
    
def show_manual():
    gui.show_manual("man/beam_prop_aniso.png",title)
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate() # because sliders may have changed
    
def calculate():
    gui.change_cursor(root,"trek")
    try:     
        epsilon = np.array([float(var_string[0].get()),float(var_string[1].get()),float(var_string[2].get())])
        theta0 = float(var_string[4].get())*np.pi
        phi0 = float(var_string[5].get())*np.pi
        theta_view = 0.3*180
        phi_view = 0.8*180
        kv = 2*np.pi*var_double[0].get() # in units of 1/w_0
        lambda0 = 1/var_double[0].get() # in units of w_0
        psi = var_double[2].get()*np.pi
        ellipticity = var_double[3].get()
        prop = var_string[3].get()
        
        if (epsilon < 1).any() or (epsilon > 12).any(): 
            gui.input_error("Tensor elements must be between 1 and 12. Re-initializing ...",reinitialize)
        elif epsilon[0] > epsilon[1] or epsilon[1] > epsilon[2]:
            gui.input_error("Tensor elements must not decrease with increasing index. Re-initializing ...",reinitialize)
        elif np.abs(epsilon[0] - epsilon[1]) + np.abs(epsilon[1] - epsilon[2]) < 0.0001:
            gui.input_error("Difference between tensor elements too small. Re-initializing ...",reinitialize)
        else:
            f.clf()
            
            # propagation along optical axis
            if prop=='oaprop':
                thetaoa = np.arcsin(np.sqrt(epsilon[2]*(epsilon[1]-epsilon[0])/(epsilon[1]*(epsilon[2]-epsilon[0]))))
                theta0 = thetaoa
                var_string[4].set(theta0/np.pi)
                phi0 = 0
                var_string[5].set(phi0/np.pi)
            
            # fix lab coordinate vectors ex,ey,ez, as Da,Db,uk
            ez = ani.uk(theta0,phi0)
            ex,ey = ani.D(ez,epsilon) 
            na,nb = ani.dr(ez,epsilon)

            Lx = 8
            Ly = Lx          
            Nx0 = 64
            Ny0 = 64
            Nx = Nx0
            Ny = Ny0
            x_calc, delta_x = np.linspace(-Lx/2,Lx/2,Nx,endpoint=False, retstep=True)
            y_calc, delta_y = np.linspace(-Ly/2,Ly/2,Ny,endpoint=False, retstep=True)
            kx_calc = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(Nx,delta_x))
            ky_calc = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(Ny,delta_y))
            
            KX,KY = np.meshgrid(kx_calc,ky_calc, indexing='xy')
            KX = KX/kv # normalize to kv as in aniso routines
            KY = KY/kv # normalize to kv as in aniso routines
            KZa,KZb = ani.KZ_mesh(KX,KY,epsilon,theta0,phi0)
            
            KZax,KZay = np.gradient(KZa)
            deltaa = np.amax(np.sqrt((KZax[KX**2+KY**2<lambda0**2]*kv/(kx_calc[1]-kx_calc[0]))**2+(KZay[KX**2+KY**2<lambda0**2]*kv/(ky_calc[1]-ky_calc[0]))**2))
            KZbx,KZby = np.gradient(KZb)
            deltab = np.amax(np.sqrt((KZbx[KX**2+KY**2<lambda0**2]*kv/(kx_calc[1]-kx_calc[0]))**2+(KZby[KX**2+KY**2<lambda0**2]*kv/(ky_calc[1]-ky_calc[0]))**2))
            deltamax = max(deltaa,deltab)
            z = var_double[1].get()/deltamax # in units of w_0
            
            Dax,Day,Daz,Dbx,Dby,Dbz = ani.D_proj_mesh(KX,KY,KZa,KZb,ex,ey,ez,epsilon)
            
            D0_x = (np.cos(psi)-1j*ellipticity*np.sin(psi))*bpm.init_2D_beam(x_calc,y_calc)
            D0_y = (np.sin(psi)+1j*ellipticity*np.cos(psi))*bpm.init_2D_beam(x_calc,y_calc)
                       
            prop_a,prop_b = bpm.init_prop_2D_aniso(KZa*kv,KZb*kv,z)
            
            D_x,D_y,D_z,FTD_x,FTD_y,FTD_z,FTDa_x,FTDa_y,FTDa_z,FTDb_x,FTDb_y,FTDb_z = bpm.propagation_2D_aniso(D0_x,D0_y,Dax,Day,Daz,Dbx,Dby,Dbz,prop_a,prop_b)
            
            FTHa_x = (KY*FTDa_z - KZa*FTDa_y)/(KX**2+KY**2+KZa**2)
            FTHa_y = (KZa*FTDa_x - KX*FTDa_z)/(KX**2+KY**2+KZa**2)
            
            FTHb_x = (KY*FTDb_z - KZb*FTDb_y)/(KX**2+KY**2+KZb**2)
            FTHb_y = (KZb*FTDb_x - KX*FTDb_z)/(KX**2+KY**2+KZb**2)
            
            FTH_x = FTHa_x + FTHb_x
            FTH_y = FTHa_y + FTHb_y
            
            E_x = D_x*np.dot(ex/epsilon,ex) + D_y*np.dot(ey/epsilon,ex) + D_z*np.dot(ez.flatten()/epsilon,ex)
            E_y = D_x*np.dot(ex/epsilon,ey) + D_y*np.dot(ey/epsilon,ey) + D_z*np.dot(ez.flatten()/epsilon,ey)
            
            H_x = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(FTH_x)))
            H_y = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(FTH_y)))             
            
            S_z = np.real(E_x*np.conjugate(H_y) - E_y*np.conjugate(H_x))
            
            D0x,x0,y0 = bpm.reshaping_2D(D0_x,x_calc,y_calc,Nx0,Ny0)
            D0y,x0,y0 = bpm.reshaping_2D(D0_y,x_calc,y_calc,Nx0,Ny0)
            Dx,x,y = bpm.reshaping_2D(D_x,x_calc,y_calc,Nx,Ny)
            Dy,x,y = bpm.reshaping_2D(D_y,x_calc,y_calc,Nx,Ny)
            
            Sz,x,y = bpm.reshaping_2D(S_z,x_calc,y_calc,Nx,Ny)
            
            vmax = max(np.amax(np.abs(D0x)),np.amax(np.abs(D0y)))
            
            a1 = f.add_subplot(231)
            im1 = a1.imshow(np.abs(D0x) ,extent=[x0[0], x0[-1], y0[0], y0[-1]] , aspect='equal', origin='lower', vmin=0, vmax=vmax, cmap='jet')
            a1.annotate(r'$|D_{0x}|/|\mathbf{D}_{0}|_{\rm max}$', xy=(0.9*x0[0],0.7*y0[-1]),horizontalalignment='left', verticalalignment='bottom', color='w', size=14)
            a1.set_ylabel(r'$y/w_0$')
            plt.colorbar(im1,location='top',shrink=0.75)
            
            a2 = f.add_subplot(232)
            im2 = a2.imshow(np.abs(D0y) ,extent=[x0[0], x0[-1], y0[0], y0[-1]] , aspect='equal', origin='lower', vmin=0, vmax=vmax, cmap='jet')
            a2.annotate(r'$|D_{0y}|/|\mathbf{D}_{0}|_{\rm max}$', xy=(0.9*x0[0],0.7*y0[-1]),horizontalalignment='left', verticalalignment='bottom', color='w', size=14)
            plt.colorbar(im2,location='top',shrink=0.75)
            
            a12 = f.add_subplot(233, projection='3d')
            a12.view_init(azim=phi_view, elev=90-theta_view)
            ani.plot_ns(a12,theta0,phi0,epsilon,'no_show','no_show',True)
            
            vmax = max(np.amax(np.abs(Dx)),np.amax(np.abs(Dy)))
            
            a3 = f.add_subplot(234)
            im3 = a3.imshow(np.abs(Dx) ,extent=[x[0], x[-1], y[0], y[-1]] , aspect='equal', origin='lower', vmin=0, vmax=vmax, cmap='jet')
            a3.annotate(r'$|D_{x}(z)|/|\mathbf{D}_{0}|_{\rm max}$', xy=(0.9*x[0],0.7*y[-1]),horizontalalignment='left', verticalalignment='bottom', color='w', size=14)
            a3.set_xlabel(r'$x/w_0$')
            a3.set_ylabel(r'$y/w_0$')
            plt.colorbar(im3,location='top',shrink=0.75)
            
            a4 = f.add_subplot(235)
            im4 = a4.imshow(np.abs(Dy) ,extent=[x[0], x[-1], y[0], y[-1]] , aspect='equal', origin='lower', vmin=0, vmax=vmax, cmap='jet')
            a4.annotate(r'$D_{y}(z)|/|\mathbf{D}_{0}|_{\rm max}$', xy=(0.9*x[0],0.7*y[-1]),horizontalalignment='left', verticalalignment='bottom', color='w', size=14)
            a4.set_xlabel(r'$x/w_0$')
            plt.colorbar(im4,location='top',shrink=0.75)

            a34 = f.add_subplot(236)
            im34 = a34.imshow(S_z/np.amax(S_z) ,extent=[x0[0], x0[-1], y0[0], y0[-1]] , aspect='equal', origin='lower', vmin=0, cmap='jet')
            a34.annotate(r'$\langle S_{z}(z)\rangle/\langle S_{z}(z)\rangle_{\rm max}$', xy=(0.9*x0[0],0.7*y0[-1]),horizontalalignment='left', verticalalignment='bottom', color='w', size=14)
            a34.set_xlabel(r'$x/w_0$')
            plt.colorbar(im34,location='top',shrink=0.75,label = r'$z=$'+str(round(z,2))+r'$w_0$')
                
#            plt.savefig('beam_prop_aniso.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()       
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[9,6])

canvas = gui.create_canvas(root,f)
canvas.draw() # for faster feedback to user on startup
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(6)
var_save = gui.create_stringvar_vector(6)
var_double = gui.create_doublevar_vector(4)

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"dielectric tensor element $\varepsilon_1=$",var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r"dielectric tensor element $\varepsilon_2=$",var_string[1],row)
row = gui.create_entry_with_latex(mainframe,r"dielectric tensor element $\varepsilon_3=$",var_string[2],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_checkbutton_with_latex(mainframe,r'propagation along optical axis','angles','oaprop',var_string[3],row)
row = gui.create_entry_with_latex(mainframe,r"azimuthal angle of propagation direction $\varphi/\pi=$",var_string[5],row)
row = gui.create_entry_with_latex(mainframe,r"polar angle of propagation direction $\theta/\pi=$",var_string[4],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_formula_with_latex(mainframe,r'$D_x =$',r'$D_0\left(\cos \Psi - \mathrm{i}\epsilon\sin \Psi \right)\exp\!\left(-\frac{x^2+y^2}{w_0^2}\right)$',row)
row = gui.create_formula_with_latex(mainframe,r'$D_y =$',r'$D_0\left(\sin \Psi + \mathrm{i}\epsilon\cos \Psi \right)\exp\!\left(-\frac{x^2+y^2}{w_0^2}\right)$',row)
row = gui.create_slider_with_latex(mainframe,r'beam width $w_0/\lambda=$',var_double[0],5,20,row,increment=.5)
row = gui.create_slider_with_latex(mainframe,r'propagation distance $z\,\delta_{\rm max}/w_0=$',var_double[1],0,2.5,row,increment=.25)
row = gui.create_slider_with_latex(mainframe,r'polarization angle $\Psi/\pi=$',var_double[2],0,0.5,row,increment=.01)
row = gui.create_slider_with_latex(mainframe,r'polarization ellipticity $\epsilon=$',var_double[3],-1,1,row,increment=.1)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)