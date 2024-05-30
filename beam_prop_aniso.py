import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import bpm_stuff as bpm
import aniso_stuff as ani

gui.set_rcParams()
root = Tk.Tk()
root.title("Beam Propagation in Homogeneous and Anisotropic Media")

def initialize():
    var_string[0].set("2")   # epsilon1
    var_string[1].set("3")   # epsilon2
    var_string[2].set("4")   # epsilon3
    var_string[3].set("D")   # show D or u
    var_double[0].set(0.3)  # theta0/pi
    var_double[1].set(0.4)  # phi0/pi    
    var_double[2].set(10)    # w_0/\lambda
    var_double[3].set(0.5)   # z/L_F  
    var_double[4].set(0)     # polarization angle
    var_double[5].set(0)     # polarization ellipticity
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate() # because sliders may have changed
    
def calculate():
    gui.change_cursor(root,"trek")
    try:     
        epsilon = np.array([float(var_string[0].get()),float(var_string[1].get()),float(var_string[2].get())])
        theta0 = var_double[0].get()*np.pi
        phi0 = var_double[1].get()*np.pi
        theta_view = 0.3*180
        phi_view = 0.8*180
        kv = 2*np.pi*var_double[2].get() # in units of 1/w_0
        lambda0 = 1/var_double[2].get() # in units of w_0
        psi = var_double[4].get()*np.pi
        ellipticity = var_double[5].get()
        show = var_string[3].get()
        
        if epsilon[0] <= 0 or epsilon[1] <= 0  or epsilon[2] <= 0: 
            gui.input_error("Tensor elements have to be positive. Re-initializing ...",reinitialize)
        else:

            f.clf()
            
            # compute thetaoa optical axis
#            thetaoa = np.arcsin(np.sqrt(epsilon[2]*(epsilon[1]-epsilon[0])/(epsilon[1]*(epsilon[2]-epsilon[0]))))
#            theta0 = thetaoa
#            var_double[0].set(theta0/np.pi)
#            phi0 = 0
#            var_double[1].set(phi0/np.pi)
            
            # fix lab coordinate vectors ex,ey,ez, as Da,Db,uk
            ez = ani.uk(theta0,phi0)
            ex,ey = ani.D(ez,epsilon) 
            na,nb = ani.dr(ez,epsilon)
            
#            x_calc,y_calc,Nx,Ny,kx_calc,ky_calc,Nkx,Nky,Nx0,Ny0 = bpm.init_2D_grid(kv*(na[0]+nb[0])/2,z)

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
            z = var_double[3].get()/deltamax # in units of w_0
            
            Dax,Day,Daz,Dbx,Dby,Dbz = ani.D_proj_mesh(KX,KY,KZa,KZb,ex,ey,ez,epsilon)
            
            D0_x = (np.cos(psi)-1j*ellipticity*np.sin(psi))*bpm.init_2D_beam(x_calc,y_calc)
            D0_y = (np.sin(psi)+1j*ellipticity*np.cos(psi))*bpm.init_2D_beam(x_calc,y_calc)
            
            FTD0_x = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(D0_x)))
            FTD0_y = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(D0_y)))

            U0_a = (Dby*FTD0_x-Dbx*FTD0_y)/(Dax*Dby-Day*Dbx)
            U0_b = (Dax*FTD0_y-Day*FTD0_x)/(Dax*Dby-Day*Dbx)
            
            u0_a = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(U0_a)))
            u0_b = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(U0_b)))
            
            prop_a,prop_b = bpm.init_prop_2D_aniso(KZa*kv,KZb*kv,z)
            
            u_a,U_a,U0_a = bpm.propagation_2D(u0_a,prop_a)
            u_b,U_b,U0_b = bpm.propagation_2D(u0_b,prop_b)
            
            FTDa_x = U_a*Dax
            FTDa_y = U_a*Day
            FTDa_z = U_a*Daz
            
            FTDb_x = U_b*Dbx
            FTDb_y = U_b*Dby
            FTDb_z = U_b*Dbz
            
            FTD_x = FTDa_x + FTDb_x
            FTD_y = FTDa_y + FTDb_y
            FTD_z = FTDa_z + FTDb_z
            
            FTHa_x = KY*FTDa_z - KZa*FTDa_y
            FTHa_y = KZa*FTDa_x - KX*FTDa_z
            
            FTHb_x = KY*FTDb_z - KZa*FTDb_y
            FTHb_y = KZa*FTDa_x - KX*FTDb_z
            
            FTH_x = FTHa_x + FTHb_x
            FTH_y = FTHa_y + FTHb_y
            
            FTE_x = FTD_x*np.dot(ex/epsilon,ex) + FTD_y*np.dot(ey/epsilon,ex) + FTD_z*np.dot(ez.flatten()/epsilon,ex)
            FTE_y = FTD_x*np.dot(ex/epsilon,ey) + FTD_y*np.dot(ey/epsilon,ey) + FTD_z*np.dot(ez.flatten()/epsilon,ey)
            
            D_x = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(FTD_x)))
            D_y = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(FTD_y)))  
            
            H_x = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(FTH_x)))
            H_y = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(FTH_y)))             
            
            E_x = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(FTE_x)))
            E_y = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(FTE_y)))  
            
            S_z = np.real(E_x*np.conjugate(H_y) - E_y*np.conjugate(H_x))
            
            u0a,x0,y0 = bpm.reshaping_2D(u0_a,x_calc,y_calc,Nx0,Ny0)
            u0b,x0,y0 = bpm.reshaping_2D(u0_b,x_calc,y_calc,Nx0,Ny0)
            ua,x,y = bpm.reshaping_2D(u_a,x_calc,y_calc,Nx,Ny)
            ub,x,y = bpm.reshaping_2D(u_b,x_calc,y_calc,Nx,Ny)
            
            D0x,x0,y0 = bpm.reshaping_2D(D0_x,x_calc,y_calc,Nx0,Ny0)
            D0y,x0,y0 = bpm.reshaping_2D(D0_y,x_calc,y_calc,Nx0,Ny0)
            Dx,x,y = bpm.reshaping_2D(D_x,x_calc,y_calc,Nx,Ny)
            Dy,x,y = bpm.reshaping_2D(D_y,x_calc,y_calc,Nx,Ny)
            
            Sz,x,y = bpm.reshaping_2D(S_z,x_calc,y_calc,Nx,Ny)
            
            vmax = max(np.amax(np.abs(D0x)),np.amax(np.abs(D0y)))
            
            a1 = f.add_subplot(231)
            if show == 'D':
                im1 = a1.imshow(np.abs(D0x) ,extent=[x0[0], x0[-1], y0[0], y0[-1]] , aspect='equal', origin='lower', vmin=0, vmax=vmax, cmap='jet')
                a1.annotate(r'$|D_{0x}|/|\mathbf{D}_{0}|_{\rm max}$', xy=(0.9*x0[0],0.8*y0[-1]),horizontalalignment='left', verticalalignment='bottom', color='w')
            else:
                im1 = a1.imshow(np.abs(u0a) ,extent=[x0[0], x0[-1], y0[0], y0[-1]] , aspect='equal', origin='lower', vmin=0, cmap='jet')
                a1.annotate(r'$|u_{0}^a|/|\mathbf{D}_{0}|_{\rm max}$', xy=(0.9*x0[0],0.8*y0[-1]),horizontalalignment='left', verticalalignment='bottom', color='w')
            a1.set_ylabel(r'$y/w_0$')
            plt.colorbar(im1,location='top',shrink=0.75)
            
            a2 = f.add_subplot(232)
            if show == 'D':
                im2 = a2.imshow(np.abs(D0y) ,extent=[x0[0], x0[-1], y0[0], y0[-1]] , aspect='equal', origin='lower', vmin=0, vmax=vmax, cmap='jet')
                a2.annotate(r'$|D_{0y}|/|\mathbf{D}_{0}|_{\rm max}$', xy=(0.9*x0[0],0.8*y0[-1]),horizontalalignment='left', verticalalignment='bottom', color='w')
            else:
                im2 = a2.imshow(np.abs(u0b) ,extent=[x0[0], x0[-1], y0[0], y0[-1]] , aspect='equal', origin='lower', vmin=0, cmap='jet')
                a2.annotate(r'$|u_{0}^b|/|\mathbf{D}_{0}|_{\rm max}$', xy=(0.9*x0[0],0.8*y0[-1]),horizontalalignment='left', verticalalignment='bottom', color='w')                
            plt.colorbar(im2,location='top',shrink=0.75)
            
            a12 = f.add_subplot(233, projection='3d')
            a12.view_init(azim=phi_view, elev=90-theta_view)
            ani.plot_ns(a12,theta0,phi0,epsilon,'no_show','no_show',True)
            
            vmax = max(np.amax(np.abs(Dx)),np.amax(np.abs(Dy)))
            
            a3 = f.add_subplot(234)
            if show == 'D':
                im3 = a3.imshow(np.abs(Dx) ,extent=[x[0], x[-1], y[0], y[-1]] , aspect='equal', origin='lower', vmin=0, vmax=vmax, cmap='jet')
                a3.annotate(r'$|D_{x}(z)|/|\mathbf{D}_{0}|_{\rm max}$', xy=(0.9*x[0],0.8*y[-1]),horizontalalignment='left', verticalalignment='bottom', color='w')
            else:
                im3 = a3.imshow(np.abs(ua) ,extent=[x[0], x[-1], y[0], y[-1]] , aspect='equal', origin='lower', vmin=0, cmap='jet')
                a3.annotate(r'$|u^a(z)|/|\mathbf{D}_{0}|_{\rm max}$', xy=(0.9*x[0],0.8*y[-1]),horizontalalignment='left', verticalalignment='bottom', color='w')
            a3.set_xlabel(r'$x/w_0$')
            a3.set_ylabel(r'$y/w_0$')
            plt.colorbar(im3,location='top',shrink=0.75)
            
            a4 = f.add_subplot(235)
            if show == 'D':
                im4 = a4.imshow(np.abs(Dy) ,extent=[x[0], x[-1], y[0], y[-1]] , aspect='equal', origin='lower', vmin=0, vmax=vmax, cmap='jet')
                a4.annotate(r'$D_{y}(z)|/|\mathbf{D}_{0}|_{\rm max}$', xy=(0.9*x[0],0.8*y[-1]),horizontalalignment='left', verticalalignment='bottom', color='w')
            else:
                im4 = a4.imshow(np.abs(ub) ,extent=[x[0], x[-1], y[0], y[-1]] , aspect='equal', origin='lower', vmin=0, cmap='jet')
                a4.annotate(r'$u^b(z)|/|\mathbf{D}_{0}|_{\rm max}$', xy=(0.9*x[0],0.8*y[-1]),horizontalalignment='left', verticalalignment='bottom', color='w')
            a4.set_xlabel(r'$x/w_0$')
            plt.colorbar(im4,location='top',shrink=0.75)

            a34 = f.add_subplot(236)
            im34 = a34.imshow(S_z/np.amax(S_z) ,extent=[x0[0], x0[-1], y0[0], y0[-1]] , aspect='equal', origin='lower', vmin=0, cmap='jet')
            a34.annotate(r'$\langle S_{z}(z)\rangle/\langle S_{z}(z)\rangle_{\rm max}$', xy=(0.9*x0[0],0.8*y0[-1]),horizontalalignment='left', verticalalignment='bottom', color='w')
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

var_string = gui.create_stringvar_vector(4)
var_save = gui.create_stringvar_vector(4)
var_double = gui.create_doublevar_vector(6)

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_1=$",var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_2=$",var_string[1],row)
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_3=$",var_string[2],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r'Azimuth of propagation direction $\varphi/\pi=$',var_double[1],0,2,row)
row = gui.create_slider_with_latex(mainframe,r'Elevation of propagation direction $\theta/\pi=$',var_double[0],0,1,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_formula_with_latex(mainframe,r'$D_x =$',r'$D_0\left(\cos \Psi - \mathrm{i}\epsilon\sin \Psi \right)\exp\!\left(-\frac{x^2+y^2}{w_0^2}\right)$',row)
row = gui.create_formula_with_latex(mainframe,r'$D_y =$',r'$D_0\left(\sin \Psi + \mathrm{i}\epsilon\cos \Psi \right)\exp\!\left(-\frac{x^2+y^2}{w_0^2}\right)$',row)
row = gui.create_slider_with_latex(mainframe,r'Beam width $w_0/\lambda=$',var_double[2],5,20,row)
row = gui.create_slider_with_latex(mainframe,r'Propagation distance $z\delta_{\rm max}/(w_0)=$',var_double[3],0,2,row)
row = gui.create_slider_with_latex(mainframe,r'Polarization angle $\Psi/\pi=$',var_double[4],0,0.5,row)
row = gui.create_slider_with_latex(mainframe,r'Polarization ellipticity $\epsilon=$',var_double[5],-1,1,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_radiobutton(mainframe,['show:','D','u'],var_string[3],2,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)