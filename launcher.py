import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("Interactive Figure Launcher")

mainframe = gui.create_mainframe(root)

python_string = gui.check_python(root)

row = 1
column = 1
row = gui.create_launch_title(mainframe,'Chapter 1',column,row)
row = gui.create_launch_button(mainframe,python_string,'epsilon.py',column,row)

row = 1
column = 2
row = gui.create_launch_title(mainframe,'Chapter 2',column,row)
row = gui.create_launch_button(mainframe,python_string,'index_ellipsoid.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'normal_surface.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'normal_surface_ua.py',column,row)

row = 1
column = 3
row = gui.create_launch_title(mainframe,'Chapter 3',column,row)
row = gui.create_launch_button(mainframe,python_string,'scalar_beam_prop.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'beam_prop_aniso.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'diffraction_fps.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'image_2f.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'image_4f.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'filter_4f.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'scalar_pulse_prop.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'gauss_bullet_prop.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'gauss_bullet_foc.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'gauss_bullet_pft_wfr.py',column,row)

row = 1
column = 4
row = gui.create_launch_title(mainframe,'Chapter 4',column,row)
row = gui.create_launch_button(mainframe,python_string,'rt_interface.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'r_interface_beam.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'r_interface_pulse.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'rt_stack_kx.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'rt_stack_lambda.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'graded_index.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'3D_Bloch.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'2D_Bloch_lambda.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'2D_Bloch_kx.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'Bloch_interface.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'Bloch_stack.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'Bloch_stack_lambda.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'FP_kx.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'FP_lambda.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'FP_image.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'FP_Airy.py',column,row)

row = 1
column = 5
row = gui.create_launch_title(mainframe,'Chapter 5',column,row)
row = gui.create_launch_button(mainframe,python_string,'film_waveguide.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'lossy_mode.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'leaky_mode.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'spp.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'spp_drude.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'strip_waveguide_fd.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'strip_waveguide_ei.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'fiber_wga.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'waveguide_array.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'prism_coupling.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'end_face_coupling.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'taper.py',column,row)

row = 1
column = 6
row = gui.create_launch_title(mainframe,'Chapter 6',column,row)
row = gui.create_launch_button(mainframe,python_string,'obe.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'3wm_sfg.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'3wm_shg.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'3wm_dfg.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'3wm.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'shg_pulse.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'soliton.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'mi.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'4wm.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'srs.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'mbe.py',column,row)

row = 1
column = 7
row = gui.create_launch_title(mainframe,'Chapter 7',column,row)
row = gui.create_launch_button(mainframe,python_string,'taylor_series.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'fourier_series.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'fft_gauss.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'bpm_parax.py',column,row)
row = gui.create_launch_button(mainframe,python_string,'bpm_nls.py',column,row)

root.geometry("+0+0")

root.mainloop()