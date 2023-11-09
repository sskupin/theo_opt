import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("interactive figure launcher")

mainframe = gui.create_mainframe(root)

row = 1
column = 1
row = gui.create_launch_title(mainframe,'Chapter 1',column,row)
row = gui.create_launch_button(mainframe,'epsilon.py',column,row)

row = 1
column = 2
row = gui.create_launch_title(mainframe,'Chapter 2',column,row)
row = gui.create_launch_button(mainframe,'normal_surface_ua.py',column,row)
row = gui.create_launch_button(mainframe,'normal_surface.py',column,row)

row = 1
column = 3
row = gui.create_launch_title(mainframe,'Chapter 3',column,row)
row = gui.create_launch_button(mainframe,'fps.py',column,row)
row = gui.create_launch_button(mainframe,'2f4f.py',column,row)
row = gui.create_launch_button(mainframe,'fps_4f_mask.py',column,row)

row = 1
column = 4
row = gui.create_launch_title(mainframe,'Chapter 4',column,row)
row = gui.create_launch_button(mainframe,'interface.py',column,row)
row = gui.create_launch_button(mainframe,'stack.py',column,row)
row = gui.create_launch_button(mainframe,'stack_lambda.py',column,row)
row = gui.create_launch_button(mainframe,'graded_index.py',column,row)
row = gui.create_launch_button(mainframe,'3D_Bloch.py',column,row)
row = gui.create_launch_button(mainframe,'2D_Bloch_lambda.py',column,row)
row = gui.create_launch_button(mainframe,'2D_Bloch_kx.py',column,row)
row = gui.create_launch_button(mainframe,'Bloch_interface.py',column,row)
row = gui.create_launch_button(mainframe,'Bloch_stack.py',column,row)
row = gui.create_launch_button(mainframe,'Bloch_stack_lambda.py',column,row)
row = gui.create_launch_button(mainframe,'FP_kx.py',column,row)
row = gui.create_launch_button(mainframe,'FP_image.py',column,row)
row = gui.create_launch_button(mainframe,'FP_Airy.py',column,row)

row = 1
column = 5
row = gui.create_launch_title(mainframe,'Chapter 5',column,row)
row = gui.create_launch_button(mainframe,'film_waveguide.py',column,row)
row = gui.create_launch_button(mainframe,'lossy_mode.py',column,row)
row = gui.create_launch_button(mainframe,'leaky_mode.py',column,row)
row = gui.create_launch_button(mainframe,'spp.py',column,row)
row = gui.create_launch_button(mainframe,'spp_drude.py',column,row)
row = gui.create_launch_button(mainframe,'strip_waveguide_fd.py',column,row)
row = gui.create_launch_button(mainframe,'strip_waveguide_ei.py',column,row)
row = gui.create_launch_button(mainframe,'fiber_wga.py',column,row)
row = gui.create_launch_button(mainframe,'waveguide_array.py',column,row)
row = gui.create_launch_button(mainframe,'end_face_coupling.py',column,row)
row = gui.create_launch_button(mainframe,'prism_coupling.py',column,row)
row = gui.create_launch_button(mainframe,'taper.py',column,row)

row = 1
column = 6
row = gui.create_launch_title(mainframe,'Chapter 6',column,row)
row = gui.create_launch_button(mainframe,'obe.py',column,row)
row = gui.create_launch_button(mainframe,'3wm_sfg.py',column,row)
row = gui.create_launch_button(mainframe,'3wm_shg.py',column,row)
row = gui.create_launch_button(mainframe,'3wm_dfg.py',column,row)
row = gui.create_launch_button(mainframe,'3wm.py',column,row)
row = gui.create_launch_button(mainframe,'shg_pulse.py',column,row)
row = gui.create_launch_button(mainframe,'soliton.py',column,row)
row = gui.create_launch_button(mainframe,'mi.py',column,row)
row = gui.create_launch_button(mainframe,'4wm.py',column,row)
row = gui.create_launch_button(mainframe,'srs.py',column,row)
row = gui.create_launch_button(mainframe,'mbe.py',column,row)

row = 1
column = 7
row = gui.create_launch_title(mainframe,'Chapter 7',column,row)
row = gui.create_launch_button(mainframe,'taylor_series.py',column,row)
row = gui.create_launch_button(mainframe,'fourier_series.py',column,row)
row = gui.create_launch_button(mainframe,'fft_gauss.py',column,row)
row = gui.create_launch_button(mainframe,'bpm_parax.py',column,row)
row = gui.create_launch_button(mainframe,'bpm_nls.py',column,row)

gui.mainloop_safe_for_mac(root)