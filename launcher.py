import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("interactive figures launcher")

mainframe = gui.create_mainframe(root)

row = 1
column = 1
row = gui.create_launch_button(mainframe,'film_waveguide.py',column,row)
row = gui.create_launch_button(mainframe,'lossy_mode.py',column,row)
row = gui.create_launch_button(mainframe,'leaky_mode.py',column,row)
row = gui.create_launch_button(mainframe,'spp.py',column,row)
row = gui.create_launch_button(mainframe,'spp_drude.py',column,row)
row = gui.create_launch_button(mainframe,'strip_waveguide_finite_difference.py',column,row)
row = gui.create_launch_button(mainframe,'strip_waveguide_effective_index.py',column,row)
row = gui.create_launch_button(mainframe,'fiber_wga.py',column,row)
row = gui.create_launch_button(mainframe,'waveguide_array.py',column,row)
row = gui.create_launch_button(mainframe,'frontface_coupling.py',column,row)
row = gui.create_launch_button(mainframe,'prism_coupler.py',column,row)

row = 1
column = 2
row = gui.create_launch_button(mainframe,'taylor_series.py',column,row)

row = 1
column = 3
row = gui.create_launch_button(mainframe,'mbe.py',column,row)

gui.mainloop_safe_for_mac(root)