import subprocess
import platform
import tkinter as Tk
from tkinter import messagebox as mbox
import numpy as np
import matplotlib as mpl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import ttk
import sympy as sp
from io import BytesIO
from PIL import Image, ImageTk
from sys import exit
import os
from IPython import get_ipython

def set_rcParams():
    mpl.rcParams['backend'] = 'tkagg'
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rc('text', usetex=True)
    if subprocess.run(["kpsewhich","cmbright.sty"], stdout=subprocess.PIPE).returncode == 0:
        mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')
    mpl.rcParams.update({'font.size': 10})
    mpl.rcParams['figure.dpi'] = 90
    
def create_canvas(root,f):
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.get_tk_widget().grid(column=1, row=1, sticky=(Tk.W, Tk.E))
    return canvas
    
def create_mainframe(root):
    mainframe = ttk.Frame(root, padding="3 3 12 12")
    mainframe.grid(column=2, row=1, sticky=(Tk.W, Tk.E))  
    style = ttk.Style(root)
    style.theme_use('clam')
    style.configure('.', background= 'lightgray')
    root.config(bg="lightgray")
    return mainframe

def latex2png(latex):
    if subprocess.run(["kpsewhich","cmbright.sty"], stdout=subprocess.PIPE).returncode == 0:
        preamble = "\\documentclass[10pt]{article}\n" \
               "\\usepackage{cmbright}\\usepackage{xcolor}\\pagestyle{empty}\\begin{document}\\pagecolor{lightgray}"
    else:
        preamble = "\\documentclass[10pt]{article}\n" \
               "\\usepackage{xcolor}\\pagestyle{empty}\\begin{document}\\pagecolor{lightgray}"
    obj = BytesIO()
    sp.preview(latex, viewer='BytesIO', output='png', outputbuffer=obj, preamble=preamble)
    obj.seek(0)
    img = Image.open(obj)
    img = img.convert("RGBA")
    datas = img.getdata()
    newData = []
    for item in datas:
        if item[0] == 191 and item[1] == 191 and item[2] == 191:
            newData.append((255, 255, 255, 0))
        else:
            newData.append(item)
    img.putdata(newData)
    labelimage = ImageTk.PhotoImage(img)
    return labelimage

def latexlabel(mainframe,latex):
    labelimage = latex2png(latex)
    label = ttk.Label(mainframe, image=labelimage)
    label.labelimage = labelimage
    return label    
    
def create_title(mainframe,text,row):
    ttk.Label(mainframe, text=text).grid(column=2, row=row, sticky=Tk.W, padx=5, pady=5)
    row=row+1
    return row

def create_launch_title(mainframe,text,column,row):
    ttk.Label(mainframe, text=text).grid(column=column, row=row, padx=5, pady=5)
    row=row+1
    return row

def create_description(mainframe,text,row):
    ttk.Label(mainframe, text=text).grid(column=1, row=row, sticky=Tk.W, padx=5, pady=5)
    row=row+1
    return row

def create_formula_with_latex(mainframe,latex1,latex2,row):
    labelimage1 = latex2png(latex1)
    label1 = ttk.Label(mainframe, image=labelimage1)
    label1.labelimage = labelimage1
    label1.grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    labelimage2 = latex2png(latex2)
    label2 = ttk.Label(mainframe, image=labelimage2)
    label2.labelimage = labelimage2
    label2.grid(column=2, row=row, sticky=Tk.W, padx=5, pady=5)
    row=row+1
    return row
    
def create_entry(mainframe,text,textvariable,row):
    ttk.Label(mainframe, text=text).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Entry(mainframe, width=7, textvariable=textvariable).grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    row=row+1
    return row

def create_entry_with_latex(mainframe,latex,textvariable,row):
    latexlabel(mainframe, latex).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Entry(mainframe, width=7, textvariable=textvariable).grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    row=row+1
    return row 

def create_double_entry(mainframe,text1,textvariable1,text2,textvariable2,row):
    Aframe1 = ttk.Frame(mainframe)
    Aframe1.grid(column=1, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    ttk.Label(Aframe1, text=text1).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Entry(Aframe1, width=7, textvariable=textvariable1).grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    Aframe2 = ttk.Frame(mainframe)
    Aframe2.grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    ttk.Label(Aframe2, text=text2).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Entry(Aframe2, width=7, textvariable=textvariable2).grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    row=row+1
    return row

def create_double_entry_with_latex(mainframe,latex1,textvariable1,latex2,textvariable2,row):
    Aframe1 = ttk.Frame(mainframe)
    Aframe1.grid(column=1, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    latexlabel(Aframe1, latex1).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Entry(Aframe1, width=7, textvariable=textvariable1).grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    Aframe2 = ttk.Frame(mainframe)
    Aframe2.grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    latexlabel(Aframe2, latex2).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Entry(Aframe2, width=7, textvariable=textvariable2).grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    row=row+1
    return row

def create_triple_entry(mainframe,text1,textvariable1,text2,textvariable2,text3,textvariable3,row):
    Aframe1 = ttk.Frame(mainframe)
    Aframe1.grid(column=1, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    ttk.Label(Aframe1, text=text1).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Entry(Aframe1, width=7, textvariable=textvariable1).grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    Aframe2 = ttk.Frame(mainframe)
    Aframe2.grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    ttk.Label(Aframe2, text=text2).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Entry(Aframe2, width=7, textvariable=textvariable2).grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    ttk.Label(Aframe2, text=text3).grid(column=3, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Entry(Aframe2, width=7, textvariable=textvariable3).grid(column=4, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    row=row+1
    return row

def create_entry_with_image(mainframe,image,textvariable,row):
    labelimage = Tk.PhotoImage(file=image)
    label = ttk.Label(mainframe, image=labelimage)
    label.labelimage = labelimage
    label.grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Entry(mainframe, width=7, textvariable=textvariable).grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    row=row+1
    return row
    
def create_spacer(mainframe,row):    
    ttk.Label(mainframe, text=" ").grid(column=1, row=row, padx=5, pady=5)
    row=row+1
    return row
    
def create_button(mainframe,text,command,row):     
    ttk.Button(mainframe, text=text, command=command).grid(column=2, row=row, sticky=Tk.W, padx=5, pady=5)
    row=row+1
    return row

def create_double_button(mainframe,text1,command1,text2,command2,row):     
    ttk.Button(mainframe, text=text1, command=command1).grid(column=1, row=row, sticky=Tk.W, padx=5, pady=5)
    ttk.Button(mainframe, text=text2, command=command2).grid(column=2, row=row, sticky=Tk.W, padx=5, pady=5)
    row=row+1
    return row

def check_python(root):
    try:
        python3_returncode = subprocess.run(['python3', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).returncode
    except:
        python3_returncode = -1
    try:
        python_returncode = subprocess.run(['python', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).returncode
        if "Python 3" not in subprocess.run(['python', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.decode('utf-8'):
            python_returncode = -1
    except:
        python_returncode = -1
    if python_returncode==0:
        python_string = "python"
    elif python3_returncode==0:
        python_string = "python3"
    else:
        mbox.showerror("Error", "Could not find Python 3")
        root.destroy()
    return python_string

def create_launch_button(mainframe,python_string,filename,column,row): 
    def command():
        subprocess.Popen([python_string,filename])
    if platform.system()=='Darwin':
        ttk.Button(mainframe, text=filename, command=command, width=17).grid(column=column, row=row, padx=5, pady=5)
    elif platform.system()=='Windows':
        ttk.Button(mainframe, text=filename, command=command, width=21).grid(column=column, row=row, padx=5, pady=5)
    else:
        ttk.Button(mainframe, text=filename, command=command, width=19).grid(column=column, row=row, padx=5, pady=5)
    row=row+1
    return row

def create_checkbutton(mainframe,text,offvalue,onvalue,variable,row):
    ttk.Checkbutton(mainframe, text=text, offvalue=offvalue, onvalue=onvalue, variable=variable).grid(column=2, row=row, sticky=Tk.W, padx=5, pady=5)
    row=row+1
    return row

def create_checkbutton_with_latex(mainframe,latex,offvalue,onvalue,variable,row):
    Aframe1 = ttk.Frame(mainframe)
    Aframe1.grid(column=1, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    latexlabel(Aframe1, latex).grid(column=2, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Checkbutton(Aframe1, text='', offvalue=offvalue, onvalue=onvalue, variable=variable).grid(column=1, row=row, sticky=Tk.W, padx=5, pady=5)
    row=row+1
    return row

def create_double_checkbutton(mainframe,text1,offvalue1,onvalue1,variable1,text2,offvalue2,onvalue2,variable2,row):
    Aframe1 = ttk.Frame(mainframe)
    Aframe1.grid(column=1, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    ttk.Checkbutton(Aframe1, text=text1, offvalue=offvalue1, onvalue=onvalue1, variable=variable1).grid(column=2, row=row, sticky=Tk.W, padx=5, pady=5)
    Aframe2 = ttk.Frame(mainframe)
    Aframe2.grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    ttk.Checkbutton(Aframe2, text=text2, offvalue=offvalue2, onvalue=onvalue2, variable=variable2).grid(column=2, row=row, sticky=Tk.W, padx=5, pady=5)
    row=row+1
    return row

def create_double_checkbutton_with_latex(mainframe,latex1,offvalue1,onvalue1,variable1,latex2,offvalue2,onvalue2,variable2,row):
    Aframe1 = ttk.Frame(mainframe)
    Aframe1.grid(column=1, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    latexlabel(Aframe1, latex1).grid(column=2, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Checkbutton(Aframe1, text='', offvalue=offvalue1, onvalue=onvalue1, variable=variable1).grid(column=1, row=row, sticky=Tk.W, padx=5, pady=5)
    Aframe2 = ttk.Frame(mainframe)
    Aframe2.grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    latexlabel(Aframe2, latex2).grid(column=2, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Checkbutton(Aframe2, text='', offvalue=offvalue2, onvalue=onvalue2, variable=variable2).grid(column=1, row=row, sticky=Tk.W, padx=5, pady=5)
    row=row+1
    return row

def create_radiobutton(mainframe,text,textvariable,N,row,optional_command='none'):
    ttk.Label(mainframe, text=text[0]).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    Vframe = ttk.Frame(mainframe)
    Vframe.grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    for index_row in range(int(np.ceil(N/2))):
        for index_col in range(np.minimum(2,N-index_row*2)):
            if optional_command=='none':
                ttk.Radiobutton(Vframe, text=text[index_row*2+index_col+1], variable=textvariable, value=text[index_row*2+index_col+1]).grid(column=2*index_col+2, row=1+index_row*2, sticky=(Tk.W, Tk.E), padx=5, pady=5)
            else:
                ttk.Radiobutton(Vframe, text=text[index_row*2+index_col+1], variable=textvariable, value=text[index_row*2+index_col+1], command=optional_command).grid(column=2*index_col+2, row=1+index_row*2, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    row=row+1
    return row  

def create_radiobutton_with_entry(mainframe,entrytext,entrytextvariable,text,textvariable,N,row,optional_command='none'):
    Aframe1 = ttk.Frame(mainframe)
    Aframe1.grid(column=1, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    ttk.Label(Aframe1, text=entrytext).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Entry(Aframe1, width=7, textvariable=entrytextvariable).grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    Vframe = ttk.Frame(mainframe)
    Vframe.grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    for index_row in range(int(np.ceil(N/2))):
        for index_col in range(np.minimum(2,N-index_row*2)):
            if optional_command=='none':
                ttk.Radiobutton(Vframe, text=text[index_row*2+index_col], variable=textvariable, value=text[index_row*2+index_col]).grid(column=2*index_col+2, row=1+index_row*2, sticky=(Tk.W, Tk.E), padx=5, pady=5)
            else:
                ttk.Radiobutton(Vframe, text=text[index_row*2+index_col], variable=textvariable, value=text[index_row*2+index_col], command=optional_command).grid(column=2*index_col+2, row=1+index_row*2, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    row=row+1
    return row  

def create_radiobutton_single_column(mainframe,text,textvariable,N,row,optional_command='none'):
    ttk.Label(mainframe, text=text[0]).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    Vframe = ttk.Frame(mainframe)
    Vframe.grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    for index_row in range(N):
        if optional_command=='none':
            ttk.Radiobutton(Vframe, text=text[index_row+1], variable=textvariable, value=text[index_row+1]).grid(column=2, row=1+index_row*2, sticky=(Tk.W, Tk.E), padx=5, pady=5)
        else:
            ttk.Radiobutton(Vframe, text=text[index_row+1], variable=textvariable, value=text[index_row+1], command=optional_command).grid(column=2, row=1+index_row*2, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    row=row+1
    return row  

def create_radiobutton_single_column_with_latex(mainframe,latex,text,textvariable,N,row,optional_command='none'):
    latexlabel(mainframe, latex[0]).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    Vframe = ttk.Frame(mainframe)
    Vframe.grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    for index_row in range(N):
        if optional_command=='none':
            latexlabel(Vframe, latex[index_row+1]).grid(column=2, row=1+index_row*2, sticky=Tk.W, padx=5, pady=5)
            ttk.Radiobutton(Vframe, text='', variable=textvariable, value=text[index_row]).grid(column=1, row=1+index_row*2, sticky=(Tk.W, Tk.E), padx=5, pady=5)
        else:
            latexlabel(Vframe, latex[index_row+1]).grid(column=2, row=1+index_row*2, sticky=Tk.W, padx=5, pady=5)
            ttk.Radiobutton(Vframe, text='', variable=textvariable, value=text[index_row], command=optional_command).grid(column=1, row=1+index_row*2, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    row=row+1
    return row  

def create_slider(mainframe,text,variable,from_value,to_value,row,optional_command='none'):
    ttk.Label(mainframe, text=text).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Label(mainframe, text=str(round(variable.get(),4))).grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    def showvalue(value):
        ttk.Label(mainframe, text=str(round(float(value),4))).grid(column=2, row=row-2, sticky=(Tk.W, Tk.E), padx=5, pady=5)
        if optional_command!='none':
            optional_command(float(value))
    ttk.Scale(mainframe, variable=variable, from_=from_value, to=to_value, command=showvalue, cursor='crosshair').grid(column=2, row=row+1, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    row=row+2
    return row  

def create_slider_with_latex(mainframe,latex,variable,from_value,to_value,row,optional_command='none'):
    latexlabel(mainframe, latex).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Label(mainframe, text=str(round(variable.get(),4))).grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    def showvalue(value):
        ttk.Label(mainframe, text=str(round(float(value),4))).grid(column=2, row=row-2, sticky=(Tk.W, Tk.E), padx=5, pady=5)
        if optional_command!='none':
            optional_command(float(value))
    ttk.Scale(mainframe, variable=variable, from_=from_value, to=to_value, command=showvalue, cursor='crosshair').grid(column=2, row=row+1, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    row=row+2
    return row  

def create_intslider_with_latex(mainframe,latex,variable,from_value,to_value,row,optional_command='none'):
    latexlabel(mainframe, latex).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Label(mainframe, text='%d' % variable.get()).grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    def showvalue(value):
        ttk.Label(mainframe, text='%d' % float(value)).grid(column=2, row=row-2, sticky=(Tk.W, Tk.E), padx=5, pady=5)
        if optional_command!='none':
            optional_command(float(value))
    ttk.Scale(mainframe, variable=variable, from_=from_value, to=to_value, command=showvalue, cursor='crosshair').grid(column=2, row=row+1, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    row=row+2
    return row  

def create_logslider_with_latex(mainframe,latex,variable,from_value,to_value,row,optional_command='none'):
    latexlabel(mainframe, latex).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Label(mainframe, text=str(round(np.exp(variable.get()),4))).grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    def showvalue(value):
        ttk.Label(mainframe, text=str(round(np.exp(float(value)),4))).grid(column=2, row=row-2, sticky=(Tk.W, Tk.E), padx=5, pady=5)
        if optional_command!='none':
            optional_command(np.exp(float(value)))
    ttk.Scale(mainframe, variable=variable, from_=np.log(from_value), to=np.log(to_value), command=showvalue, cursor='crosshair').grid(column=2, row=row+1, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    row=row+2
    return row  
    
def create_entry_matrix_3x3_sym(mainframe,text,textvariable,row): 
    ttk.Label(mainframe, text=text).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    Aframe = ttk.Frame(mainframe)
    Aframe.grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    ttk.Entry(Aframe, width=3, textvariable=textvariable[0]).grid(column=1, row=1, sticky=(Tk.W, Tk.E))
    ttk.Entry(Aframe, width=3, textvariable=textvariable[1]).grid(column=2, row=1, sticky=(Tk.W, Tk.E))
    ttk.Entry(Aframe, width=3, textvariable=textvariable[2]).grid(column=3, row=1, sticky=(Tk.W, Tk.E))
    ttk.Label(Aframe, textvariable=textvariable[3]).grid(column=1, row=2, sticky=(Tk.W, Tk.E))
    ttk.Entry(Aframe, width=3, textvariable=textvariable[4]).grid(column=2, row=2, sticky=(Tk.W, Tk.E))
    ttk.Entry(Aframe, width=3, textvariable=textvariable[5]).grid(column=3, row=2, sticky=(Tk.W, Tk.E))
    ttk.Label(Aframe, textvariable=textvariable[6]).grid(column=1, row=3, sticky=(Tk.W, Tk.E))
    ttk.Label(Aframe, textvariable=textvariable[7]).grid(column=2, row=3, sticky=(Tk.W, Tk.E))
    ttk.Entry(Aframe, width=3, textvariable=textvariable[8]).grid(column=3, row=3, sticky=(Tk.W, Tk.E))
    for child in Aframe.winfo_children(): child.grid_configure(padx=2, pady=2)   
    row=row+1
    return row

def create_label(mainframe,text,textvariable,row):
    ttk.Label(mainframe, text=text).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Label(mainframe, textvariable=textvariable).grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    row=row+1
    return row 

def create_label_with_image(mainframe,image,textvariable,row):
    labelimage = Tk.PhotoImage(file=image)
    label = ttk.Label(mainframe, image=labelimage)
    label.labelimage = labelimage
    label.grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Label(mainframe, textvariable=textvariable).grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    row=row+1
    return row 

def create_label_with_latex(mainframe,latex,textvariable,row):
    latexlabel(mainframe, latex).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    ttk.Label(mainframe, textvariable=textvariable).grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    row=row+1
    return row 

def create_label_vector_3(mainframe,text,textvariable,row):
    ttk.Label(mainframe, text=text).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    Vframe = ttk.Frame(mainframe)
    Vframe.grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    ttk.Label(Vframe, text="(").grid(column=1, row=1, sticky=(Tk.W, Tk.E))
    ttk.Label(Vframe, textvariable=textvariable[0]).grid(column=2, row=1, sticky=(Tk.W, Tk.E))
    ttk.Label(Vframe, text=",").grid(column=3, row=1, sticky=(Tk.W, Tk.E))
    ttk.Label(Vframe, textvariable=textvariable[1]).grid(column=4, row=1, sticky=(Tk.W, Tk.E))
    ttk.Label(Vframe, text=",").grid(column=5, row=1, sticky=(Tk.W, Tk.E))
    ttk.Label(Vframe, textvariable=textvariable[2]).grid(column=6, row=1, sticky=(Tk.W, Tk.E))
    ttk.Label(Vframe, text=")").grid(column=7, row=1, sticky=(Tk.W, Tk.E))
    row=row+1
    return row   
    
def create_label_vector(mainframe,text,textvariable,N,row):
    ttk.Label(mainframe, text=text).grid(column=1, row=row, sticky=Tk.E, padx=5, pady=5)
    Vframe = ttk.Frame(mainframe)
    Vframe.grid(column=2, row=row, sticky=(Tk.W, Tk.E), padx=5, pady=5)
    ttk.Label(Vframe, text="(").grid(column=1, row=1, sticky=(Tk.W, Tk.E))
    for index in range(N):
        ttk.Label(Vframe, textvariable=textvariable[index]).grid(column=2*index+2, row=1, sticky=(Tk.W, Tk.E))
        if index < N-1: ttk.Label(Vframe, text=",").grid(column=2*index+3, row=1, sticky=(Tk.W, Tk.E))
    ttk.Label(Vframe, text=")").grid(column=2*N+1, row=1, sticky=(Tk.W, Tk.E))
    row=row+1
    return row  

def change_cursor(root,cursor):
    root.config(cursor=cursor)
    root.update()
    
def input_error(message,optional_command='none'):
    mbox.showerror("Error", message)
    if optional_command!='none':
        optional_command()
    pass

def input_warning(message,optional_command='none'):
    mbox.showwarning("Warning", message)
    if optional_command!='none':
        optional_command()
    pass

def create_stringvar_vector(N):
    var_string = []
    for index in range(N): var_string.append(Tk.StringVar())
    return var_string

def create_doublevar_vector(N):
    var_double = []
    for index in range(N): var_double.append(Tk.DoubleVar())
    return var_double

def copy_stringvar_vector(var1_string,var2_string):
    for index in range(len(var1_string)): var2_string[index].set(var1_string[index].get())

def mainloop_safe_for_mac(root):
    while True:
        try:
            if get_ipython() is None or 'SPY_PYTHONPATH' not in os.environ:
                root.protocol("WM_DELETE_WINDOW", exit)
            root.mainloop()
            break
        except UnicodeDecodeError:
            pass

def read_image(filename):
    img = ImageTk.PhotoImage(Image.open(filename))
    return img
    
def show_image(root,title,img):
    root.title(title)
    dimensions = "%dx%d" % (img.width()+10, img.height()+10)
    root.geometry(dimensions)
    frame = Tk.Frame(root)
    frame.pack()
    frame.place(anchor='center', relx=0.5, rely=0.5)
    label = Tk.Label(frame, image = img)
    label.configure(background= 'white')
    label.pack()

def show_manual(filename,title):
    top = Tk.Toplevel()
    style = ttk.Style(top)
    style.theme_use('clam')
    top.config(bg="white")
    img = read_image(filename)
    show_image(top,title,img)
    top.mainloop()  