This repository contains supplemental interactive figures to the textbook 
[***Lectures in Theoretical Optics***](https://link.springer.com/book/9783662726334) 
by Falk Lederer and Stefan Skupin. The figures are written in [python 3](https://www.python.org/), and require 
[os](https://docs.python.org/3/library/os.html),
[matplotlib](https://matplotlib.org/),
[numpy](https://numpy.org/),
[scipy](https://scipy.org/),
[sys](https://docs.python.org/3/library/sys.html),
[tkinter](https://docs.python.org/3/library/tkinter.html), and a working [LaTeX](https://www.latex-project.org/) installation. We commend [spyder](https://www.spyder-ide.org/) to run and interact with the figures.

## Installation

**Linux**

Install the required packages through the package manager of your Linux distribution. Most distributions come with python and latex installed by default. For example, to install tkinter on [Ubuntu](https://ubuntu.com/), use `sudo apt update` then `sudo apt install python3-tk`. Spyder is available as a Ubuntu package as well, run `sudo apt install spyder`.

**MacOS**

We recommend to install the [spyder](https://www.spyder-ide.org/) standalone app for macOS and use the internal python interpreter `Spyder 6 --> Preferences... --> Python interpreter --> internal (same used by Spyder)`. It brings all necessary python dependencies by default. We also recommend to install [MacTex](https://www.tug.org/mactex/), which includes a complete TeX system with LaTeX. To make spyder aware of the LateX installation, you have to add its path in `Spyder 6 --> Preferences... --> IPython console --> Startup --> Run code --> Lines:`. In case of a standard MacTex installation, you have to add `import os; os.environment[PATH] += ':/Library/Tex/texbin'`.

**Windows**

We recommend to install the [spyder](https://www.spyder-ide.org/) standalone app for Windows and use the internal python interpreter `Tools --> Preferences... --> Python interpreter --> internal (same used by Spyder)`. It brings all necessary python dependencies by default. We also recommend to install [MikTex](https://miktex.org/), which includes a complete TeX system with LaTeX.

## Usage

## F.A.Q.