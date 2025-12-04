import numpy as np
import scipy.optimize as spo
import scipy.signal as sps
import gui_stuff as gui

def nnarctan(argument): # define non-negative arctan function 
    nnarctan = np.arctan(argument)
    if nnarctan < 0 and np.imag(argument) == 0:
        nnarctan = nnarctan + np.pi
        
    return nnarctan

def number_of_modes(epsilon_f,epsilon_s,epsilon_c,d,polarization): 
    number_of_modes = np.int_(0)
    if epsilon_f > 0:
        if epsilon_f > epsilon_s:
            a = (epsilon_s-epsilon_c)/(epsilon_f-epsilon_s) 
            fs = np.sqrt(epsilon_f-epsilon_s)
            if polarization == 'TE':
                if epsilon_s < 0:
                    if np.sqrt(epsilon_c*epsilon_s) == epsilon_f:
                        number_of_modes = np.ceil(2*d*np.sqrt(epsilon_f)-1/2)
                    else:
                        number_of_modes = np.ceil(2*d*np.sqrt(epsilon_f)-nnarctan(np.sqrt(epsilon_f)*(np.sqrt(-epsilon_c)+np.sqrt(-epsilon_s))/(epsilon_f-np.sqrt(epsilon_c*epsilon_s)))/np.pi)
                else:
                    number_of_modes = np.ceil(2*d*fs-nnarctan(np.sqrt(a))/np.pi)
            elif polarization == 'TM':
                if epsilon_s < 0:
                    if np.sqrt(epsilon_c*epsilon_s) == epsilon_f:
                        number_of_modes = np.ceil(2*d*np.sqrt(epsilon_f)-1/2)
                    else:
                        number_of_modes = np.ceil(2*d*np.sqrt(epsilon_f)-nnarctan(epsilon_f*np.sqrt(epsilon_f)*(epsilon_s*np.sqrt(-epsilon_c)+epsilon_c*np.sqrt(-epsilon_s))/(epsilon_s*epsilon_c*epsilon_f-epsilon_f**2*np.sqrt(epsilon_c*epsilon_s)))/np.pi)
                    number_of_modes = np.maximum(1.,number_of_modes) # look for mode with neff(d) not monotoneous increasing
                    if -epsilon_s > epsilon_f:
                        number_of_modes = number_of_modes + 1 # mode which converges to film-substrate SPP for d -> infinity
                else:
                    number_of_modes = np.ceil(2*d*fs-nnarctan(epsilon_f/epsilon_c*np.sqrt(a))/np.pi)
                    if -epsilon_c > epsilon_s:
                        number_of_modes = number_of_modes + 1 # mode which converges to SPP for d -> 0
    elif polarization == 'TM':
        if -epsilon_f > epsilon_s: 
            number_of_modes = number_of_modes + 1 # mode which converges to film-substrate SPP for d -> infinity
        if -epsilon_f > epsilon_c and epsilon_f*epsilon_c/(epsilon_f+epsilon_c) > epsilon_s:
            number_of_modes = number_of_modes + 1 # mode which converges to film-claddinge SPP for d -> infinity
        if -epsilon_f <= epsilon_c:
            number_of_modes = number_of_modes + 1 # mode which starts at sqrt(epsilon_s) 
            
    return number_of_modes.astype(int)

def invdisp(epsilon_f,epsilon_s,epsilon_c,n_eff,d,polarization): # compute inverse dispersion relation for fundamental mode mu=0
    kx = np.sqrt(epsilon_f-n_eff**2+0j)
    gs = np.real(np.sqrt(n_eff**2-epsilon_s+0j))
    gc = np.real(np.sqrt(n_eff**2-epsilon_c+0j))
    if polarization == 'TE':
        if kx**2 == gc*gs:
            invdisp = np.pi/2
        else:
            invdisp = nnarctan((kx*gc+kx*gs)/(kx**2-gc*gs))
    elif polarization == 'TM': 
        if epsilon_s*epsilon_c*kx**2 == epsilon_f**2*gc*gs:
            invdisp = np.pi/2
        else:
            argument = (epsilon_f*epsilon_s*kx*gc+epsilon_f*epsilon_c*kx*gs)/(epsilon_s*epsilon_c*kx**2-epsilon_f**2*gc*gs)
            if argument == 1j or (-epsilon_c > epsilon_f and n_eff == np.real(np.sqrt(epsilon_f*epsilon_c/(epsilon_f+epsilon_c)+0j))) \
                              or (-epsilon_s > epsilon_f and n_eff == np.real(np.sqrt(epsilon_f*epsilon_s/(epsilon_f+epsilon_s)+0j))):
                invdisp = 4*d*np.pi*kx # get SPP for d -> infinity
            elif epsilon_s > 0 and -epsilon_c > epsilon_s:
                invdisp = np.arctan(argument)
            else:
                invdisp = nnarctan(argument)
    if epsilon_f == n_eff**2:
        invdisp = np.maximum(0.,(epsilon_f*epsilon_s*gc+epsilon_f*epsilon_c*gs)/(-epsilon_f**2*gc*gs*2*np.pi)) # get kx=0
    else:
        invdisp = invdisp/(2*np.pi*kx)
    if np.imag(invdisp) != 0: # sanity check
        invdisp = np.nan
        
    return invdisp

def n_eff_d(epsilon_f,epsilon_s,epsilon_c,n_eff,d,mu,polarization): # compute n_eff at d for mode with index mu
    def func(n_eff):
        if mu == 0:
            func = np.real(invdisp(epsilon_f,epsilon_s,epsilon_c,n_eff,d,polarization)-d)
        else:
            func = np.real(invdisp(epsilon_f,epsilon_s,epsilon_c,n_eff,d,polarization)-d+mu/(2*np.sqrt(epsilon_f-n_eff**2+0j)))
        return func
    vfunc = np.vectorize(func)
    argmin = sps.argrelextrema(vfunc(n_eff), np.less)
    argmax = sps.argrelextrema(vfunc(n_eff), np.greater)
    argext = np.sort(np.concatenate((np.int_([0]),argmin[0],argmax[0],np.int_([n_eff.size-1]))))
    if argext.size > 4:
        gui.input_warning("DR too complicated to plot all modes...")
    n_eff_select = -np.ones(3)
    for index in range(argext.size-1):
        if func(n_eff[argext[index]])*func(n_eff[argext[index+1]]) < 0:
           n_eff_select[index] = spo.bisect(func, n_eff[argext[index]], n_eff[argext[index+1]])
    return np.delete(n_eff_select, np.where( n_eff_select < 0 ))

def mode_profile(epsilon_f,epsilon_s,epsilon_c,n_eff,d,x,polarization): # computing mode profile
    u = np.sqrt(epsilon_f-n_eff**2+0j)*2*np.pi*d
    v = np.real(np.sqrt(n_eff**2-epsilon_s+0j))*2*np.pi*d
    w = np.real(np.sqrt(n_eff**2-epsilon_c+0j))*2*np.pi*d
    if np.abs(np.imag(u)) > 10:
        gui.input_warning("Mode profile not reliable with present numerical precision. Plotting DR only...")
        F=np.ones_like(x)*np.nan
        Gx=np.ones_like(x)*np.nan
        Gz=np.ones_like(x)*np.nan
    else:
        if polarization == 'TE':
            Fs = np.exp(v*(x[(x<=-0.5)]+0.5))
            Ff = np.cos(u*(x[(x>-0.5) & (x<0.5)]+0.5))+v/u*np.sin(u*(x[(x>-0.5) & (x<0.5)]+0.5))
            Fc = (np.cos(u)+v/u*np.sin(u))*np.exp(-w*(x[(x>=0.5)]-0.5))
            F = np.real(np.concatenate((Fs,Ff,Fc)))
            Gx = F
            Gzs = Fs
            Gzf = -u/v*np.sin(u*(x[(x>-0.5) & (x<0.5)]+0.5))+np.cos(u*(x[(x>-0.5) & (x<0.5)]+0.5))
            Gzc = -w/v*Fc
            Gz = np.real(np.concatenate((Gzs,Gzf,Gzc)))
        elif polarization == 'TM':    
            Fs = np.exp(v*(x[(x<=-0.5)]+0.5))
            Ff = np.cos(u*(x[(x>-0.5) & (x<0.5)]+0.5))+epsilon_f/epsilon_s*v/u*np.sin(u*(x[(x>-0.5) & (x<0.5)]+0.5))
            Fc = (np.cos(u)+epsilon_f/epsilon_s*v/u*np.sin(u))*np.exp(-w*(x[(x>=0.5)]-0.5))
            F = np.real(np.concatenate((Fs,Ff,Fc)))
            Gxs = Fs
            Gxf = Ff*epsilon_s/epsilon_f
            Gxc = Fc*epsilon_s/epsilon_c
            Gx = np.real(np.concatenate((Gxs,Gxf,Gxc))) 
            Gzs = Fs
            Gzf = -epsilon_s/epsilon_f*u/v*np.sin(u*(x[(x>-0.5) & (x<0.5)]+0.5))+np.cos(u*(x[(x>-0.5) & (x<0.5)]+0.5))
            Gzc = -epsilon_s/epsilon_c*w/v*Fc
            Gz = np.real(np.concatenate((Gzs,Gzf,Gzc))) 
        
    return F,Gx,Gz