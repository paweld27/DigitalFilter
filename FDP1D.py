
import pandas as pd
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt


fs = 100
dt = 1.0/fs
fMax = 1000.0
nN = 10000

fg = 10
T = 1/(2*np.pi*fg)
alfa = dt/(T + dt)
alfa1 = 1 - alfa

db3 = 0.7071

freq = np.linspace(0, fMax-1/nN, nN) 
first_col_time = True
title = 'LPF1'

b0 = [1.]             #    1      >     1
a0 = [T, 1.]          #  1 + T*s  >  T*s + 1
w0, h0 = signal.freqs(b0, a0, worN=freq)



def make_PY(*args):
    
    b = [alfa]                 #        alfa
    a = [1., -(1 - alfa)]      #   1 -(1 - alfa)*z^-1
    
    #    b, a = signal.butter(1, fg, 'low', fs=fs)
    
    b, a = signal.bilinear(b0, a0, fs=fs)
    
    w, h = signal.freqz(b, a, fs=fs, worN=freq)
    
    fileDAT = pd.DataFrame(freq, columns=['Freq'])

    Lin = list(map(lambda x: abs(1/(1+1j*2*np.pi*x*fMax/nN*T)), range(1,nN+1)))

    Dig = list(map(lambda x: abs(alfa/(1-alfa1*np.exp(-1j*2*np.pi*x*fMax/nN*dt))), range(1,nN+1)))

    fileDAT['Digital'] = abs(h)
    fileDAT['Analog'] = Lin

#    fileDAT['Err'] = 10*(fileDAT['Dig'] - fileDAT['Lin']) / fileDAT['Lin']

    ncols = 1

    return fileDAT, title, first_col_time, fs, ncols




if __name__ == "__main__":
    fileDAT, _, _, _, _ = make_PY()

    nerr = 1
    colors=('#C3063C', '#0000FF', 'green')
   
    print(f'fs = {fs}')
    print(f'fg = {fg}')
    print(f'T  = {T}')
    print(f'alfa=dt/(T+dt)= {alfa}')
    print(f'(1-alfa) = {alfa1}')
    
    try:
        idx1 = np.nonzero(fileDAT['Digital'] <= db3)[0][0]
    except:
        nerr = 0
        fgc = 0.001
    if(nerr):
      fgc = fileDAT['Freq'].iloc[idx1]
      print(f'fgc= {fgc}')

    """    
    kanals = fileDAT.columns
    kan_len = len(kanals)-1
    fig, ax = plt.subplots()
    fig.set_figwidth(8) # inches
    fig.set_figheight(4)

    for i in range(kan_len):
        fileDAT.plot(kind='line', x='Freq', y=kanals[i+1], ax=ax, c=colors[i])

    plt.legend(bbox_to_anchor=(1.1, 0.8),loc='center right', draggable=True)
    plt.xlabel('Freq [Hz]')
    plt.ylabel("[V]")
    plt.ylim([0, 1.1])
    plt.xlim([0.1, 1000])
    plt.xscale("log")
    ax.grid(which='both', axis='x')
    ax.axhline(y=db3, ls='--')
    ax.axvline(x=fgc, ls='--', c=colors[0], lw=1)
    ax.axvline(x=fg,  ls='--', c=colors[1], lw=1)

    plt.grid()
    plt.minorticks_on()
    fig_mgr = plt.get_current_fig_manager()
    tb = fig_mgr.toolbar
    tb.zoom()
    """

    figr, axr = plt.subplots()

    """
    instead of substituting s-> 1j*w or z->exp(1j*w)
    in the filter equations, you can utilize the signal function
    from scipy module to directly use the filter coefficients themselves
    """
    
    b = [alfa]            #            alfa
    a = [1., -alfa1]      #   1 -  (1 - alfa)*z^-1
#    b, a = signal.butter(1, fg, 'low', fs=fs)
    b, a = signal.bilinear(b0, a0, fs=fs)
    w, h = signal.freqz(b, a, fs=fs, worN=freq)
    
    axr.semilogx(w , 20 * np.log10(np.maximum(abs(h), 1e-5)), c=colors[0],
                 label='Digital')
    
    axr.semilogx(w0 / (2*np.pi), 20 * np.log10(np.maximum(abs(h0), 1e-5)), c=colors[1],
                 label='Analog')
    
    plt.legend(bbox_to_anchor=(0.01, 0.85),loc='center left', draggable=True)
    axr.axhline(y=-3., ls='--')
    axr.axvline(x=fgc, ls='--', c=colors[0], lw=1)
    axr.axvline(x=fg,  ls='--', c=colors[1], lw=1)
    axr.set_title('Lowpass filter frequency response')
    axr.set_xlabel('Frequency [Hz]')
    axr.set_ylabel('Amplitude [dB]')
    axr.axis((0.1, 100, -15, 1))
    axr.grid(which='both', axis='both')    

    axp = axr.twinx()

    phase0 = np.unwrap(np.angle(h0))*180./np.pi
    axp.semilogx(w0 / (2*np.pi), phase0, ls='--', c=colors[1])

    phase = np.unwrap(np.angle(h))*180./np.pi
    axp.semilogx(w , phase, ls='--', c=colors[0])
            
    axp.set_ylabel('Phase [--deg]', color='k')
    
    axp.axvline(x=fgc, ls='--', c=colors[0], lw=1)
    axp.axvline(x=fg,  ls='--', c=colors[1], lw=1)
    axp.grid(axis='x')
    axp.axis('tight')
    axp.set_ylim([-45, 45])
    nticks = 8

    print(f'b = {b}')
    print(f'a = {a}')

    plt.show()
#sys.exit("S T O P")


