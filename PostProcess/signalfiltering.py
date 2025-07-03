import scipy.signal as sp
import numpy as np
import matplotlib.pyplot as plt
import struct

#Kode fra David Kristiansen - Ingen kode herfra brukt i master, men logikk brukt til å lage egne funksjoner


def testsignal(Nsamples, f1, dt):

     # Create a noisy time signal
     a1 = 0.4   # First harmonic amplitude
     a2 = 0.08  # Second harmonic amplitude
     a3 = 0.005 # Third harmonic amplitude
     a4 = 0.0003 # Fourth harmonic amplitude
     a5 = 0.00002 # Fifth harmonic amplitude

     omg1 = 2.0*np.pi*f1 # [rad/s] Signal base circular frequency
     tvec = np.arange(0, (Nsamples)*dt, dt)

     # Create white noise signal (normal distributed random variable)
     noise = a2 * np.random.normal(0, 1, len(tvec))

     xsig = a1*np.sin(omg1*tvec) + a2*np.sin(2*omg1*tvec+0.03) + a3*np.sin(3*omg1*tvec-0.8) + a4*np.sin(4*omg1*tvec-0.4) \
          + a5*np.sin(5*omg1*tvec) + noise

     return xsig, tvec

def moving_average(xsig, w):
     # Moving average of (window) length w
     return np.convolve(xsig, np.ones(w), 'same') / w

def bpass_butter(xsig, dt, flow, fhigh, N=3):

     fs = 1.0/dt
     fnyq = 0.5 * fs

     # Compute filter coefficients
     [b, a] = sp.butter(N, [flow/fnyq, fhigh/fnyq], btype='band', output='ba')  # Band-pass filter

     # Filter signal
     fx_bp = sp.filtfilt(b, a, xsig)

     # Compute frequency filter response amplitude
     
     wn, h = sp.freqz(b, a, worN=2000)  # Compute filter frequency response amplitude


     return fx_bp, (wn, h)

def lpass_butter(xsig, dt, fcut, N=3):

     fs = 1.0/dt
     fnyq = 0.5 * fs

     # Compute filter coefficients
     [b, a] = sp.butter(N, fcut/fnyq, btype='low', output='ba')  # Band-pass filter

     # Filter signal
     fx_bp = sp.filtfilt(b, a, xsig, axis=0)

     # Compute frequency filter response amplitude
     wn, h = sp.freqz(b, a, worN=2000)  # Compute filter frequency response amplitude


     return fx_bp, (wn, h)

def taperWindow(Nsamples, samplingFrequency, rampduration):

     iramp = int(np.round(rampduration*samplingFrequency))

     svec = np.linspace(0, np.pi/2, iramp)
     lramp = np.sin(svec)**2
     rramp = lramp[::-1]
     oneband = np.ones(np.max([Nsamples-2*iramp, 1])) # Unity pass-band
     winfunc = np.concatenate((lramp, oneband, rramp)).reshape((Nsamples, 1)) # Combined window function

     return winfunc

def bpass_gaussian(xsig, dt, flow, fhigh):

     xsig_fft = np.fft.fft(xsig)
     fvec = np.fft.fftfreq(np.shape(xsig)[0], d=dt)
     Nsamples = len(fvec)

     # Create Gaussian window filter:
     M = 50  # Total number of samples in Gaussian ramp function
     gausswin = sp.windows.gaussian(M, std=7)
     lramp = np.array(gausswin[:M//2])  # Left ramp function (left half of Gaussian)
     rramp = np.array(gausswin[M//2:])  # Right ramp function (right half of Gaussian)
     i1 = np.argwhere(np.abs(fvec[:Nsamples//2]) > flow)[0][0]   # Index of low cut frequency
     i2 = np.argwhere(np.abs(fvec[:Nsamples//2]) < fhigh)[-1][0] # Index of high cut frequency
     lwin = np.zeros(np.amax([i1-M//2, 0]))     # Left zero-padding
     rwin = np.zeros(np.min((Nsamples//2 - (i2+M//2), Nsamples//2-1))) # Right zero-padding
     oneband = np.ones(np.max([Nsamples//2-len(lwin)-len(rwin)-len(lramp)-len(rramp), 1])) # Unity pass-band
     tmp = np.concatenate((lwin, lramp, oneband, rramp, rwin)) # Combined window function
     if 2*len(tmp) == len(xsig_fft) - 1:
          winfunc = np.concatenate((tmp, [0], tmp[::-1]))
     else:
          winfunc = np.concatenate((tmp, tmp[::-1]))


     # Filter signal by multiplication with window function in frequency domain
     # Inverse FFT to obtain filtered time-series
     #??????????????
     
     # if len(winfunc) > len(xsig_fft):
     #      winfunc = winfunc[:len(xsig_fft)]
     # elif len(winfunc) < len(xsig_fft):
     #      winfunc = np.pad(winfunc, (0, len(xsig_fft) - len(winfunc)), mode='constant')

     #???????????????
     xsig_bp = np.fft.ifft(xsig_fft * winfunc).real

     return xsig_bp, (fvec, winfunc)

def bpass_gaussian_v2(xsig, dt, flow, fhigh, axis=0, plotfilter=False):

     xsig_fft = np.fft.fft(xsig, axis=axis)
     fvec = np.fft.fftfreq(np.shape(xsig)[axis], d=dt) # Return frequency vector with positive frequencies for the first N//2 elements
     Nsamples = len(fvec)

     winfunc = np.ones_like(fvec[:Nsamples//2])


     # Create Gaussian window filter:
     M = 50  # Total number of samples in Gaussian ramp function
     gausswin = sp.windows.gaussian(M, std=7)
     lramp = np.array(gausswin[:M//2])  # Left ramp function (left half of Gaussian)
     rramp = np.array(gausswin[M//2:])  # Right ramp function (right half of Gaussian)
     winfunc[np.where(np.abs(fvec[:Nsamples//2]) < flow)] = 0 #+np.argwhere(np.abs(fvec[:Nsamples//2]) < fhigh)
     winfunc[np.where(np.abs(fvec[:Nsamples//2]) > fhigh)] = 0 #+np.argwhere(np.abs(fvec[:Nsamples//2]) < fhigh)
     i1 = np.argwhere(np.diff(winfunc) > 0)[0][0]
     i2 = np.argwhere(np.diff(winfunc) < 0)[0][0]+1
     winfunc[i1:np.amin([i1+M//2, len(winfunc)])] *= lramp
     winfunc[np.amax([i2-M//2, 0]):i2] *= rramp

     if len(winfunc)*2 == len(fvec):
          gaussian_window = np.concatenate((winfunc, winfunc[::-1])).reshape(-1, 1)
     else:
          gaussian_window = np.concatenate((winfunc, [0], winfunc[::-1])).reshape(-1, 1)

     if plotfilter:
         plt.figure()
         plt.plot(fvec, np.abs(xsig_fft), color='b', label='FFT signal')
         plt.plot(fvec, gaussian_window, color='r', label='Filter')
         plt.plot(fvec, np.abs(xsig_fft * gaussian_window), color='g', label='Filtered')
         plt.ylim(ymin=0, ymax=5)
         plt.xlabel('Frequency [Hz]')
         plt.show()

     A = xsig_fft * gaussian_window
     print(np.shape(xsig_fft), np.shape(A))
     # Filter signal by multiplication with window function in frequency domain
     # Inverse FFT to obtain filtered time-series
     xsig_bp = np.fft.ifft(xsig_fft * gaussian_window).real

     return xsig_bp, (fvec, gaussian_window)





if __name__ == "__main__":

    # Load binary file
    filename = r"C:\Users\HP\OneDrive - NTNU\Desktop\Master\Code\PostProcess\signalfiltering.py"
    # Adjust reading format based on your data structure
    with open(filename, "rb") as f:
        raw_data = f.read()

    # Convert binary data to a NumPy array
    data_format = "f"  # Assuming 32-bit float values, modify if needed
    fx = np.array(struct.unpack(f"{len(raw_data) // 4}{data_format}", raw_data))

    # Set up time vector (modify dt based on actual sampling rate)
    dt = 0.005  # Adjust this if you know the real sampling time
    tvec = np.arange(0, len(fx) * dt, dt)

    # Define filter frequencies
    f1 = 0.96  # Adjust based on your signal's base frequency
    flow = f1 * 0.85
    fhigh = f1 * 1.15

    # Apply band-pass filters
    fx_bp_butter, _ = bpass_butter(fx, dt, flow, fhigh)
    fx_bp_gauss, _ = bpass_gaussian(fx, dt, flow, fhigh)

    # Plot results
    plt.figure()
    plt.plot(tvec, fx, label="Original Signal")
    plt.plot(tvec, fx_bp_butter, label="Butterworth Filtered")
    plt.plot(tvec, fx_bp_gauss, label="Gaussian Filtered")
    plt.xlabel("Time [s]")
    plt.ylabel("Signal Amplitude")
    plt.legend()
    plt.show()

    
    
    





