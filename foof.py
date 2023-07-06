import h5py
from fooof import FOOOF
from fooof.sim.gen import gen_aperiodic
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch,decimate
from scipy.signal.windows import hann as hanning
#%matplotlib inline

CASE = 1 # 1 or 2

# Load data
ecp_h5_location = 'outputECP/ecp.h5' #'tyler_case2_ecp.h5' if CASE == 2 else 'tyler_case1_ecp.h5'
ecp_channel = 0
f = h5py.File(ecp_h5_location)
data_raw = np.array(f['ecp']['data'])
ecp = data_raw.T[ecp_channel] #flip verts and grab channel 0

# Format data
dt = 0.05
steps_per_ms = 1/dt
skip_seconds = 5
skip_ms = skip_seconds*1000
skip_n = int(skip_ms * steps_per_ms)
end_ms = 15000
downsample = int(1.0/dt)
nfft = 1024
fs = 1000
noverlap = 0

data = ecp[skip_n:]
lfp_d = decimate(data,downsample) #downsample the data to fit ms (steps used 20=1/.05 step)
#win = hanning(nfft, True)
#f,pxx = welch(lfp_d,fs,window=win,noverlap=noverlap,nfft=nfft)
f,pxx = welch(lfp_d,fs=1000,nfft=1024)

# FOOOF
freqs,spectrum = np.array(f),np.array(pxx)
fm = FOOOF(aperiodic_mode='knee')
fm.fit(freqs, spectrum, [1,150])
ap_fit = fm._ap_fit
residual_spec = spectrum[0:152] - 10**ap_fit

# Plot
#plt.plot([i for i in range(len(residual_spec))],residual_spec)
plt.plot(freqs[:len(residual_spec)], residual_spec)
#plt.plot(10**ap_fit)
#plt.plot([i for i in range(4,13)],residual_spec[4:13]) # Only theta range
plt.xlabel("Hz")
plt.ylabel("Power")
plt.title('residual from spectrum - aperiodic component')
plt.show()
theta = residual_spec[4:13]

psd_power = max(theta)
print(f"PSD Theta Power for case {CASE}: {psd_power}")
