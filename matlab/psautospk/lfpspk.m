function [data] = lfpspk(h5path, network, nfft, noverlap)
% LFP takes an h5 ecp file output from bmtk and a channel number
% Plots severla plots for frequency power USING PSAUTOSPK.m
% Example usage:
% lfpspk('../../vpsi_inh_spikes.h5','vpsi_inh');


timestamp_path = join(['/spikes/',network,'/timestamps']);
timestamps = h5read(h5path,timestamp_path);

% assuming LFP is stored in an array called LFP
fs = 10000; % sampling frequency in Hz
dt = 1; % delta time
sp = 1/dt; %samples per
p = fs*sp; %total bins where spike is possible (total timesteps)
spk = zeros(1,15000); %array of zeros, one for each timestep
s = int64(timestamps); %change timestamps of spikes to indices
spk(s+1) = 1; %set spike times

if ( nargin == 2 )
  nfft = 1024;
  noverlap = 512;
end;

if ( nargin == 3 )
  noverlap = 512;
end;


%psautospk(spk,dt,nfft)
window = hanning(1024);
[f,Pxxn,tvect,Cxx] = psautospk(spk,dt,nfft,window,noverlap);
plot(f,Pxxn)
xlim([0 25])
end