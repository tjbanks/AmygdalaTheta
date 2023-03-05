function [f,Pxxn,tvect,Cxx] = psautospk(spk,tstep,nfft,window,noverlap)
% psautospk: function which returns the autocorrelation function and power
% spectral density of a neuronal spike train. If called without output
% parameters, the function plots the autocorrelation and the power
% spectrum of the spike train. 
%  
%   [f,Pxxn,tvect,Cxx] = psautospk(spk,tstep)
%
%   where
%       spk = spike train
%	tstep = sampling time step (in msec)
%   
%   The function may also be called as follows:
%
%   [f,Pxxn,tvect,Cxx] = psautospk(spk,tstep,nfft,window,noverlap,dflag)
%
%   where the additional parameters control the power spectrum calculation
%   (see matlab psd function) and replace the default values:
%
%       nfft = number of points used for a single fft operation (default: 2048)
%       window = window function (default: bartlett(nfft))
%       noverlap = number of overlapping points per segment (default: 1024)
%       dflag = detrending option flag (default:'none')
%
%   The return parameters are:
%
%   	f = frequency samples
%	Pxxn = power spectral density at the frequency samples
%	tvect = time domain samples
%	Cxx = autocorrelation function at the time domain samples
%

if ( (nargin ~= 5) & (nargin ~= 2) )
  disp(' ');
  disp('usage1: psautospk(spk,tstep) ');
  disp(' ');
  disp('usage2: psautospk(spk,tstep,nfft,window,noverlap) ');
  disp(' ');
  disp('       for more information type "help psautospk" in the main');
  disp('       matlab window');
  disp(' ');
  return;
end;


%The parameters are setup with the following simulations in mind:
%a sampling rate of 0.5 msec and a 
%simulation time of 135 sec so that a reliable estimate of the power
%spectrum can usually be obtained at a resolution of approx. 1Hz. 
if ( nargin == 2 )
  nfft = 1024;
  window = bartlett(nfft);
  noverlap = 512;
end;

%computes the sampling frequency in Hz
tstep_s = tstep*1e-3;  %converts to sec
Fs = 1/tstep_s; %in Hz

%computes and subtracts the mean firing rate
spk = spk(:); %convertes to column vector if necessary
spk = spk*Fs; %converts to units of spikes/sec
l_spk = length(spk);
s_spk = sum(spk);
m_spk = s_spk/l_spk;
spk = spk - m_spk;

%[Pxx,f] = psd(spk,nfft,Fs,window,noverlap,dflag);
spk = lowpass(spk,500,Fs);
[Pxx,f] = pwelch(spk,window ,noverlap, nfft,Fs);

%converts to units of (spk/Hz)^2
Pxxn = Pxx * tstep_s;

%prepares the data to compute the autocorrelation
Pxxx = zeros(nfft,1);
Pxxx(1:nfft/2+1,1) = Pxx(1:nfft/2+1,1);
for k = 2:nfft/2
  Pxxx(nfft+2-k,1) = Pxx(k,1);
end;

%computes the autocorrelation function
Cxxx = fft(Pxxx,nfft);
%normalizes to get the usual definition of autocorrelation
Cxxx = Cxxx/nfft;

tvect = -(nfft/2)*tstep:tstep:(nfft/2)*tstep;
Cxx = zeros(nfft+1,1);
for k = 1:nfft/2
  Cxx(k,1) = real(Cxxx(nfft/2 + k,1));
end;