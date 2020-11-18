function [data] = lfp(h5path, channel,nfft,noverlap)
% LFP takes an h5 ecp file output from bmtk and a channel number
% Plots severla plots for frequency power
% Example usage:
% lfp('../100_output/ecp.h5',1);

fs = 10000; % sampling frequency in Hz

if ( nargin < 2 )
  channel = 1;
end;
if ( nargin < 3 )
  nfft = 1024/2;
end;

data = h5read(h5path,'/ecp/data');
LFP = data(channel,:);
%save('lfp1200.mat','LFP');

% assuming LFP is stored in an array called LFP

%LFP = lowpass(LFP,500,fs);
%[pxx,f] = pwelch(LFP,bartlett(nfft),0,nfft,fs);
[pxx,f] = pwelch(LFP,nfft,0,nfft,fs);
plot(f,pxx);
%xlim([0,100]);
set(gca,'yscale','log');
set(gca,'xscale','log');

%window = 1*floor(size(LFP,2)/7.8);
%window = size(LFP,2)
%F = 0:1:50; % Limit frequencies from 1 to 50 Hz
%spectrogram(LFP,window,window-2,F,fs,'yaxis');
%hold on;
%colormap(jet)

end


