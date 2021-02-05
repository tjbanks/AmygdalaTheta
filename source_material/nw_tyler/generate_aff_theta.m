clear all;close all;
clc;

rng(1)
F=8; %%%%%Hz  define rythmicity f(Hz)-main
FR=1; %20 %%%%control mean FR of input(Hz)
cell_num=1000;%%%%% number of cells needed to generate, including PFC,ER and PR,1800 as original
t_len=50000; %defind time length to generate(ms)
sigma_FR=120;%0.2; %%Hz, to control the variability among cells' FR
sigma_spk=40;%5; %%ms, to control the variability within one cell's spk time

mu = log((FR^2)/sqrt(sigma_FR+FR^2));
sigma = sqrt(log(sigma_FR/(FR^2)+1));

FR_cell = lognrnd(mu,sigma,[cell_num,1]); %%%randomize FR
spk_num_cell=round(t_len*FR_cell/1000);
figure (1)
hist(FR_cell);set(gca,'xscale','log');xlabel('Hz'); title('cell FR dis.');


%%%%generate nT(interval) based on lognormal distribution
T_main=1000/F; %%%%%Hz  define main rythmicity interval(ms)
T_rand=T_main;
sigma_T=1*7; %%to control the variability among interval(ms)

mu = log((T_main^2)/sqrt(sigma_T+T_main^2));
sigma = sqrt(log(sigma_T/(T_main^2)+1));
k=0;
while T_rand(end)<=t_len;
    k=k+1
    T_temp=[];T_temp=T_rand(end)+lognrnd(mu,sigma);
    T_rand=[T_rand,T_temp];
end
figure (345)
hist(diff(T_rand),[0:0.01:100]);set(gca,'xscale','log');xlabel('ms');title('interval distr. between T(ms)');

for j=1:cell_num;
    clear spk_int
    spk_int=datasample([1:numel(T_rand)],spk_num_cell(j));
    spiketime_cell{j,1}=T_rand(spk_int)+normrnd(0,sigma_spk,[1,numel(spk_int)]); %%%to store spk time(ms) per cell
    
    spiketime_cell{j,1}=unique(round(spiketime_cell{j,1}(1,spiketime_cell{j,1}>0)*1e2)/1e2);%%%%make sure event diff is greater thn 0.01
    spiketime_cell{j,1}=sort(spiketime_cell{j,1});
end

%%%plot generated spike raster
%figure (2)
% for i=1:cell_num;
% hold on;plot(spiketime_cell{i,1},i*ones(1,numel(spiketime_cell{i,1})),'red.')
% end;
% axis tight;

%%%%plot moving-windowed firing rate     
spiketimes=[spiketime_cell{:}]';
spiketimes=spiketimes(spiketimes<=t_len&spiketimes>0.5);
binSize = 20; %%ms
timeMax = max(spiketimes);
timeMax = timeMax - binSize;
binSpikeCount = [];

% tt_low=[0:1:tstop]-(binSize/2);
tt_up=[0:1:t_len]+(binSize/2);
N_up = histc(spiketimes, tt_up);

movingSum = conv(N_up, ones(1, binSize));% N_up = histc(spiketimes, tt_up);

%%%%to plot P and I group rate sepertated.
% spiketimes_P=[spiketime_cell{1:27*900}]';
% spiketimes_I=[spiketime_cell{27*900+1:27*1000}]';
% movingSum_P=conv(histc(spiketimes_P, tt_up),ones(1, binSize));
% movingSum_I=conv(histc(spiketimes_I, tt_up),ones(1, binSize));
% figure(789)
% plot(movingSum_P/(27*900));hold on;plot(movingSum_I/(27*100),'r');legend('norm. group rate for p','norm. group rate for i');

figure (12222)
% psd(movingSum);
nfft=1024/2;fs=1000;
figure;plot(movingSum)
[pxx,f] = pwelch(movingSum,nfft,0,nfft,fs);

figure;plot(f,pxx); set(gca, 'YScale', 'log');set(gca, 'XScale', 'log'); 

% rSPT=round(spiketimes*1);
% rSPT_sort=sort(rSPT);
% rSPT_unique=unique(rSPT_sort);
% [n, bin] = histc(rSPT_sort, rSPT_unique);
% ev_P=zeros(round(max(rSPT_unique)),1);
% ev_P(rSPT_unique)=n;
% 
%      %ev_P(1:exclusive_time) = [];
% windowLen = 1024; 
% ev_P=ev_P-mean(ev_P);
% [f,Pxxn,tvect,Cxx] = psautospk(ev_P, 1, windowLen, bartlett(windowLen), windowLen/2, 'none') ;
% figure(1111)
%      plot(f,Pxxn,'blue', 'LineWidth', 2);
% set(gca,'yscale','log'); set(gca,'xscale','log');legend('PSD for PN firing rate');xlabel('Hz');ylabel('PSD');

%   for tt=0:1:tstop %timeMax
%    binSpikeTimes = spiketimes((spiketimes(:) >= tt-(binSize/2) & spiketimes(:) < tt+ (binSize/2) ) );
%    binSpikeCount(end+1) = length(binSpikeTimes);
% end
% firingRate= (1000*  binSpikeCount / binSize); % / length(nID) ; 
figure (1212123)
% firingRate=firingRate';
%subplot(2,1,1)
% plot(firingRate,'r');% plot(firingRate)
hold on; plot(movingSum*1000/binSize,'r');

dname=''

fid=fopen(fullfile(dname,'spikesmatrix_op_ryt'),'w');  %%%%for saving all pre (both INH and EXC) cells for PNs&ITNs 
for i=1:cell_num;
clear fmt; fmt=[repmat('%.2f\t',1,numel(spiketime_cell{i,:})) '\n'];   %%%%format
spikesnum(i,1) = numel((spiketime_cell{i,:}));   %%%%find spk time, unit is ms
fprintf(fid,fmt,spiketime_cell{i,:});
end

ind_num=cumsum(spikesnum);ind_num=[0;ind_num];
fid_ind=fopen(fullfile(dname,'spikes_ind_ryt'),'w'); %%%%%for saving index for segerate fid file for each neuron
%fclose(fid_ind)
fmt=[repmat('%d\t',1,2) '\n'];
for i=1:cell_num;
    i
 temp=[];temp=[ind_num(i),ind_num(i+1)-1];
 fprintf(fid_ind,fmt,temp);
end

