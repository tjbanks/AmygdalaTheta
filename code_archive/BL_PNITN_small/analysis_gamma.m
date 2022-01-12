clc;close all;clear all;
%addpath('F:\Feng_gammNSG')
data = load('data');
% data=data(data(:,1)<=5000,:);% data = load('F:\Feng_gammNSG\competition\small network\NEW noise\smallnetwork-noise52-41-results.tar\smallnetwork-noise52-41-results\BLA_for_NSG_extrinsic_invivo_QW_disdelay_multielec_forNSG\data');
lfp = load('./LFPs/LFP_elec_combine');
tstop=length(lfp)
burst_len=1000;rest_len=500;
burst_start=[1:burst_len+rest_len:tstop];
burst_stop=[1+burst_len:burst_len+rest_len:tstop];

burstflag=3; %%%%1 means burst period. 2 means resting peroid, 3 means all

lfp_burst=[];lfp_rest=[];
if burstflag==1;
    for i=1:numel(burst_stop)
    lfp_burst=[lfp_burst;lfp(burst_start(i):burst_stop(i))]
    end
    lfp_test=lfp_burst;
elseif burstflag==2;
     for i=1:numel(burst_stop)-1
    lfp_rest=[lfp_rest;lfp(burst_stop(i):burst_start(i+1))]
     end
     lfp_test=lfp_rest;
elseif burstflag==3;
     lfp_test=lfp;
end


% lfp=lfp(1:5000);
% lfp=load('F:\Feng_gammNSG\optimize&simplify\gap verify\full_brigdges\5000r\LFPs\LFP_elec_0');lfp=lfp(1000:6000,:); %%%reference
nfft=1024;fs=1000;
figure(666);plot(lfp_test);
windowLen = 1024;

[pxx,f] = pwelch(lfp_test,nfft,0,nfft,fs);

figure(1213)
plot(f,pxx,'LineWidth', 2);
 set(gca, 'YScale', 'log');
 set(gca, 'XScale', 'log');axis tight;

     
%[f,Pxxn,tvect,Cxx] = psautospk(lfp_test, 1, windowLen, bartlett(windowLen), windowLen/2, 'none') ;
%figure (1212)
%plot(f,Pxxn,'b', 'LineWidth', 2);
%set(gca,'xscale','log');set(gca,'yscale','log'); legend('PSD of LFP');xlabel('Hz');ylabel('PSD');
%axis tight

start_time=1;
stop_time=max(data(:,1));     



%load LP.mat;
upscale=1;
num=900*upscale; %define the cell number used for plots
inter_num=100*upscale;

all_num=num+inter_num;

principal_cell = [1:num]' - 1;
Inter_cell = [num+1:inter_num+num]' - 1;


figure (12)
data_P=data(find(data(:,2)<num),:);
data_I=data(find(data(:,2)>=num&data(:,2)<num+inter_num),:);
plot(data_P(:,1),data_P(:,2)+1,'blue.')
hold on; plot(data_I(:,1),data_I(:,2)+1,'r.'); axis tight;  

data_analysis=data(data(:,1)>=start_time&data(:,1)<=stop_time,:);
spikes_sort=sortrows(data_analysis,2);
[n, bin] = histc(spikes_sort(:,2), unique(spikes_sort(:,2)));
n_cum=cumsum(n);
cell_freq=zeros(all_num,1);
cell_freq(spikes_sort(n_cum(1),2)+1,1)=n(1)/(stop_time-start_time)*1e3;
for i=2:length(n);                  
cell_freq(spikes_sort(n_cum(i),2)+1,1)=n(i)/(stop_time-start_time)*1e3;
end;
plastic_cell_freq_mean=mean(cell_freq(1:num));
plastic_cell_freq_std=std(cell_freq(1:num));
% plastic_cell_A_freq=cell_freq(1:TypeA);
% plastic_cell_C_freq=cell_freq(TypeA+1:num);
plastic_cell_freq=cell_freq(1:num);

ITN_cell_freq_mean=mean(cell_freq(num+1:num+inter_num));
ITN_cell_freq_std=std(cell_freq(num+1:num+inter_num));
ITN_cell_freq=cell_freq(num+1:num+inter_num);

figure (1)
A=[0:0.05:100];
[nb_A,xb_A]=hist(plastic_cell_freq,A);
[X,I]=find(nb_A>0);   %%%delete any zero histogram bar
bh_A=bar(xb_A(:,I),nb_A(:,I)/numel(plastic_cell_freq),'EdgeColor','none','BarWidth',1);
set(bh_A,'facecolor',[0.4 0.4 0.8]);

% nb_A_plot=xb_A(:,I)';
% xb_A_plot=nb_A(:,I)';

% set(gca,'xscale','log'); 
% xlim([10^(-2) 10^2]);

% [nb_C,xb_C]=hist(plastic_cell_C_freq,A);
% [X,I]=find(nb_C>0);   %%%delete any zero histogram bar
% 
% figure (1)
% hold on; bh_C=bar(xb_C(:,I),nb_C(:,I)/numel(plastic_cell_C_freq),'EdgeColor','none','BarWidth',1);
% set(bh_C,'facecolor',[1 0.4 0.4]);

A=[0:3:100];
[nb,xb]=hist(ITN_cell_freq,A);
[X,I]=find(nb>0);
hold on; bh=bar(xb(:,I),nb(:,I)/numel(ITN_cell_freq),'EdgeColor','none','BarWidth',1);
set(bh,'facecolor',[0.8 0.4 0.4]);
     
set(gca,'xscale','log');xlabel('Hz'); 
xlim([10^(-2) 10^(2)]);
x=sprintf('Spiking freq histogram for PNs is %3.2fHz+/-%3.2fHz;PVs is %3.2fHz+/-%3.2fHz; ',plastic_cell_freq_mean,plastic_cell_freq_std,ITN_cell_freq_mean,ITN_cell_freq_std);
title(x);ylabel('Percentage(#)');
legend('PN','ITN')

spikefreqinfo=cell(2,3);
spikefreqinfo(1,1:2)={'PNs','PVs'};
%%%%fill in CPs
cell_temp=[];cell_temp=cell(num+1,2);
cell_temp(1,1:2)={'Freq(Hz)','IDs'};
cell_temp(2:end,:)=num2cell([plastic_cell_freq,[1:num]']);
spikefreqinfo{2,1}=cell_temp;

%%%%fill in PVs

cell_temp=[];cell_temp=cell(inter_num+1,2);
cell_temp(1,1:2)={'Freq(Hz)','IDs'};
cell_temp(2:end,:)=num2cell([cell_freq(num+1:num+inter_num),[num+1:num+inter_num]']);
spikefreqinfo{2,2}=cell_temp;


save('spikefreqinfo_shortburst','spikefreqinfo');

