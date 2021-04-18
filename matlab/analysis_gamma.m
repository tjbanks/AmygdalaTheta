function [data] = analysis(ecph5path,spikeh5path)
%analysis('../outputECP/ecp.h5','../outputECP/spikes.h5')
%clear all;
%close all;
clc;
channel = 1;
skip_n = 50000;
data = h5read(ecph5path,'/ecp/data');
lfp = data(channel,:);
lfp = lfp(skip_n:end);

%lfp=load ('ecp.mat');
lfp_d = downsample(lfp,20);%x1000 mV to V fix
nfft=1024;fs=1000;
figure(1);plot(lfp_d*1e3);
[pxx,f] = pwelch(lfp_d,nfft,0,nfft,fs);

 figure(2);plot(f,pxx*1e3);
 set(gca, 'YScale', 'log');
 set(gca, 'XScale', 'log');
 %return
%load ('spikes.mat');
timestamps = h5read(spikeh5path,'/spikes/BLA/timestamps');
node_ids = h5read(spikeh5path,'/spikes/BLA/node_ids');
node_ids=double(node_ids);


data=([timestamps,node_ids]);
data = data(data>skip_n/10,:)
start_time=1;
stop_time=max(data(:,1));  


%load LP.mat;
upscale=1;
num=900*upscale; %define the cell number used for plots
inter_num=100*upscale;
all_num=num+inter_num;
principal_cell = [1:num]' - 1;

%load NP_cell.mat;
Inter_cell = [num+1:all_num]' - 1;
TypeA=640*upscale;TypeC=260*upscale;

principal_cell_A = [1:TypeA]' - 1;
principal_cell_C = [TypeA+1:num]' - 1;

figure (12)
data_P=data(find(data(:,2)<num),:);
data_I=data(find(data(:,2)>=num),:);
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
plastic_cell_A_freq=cell_freq(1:TypeA);
plastic_cell_C_freq=cell_freq(TypeA+1:num);
plastic_cell_freq=cell_freq(1:num);

Nplastic_cell_freq_mean=mean(cell_freq(num+1:end));
Nplastic_cell_freq_std=std(cell_freq(num+1:end));
Nplastic_cell_freq=cell_freq(num+1:end);

% GG_list = ['PN_freq','.txt'];
% dlmwrite(GG_list,plastic_cell_freq,'delimiter','\t','precision', '%f');
% 
% GG_list = ['ITN_freq','.txt'];
% dlmwrite(GG_list,Nplastic_cell_freq,'delimiter','\t','precision', '%f');
figure (10)
%subplot (2,1,1)
A=[0:0.05:100];
[nb_A,xb_A]=hist(plastic_cell_A_freq,A);
[X,I]=find(nb_A>0);   %%%delete any zero histogram bar
bh_A=bar(xb_A(:,I),nb_A(:,I)/numel(plastic_cell_A_freq),'EdgeColor','none','BarWidth',1);
set(bh_A,'facecolor',[1 0.8 0.4]);

% nb_A_plot=xb_A(:,I)';
% xb_A_plot=nb_A(:,I)';

% set(gca,'xscale','log'); 
% xlim([10^(-2) 10^2]);

[nb_C,xb_C]=hist(plastic_cell_C_freq,A);
[X,I]=find(nb_C>0);   %%%delete any zero histogram bar

figure (10)
hold on; bh_C=bar(xb_C(:,I),nb_C(:,I)/numel(plastic_cell_C_freq),'EdgeColor','none','BarWidth',1);
set(bh_C,'facecolor',[1 0.4 0.4]);

A=[0:3:100];
[nb,xb]=hist(Nplastic_cell_freq,A);
[X,I]=find(nb>0);
hold on; bh=bar(xb(:,I),nb(:,I)/numel(Nplastic_cell_freq),'EdgeColor','none','BarWidth',1);
set(bh,'facecolor',[0.6 0.6 0.6]);
nb_I_plot=xb(:,I)';
xb_I_plot=nb(:,I)';
     
set(gca,'xscale','log');xlabel('Hz'); 
xlim([10^(-2) 10^(2)]);
x=sprintf('Spiking freq histogram for( PNs is %3.2fHz+/-%3.2fHz);ITNs is %3.2fHz+/-%3.2fHz) ',plastic_cell_freq_mean,plastic_cell_freq_std,Nplastic_cell_freq_mean,Nplastic_cell_freq_std);
title(x);ylabel('Percentage(#)');
legend('PNa','PNc','FSI')

    
