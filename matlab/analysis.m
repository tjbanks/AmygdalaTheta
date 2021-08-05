function [data] = analysis(ecph5path,spikeh5path)
%analysis('../outputECP/ecp.h5','../outputECP/spikes.h5')
%clear all;
close all;
clc;
channel = 1;
dt = 0.05;
steps_per_ms = 1/dt;
skip_seconds = 5;
skip_ms = skip_seconds*1000
skip_n = skip_ms * steps_per_ms;%50000;
end_ms = 15000
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
data = data(data>skip_ms,:)
start_time=1+skip_ms;
stop_time=end_ms;%max(data(:,1));  


%load LP.mat;
upscale=1;
TypeA_num=569*upscale;
TypeC_num=231*upscale;
num_pyr=TypeA_num+TypeC_num; %define the cell number used for plots
inter_num=93*upscale;
som_num=51*upscale;
cr_num=56*upscale;

all_num=num_pyr+inter_num+som_num+cr_num;

%principal_cell = [1:num_pyr]' - 1;

%load NP_cell.mat;

principal_cell_A = [1:TypeA_num]' - 1;
principal_cell_C = [TypeA_num+1:num_pyr]' - 1;
Inter_cell = [num_pyr+1:num_pyr+inter_num]' - 1;
Som_cell = [num_pyr+inter_num+1:num_pyr+inter_num+som_num]' - 1;
Cr_cell = [num_pyr+inter_num+som_num+1:all_num]' - 1;

figure (12)
data_P=data(find(data(:,2)<=num_pyr),:);
data_I=data(find(data(:,2)>num_pyr & data(:,2)<=num_pyr+inter_num),:);
data_SOM=data(find(data(:,2)>num_pyr+inter_num & data(:,2)<=num_pyr+inter_num+som_num),:);
data_CR=data(find(data(:,2)>num_pyr+inter_num+som_num & data(:,2)<=all_num),:);

plot(data_P(:,1),data_P(:,2)+1,'blue.')
hold on; 
plot(data_I(:,1),data_I(:,2)+1,'r.');
hold on; 
plot(data_SOM(:,1),data_SOM(:,2)+1,'green.');
hold on; 
plot(data_CR(:,1),data_CR(:,2)+1,'m.');
axis tight;  
legend('PN','FSI','SOM','CR');


data_analysis=data(data(:,1)>=start_time&data(:,1)<=stop_time,:);
spikes_sort=sortrows(data_analysis,2);
[n, bin] = histcounts(spikes_sort(:,2), unique(spikes_sort(:,2)));
n_cum=cumsum(n);
cell_freq=zeros(all_num,1);
cell_freq(spikes_sort(n_cum(1),2)+1,1)=n(1)/(stop_time-start_time)*1e3;
for i=2:length(n);                  
cell_freq(spikes_sort(n_cum(i),2)+1,1)=n(i)/(stop_time-start_time)*1e3;
end;
plastic_cell_freq_mean=mean(cell_freq(1:num_pyr));
plastic_cell_freq_std=std(cell_freq(1:num_pyr));
plastic_cell_A_freq=cell_freq(1:TypeA_num);
plastic_cell_C_freq=cell_freq(TypeA_num+1:num_pyr);
plastic_cell_freq=cell_freq(1:num_pyr);

Inter_plastic_cell_freq_mean=mean(cell_freq(num_pyr+1:num_pyr+inter_num));
Inter_plastic_cell_freq_std=std(cell_freq(num_pyr+1:num_pyr+inter_num));
Inter_plastic_cell_freq=cell_freq(num_pyr+1:num_pyr+inter_num);

SOM_start = num_pyr+inter_num+1;
SOM_end = num_pyr+inter_num+som_num;
SOM_plastic_cell_freq_mean=mean(cell_freq(SOM_start:SOM_end));
SOM_plastic_cell_freq_std=std(cell_freq(SOM_start:SOM_end));
SOM_plastic_cell_freq=cell_freq(SOM_start:SOM_end);

CR_start = num_pyr+inter_num+som_num+1;
CR_end = num_pyr+inter_num+som_num+cr_num;
CR_plastic_cell_freq_mean=mean(cell_freq(CR_start:CR_end));
CR_plastic_cell_freq_std=std(cell_freq(CR_start:CR_end));
CR_plastic_cell_freq=cell_freq(CR_start:CR_end);

% GG_list = ['PN_freq','.txt'];
% dlmwrite(GG_list,plastic_cell_freq,'delimiter','\t','precision', '%f');
% 
% GG_list = ['ITN_freq','.txt'];
% dlmwrite(GG_list,Inter_plastic_cell_freq,'delimiter','\t','precision', '%f');
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
hold on; 
bh_C=bar(xb_C(:,I),nb_C(:,I)/numel(plastic_cell_C_freq),'EdgeColor','none','BarWidth',1);
set(bh_C,'facecolor',[1 0.4 0.4]);

A=[0:3:100];
[nb,xb]=hist(Inter_plastic_cell_freq,A);
[X,I]=find(nb>0);
hold on; 
bh=bar(xb(:,I),nb(:,I)/numel(Inter_plastic_cell_freq),'EdgeColor','none','BarWidth',1);
set(bh,'facecolor',[0.6 0.6 0.6]);
nb_I_plot=xb(:,I)';
xb_I_plot=nb(:,I)';

%SOM
A=[0:3:100];
[nb,xb]=hist(SOM_plastic_cell_freq,A);
[X,I]=find(nb>0);
hold on; 
bh=bar(xb(:,I),nb(:,I)/numel(SOM_plastic_cell_freq),'EdgeColor','none','BarWidth',1);
set(bh,'facecolor',[0.2 0.6 0.2]);
nb_I_plot=xb(:,I)';
xb_I_plot=nb(:,I)';

%CR
A=[0:3:100];
[nb,xb]=hist(CR_plastic_cell_freq,A);
[X,I]=find(nb>0);
hold on; 
bh=bar(xb(:,I),nb(:,I)/numel(CR_plastic_cell_freq),'EdgeColor','none','BarWidth',1);
set(bh,'facecolor',[0.8 0.2 0.8]);
nb_I_plot=xb(:,I)';
xb_I_plot=nb(:,I)';
     
set(gca,'xscale','log');xlabel('Hz'); 
xlim([10^(-2) 10^(2)]);
title_p=sprintf('Spiking freq histogram for (PNs is %3.2fHz+/-%3.2fHz))',plastic_cell_freq_mean,plastic_cell_freq_std);
title_i=sprintf('ITNs is %3.2fHz+/-%3.2fHz)',Inter_plastic_cell_freq_mean,Inter_plastic_cell_freq_std);
title_s=sprintf('SOMSs is %3.2fHz+/-%3.2fHz) ',SOM_plastic_cell_freq_mean,SOM_plastic_cell_freq_std);
title_c=sprintf('CRs is %3.2fHz+/-%3.2fHz) ',CR_plastic_cell_freq_mean,CR_plastic_cell_freq_std);
x={[title_p],[title_i],[title_s],[title_c]}

title(x);ylabel('Percentage(#)');
legend('PNa','PNc','FSI','SOM','CR')

    
