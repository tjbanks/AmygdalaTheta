function [data] = ipsc_analysis(syn_report_h5path)
%ipsc_analysis('../outputECP/syn_report.h5')
%clear all;
close all;
clc;
skip_n = 50000;
data_raw = h5read(syn_report_h5path,'/report/BLA/data');
data = sum(data_raw);
lfp = data(skip_n:end);

lfp_d = downsample(lfp,20);%x1000 mV to V fix
nfft=1024;fs=1000;
figure(1);plot(lfp_d*1e3);
[pxx,f] = pwelch(lfp_d,nfft,0,nfft,fs);

 figure(2);plot(f,pxx*1e3);
 set(gca, 'YScale', 'log');
 set(gca, 'XScale', 'log');
 %return
 
 %Check on voltage traces to be sure it is clamped
 %voltage_raw = h5read('../outputECP/cell_vars.h5','/report/BLA/data');
 %plot(voltage_raw(1,:))