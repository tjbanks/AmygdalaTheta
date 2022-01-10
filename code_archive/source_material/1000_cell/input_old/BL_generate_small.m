clear all;close all;clc;
seed= 12345789; %(random1) % control the sequencies of generation the random values.
%seed= 3344;%(random2)
%seed=271178;%random5 5801 %random4 96045 %random3 10858;
rng(seed);

tic

%INSERT INTERNEURONS EVENLY TO FILL VOLUME, 120 SPOTS CURRENTLY
AP=0.6;%%1.4;
LAT=0.6;%1.4;
DV=0.6;%1.4;
AP_I=0.6;%1.4;   %%%%inner dimension
LAT_I=0.6;%1.4;
DV_I=0.6;%1.4;


% SET NUERON DENSITY.
bla_volume = 610;       % mm^3 from human being Rubinow et al 2016
nrn_num = 12.68e6;  % from human being Rubinow et al 2016
dens = nrn_num/bla_volume;  % cell density mm^-3.
avg_dist = round((bla_volume*1e9/nrn_num)^(1/3));   % grid length = 36 um.
dmin = 0.025;   % minimum distance between soma.


upscale=1;
PCELL=900*upscale;

INTCELL=100*upscale;   %%according to Drew, only 10% is FS.
NCELL=PCELL+INTCELL;

TypeA=640*upscale;Anone=110*upscale;ADA=190*upscale;ANE=120*upscale; ADANE=220*upscale;    %%70,110,79,141  %%%%%according to Drew data, A took up ~65% of Pyr
TypeB=0;%TypeB=240; Bnone=40;BDA=65;BNE=48; BDANE=87;
TypeC=260*upscale; Cnone=57*upscale;CDA=64*upscale;CNE=50*upscale; CDANE=89*upscale;    %%27,44,32,57

% volume = NCELL/dens;
% scale = (volume/(AP*LAT*DV))^(1/3);
% scale = 1; %%%%for large network, no need to scale
% AP = AP*scale;
% LAT = LAT*scale;
% DV = DV*scale;


n=0;
m=0;
xI=(0.05:0.05:AP_I)';
yI=(0.05:0.05:LAT_I)';
zI=(0.05:0.05:DV_I)';
loca=[];loc=[];
for z=dmin:dmin:DV;
    for y=dmin:dmin:LAT;              % 0.05 seems the minimal we could make, if less than this, would be too small for placing one neuron
        for x=dmin:dmin:AP;
            if any(abs(xI-x)<=dmin/2) && any(abs(yI-y)<=dmin/2) && any(abs(zI-z)<=dmin/2)
            %if ((any(uint16(xI.*100)==uint16(x.*100)))&&(any(uint16(yI.*100)==uint16(y.*100)))&&(any(uint16(zI.*100)==uint16(z.*100))))
                m=m+1;
                loca(m,1)=x;           % for put interneuron
                loca(m,2)=y;
                loca(m,3)=z;
            else
                n=n+1;
                loc(n,1)=x;              % for put primary cell
                loc(n,2)=y;
                loc(n,3)=z;
            end
        end
    end
end

p = randperm(n,PCELL);
p=p';
NM=zeros(NCELL,1);

location=zeros(NCELL,3);
type=zeros(NCELL,1);
%put in TypeA
% for i=1:Anone;
%     location(i,1)=loc(p(i),1);
%     location(i,2)=loc(p(i),2);
%     location(i,3)=loc(p(i),3);
%     location(i,4)=1;
%     type(i,1)=1;
%     NM(i,1)=0;      %%%none type
% end
index=[1:Anone];
location(index,1:3)=loc(p(index),1:3);
type(index,1)=1;
NM(index,1)=0;



% for i=Anone+1:Anone+ADA;
%     location(i,1)=loc(p(i),1);
%     location(i,2)=loc(p(i),2);
%     location(i,3)=loc(p(i),3);
%     location(i,4)=1;
%     type(i,1)=1;
%     NM(i,1)=1;     %%%DA type
% end
index=[Anone+1:Anone+ADA];
location(index,1:3)=loc(p(index),1:3);
type(index,1)=1;
NM(index,1)=1;

% for i=Anone+ADA+1:Anone+ADA+ANE;
%     location(i,1)=loc(p(i),1);
%     location(i,2)=loc(p(i),2);
%     location(i,3)=loc(p(i),3);
%     location(i,4)=1;
%     type(i,1)=1;
%     NM(i,1)=2;   %%%NE type
% end

index=[Anone+ADA+1:Anone+ADA+ANE];
location(index,1:3)=loc(p(index),1:3);
type(index,1)=1;
NM(index,1)=2;


% for i=Anone+ADA+ANE+1:TypeA;
%     location(i,1)=loc(p(i),1);
%     location(i,2)=loc(p(i),2);
%     location(i,3)=loc(p(i),3);
%     location(i,4)=1;
%     type(i,1)=1;
%     NM(i,1)=3; %%%DANE type
% end

index=[Anone+ADA+ANE+1:TypeA];
location(index,1:3)=loc(p(index),1:3);
type(index,1)=1;
NM(index,1)=3;

% %put in TypeB
% for i=TypeA+1:TypeA+Bnone;
% location(i,1)=loc(p(i),1);
% location(i,2)=loc(p(i),2);
% location(i,3)=loc(p(i),3);
% location(i,4)=2;
% type(i,1)=6;
% NM(i,1)=0;
% end
%
% for i=TypeA+Bnone+1:TypeA+Bnone+BDA;
% location(i,1)=loc(p(i),1);
% location(i,2)=loc(p(i),2);
% location(i,3)=loc(p(i),3);
% location(i,4)=2;
% type(i,1)=6;
% NM(i,1)=1;
% end
%
% for i=TypeA+Bnone+BDA+1:TypeA+Bnone+BDA+BNE;
% location(i,1)=loc(p(i),1);
% location(i,2)=loc(p(i),2);
% location(i,3)=loc(p(i),3);
% location(i,4)=2;
% type(i,1)=6;
% NM(i,1)=2;
% end
%
% for i=TypeA+Bnone+BDA+BNE+1:TypeA+TypeB;
% location(i,1)=loc(p(i),1);
% location(i,2)=loc(p(i),2);
% location(i,3)=loc(p(i),3);
% location(i,4)=2;
% type(i,1)=6;
% NM(i,1)=3;
% end

%put in TypeC
% for i=TypeA+TypeB+1:TypeA+TypeB+Cnone;
%     location(i,1)=loc(p(i),1);
%     location(i,2)=loc(p(i),2);
%     location(i,3)=loc(p(i),3);
%     location(i,4)=3;
%     type(i,1)=10;
%     NM(i,1)=0;
% end

index=[TypeA+TypeB+1:TypeA+TypeB+Cnone];
location(index,1:3)=loc(p(index),1:3);



type(index,1)=10;
NM(index,1)=0;











% for i=TypeA+TypeB+Cnone+1:TypeA+TypeB+Cnone+CDA;
%     location(i,1)=loc(p(i),1);
%     location(i,2)=loc(p(i),2);
%     location(i,3)=loc(p(i),3);
%     location(i,4)=3;
%     type(i,1)=10;
%     NM(i,1)=1;
% end
index=[TypeA+TypeB+Cnone+1:TypeA+TypeB+Cnone+CDA];
location(index,1:3)=loc(p(index),1:3);
type(index,1)=10;
NM(index,1)=1;


% for i=TypeA+TypeB+Cnone+CDA+1:TypeA+TypeB+Cnone+CDA+CNE;
%     location(i,1)=loc(p(i),1);
%     location(i,2)=loc(p(i),2);
%     location(i,3)=loc(p(i),3);
%     location(i,4)=3;
%     type(i,1)=10;
%     NM(i,1)=2;
% end
index=[TypeA+TypeB+Cnone+CDA+1:TypeA+TypeB+Cnone+CDA+CNE];
location(index,1:3)=loc(p(index),1:3);
type(index,1)=10;
NM(index,1)=2;

% for i=TypeA+TypeB+Cnone+CDA+CNE+1:TypeA+TypeB+TypeC;
%     location(i,1)=loc(p(i),1);
%     location(i,2)=loc(p(i),2);
%     location(i,3)=loc(p(i),3);
%     location(i,4)=3;
%     type(i,1)=10;
%     NM(i,1)=3;
% end
index=[TypeA+TypeB+Cnone+CDA+CNE+1:TypeA+TypeB+TypeC];
location(index,1:3)=loc(p(index),1:3);
type(index,1)=10;
NM(index,1)=3;


%put in Interneuorn
p=randperm(m,INTCELL);
p=p';

% for i=1:INTCELL;
%     location(i+PCELL,1)=loca(p(i),1);
%     location(i+PCELL,2)=loca(p(i),2);
%     location(i+PCELL,3)=loca(p(i),3);
%     location(i+PCELL,4)=4;
%     type(i+PCELL,1)=100;
% end

index=[1:INTCELL];
location(index+PCELL,1:3)=loca(p(index),1:3);
type(index+PCELL,1)=100;
NM(index+PCELL,1)=0;

%save('location','location');
%%%%select location on the edges%%%%%%
% loc_edge(:,1:3)=location(find(location(:,1)>2.2|location(:,1)<0.3|location(:,2)>0.7|location(:,2)<0.3),:);
% 
% loc_edge_ID(:,1)=find(location(:,1)>2.2|location(:,1)<0.3|location(:,2)>0.7|location(:,2)<0.3);


figure
x_p=location(1:PCELL,1);y_p=location(1:PCELL,2);z_p=location(1:PCELL,3);
x_I=location(PCELL+1:NCELL,1);y_I=location(PCELL+1:NCELL,2);z_I=location(PCELL+1:NCELL,3);
scatter3(x_p,y_p,z_p,'red','filled');hold on;
scatter3(x_I,y_I,z_I,'filled','blue');hold on;


% figure 
% scatter3(loc_edge(:,1),loc_edge(:,2),loc_edge(:,3),'blue','filled');hold on;

% a=[PCELL+1:NCELL]';b=num2str(a);c=cellstr(b);
% dx=0.001;dy=0.001;dz=0.001;
% text(x_I+dx,y_I+dy,z_I+dz,c);

%%%calculate inter-soma distances
x_dis=pdist(location(:,1),@(x,y) x-y);
x_dis=x_dis.^2;
y_dis=pdist(location(:,2),@(x,y) x-y);
y_dis=y_dis.^2;
z_dis=pdist(location(:,3),@(x,y) x-y);
z_dis=z_dis.^2;

 dis=sqrt(x_dis+y_dis+z_dis);

clear -regexp _dis;
% clear -regexp loc;
dis_matrix=zeros(NCELL,NCELL);

for i=1:NCELL-1;
    ind_start=0.5*(2*NCELL-i)*(i-1)+1;
    ind_end=0.5*(2*NCELL-1-i)*(i);
    dis_matrix(i,i+1:end)=dis(ind_start:ind_end);
%     dis_matrix(i+1:end,i)=dis(ind_start:ind_end);
   
end

clear dis;
% tic
% for i=1:PCELL;
%     cellID=find(dis_matrix(i,1:PCELL)<=0.050);%%%within 50um
%     P2P_dis_50(i,1:length(cellID))=cellID; 
%     cellID=find(dis_matrix(i,1:PCELL)<=0.100&dis_matrix(i,1:PCELL)>0.050);%%%within 50um-100um
%     P2P_dis_100(i,1:length(cellID))=cellID;%%%within 100um-200um
%     cellID=find(dis_matrix(i,1:PCELL)<=0.200&dis_matrix(i,1:PCELL)>0.100);
%     P2P_dis_200(i,1:length(cellID))=cellID; 
% end
% toc

slice_dimension=0.3*10e20;%%%%%slice dimension is defined according to Wooduff 2007
%%%%for I-P pair%%%%%%%
[a_300,b_300]=find((dis_matrix(1:PCELL,PCELL+1:NCELL)-slice_dimension)<=0.0001&dis_matrix(1:PCELL,PCELL+1:NCELL)>0);
IP_dis_300_pair=[a_300,b_300+PCELL];  %%%%IP pair within 300um.
IP_dis_300_pair=sortrows(IP_dis_300_pair,2);
[n, bin] = histc((IP_dis_300_pair(:,2)), unique(IP_dis_300_pair(:,2)));
n_cum=cumsum(n);IP_dis_300=zeros(INTCELL,max(n));
IP_dis_300(IP_dis_300_pair(n_cum(1),2)-PCELL,1:n(1))=IP_dis_300_pair(0+1:n_cum(1),1);

for i=2:length(n);
    IP_dis_300(IP_dis_300_pair(n_cum(i),2)-PCELL,1:n(i))=IP_dis_300_pair(n_cum(i-1)+1:n_cum(i),1);
end   %%%%IP_dis_300 stores for each I cells, corresponding P cells within 300um to connect possibly.

% check each cells connectivity
% for i=1:INTCELL
% I2P(i)=numel(find(IP_dis_300(i,:)>0));
% end


IP_dis_300_pair=sortrows(IP_dis_300_pair,1);
[n, bin] = histc((IP_dis_300_pair(:,1)), unique(IP_dis_300_pair(:,1)));
n_cum=cumsum(n);PI_dis_300=zeros(PCELL,max(n));
PI_dis_300(IP_dis_300_pair(n_cum(1),1),1:n(1))=IP_dis_300_pair(0+1:n_cum(1),2);

for i=2:length(n);
    PI_dis_300(IP_dis_300_pair(n_cum(i),1),1:n(i))=IP_dis_300_pair(n_cum(i-1)+1:n_cum(i),2);
end   %%%%PI_dis_300 stores for each P cells, corresponding I cells within 300um to connect possibly.
save('PI_dis_300','PI_dis_300');
clear PI_dis_300;

% check each cells connectivity
% for i=1:PCELL
% P2I(i)=numel(find(PI_dis_300(i,:)>0));
% end
 clear a_300; clear a_300; clear IP_dis_300_pair;

%%%%filled in symmetric information for holistic II and PP pairs
[n,m]=size(dis_matrix);
dis_matrix_full=dis_matrix'+dis_matrix;
dis_matrix_full(1:n+1:end)=diag(dis_matrix);
dis_matrix=dis_matrix_full;
clear dis_matrix_full;


%%%%for I-I pair%%%%%%%
[a_300,b_300]=find((dis_matrix(PCELL+1:NCELL,PCELL+1:NCELL)-slice_dimension)<=0.0001&dis_matrix(PCELL+1:NCELL,PCELL+1:NCELL)>0);
I2I_dis_300_pair=[a_300+PCELL,b_300+PCELL]; %%%%II pair within 300um.

I2I_dis_300_pair=sortrows(I2I_dis_300_pair,2);
[n, bin] = histc((I2I_dis_300_pair(:,2)), unique(I2I_dis_300_pair(:,2)));
n_cum=cumsum(n);II_dis_300=zeros(INTCELL,max(n));
II_dis_300(I2I_dis_300_pair(n_cum(1),2)-PCELL,1:n(1))=I2I_dis_300_pair(0+1:n_cum(1),1);

for i=2:length(n);
    II_dis_300(I2I_dis_300_pair(n_cum(i),2)-PCELL,1:n(i))=I2I_dis_300_pair(n_cum(i-1)+1:n_cum(i),1);
end   %%%%II_dis_300 stores for each I cells, corresponding I cells within 300um to connect possibly.
clear a_300; clear b_300; clear I2I_dis_300_pair;

%%%%for P-P pair%%%%%%%    PP connectivity is defined according to Marios' email
%%%within 50um     
P2P_slice_offset=0.00;
[a_50,b_50]=find((dis_matrix(1:PCELL,1:PCELL)-(0.05+P2P_slice_offset))<=0.0001&dis_matrix(1:PCELL,1:PCELL)>0);%%%within 50um
P2P_dis_50_pair=[a_50,b_50];
P2P_dis_50_pair=sortrows(P2P_dis_50_pair,2);
[n, bin] = histc((P2P_dis_50_pair(:,2)), unique(P2P_dis_50_pair(:,2)));
n_cum=cumsum(n);PP_dis_50=zeros(PCELL,max(n));
PP_dis_50(P2P_dis_50_pair(n_cum(1),2),1:n(1))=P2P_dis_50_pair(0+1:n_cum(1),1);

for i=2:length(n);
    PP_dis_50(P2P_dis_50_pair(n_cum(i),2),1:n(i))=P2P_dis_50_pair(n_cum(i-1)+1:n_cum(i),1);
end   %%%%PP_dis_50 stores for each P cells, corresponding P cells within 50um to connect possibly.
clear a_50; clear b_50; clear P2P_dis_50_pair;

%%%within 100um    
[a_100,b_100]=find((dis_matrix(1:PCELL,1:PCELL)-(0.1+P2P_slice_offset))<=0.0001&dis_matrix(1:PCELL,1:PCELL)>(0.05+P2P_slice_offset));%%%within 50um-100um
P2P_dis_100_pair=[a_100,b_100];
P2P_dis_100_pair=sortrows(P2P_dis_100_pair,2);
[n, bin] = histc((P2P_dis_100_pair(:,2)), unique(P2P_dis_100_pair(:,2)));
n_cum=cumsum(n);PP_dis_100=zeros(PCELL,max(n));
PP_dis_100(P2P_dis_100_pair(n_cum(1),2),1:n(1))=P2P_dis_100_pair(0+1:n_cum(1),1);

for i=2:length(n);
    PP_dis_100(P2P_dis_100_pair(n_cum(i),2),1:n(i))=P2P_dis_100_pair(n_cum(i-1)+1:n_cum(i),1);
end   %%%%PP_dis_100 stores for each P cells, corresponding P cells within 100um to connect possibly.
    clear a_100; clear b_100; clear P2P_dis_100_pair;
 
%%%within 200um 
[a_200,b_200]=find((dis_matrix(1:PCELL,1:PCELL)-(0.2+P2P_slice_offset))<=0.0001&dis_matrix(1:PCELL,1:PCELL)>(0.1+P2P_slice_offset));%%%within 100um-200um
P2P_dis_200_pair=[a_200,b_200];
P2P_dis_200_pair=sortrows(P2P_dis_200_pair,2);
[n, bin] = histc((P2P_dis_200_pair(:,2)), unique(P2P_dis_200_pair(:,2)));
n_cum=cumsum(n);PP_dis_200=zeros(PCELL,max(n));
PP_dis_200(P2P_dis_200_pair(n_cum(1),2),1:n(1))=P2P_dis_200_pair(0+1:n_cum(1),1);

for i=2:length(n);
    PP_dis_200(P2P_dis_200_pair(n_cum(i),2),1:n(i))=P2P_dis_200_pair(n_cum(i-1)+1:n_cum(i),1);
end   %%%%PP_dis_200 stores for each P cells, corresponding P cells within 200um to connect possibly.
    clear a_200; clear b_200; clear P2P_dis_200_pair; 


 %%%within 300um     
P2P_slice_offset=0.00;
[a_300,b_300]=find((dis_matrix(1:PCELL,1:PCELL)-(0.3+P2P_slice_offset))<=0.0001&dis_matrix(1:PCELL,1:PCELL)>(0.2+P2P_slice_offset));%%%within 300um
P2P_dis_300_pair=[a_300,b_300];
P2P_dis_300_pair=sortrows(P2P_dis_300_pair,2);
[n, bin] = histc((P2P_dis_300_pair(:,2)), unique(P2P_dis_300_pair(:,2)));
n_cum=cumsum(n);PP_dis_300=zeros(PCELL,max(n));
PP_dis_300(P2P_dis_300_pair(n_cum(1),2),1:n(1))=P2P_dis_300_pair(0+1:n_cum(1),1);

for i=2:length(n);
    PP_dis_300(P2P_dis_300_pair(n_cum(i),2),1:n(i))=P2P_dis_300_pair(n_cum(i-1)+1:n_cum(i),1);
end   %%%%PP_dis_300 stores for each P cells, corresponding P cells within 300um to connect possibly.
clear a_300; clear b_300; clear P2P_dis_300_pair;   


 %%%within 400um     
P2P_slice_offset=0.00;
[a_400,b_400]=find((dis_matrix(1:PCELL,1:PCELL)-(0.4+P2P_slice_offset))<=0.0001&dis_matrix(1:PCELL,1:PCELL)>(0.3+P2P_slice_offset));%%%within 400um
P2P_dis_400_pair=[a_400,b_400];
P2P_dis_400_pair=sortrows(P2P_dis_400_pair,2);
[n, bin] = histc((P2P_dis_400_pair(:,2)), unique(P2P_dis_400_pair(:,2)));
n_cum=cumsum(n);PP_dis_400=zeros(PCELL,max(n));
PP_dis_400(P2P_dis_400_pair(n_cum(1),2),1:n(1))=P2P_dis_400_pair(0+1:n_cum(1),1);

for i=2:length(n);
    PP_dis_400(P2P_dis_400_pair(n_cum(i),2),1:n(i))=P2P_dis_400_pair(n_cum(i-1)+1:n_cum(i),1);
end   %%%%PP_dis_400 stores for each P cells, corresponding P cells within 400um to connect possibly.
clear a_400; clear b_400; clear P2P_dis_400_pair;   

 %%%within 500um     
P2P_slice_offset=0.00;
[a_500,b_500]=find((dis_matrix(1:PCELL,1:PCELL)-(0.5+P2P_slice_offset))<=0.0001&dis_matrix(1:PCELL,1:PCELL)>(0.4+P2P_slice_offset));%%%within 500um
P2P_dis_500_pair=[a_500,b_500];
P2P_dis_500_pair=sortrows(P2P_dis_500_pair,2);
[n, bin] = histc((P2P_dis_500_pair(:,2)), unique(P2P_dis_500_pair(:,2)));
n_cum=cumsum(n);PP_dis_500=zeros(PCELL,max(n));
PP_dis_500(P2P_dis_500_pair(n_cum(1),2),1:n(1))=P2P_dis_500_pair(0+1:n_cum(1),1);

for i=2:length(n);
    PP_dis_500(P2P_dis_500_pair(n_cum(i),2),1:n(i))=P2P_dis_500_pair(n_cum(i-1)+1:n_cum(i),1);
end   %%%%PP_dis_500 stores for each P cells, corresponding P cells within 50um to connect possibly.
clear a_500; clear b_500; clear P2P_dis_500_pair;   


 %%%within all the rest range     
P2P_slice_offset=0.00;
[a_600,b_600]=find((dis_matrix(1:PCELL,1:PCELL)-(0.6*10e20+P2P_slice_offset))<=0.0001&dis_matrix(1:PCELL,1:PCELL)>(0.5+P2P_slice_offset));%%%within 500um
P2P_dis_600_pair=[a_600,b_600];
P2P_dis_600_pair=sortrows(P2P_dis_600_pair,2);
[n, bin] = histc((P2P_dis_600_pair(:,2)), unique(P2P_dis_600_pair(:,2)));
n_cum=cumsum(n);PP_dis_600=zeros(PCELL,max(n));
PP_dis_600(P2P_dis_600_pair(n_cum(1),2),1:n(1))=P2P_dis_600_pair(0+1:n_cum(1),1);

for i=2:length(n);
    PP_dis_600(P2P_dis_600_pair(n_cum(i),2),1:n(i))=P2P_dis_600_pair(n_cum(i-1)+1:n_cum(i),1);
end   %%%%PP_dis_600 stores for each P cells, corresponding P cells within 50um to connect possibly.
clear a_600; clear b_600; clear P2P_dis_600_pair;    clear dis_matrix; clear n; clear n_cum; clear bin;
    

    
%%%%make connection%%%%%%

%%%make P2P first
P2P_conn_scale=1; P2P_conn_scale_local=1;   
P2P_conn_50=0.02*P2P_conn_scale_local;  %%%connectivity within 50um
P2P_conn_100=0.02*P2P_conn_scale_local;  %%%connectivity within 50-100um
P2P_conn_200=0.02*P2P_conn_scale_local;  %%%connectivity within 100-200um
P2P_conn_300=0.02*P2P_conn_scale;  %%%connectivity within 200-300um
P2P_conn_400=0.02*P2P_conn_scale;  %%%connectivity within 300-400um
P2P_conn_500=0.02*P2P_conn_scale;  %%%connectivity within 400-500um
P2P_conn_600=0.02*P2P_conn_scale;  %%%connectivity within 500-600um

%%%%make P2P within 50um
index_nonezero=find(PP_dis_50);
index_nonezero = index_nonezero(randperm(length(index_nonezero)));  %%%sort found index randomly
index=index_nonezero(1:floor(length(index_nonezero)*P2P_conn_50));
[PP_50_pre col]=ind2sub(size(PP_dis_50),index);
active_syn_P2P_pair_50=[PP_50_pre,PP_dis_50(index)]; clear PP_50_pre;clear PP_50_pre; clear PP_dis_50;

%%%%make P2P within 100um
index_nonezero=find(PP_dis_100);
index_nonezero = index_nonezero(randperm(length(index_nonezero)));  %%%sort found index randomly
index=index_nonezero(1:floor(length(index_nonezero)*P2P_conn_100));
[PP_100_pre col]=ind2sub(size(PP_dis_100),index);
active_syn_P2P_pair_100=[PP_100_pre,PP_dis_100(index)]; clear PP_100_pre;clear PP_100_pre;clear PP_dis_100;

%%%%make P2P within 200um
index_nonezero=find(PP_dis_200);
index_nonezero = index_nonezero(randperm(length(index_nonezero)));  %%%sort found index randomly
index=index_nonezero(1:floor(length(index_nonezero)*P2P_conn_200));
[PP_200_pre col]=ind2sub(size(PP_dis_200),index);
active_syn_P2P_pair_200=[PP_200_pre,PP_dis_200(index)]; clear PP_200_pre; clear col;clear PP_dis_200;

%%%%make P2P within 300um
index_nonezero=find(PP_dis_300);
index_nonezero = index_nonezero(randperm(length(index_nonezero)));  %%%sort found index randomly
index=index_nonezero(1:floor(length(index_nonezero)*P2P_conn_300));
[PP_300_pre col]=ind2sub(size(PP_dis_300),index);
active_syn_P2P_pair_300=[PP_300_pre,PP_dis_300(index)]; clear PP_300_pre; clear col;clear PP_dis_300;

%%%%make P2P within 400um
index_nonezero=find(PP_dis_400);
index_nonezero = index_nonezero(randperm(length(index_nonezero)));  %%%sort found index randomly
index=index_nonezero(1:floor(length(index_nonezero)*P2P_conn_400));
[PP_400_pre col]=ind2sub(size(PP_dis_400),index);
active_syn_P2P_pair_400=[PP_400_pre,PP_dis_400(index)]; clear PP_400_pre; clear col;clear PP_dis_400;

%%%%make P2P within 500um
index_nonezero=find(PP_dis_500);
index_nonezero = index_nonezero(randperm(length(index_nonezero)));  %%%sort found index randomly
index=index_nonezero(1:floor(length(index_nonezero)*P2P_conn_500));
[PP_500_pre col]=ind2sub(size(PP_dis_500),index);
active_syn_P2P_pair_500=[PP_500_pre,PP_dis_500(index)]; clear PP_500_pre; clear col;clear PP_dis_500;

%%%%make P2P within 600um
index_nonezero=find(PP_dis_600);
index_nonezero = index_nonezero(randperm(length(index_nonezero)));  %%%sort found index randomly
index=index_nonezero(1:floor(length(index_nonezero)*P2P_conn_600));
[PP_600_pre col]=ind2sub(size(PP_dis_600),index);
active_syn_P2P_pair_600=[PP_600_pre,PP_dis_600(index)]; clear PP_600_pre; clear col;clear PP_dis_600;


active_syn_P2P_pair=[active_syn_P2P_pair_50;active_syn_P2P_pair_100;active_syn_P2P_pair_200;active_syn_P2P_pair_300;...
    active_syn_P2P_pair_400;active_syn_P2P_pair_500;active_syn_P2P_pair_600];
clear active_syn_P2P_pair_50;clear active_syn_P2P_pair_100;clear active_syn_P2P_pair_200;
clear active_syn_P2P_pair_300;clear active_syn_P2P_pair_400;clear active_syn_P2P_pair_500;clear active_syn_P2P_pair_600;

active_syn_P2P_pair=sortrows(active_syn_P2P_pair,1);

[n, bin] = histc((active_syn_P2P_pair(:,1)), unique(active_syn_P2P_pair(:,1)));
n_cum=cumsum(n);active_syn_P2P=zeros(PCELL,max(n));
active_syn_P2P(active_syn_P2P_pair(n_cum(1),1),1:n(1))=active_syn_P2P_pair(0+1:n_cum(1),2);

for i=2:length(n);
    active_syn_P2P(active_syn_P2P_pair(n_cum(i),1),1:n(i))=active_syn_P2P_pair(n_cum(i-1)+1:n_cum(i),2);
end   %%%%active_syn_P2P stores active P2P, row and col corresponding to Pre and Post;
    clear active_syn_P2P_pair;
    
    
    
%%%%check minimum dist%%%%%
% min_dist=Inf;
% for p1 = 1:length(location);
%     for p2 = 1:length(location);
%         d = pdist([location(p1,1:3);location(p2,1:3)]);
%         if d < min_dist && p1 ~= p2
%             min_p1 = p1;
%             min_p2 = p2;
%             min_dist = d;
%         end
%     end
% end

%%%%check distance info %%%%%
% for p1 = 1:length(location);
%     for p2 = p1+1:length(location);
% %         if  p1 ~= p2
%         dis(p1,p2) = pdist([location(p1,1:3);location(p2,1:3)]);
%             %min_dist = d;
% %         else
%             
% %         end
%     end
% end

% netsheet=fopen('netsheet','w');
%
% for i=1:size(type,1);
% fprintf(netsheet, '%f\n', type(i,1));
% end
% for i=1:size(location,1)
% fprintf(netsheet, '%f\n', location(i,1));
% end
% for i=1:size(location,1)
% fprintf(netsheet, '%f\n', location(i,2));
% end
% for i=1:size(location,1)
% fprintf(netsheet, '%f\n', location(i,3));
% end
% for i=1:size(layer,1)
%   fprintf(netsheet, '%f\n', layer(i,1));
% end
%  fclose(netsheet);

%%%%------SETTING UP WORKSHEET OF CONNECTIONS HERE------------
%%%%--- to form GAP JUNCTION among neigbouring ITNs first%%%
n=0;m=0;r=0;
% GAP_matrix=zeros(INTCELL,INTCELL);
% for i=1+PCELL:PCELL+INTCELL-1;   %%%source ID
%     for j=i+1:PCELL+INTCELL;   %%%target ID  no replicate
%         dist = pdist([location(i,1:3);location(j,1:3)]);
%         if 1 % dist<=0.32*scale;  %%%%%%limiting the GP within 300um slices. ***used to be0.3*1000000
%             r=r+1;
%             if unifrnd(0,100) <= 8                %8;
%                 GAP_matrix(i-PCELL,j-PCELL)=1;
%                 GAP_matrix(j-PCELL,i-PCELL)=1;
%                 n=n+1;
%                 GAPID(n,1)=i;      %%%%%%to store GAP pairs;no replicate
%                 GAPID(n,2)=j;
%             else
%                 m=m+1;
%                 GAPnoID(m,1)=i;    %%%%%%to store noGAP pairs;no replicate
%                 GAPnoID(m,2)=j;
%             end
%         else
%             m=m+1;
%             GAPnoID(m,1)=i;    %%%%%%to store noGAP pairs;no replicate
%             GAPnoID(m,2)=j;
%         end
%     end
% end


conn_factor=1;   %%%%consider the size-effect

GAP_conn=0.08*conn_factor;
active_syn_GAP1=zeros(INTCELL,ceil(INTCELL*GAP_conn));
active_syn_GAP2=zeros(INTCELL,ceil(INTCELL*GAP_conn));
post_nouni=zeros(INTCELL,INTCELL-ceil(INTCELL*GAP_conn));
% for i=1:INTCELL-1;
% post_cadnidate=[i+1:INTCELL]+PCELL;
% post_temp=randperm(numel(post_cadnidate),ceil(numel(post_cadnidate)*GAP_conn));
% post=post_cadnidate(post_temp);
% post_nouni_temp=setdiff(post_cadnidate,post);
% post_nouni(i,1:length(post_nouni_temp))=post_nouni_temp;
% active_syn_GAP1(i,1:length(post(post~=i)))=post(post~=i);
% end
for i=1:INTCELL-1;                 %%%%only select within 300um
  index=find(II_dis_300(i,:)>i+PCELL+1);
  index = index(randperm(length(index)));  %%%sort found index randomly
  index_gap=index(1:floor(length(index)*GAP_conn));  %%%randomly select 8%
  active_syn_GAP1(i,1:length(index_gap))=II_dis_300(i,index_gap);
  
  index_nouni=index(floor(length(index)*GAP_conn)+1:end);
  post_nouni(i,1:length(index_nouni))=II_dis_300(i,index_nouni);  %%%cells are not coupled within 300um
end
[row,col]=ind2sub(size(active_syn_GAP1),find(active_syn_GAP1>0));
colmax=max(col);
active_syn_GAP1(:,colmax+1:end)=[]; clear row; clear col;  %%%%tailor reduntant matrix;

[row,col]=ind2sub(size(post_nouni),find(post_nouni>0));
colmax=max(col);
post_nouni(:,colmax+1:end)=[]; clear row; clear col;  %%%%tailor reduntant matrix;
  




for i=1+PCELL:INTCELL+PCELL;                                       %%%%reverse pre and post
 [a,b]=find(active_syn_GAP1==i);
 active_syn_GAP2(i-PCELL,1:length(a))=a+PCELL;
end
[row,col]=ind2sub(size(active_syn_GAP2),find(active_syn_GAP2>0));
colmax=max(col);
active_syn_GAP2(:,colmax+1:end)=[]; clear row; clear col;  %%%%tailor reduntant matrix;

active_syn_GAP=[active_syn_GAP1,active_syn_GAP2]; %%%%active_syn_GAP1 and active_syn_GAP2 should be symmetric
 clear active_syn_GAP2;
%   GAPID=[];
%  for i=1:size(active_syn_GAP1,1)                    %%%%to form similar GAPID matrix as in original code
%      nn=numel(find(active_syn_GAP1(i,:)>0));
%      GAPID_temp=[repmat(i+PCELL,nn,1),active_syn_GAP1(i,1:nn)'];
%      GAPID=[GAPID;GAPID_temp];
%  end
 
 
%  GAPnoID=[];
%  for i=1:size(post_nouni,1)
%      nn=numel(find(post_nouni(i,:)>0));
%      GAPnoID_temp=[repmat(i+PCELL,nn,1),post_nouni(i,1:nn)'];
%      GAPnoID=[GAPnoID;GAPnoID_temp];
%  end
 

%%%%%select candidate pairs to form synaptic connections%%%%%%
% GAPpairnum=size(GAPID,1);
% GAPnopairnum=size(GAPnoID,1);


couple_uni_p=0.5*conn_factor;couple_bi_p=0.5*conn_factor;uncouple_uni_p=0.19*2*conn_factor;uncouple_bi_p=0.03*2*conn_factor;

% p=randperm(GAPpairnum,floor(GAPpairnum*couple_uni_p));
% p=p';%couple_uni_id=GAPID(p,:);p1=p;
% k=1;h=1;
% for i=1:GAPpairnum;
%     if any(p==i);
%         couple_uni_id(k,1:2)=GAPID(i,:);
%         couple_uni_id(k,3)=i;
%         k=k+1;
%     elseif unifrnd(0,100) <= couple_bi_p*100
%         couple_bi_id(h,1:2)=GAPID(i,:);
%         h=h+1;
%     end
% end

%%%%%coupl pairs
%%%%%select for couple-uni pair
aa=find(active_syn_GAP1);           %%%%active_syn_GAP1 already consider 300um
p=randperm(numel(aa),floor(numel(aa)*couple_uni_p));
p=p';%couple_uni_id=GAPID(p,:);p1=p;

couple_uni_post=active_syn_GAP1(aa(p));

[couple_uni_pre col]=ind2sub(size(active_syn_GAP1),aa(p));
couple_uni_pre=couple_uni_pre+PCELL;
couple_uni_id=[couple_uni_pre,couple_uni_post];


%%%%%select rest for couple-bi pair
aa_nopick=setdiff(aa,aa(p));
%aa_nopick=aa_nopick';
p=randperm(numel(aa_nopick),floor(numel(aa_nopick)*couple_bi_p));
p=p';%couple_uni_id=GAPID(p,:);p1=p;

couple_bi_post=active_syn_GAP1(aa_nopick(p));

[couple_bi_pre col]=ind2sub(size(active_syn_GAP1),aa_nopick(p));
couple_bi_pre=couple_bi_pre+PCELL;
couple_bi_id=[couple_bi_pre,couple_bi_post];

%%%%%uncoupl pairs
%%%%%select for uncouple-uni pair
aa=find(post_nouni);
p=randperm(numel(aa),floor(numel(aa)*uncouple_uni_p));
p=p';%couple_uni_id=GAPID(p,:);p1=p;

uncouple_uni_post=post_nouni(aa(p));

[uncouple_uni_pre col]=ind2sub(size(post_nouni),aa(p));
uncouple_uni_pre=uncouple_uni_pre+PCELL;
uncouple_uni_id=[uncouple_uni_pre,uncouple_uni_post];


%%%%%select rest for uncouple-bi pair
aa_nopick=setdiff(aa,aa(p));
%aa_nopick=aa_nopick';
p=randperm(numel(aa_nopick),floor(numel(aa_nopick)*uncouple_bi_p));
p=p';%couple_uni_id=GAPID(p,:);p1=p;

uncouple_bi_post=post_nouni(aa_nopick(p));

[uncouple_bi_pre col]=ind2sub(size(post_nouni),aa_nopick(p));
uncouple_bi_pre=uncouple_bi_pre+PCELL;
uncouple_bi_id=[uncouple_bi_pre,uncouple_bi_post];






% p=randperm(GAPnopairnum,floor(GAPnopairnum*uncouple_uni_p));
% p=p';%uncouple_uni_id=GAPnoID(p,:);p3=p;
% 
% k=1;h=1;
% for i=1:GAPnopairnum;
%     if any(p==i);
%         uncouple_uni_id(k,1:2)=GAPnoID(i,:);
%         uncouple_uni_id(k,3)=i;
%         k=k+1;
%     elseif unifrnd(0,100) <= uncouple_bi_p*100
%         uncouple_bi_id(h,1:2)=GAPnoID(i,:);
%         h=h+1;
%     end
% end

% p=randperm(GAPnopairnum,floor(GAPnopairnum*uncouple_bi_p));
% p=p';%uncouple_bi_id=GAPnoID(p,:);p4=p;

%%%%--- to form I-I synaptically among ITNs second%%%

%syn_matrix=zeros(INTCELL+PCELL,INTCELL+PCELL);
% active_syn_I2I=zeros(INTCELL,INTCELL*couple_uni_p);  %%%%consider maximum connectivity to define matrix size

temp=zeros(length(couple_uni_id),2);%% to store pre and post
for i=1:length(couple_uni_id);   %%%assign couple_uni
    if round(rand)==0;    %%%% distribute pre and post randomly
        %syn_matrix(couple_uni_id(i,1),couple_uni_id(i,2))=1;
        temp(i,1:2)=[couple_uni_id(i,1),couple_uni_id(i,2)];
    else
        %syn_matrix(couple_uni_id(i,2),couple_uni_id(i,1))=1;
        temp(i,1:2)=[couple_uni_id(i,2),couple_uni_id(i,1)];
    end   
end
temp_I2I_couple_uni=sortrows(temp,1);
% temp_unique=unique(temp(:,1));   %%%%get rid of repetitive pre
% for i=1:length(temp_unique);    
%     active_syn_I2I(temp_unique(i,1)-PCELL,1:length(temp(temp(:,1)==temp_unique(i,1),2)))=temp(temp(:,1)==temp_unique(i,1),2);
% end


temp=zeros(length(uncouple_uni_id),2);%% to store pre and post
for i=1:length(uncouple_uni_id);   %%%assign uncouple_uni
    if round(rand)==0;    %%%% distribute pre and post randomly
        %syn_matrix(uncouple_uni_id(i,1),uncouple_uni_id(i,2))=1;
         temp(i,1:2)=[uncouple_uni_id(i,1),uncouple_uni_id(i,2)];
    else
        %syn_matrix(uncouple_uni_id(i,2),uncouple_uni_id(i,1))=1;
         temp(i,1:2)=[uncouple_uni_id(i,2),uncouple_uni_id(i,1)];
    end
end

temp_I2I_uncouple_uni=sortrows(temp,1);
% temp_unique=unique(temp(:,1));   %%%%get rid of repetitive pre
% for i=1:length(temp_unique);    %%%%%keep filling
%     start_index=min(find(active_syn_I2I(temp_unique(i,1)-PCELL,:)==0));
%     end_index=start_index+length(temp(temp(:,1)==temp_unique(i,1),2))-1;
%     active_syn_I2I(temp_unique(i,1)-PCELL,start_index:end_index)=temp(temp(:,1)==temp_unique(i,1),2);
% end


temp_I2I_couple_bi=sortrows([couple_bi_id;[couple_bi_id(:,2),couple_bi_id(:,1)]],1);
%for i=1:length(couple_bi_id);   %%%assign couple_bi
    
    %syn_matrix(couple_bi_id(i,1),couple_bi_id(i,2))=1;
    %syn_matrix(couple_bi_id(i,2),couple_bi_id(i,1))=1;
%end


%for i=1:length(uncouple_bi_id);   %%%assign uncouple_bi
    %syn_matrix(uncouple_bi_id(i,1),uncouple_bi_id(i,2))=1;
    %syn_matrix(uncouple_bi_id(i,2),uncouple_bi_id(i,1))=1;
%end

temp_I2I_uncouple_bi=sortrows([uncouple_bi_id;[uncouple_bi_id(:,2),uncouple_bi_id(:,1)]],1);

temp_allI2I=[temp_I2I_couple_uni;temp_I2I_uncouple_uni;temp_I2I_couple_bi;temp_I2I_uncouple_bi];
temp_allI2I=sortrows(temp_allI2I,1);

% temp_allI2I_unique=unique(temp_allI2I(:,1));   %%%%get rid of repetitive pre
% for i=1:length(temp_allI2I_unique);    
%     start_index=1;
%     end_index=length(temp_allI2I(temp_allI2I(:,1)==temp_allI2I_unique(i,1)));
%      active_syn_I2I(temp_allI2I_unique(i,1)-PCELL,start_index:end_index)=temp_allI2I(temp_allI2I(:,1)==temp_allI2I_unique(i,1),2);
% end

[n, bin] = histc((temp_allI2I(:,1)), unique(temp_allI2I(:,1)));
n_cum=cumsum(n);active_syn_I2I=zeros(INTCELL,max(n));
active_syn_I2I(temp_allI2I(n_cum(1),1)-PCELL,1:n(1))=temp_allI2I(0+1:n_cum(1),2);

for i=2:length(n);
    active_syn_I2I(temp_allI2I(n_cum(i),1)-PCELL,1:n(i))=temp_allI2I(n_cum(i-1)+1:n_cum(i),2);
end   %%%%PP_dis_100 stores for each P cells, corresponding P cells within 100um to connect possibly.
    clear n; clear n_cum; clear bin; clear temp_allI2I; 
    clear -regexp temp_I2I; clear -regexp uncouple_;
    clear -regexp couple_;

%%%%check GJ connectivity%%%%%
%GAP_connectivity = 2*100*GAPpairnum/(INTCELL)/(INTCELL-1);

%I2I_overall_connection = 100*length(temp_allI2I)/INTCELL/(INTCELL-1);




%%%%--- to form I-P synaptically third%%%
k=1;m=1;I2P_uni_conn=0.34*conn_factor;P2I_uni_conn=0.18*conn_factor;I2P_bi_conn=0.3*conn_factor;
active_syn_I2P=zeros(INTCELL,ceil(PCELL*(I2P_uni_conn))); 
I2P_nouni=zeros(INTCELL,ceil(PCELL*(1-I2P_uni_conn))); 

%%%%I-P uni%%%%%
for i=1:INTCELL;                 %%%%only select within 300um
  index=find(IP_dis_300(i,:));
  index = index(randperm(length(index)));  %%%sort found index randomly
  index_IP=index(1:floor(length(index)*I2P_uni_conn));  %%%randomly select to connect
  active_syn_I2P(i,1:length(index_IP))=IP_dis_300(i,index_IP);
  
  index_noIP=index(floor(length(index)*I2P_uni_conn)+1:end);
  I2P_nouni(i,1:length(index_noIP))=IP_dis_300(i,index_noIP);  %%%cells are not connected within 300um
end
clear IP_dis_300;

[row,col]=ind2sub(size(active_syn_I2P),find(active_syn_I2P>0));
colmax=max(col);
active_syn_I2P(:,colmax+1:end)=[]; clear row; clear col;  %%%%tailor reduntant matrix;

[row,col]=ind2sub(size(I2P_nouni),find(I2P_nouni>0));
colmax=max(col);
I2P_nouni(:,colmax+1:end)=[]; clear row; clear col;  %%%%tailor reduntant matrix;




%%%%I-P uni%%%%%
% active_syn_I2P=zeros(INTCELL,ceil(PCELL*(I2P_uni_conn))); 
% for i=1+PCELL:PCELL+INTCELL;   %%%source ID      %%%uni I2P
%     post=randperm(PCELL,floor(PCELL*I2P_uni_conn));
%     active_syn_I2P(i-PCELL,1:length(post))=post;
%     post_nouni=setdiff([1:PCELL],post);
%     I2P_nouni(i-PCELL,1:length(post_nouni))=post_nouni;
% end
% 
% k=1;m=1;

%%%%P-I uni%%%%%
P2I_nouni=zeros(PCELL,ceil(INTCELL*(1-I2P_uni_conn*1))); %%%%
%active_syn_P2I=zeros(PCELL,ceil(INTCELL*(P2I_uni_conn))); 
active_syn_P2I_reversed=zeros(INTCELL,floor(size(I2P_nouni,2)*P2I_uni_conn));


%%%%%select for candidates of P2I within 300um
aa=find(I2P_nouni);
p=randperm(numel(aa),floor(numel(aa)*P2I_uni_conn));
p=p';%couple_uni_id=GAPID(p,:);p1=p;
p_nopick=setdiff([1:numel(aa)],p');

P2I_pre=I2P_nouni(aa(p));

[P2I_post col]=ind2sub(size(I2P_nouni),aa(p));
P2I_post=P2I_post+PCELL;
P2I_active_pair=[P2I_pre,P2I_post];
%P2I_active_pair=sortrows(P2I_active_pair,1)


for i=1:INTCELL;                                       
 [a,b]=find(P2I_active_pair(:,2)==i+PCELL);
 active_syn_P2I_reversed(i,1:length(a))=P2I_active_pair(a,1); % CAUTIONS: PRE POST is reversed
end

[row,col]=ind2sub(size(active_syn_P2I_reversed),find(active_syn_P2I_reversed>0));
colmax=max(col);
active_syn_P2I_reversed(:,colmax+1:end)=[]; clear row; clear col;  %%%%tailor reduntant matrix;

%%%%find P2I-bi-connection
P2I_I2P_no_pre=I2P_nouni(aa(p_nopick));
[P2I_I2P_post col]=ind2sub(size(I2P_nouni),aa(p_nopick));
P2I_I2P_post=P2I_I2P_post+PCELL;
P2I_I2P_inactive_pair=[P2I_I2P_no_pre,P2I_I2P_post];

clear p_nopick;

aa=size(P2I_I2P_inactive_pair,1);
p=randperm(aa,floor(aa*I2P_bi_conn));
p=p';%couple_uni_id=GAPID(p,:);p1=p;

P2I_I2P_bi_pair=P2I_I2P_inactive_pair(p,:);
%%%%add to previous I2P matrix
for i=1:INTCELL;                                       
 [a,b]=find(P2I_I2P_bi_pair(:,2)==i+PCELL);
 active_syn_I2P_bi(i,1:length(a))=P2I_I2P_bi_pair(a,1);%+PCELL;
end


%%%combine I<->P matrix%%%
active_syn_I2P=[active_syn_I2P,active_syn_I2P_bi];
active_syn_P2I_reversed=[active_syn_P2I_reversed,active_syn_I2P_bi];
clear active_syn_I2P_bi;


% %%%%to form matrix for candidates of P2I
% for i=1:PCELL;        %%%uni P2I
%     [a,b]=find(I2P_nouni==i);
%     I_temp=a';
%     %LENGTH_a(i,1)=length(a);
%     P2I_nouni(i,1:length(a))=I_temp+PCELL;
%     post_temp=randperm(length(a),floor(length(a)*P2I_uni_conn)); 
%     post=P2I_nouni(i,post_temp);
%     active_syn_P2I(i,1:length(post))=post;
%     notpost=setdiff(P2I_nouni(i,:),post);
%     notpost=notpost(notpost>0);
%     P2I_I2P_nouni(i,1:length(notpost))=notpost; 
%     
% %     post_temp=randperm(size(I2P_nouni,2),floor(size(I2P_nouni,2)*P2I_uni_conn));    %%%%select pre-P
% %     pre=I2P_nouni(i,pre_temp);
% %     pre_P2II2Pnouni=setdiff(I2P_nouni(i,:),pre);
% %     P2I_I2P_nouni(i,1:length(pre_P2II2Pnouni))=pre_P2II2Pnouni;
% %     active_syn_P2I_reversed(i,:)=pre;
% end
% 
% % for i=1:INTCELL;        %%%uni P2I 
% %     pre_temp=randperm(size(I2P_nouni,2),floor(size(I2P_nouni,2)*P2I_uni_conn));    %%%%select pre-P
% %     pre=I2P_nouni(i,pre_temp);
% %     pre_P2II2Pnouni=setdiff(I2P_nouni(i,:),pre);
% %     P2I_I2P_nouni(i,1:length(pre_P2II2Pnouni))=pre_P2II2Pnouni;
% %     active_syn_P2I_reversed(i,:)=pre;
% % end
% % 
% % for i=1:PCELL;                                       %%%%reverse matrix
% %  [a,b]=find(active_syn_P2I_reversed==i);
% %  active_syn_P2I(i,1:length(a))=a+PCELL;
% % end
% %  clear active_syn_P2I_reversed;  %%%% clear unintended matrix
%  
%  
%  
%  
% for i=1:size(P2I_I2P_nouni,1);        %%%bi P2I-I2P 
%     P2I_I2P_nouni_temp=P2I_I2P_nouni(i,:);
%     P2I_I2P_nouni_nozero_temp=P2I_I2P_nouni_temp(P2I_I2P_nouni_temp>0);
%     bi_temp=randperm(size(P2I_I2P_nouni_nozero_temp,2),floor(size(P2I_I2P_nouni_nozero_temp,2)*I2P_bi_conn));    %%%%select pre-P
%     
%     bi=P2I_I2P_nouni_nozero_temp(bi_temp);
%     P2I_I2P_bi(i,1:length(bi))=bi;
%     %temp=diff(cumsum(active_syn_P2I,2),1,2);
%    %[a,b]=find(active_syn_P2I'==0);
%                                            %%%%to fill in P2I matrix
% %     if unifrnd(0,100) <= 30   %%%25/(160-55-19)=29%
% %         syn_matrix(P2I_I2P_nouni(i,2),P2I_I2P_nouni(i,1))=1;
% %         syn_matrix(P2I_I2P_nouni(i,1),P2I_I2P_nouni(i,2))=1;
% %     end
% end
% for i=1:INTCELL;                                       %%%%reverse matrix
%  [a,b]=find(P2I_I2P_bi==i+PCELL);
%  active_syn_I2P_bi(i,1:length(a))=a;%+PCELL;
% end
% active_syn_P2I_bi=P2I_I2P_bi;
% 
% %%%combine I<->P matrix%%%
% active_syn_I2P=[active_syn_I2P,active_syn_I2P_bi];
% active_syn_P2I=[active_syn_P2I,active_syn_P2I_bi];


%%%%check both I-P connections simultaneously%%%%%%%
% l=0;n=0;m=0;o=0;
% for i=1:PCELL;   %%%source ID
%     for j=PCELL+1:PCELL+INTCELL;   %%%target ID  no replicate
%         if syn_matrix(i,j)&&~syn_matrix(j,i)        %%%uniP-I
%             l=l+1;
%         elseif ~syn_matrix(i,j)&&syn_matrix(j,i)    %%%uniI-P
%             n=n+1;
%         elseif syn_matrix(i,j)&&syn_matrix(j,i)     %%%biI-P,P-I
%             m=m+1;
%         else
%             o=o+1;                                  %%%no connect
%         end
%     end
% end
% 
% uniI2P_p=100*n/(PCELL*INTCELL);
% uniP2I_p=100*l/(PCELL*INTCELL);
% biIP_p=100*m/(PCELL*INTCELL);
% noIP_p=100*o/(PCELL*INTCELL);
% overallI2P_p=100-noIP_p;

%%%%--- to form P-P synaptically third%%%
% for i=1:PCELL;   %%%source ID      %%%uni
%     for j=1:PCELL;   %%%target ID
%         if i~=j;
%             if unifrnd(0,100) <= 2
%                 syn_matrix(i,j)=1;
%             end
%         end
%     end
% end





% % %%%%--- to disconnect I-P %%%
% for i=1+PCELL:PCELL+INTCELL;   %%%source ID      %%%uni
%     for j=1:PCELL;   %%%target ID
%         if i~=j;
%             if syn_matrix(i,j) > 0
%                 syn_matrix(i,j)=0;
%             end
%         end
%     end
% end

% %%%%--- to disconnect P-I %%%
% for i=1:PCELL;   %%%source ID      %%%uni
%     for j=1+PCELL:PCELL+INTCELL;   %%%target ID
%         if i~=j;
%             if syn_matrix(i,j) > 0
%                 syn_matrix(i,j)=0;
%             end
%         end
%     end
% end

% %%%%--- to disconnect I-I %%%
% for i=1+PCELL:PCELL+INTCELL;   %%%source ID      %%%uni
%     for j=1+PCELL:PCELL+INTCELL;   %%%target ID
%         if i~=j;
%             if syn_matrix(i,j) > 0
%                 syn_matrix(i,j)=0;
%             end
%         end
%     end
% end

%%%%% transfer ID to be used in NEURON %%%%
clear -regexp _pair; clear -regexp _nopick; clear p;
clear -regexp P2I_I2P_; clear -regexp _nouni; clear P2I_post; clear P2I_pre; clear -regexp _dis_; 
clear -regexp index; clear -regexp temp;
active_syn_GAP=active_syn_GAP-1;
active_syn_P2P=active_syn_P2P-1;
active_syn_I2I=active_syn_I2I-1;
active_syn_I2P=active_syn_I2P-1;
active_syn_P2I_reversed=active_syn_P2I_reversed-1;

%%%%%%%%%%Intrinsic weights%%%%%%%%%%%%
wE2E_mean = 5; wE2E_std = 2;
mu_E2E = log((wE2E_mean^2)/sqrt(wE2E_std+wE2E_mean^2));  %%%convert to
sigma_E2E = sqrt(log(wE2E_std/(wE2E_mean^2)+1));

wE2I_mean= 3;wE2I_std = 2;
mu_E2I = log((wE2I_mean^2)/sqrt(wE2I_std+wE2I_mean^2));  %%%convert to
sigma_E2I = sqrt(log(wE2I_std/(wE2I_mean^2)+1));

wI2E_mean= 5.0;wI2E_std = 2;
mu_I2E = log((wI2E_mean^2)/sqrt(wI2E_std+wI2E_mean^2));  %%%convert to
sigma_I2E = sqrt(log(wI2E_std/(wI2E_mean^2)+1));

wI2I_mean= 20;wI2I_std = 10;
mu_I2I = log((wI2I_mean^2)/sqrt(wI2I_std+wI2I_mean^2));  %%%convert to
sigma_I2I = sqrt(log(wI2I_std/(wI2I_mean^2)+1));

% weight=zeros(NCELL,NCELL);
% delays=zeros(NCELL,NCELL);

delay_min=0.5; %1
delay_max=5; %5  %%%define delay ranges

% generate random weights/delays in par with syn matrix
%%%%E2E
weight_PP=lognrnd(mu_E2E,sigma_E2E,[numel(find(active_syn_P2P>=0)),1]);
active_syn_P2P_weight=active_syn_P2P;
active_syn_P2P_weight(active_syn_P2P_weight>=0)=weight_PP;

delay=delay_min + (delay_max-delay_min)*rand(numel(find(active_syn_P2P>=0)),1);
active_syn_P2P_delay=active_syn_P2P;
active_syn_P2P_delay(active_syn_P2P_delay>=0)=delay;

%%%%I2I
weight_II=lognrnd(mu_I2I,sigma_I2I,[numel(find(active_syn_I2I>=0)),1]);
active_syn_I2I_weight=active_syn_I2I;
active_syn_I2I_weight(active_syn_I2I_weight>=0)=weight_II;

delay=delay_min + (delay_max-delay_min)*rand(numel(find(active_syn_I2I>=0)),1);
active_syn_I2I_delay=active_syn_I2I;
active_syn_I2I_delay(active_syn_I2I_delay>=0)=delay;

%%%%I2P
weight_IP=lognrnd(mu_I2E,sigma_I2E,[numel(find(active_syn_I2P>=0)),1]);
active_syn_I2P_weight=active_syn_I2P;
active_syn_I2P_weight(active_syn_I2P_weight>=0)=weight_IP;

delay=delay_min + (delay_max-delay_min)*rand(numel(find(active_syn_I2P>=0)),1);
active_syn_I2P_delay=active_syn_I2P;
active_syn_I2P_delay(active_syn_I2P_delay>=0)=delay;

%%%%P2I
weight_PI=lognrnd(mu_E2I,sigma_E2I,[numel(find(active_syn_P2I_reversed>=0)),1]);
active_syn_P2I_weight=active_syn_P2I_reversed;
active_syn_P2I_weight(active_syn_P2I_weight>=0)=weight_PI;

delay=delay_min + (delay_max-delay_min)*rand(numel(find(active_syn_P2I_reversed>=0)),1);
active_syn_P2I_delay=active_syn_P2I_reversed;
active_syn_P2I_delay(active_syn_P2I_delay>=0)=delay;


% for i=1:PCELL;   %%%source ID      %%%uni
%     for j=1:PCELL;   %%%target ID
%         if i~=j;
%             if syn_matrix(i,j)
%                 syn_matrix(i,j)=1;
%                 weight(i,j)=lognrnd(mu_E2E,sigma_E2E);
%                 %delays(i,j)=randi([delay_min,delay_max]);
%                 delays(i,j)=delay_min + (delay_max-delay_min)*rand(1,1);
%             end
%         end
%     end
% end
% 
% 
% for i=1:PCELL;   %%%source ID      %%%uni
%     for j=PCELL+1:NCELL;   %%%target ID
%         if i~=j;
%             if syn_matrix(i,j)
%                 syn_matrix(i,j)=2;
%                 weight(i,j)=lognrnd(mu_E2I,sigma_E2I);
%                 %delays(i,j)=randi([delay_min,delay_max]);
%                 delays(i,j)=delay_min + (delay_max-delay_min)*rand(1,1);
%             end
%         end
%     end
% end
% 
% 
% for i=PCELL+1:NCELL;   %%%source ID      %%%uni
%     for j=1:PCELL;   %%%target ID
%         if i~=j;
%             if syn_matrix(i,j)
%                 syn_matrix(i,j)=3;
%                 weight(i,j)=lognrnd(mu_I2E,sigma_I2E);
%                 %delays(i,j)=randi([delay_min,delay_max]);
%                 delays(i,j)=delay_min + (delay_max-delay_min)*rand(1,1);
%             end
%         end
%     end
% end
% 
% %values=[];n=1;
% for i=PCELL+1:NCELL;   %%%source ID      %%%uni
%     for j=PCELL+1:NCELL;   %%%target ID
%         if i~=j;
%             if syn_matrix(i,j)
%                 syn_matrix(i,j)=4;
%                 weight(i,j)=lognrnd(mu_I2I,sigma_I2I);
%                 %delays(i,j)=randi([delay_min,delay_max]);
%                 delays(i,j)=delay_min + (delay_max-delay_min)*rand(1,1);
%                 %values(n,1)=weight(i,j);
%                 %n=n+1;
%             end
%         end
%     end
% end



%%%%%%%oritation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=1000000;                       %%%method is from website http://mathworld.wolfram.com/SpherePointPicking.html%%%%
X=randn(N,1);
X2=X.^2;
Y=randn(N,1);
Y2=Y.^2;
Z=randn(N,1);
Z2=Z.^2;
square=1./(sqrt(X2+Y2+Z2));
[X_n]=square.*X;
[Y_n]=square.*Y;
[Z_n]=square.*Z;
randompo=[X_n,Y_n,Z_n];
nooritation = repmat([0 0 1],NCELL,1); 

P=1; %%%%portions of cells that are oritated. 0.3 is the default one

p_oritation=nooritation;
oritation_seq = randperm(N,NCELL*P);
oritation_seq=oritation_seq';
randompo=randompo(oritation_seq,:);
rand_cell_oritated=randperm(NCELL,NCELL*P);
rand_cell_oritated=rand_cell_oritated';
p_oritation(rand_cell_oritated,:)=randompo;


%%%%%%save various files%%%%%%%%%

fid=fopen('oritation.txt','w');

matrix=p_oritation;    %mm saved%*1000; %%%%%convert to um, for NEURON use.
[m,n]=size(matrix);
for i=1:m
    for j=1:n
        if j==n
            fprintf(fid,'%g\n',matrix(i,j));
        else
            fprintf(fid,'%g\t',matrix(i,j));
        end
    end
end
fclose(fid);


fid=fopen('location.txt','w');

matrix=location;    %mm saved%*1000; %%%%%convert to um, for NEURON use.
[m,n]=size(matrix);
for i=1:m
    for j=1:n
        if j==n
            fprintf(fid,'%g\n',matrix(i,j));
        else
            fprintf(fid,'%g\t',matrix(i,j));
        end
    end
end
fclose(fid);

%GG_type = ['location','.txt'];
%dlmwrite(GG_type,location,'delimiter','\t','precision', '%.7f');

% fid=fopen('Cell_type.txt','w');
% matrix=type;
% [m,n]=size(matrix);
% for i=1:m
%     for j=1:n
%         if j==n
%             fprintf(fid,'%g\n',matrix(i,j));
%         else
%             fprintf(fid,'%g\t',matrix(i,j));
%         end
%     end
% end
% fclose(fid);

Cell_list=[0:NCELL-1]';  %%%%generate cell_list for being used in NEURON
GG_list = ['Cell_list','.txt'];
dlmwrite(GG_list,Cell_list,'delimiter','\t','precision', '%d');
GG_type = ['Cell_type','.txt'];
dlmwrite(GG_type,type,'delimiter','\t','precision', '%d');

% fid=fopen('NM.txt','w');
% matrix=NM;
% [m,n]=size(matrix);
% for i=1:m
%     for j=1:n
%         if j==n
%             fprintf(fid,'%g\n',matrix(i,j));
%         else
%             fprintf(fid,'%g\t',matrix(i,j));
%         end
%     end
% end
% fclose(fid);

GG_NM = ['NM','.txt'];
dlmwrite(GG_NM,NM,'delimiter','\t','precision', '%d');

% fid=fopen('gj_matrix.txt','w');
% matrix=GAP_matrix;
% [m,n]=size(matrix);
% for i=1:m
%     for j=1:n
%         if j==n
%             fprintf(fid,'%g\n',matrix(i,j));
%         else
%             fprintf(fid,'%g\t',matrix(i,j));
%         end
%     end
% end
% fclose(fid);

% GG0 = ['gj_matrix','.txt'];
% dlmwrite(GG0,GAP_matrix,'delimiter','\t','precision', '%d');

%%%GAP
% GGG0 = ['active_syn_GAP','.txt'];
% dlmwrite(GGG0,active_syn_GAP,'delimiter','\t','precision', '%d');

%save('active_syn_GAP','active_syn_GAP','-ascii')

GAP_size=[size(active_syn_GAP,1),size(active_syn_GAP,2)];
% GGG1 = ['GAP_size','.txt'];
% dlmwrite(GGG1,GAP_size,'delimiter','\t','precision', '%d');
%save('GAP_size','GAP_size','-ascii')

%%%%EE
% aa = ['active_syn_P2P','.txt'];
% dlmwrite(aa,active_syn_P2P,'delimiter','\t','precision', '%d');
%save('active_syn_P2P','active_syn_P2P','-ascii')

% bb = ['active_syn_P2P_weight','.txt'];
% dlmwrite(bb,active_syn_P2P_weight,'delimiter','\t','precision', '%.1f');
% save('active_syn_P2P_weight','active_syn_P2P_weight','-ascii')

% cc = ['active_syn_P2P_delay','.txt'];
% dlmwrite(cc,active_syn_P2P_delay,'delimiter','\t','precision', '%.1f');
% save('active_syn_P2P_delay','active_syn_P2P_delay','-ascii')

%PP_size=[size(active_syn_P2P,1),size(active_syn_P2P,2)];
% aa1 = ['PP_size','.txt'];
% dlmwrite(aa1,PP_size,'delimiter','\t','precision', '%d');
%save('PP_size','PP_size','-ascii')

%%%%II
% dd = ['active_syn_I2I','.txt'];
% dlmwrite(dd,active_syn_I2I,'delimiter','\t','precision', '%d');
%save('active_syn_I2I','active_syn_I2I','-ascii')

% ee = ['active_syn_I2I_weight','.txt'];
% dlmwrite(ee,active_syn_I2I_weight,'delimiter','\t','precision', '%.1f');
% save('active_syn_I2I_weight','active_syn_I2I_weight','-ascii')

% ff = ['active_syn_I2I_delay','.txt'];
% dlmwrite(ff,active_syn_I2I_delay,'delimiter','\t','precision', '%.1f');
% save('active_syn_I2I_delay','active_syn_I2I_delay','-ascii')

%II_size=[size(active_syn_I2I,1),size(active_syn_I2I,2)];
% dd1 = ['II_size','.txt'];
% dlmwrite(dd1,II_size,'delimiter','\t','precision', '%d');
%save('II_size','II_size','-ascii')

%%%%IP
% gg = ['active_syn_I2P','.txt'];
% dlmwrite(gg,active_syn_I2P,'delimiter','\t','precision', '%d');
%save('active_syn_I2P','active_syn_I2P','-ascii')

% hh = ['active_syn_I2P_weight','.txt'];
% dlmwrite(hh,active_syn_I2P_weight,'delimiter','\t','precision', '%.1f');
% save('active_syn_I2P_weight','active_syn_I2P_weight','-ascii')

% ii = ['active_syn_I2P_delay','.txt'];
% dlmwrite(ii,active_syn_I2P_delay,'delimiter','\t','precision', '%.1f');
% save('active_syn_I2P_delay','active_syn_I2P_delay','-ascii')

 %IP_size=[size(active_syn_I2P,1),size(active_syn_I2P,2)];
% gg1 = ['IP_size','.txt'];
% dlmwrite(gg1,IP_size,'delimiter','\t','precision', '%d');
%save('IP_size','IP_size','-ascii')

%%%%PI
% jj = ['active_syn_P2I_reversed','.txt'];
% dlmwrite(jj,active_syn_P2I_reversed,'delimiter','\t','precision', '%d');
%save('active_syn_P2I_reversed','active_syn_P2I_reversed','-ascii')

% kk = ['active_syn_P2I_weight','.txt'];
% dlmwrite(kk,active_syn_P2I_weight,'delimiter','\t','precision', '%.1f');
% save('active_syn_P2I_weight','active_syn_P2I_weight','-ascii')

% ll = ['active_syn_P2I_delay','.txt'];
% dlmwrite(ll,active_syn_P2I_delay,'delimiter','\t','precision', '%.1f');
%save('active_syn_P2I_delay','active_syn_P2I_delay','-ascii')

%PI_size=[size(active_syn_P2I_reversed,1),size(active_syn_P2I_reversed,2)];
% jj1 = ['PI_size','.txt'];
% dlmwrite(jj1,PI_size,'delimiter','\t','precision', '%d');
%save('PI_size','PI_size','-ascii')
toc




%%%%%%%%for optimize%%%%%%%%
mm=1;
fid=fopen('active_syn_GAP_op','w');
active_syn_GAP1=active_syn_GAP1-1;    %%%just use half of GAP matrix, because it is symmatric
for i=1:size(active_syn_GAP1,1)
    temp=[];temp=active_syn_GAP1(i,:);
    temp(temp<0)=[];
    if length(temp)>0
fmt=[repmat('%d\t',1,length(temp)) '\n'];   %%%%format
fprintf(fid,fmt,temp);
        pregapnum(i,1)=length(temp);
        preid=i+24300-1;
        for jj=1:length(temp);
        pairtest(mm,1:2)=[2*mm-1,2*mm];
        mm=mm+1;
        end
    else 
        pregapnum(i,1)=0;
    end
end



ind_num=[];ind_num=cumsum(pregapnum);ind_num=[0;ind_num];
fid_ind=fopen('active_GAP_ind','w'); %%%%%for saving index for segerate fid file for each neuron
fmt=[repmat('%d\t',1,2) '\n'];
for i=1:INTCELL;
    
        temp=[];temp=[ind_num(i),ind_num(i+1)-1];
 fprintf(fid_ind,fmt,temp);
end

GG_gapid = ['GAPid'];
dlmwrite(GG_gapid,ind_num,'delimiter','\t','precision', '%d');



fid=fopen('active_syn_op','w');  %%%%for saving all pre (both INH and EXC) cells for PNs&ITNs 
for i=1:PCELL;
    i
    %temp=[];temp=active_syn_P2P(i,:);
    %temp(temp<=0)=[];
%fmt=[repmat('%d\t',1,length(temp)) '\n'];   %%%%format
%fprintf(fid,fmt,temp);

clear row_P2P; clear col; [row_P2P,col]=find(active_syn_P2P==i-1);  %%%row represents pre ID
clear row_I2P; clear col; [row_I2P,col]=find(active_syn_I2P==i-1); row_I2P=row_I2P+PCELL;    %%%transfer to ITN ID


if length(row_P2P)>=1||length(row_I2P)>=1;
clear row_com; row_com=[row_P2P;row_I2P];
row_com=sort(row_com)-1; %%%covert to neuron ID
%row_com=[i-1;numel(row_com);row_com]; %first num is the post ID, this is to avoid cell who receives nothing
                       %second num is indicating how many pre cells are connected to this post cell, for being used in
                       %NEURON
clear fmt; fmt=[repmat('%d\t',1,length(row_com)) '\n'];   %%%%format
fprintf(fid,fmt,row_com);

prenum(i,1)=length(row_com);  %%%to store pre num for each post


else
prenum(i,1)=0;  %%%to store pre num for each post

    
end
%[row,col]=ind2sub(size(active_syn_P2I_reversed),find(active_syn_P2I_reversed>0));

end


%fid=fopen('active_syn_I_op','w');%%%%for saving all pre (both INH and EXC) cells for ITNs 

for i=1:INTCELL;
    i
    
%     temp=[];temp=active_syn_I2I(i,:);
%     temp(temp<0)=[];
clear row_I2I; clear col; [row_I2I,col]=find(active_syn_I2I==i+PCELL-1); row_I2I=row_I2I+PCELL;%%%row represents pre ID
clear row_P2I; row_P2I=active_syn_P2I_reversed(i,:);row_P2I(row_P2I<0)=[];row_P2I=row_P2I'+1;

% fmt=[repmat('%d\t',1,length(temp)) '\n'];   %%%%format
% fprintf(fid,fmt,temp);

if length(row_P2I)>=1||length(row_I2I)>=1;
%P_prenum(i,1)=length(row_P2P)+length(row_I2P);  %%%to store pre num for each post
clear row_com; row_com=[row_P2I;row_I2I];
row_com=sort(row_com)-1;   %%%covert to neuron ID
%row_com=[i+PCELL-1;numel(row_com);row_com]; %first num is the post ID, this is to avoid cell who receives nothing
                       %second num is indicating how many pre cells are connected to this post cell, for being used in
                       %NEURON
clear fmt; fmt=[repmat('%d\t',1,length(row_com)) '\n'];   %%%%format
fprintf(fid,fmt,row_com);
prenum(i+PCELL,1)=length(row_com);  %%%to store pre num for each post
else
  prenum(i+PCELL,1)=0;  %%%to store pre num for each post  
end
end
ind_num=cumsum(prenum);ind_num=[0;ind_num];
fid_ind=fopen('active_syn_ind','w'); %%%%%for saving index for segerate fid file for each neuron
fmt=[repmat('%d\t',1,3) '\n'];
for i=1:INTCELL+PCELL;
    i
        temp=[];temp=[i-1,ind_num(i),ind_num(i+1)-1];
 fprintf(fid_ind,fmt,temp);
end

% fid=fopen('active_syn_I2P_op','w');
% 
% for i=1:size(active_syn_I2P,1)
%     temp=[];temp=active_syn_I2P(i,:);
%     temp(temp<=0)=[];
% fmt=[repmat('%d\t',1,length(temp)) '\n'];   %%%%format
% fprintf(fid,fmt,temp);
% end

% fid=fopen('active_syn_P2I_reversed_op','w');
% 
% for i=1:size(active_syn_P2I_reversed,1)
%     temp=[];temp=active_syn_P2I_reversed(i,:);
%     temp(temp<=0)=[];
% fmt=[repmat('%d\t',1,length(temp)) '\n'];   %%%%format
% fprintf(fid,fmt,temp);
% end
