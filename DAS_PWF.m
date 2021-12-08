clear all; close all
addpath /home/gilberto/MatlabCodes/
% addpath /home/gilberto/MatlabCodes/PSpicker/

% Simple code for performing Plane Wave Fitting of array data
% using a sliding window. Code used for analyzing the DAS recording
% of the Mw 4.3 Hawtorne earthquake discussed in Piana Agostinetti et al. 
% (2021); the final output is Figure 12 in that paper. 
%
% Data are also provided: they are sac files cointaining stacked
% DAS channels; Station coordinates are in Stla, Stlo header fields.
%
% Additional functions (posted separately):
% Readsac -> to read sac data files
% wgs2utm -> to convert (lat,lon) to UTM coordinates
% G. Saccorotti, INGV. 8 December 2021
%

set(0,'defaultaxeslinewidth',1.)
set(0,'defaultaxesfontsize',10)
set(0,'defaultlinelinewidth',1.)

%%% Directory where to seek for seismograms
% 
datadir='../DATA/'

%%% Processing parameters: 

f1=0.5; f2=2;   % corner frequency of filter (not used: 
                % seismograms are already filtered); 
wls=4;          % length of sliding window (s);
ct=0.85         % Correlation threshold
iref=59;        % reference channel (just for plotting)
wshift=8        % window shift is wls/wshift


D=dir([datadir '*.sac']);
nsta=length(D);
for k=1:nsta
    A=Readsac([D(k).folder '/' D(k).name]);
    stla(k)=A.STLA;
    stlo(k)=A.STLO;
    S(:,k)=A.DATA1;
    dt=A.DELTA;
end

[np,~]=size(S); T=0:dt:(np-1)*dt;
[xj,yj]=wgs2utm(stla,stlo);
xp=xj./1000; yp=yj./1000;


wl=round(wls/dt); 
sig=2*dt;   % This is the standard deviation of delay 
            % time estimates, assumed equal to twice 
            % the sampling interval
          
% Begins PWF estimates over sliding windows

ipwf=0;H=[];JH=[];

% Clean up some variables
k=0; azpwf=[];raypwf=[];ccmax=[];TPWF=[]; TPWFG=[];


for i1=1:wl/8:np-wl
    i2=i1+wl-1;
    k=k+1;
    
disp(['PWF processing time window: ',num2str(k)])
ip=0; W=[]; G=[]; ig=0; delay=[];
    
    for ii1=1:nsta-1
    y1=S(i1:i2,ii1); y1=y1.*tukeywin(numel(y1),0.05);
        for ii2=ii1+1:nsta
        ip=ip+1;
        y2=S(i1:i2,ii2); y1=y1.*tukeywin(numel(y2),0.05);
        [xcj,lags]=xcorr(y1,y2,'coef');          
        [xcmax(ip),im]=max(xcj);

        % Only keep delayu time estimates associated with 
        % max(CC) > correlation threshold
        
        if xcmax(ip)>ct
                ig=ig+1;
                delay(ig)=lags(im)*dt;
                dx=xp(ii2)-xp(ii1);
                dy=yp(ii2)-yp(ii1);
                % Kernel matrix (differential copordinates)
                G(ig,:)=[dx dy];
                % Weigth matrix
                W(ig,ig)=xcmax(ip)./(1-xcmax(ip));
            end
        end
    end
    
    JH=histogram(xcmax,[0:0.01:1]);
    H(:,k)=JH.Values;    
    [nrg(k),ncg]=size(G);
    
    % Time vector
    TPWFG(k)=mean([T(i1),T(i2)]);
    
    % Fisher's z-transform for obtaining mean correlation coeff.
    z=0.5*log((1+xcmax)./(1-xcmax)); 
%     ccmaxg(k)=median(xcmax);  
    ccmaxg(k)=tanh(mean(z));  

    % if at least three delay times estimates are present, proceed with
    % slownmess measurement
if nrg(k)>2
    ipwf=ipwf+1;
    GW=W*G; dw=W*delay';
    mpwf=inv(GW'*GW)*GW'*dw;

    % Fisher's z-transform for obtaining mean correlation coeff.
    z=0.5*log((1+xcmax)./(1-xcmax)); 
    ccmax(ipwf)=tanh(mean(z));  
  
%     ccmax(ipwf)=median(xcmax);      

    res(ipwf)=(dw-GW*mpwf)'*(dw-GW*mpwf);       % residuals
    covar=(sig.^2)*inv(G'*G);                   % covar. matrix
    
    slo(ipwf,1:2)=mpwf;
    % 67% and 95% error bounds on slowness
    conf95(ipwf,1:2)=1.96*sqrt(diag(covar));    
    conf67(ipwf,1:2)=sqrt(diag(covar));
    % Propagation azimuth
    azpwf(ipwf)=rad2deg(atan2(mpwf(1),mpwf(2)));
    % Ray Parameter
    raypwf(ipwf)=sqrt(mpwf(1)^2+mpwf(2)^2);
    % Time vector
    TPWF(ipwf)=mean([T(i1),T(i2)]);
    NPWF(ipwf)=nrg(k);
    
% Upper and lower error bounds on azimuth
azpwfe1(ipwf)=rad2deg(atan2(mpwf(1)+conf95(ipwf,1),mpwf(2)-conf95(ipwf,2)));
azpwfe2(ipwf)=rad2deg(atan2(mpwf(1)-conf95(ipwf,1),mpwf(2)+conf95(ipwf,2)));

% Upper and lower error bounds on ray par
raypwfe1(ipwf)=sqrt((mpwf(1)+conf95(ipwf,1))^2+(mpwf(2)+conf95(ipwf,2))^2);
raypwfe2(ipwf)=sqrt((mpwf(1)-conf95(ipwf,1))^2+(mpwf(2)-conf95(ipwf,2))^2);
        
end
end

% 
save([datadir 'PWF.data.mat'],'dt','np','S','ccmax','TPWFG','ccmaxg',...
'TPWF','slo','conf95','conf67');
 
 
% 
% %%%%%%% Plotting results
% 

azpwf(azpwf<0)=azpwf(azpwf<0)+360;
azpwfe1(azpwfe1<0)=azpwfe1(azpwfe1<0)+360;
azpwfe2(azpwfe2<0)=azpwfe2(azpwfe2<0)+360;


lat0=38.479; lon0=-118.366;     %epicentral coordinates
T0=datenum(2016,3,21,7,37,10);

% theoretical propagation azimuth
azteo=azimuth(lat0,lon0,mean(stla),mean(stlo));


figure(6)
colormap(gray)
subplot(411)
tax=[0:dt:(np-1).*dt];
plot(tax,S(:,iref),'-k')
ylim([-11 11])
ylabel('n strain s^{-1}')
set(gca,'XTickLabel',[])

scale=10.*ccmax./max(ccmax);

subplot(412)
plot(TPWFG,ccmaxg,'ok','MarkerFaceColor',[.5 .5 .5])
ylabel('CC_{med}'); 
ylim([0.2 1]); xlim([0 tax(end)])
grid on
set(gca,'XTickLabel',[])

subplot(413)
scatter(TPWF,azpwf,(ccmax.*scale).^2,[.5 .5 .5],'o','filled'); hold on
for k=1:numel(TPWF)
    plot([TPWF(k) TPWF(k)],[azpwfe1(k) azpwfe2(k)],'-k'); 
end
plot(tax,repmat(azteo,length(tax),1),'--k')
box on
ylim([-10,370]); xlim([0 tax(end)])
ylabel('Azim.(°)');set(gca,'XTickLabel',[]); 

subplot(414)
scatter(TPWF,1./raypwf,(ccmax.*scale).^2,[.5 .5 .5],'o','filled'); hold on
for k=1:numel(TPWF)
    plot([TPWF(k) TPWF(k)],[1./raypwfe1(k) 1./raypwfe2(k)],'-k'); 
end

box on; hold on
plot([0 tax(end)],[4 4],'--k')
plot([0 tax(end)],[6 6],'--k')
ylim([0,10]); xlim([0 tax(end)])
ylabel('V_{app} (km s^{-1})')
xlabel('TIME (s)')
set(gcf,'Position',[740  626 560 800])        


% Uncomment for saving the figure
% saveas(gcf,'../FIGURES/fig12.eps','epsc')
















