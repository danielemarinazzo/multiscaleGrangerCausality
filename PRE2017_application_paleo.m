%% Multiscale GC analysis of Temperature and CO2 data - paleolithic climate data
% application of of Faes et al. PRE 2017
clc; clear; close all;

interpola='y';

% multiscale parameter
tauv=(1:1:55)';
whichfilter='F'; % 'A' for averaging, 'F' for FIR
ncoeff=6; % if FIR, set the number of coeff

% parameters for the analysis
pmax=20; %This is pmax akaike etc
p_crit='b'; % 'a'=Akaike, 'b'=BIC, 'c' p=pimp
pimp=5;

% IAAFT surrogates
numsurro=5; %no. of surrogates, 100 in the paper
lo=5; hi=95; %percentiles of surro distribution

%% read, filter
data0=load('data800.mat');
data1=data0.data;

times=2000-data0.age_inv; %times=times-min(times);
times(565)=[]; data1(565,:)=[]; % non-strict monotonicity of times!

figure;
plotyy(times,data1(:,1),times,data1(:,2))

deltaT=(max(times)-min(times))/(length(times)-1);
timesU=(min(times):deltaT:max(times))';

if interpola=='y'
    data(:,1)=interp1(times,data1(:,1),timesU,'spline');
    data(:,2)=interp1(times,data1(:,2),timesU,'spline');
    times=timesU;
else
    data=data1;
end

temp=data(:,1); co2=data(:,2);

temp_n=(temp-mean(temp))./std(temp);
co2_n=(co2-mean(co2))./std(co2);

Y=[temp_n co2_n]';
[M,N]=size(Y);

out.data1=data1; %original data
out.data=data; %interpolated data
out.Y=Y; %analyzed data (normalized)

figure;
subplot(2,1,1); plot(timesU,temp_n,'k','linewidth',1.5);title('tempo');hold on; 
xlim([min(timesU) max(timesU)]);xlabel('years'); 
subplot(2,1,2);plot(timesU,co2_n,'r','linewidth',1.5); title('CO2'); zoom xon;
xlim([min(timesU) max(timesU)]);xlabel('years'); 


%% model identification
nscales=length(tauv);

% model order selection
[p_aic,p_bic,aic,bic] = eMVAR_mos_idMVAR(Y,pmax,0); %model order selection from eMVAR toolbox 
switch p_crit
    case 'a'
        p=p_aic;
    case 'b'
        p=p_bic;
    case 'c'
        p=pimp;
end

[Am,Su]=eMVAR_idMVAR(Y,p,0); %model identification from eMVAR toolbox 
E=eye(M*p);AA=[Am;E(1:end-M,:)];lambda=eig(AA);lambdamaxo=max(abs(lambda));
if lambdamaxo>=1,
    warning('Non-stable VAR process');
end

%% surrogate data
Ys=nan*ones(size(Y,1),size(Y,2),numsurro);
for ns=1:numsurro
    for m=1:M
        Ys(m,:,ns)=(surriaafft(Y(m,:)'))';
    end
end

%% MULTISCALE ANALYSIS
for s=1:nscales
    tau=tauv(s);
    clc; disp(['scale ' int2str(s) ' of ' int2str(nscales)]);
    
    % GC on original data at scale tau
    GC(:,:,s) = msgc(Am,Su,tau,ncoeff,whichfilter);  
    
    %%%% surrogate data - multiscale analysis
    for ns=1:numsurro
        % model order selection and identification
        [p_aic,p_bic,aic,bic] = eMVAR_mos_idMVAR(Ys(:,:,ns),pmax,0);
        switch p_crit
            case 'a', ps=p_aic;
            case 'b', ps=p_bic;
            case 'c', ps=pimp;
        end
        [Ams,Sus]=eMVAR_idMVAR(Ys(:,:,ns),ps,0);
        
        GCtmp=msgc(Ams,Sus,tau,ncoeff,whichfilter);
        GC_s(:,:,s,ns)=GCtmp;
        
    end
    
    
end

%%%% surogate distributions
for s=1:nscales
    eGCdws_d12(:,s)=prctile(squeeze(GC_s(1,2,s,:)),[lo 50 hi])';
    eGCdws_d21(:,s)=prctile(squeeze(GC_s(2,1,s,:)),[lo 50 hi])';
end


%% plots
GC12=squeeze(GC(1,2,:));
GC21=squeeze(GC(2,1,:));
scale=deltaT*tauv;

disp(['model order p=' int2str(p)]);

ymax=max([GC12; GC21]);
figure(3);clf;
subplot(1,2,1);
plot(scale,GC12,'k.-');
hold on; plot(scale,eGCdws_d12,':','color',[0.5 0.5 0.5 ],'linewidth',2);
title(['data800.mat, interpola=' interpola ', GC_{temp \leftarrow CO2}']);
legend('original', 'IAAFT surrogates');
ylim([0 1.1*ymax])
xlim([min(scale) max(scale)]);

subplot(1,2,2);
plot(scale,GC21,'r.-');
hold on; plot(scale,eGCdws_d21,':','color',[1 0.75 0.25],'linewidth',2);
title(['data800.mat, interpola=' interpola ', GC_{CO2 \leftarrow Temp}']);
legend('original', 'IAAFT surrogates');
ylim([0 1.1*ymax])
xlim([min(scale) max(scale)]);

out.GC=[tauv scale GC12 GC21];
out.surroGC=GC_s;
out.surroGC12_perctiles=eGCdws_d12;
out.surroGC21_perctiles=eGCdws_d21;





