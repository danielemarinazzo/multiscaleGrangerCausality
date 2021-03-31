%% Multiscale GC analysis -test for linear VAR processes
% analysis at varying scale tau - simulation A of Faes et al. PRE 2017
clear; close all; clc;

%% Parameters
tauv=(1:12)'; 
whichfilter='F'; % 'A' for averaging, 'F' for FIR
ncoeff=6; % if FIR, set the number of coeff

whichsimu=1; % 1 for unidirectional coupling (Fig.2), 2 for bidirectional coupling (Fig.2)
numsimu=5; %number of realizations
N=500; % length of realizations
pcrit='bic'; % 'aic', 'bic', or number for fixed model order
pmax=12; % max model order scanned by AIC/BIC

lo=25; hi=75; %percentiles (for visualization)

nscales=length(tauv);

%% Simulation parameters
M=2;
switch whichsimu
    case 1 %%%%% 1) Unidirectional
        a1=0.5; b1=1; %auto-coupling and lag for y1
        a2=0; b2=2; %auto-coupling and lag for y2
        c1=0; d1=3; %coupling and lag from y2 to y1
        c2=0.5; d2=2; %coupling and lag from y1 to y2

    case 2%%%%% 2) Other bidirectional
        a1=0.25; b1=1; %auto-coupling and lag for y1
        a2=0.25; b2=1; %auto-coupling and lag for y2
        c1=0.75; d1=2; %coupling and lag from y2 to y1
        c2=0.5; d2=7; %coupling and lag from y1 to y2
end

Su=1*eye(M);
p=max([b1*(a1~=0) b2*(a2~=0) d1*(c1~=0) d2*(c2~=0)]);

% VAR parameters
Ak=zeros(M,M,p);
if a1~=0, Ak(1,1,b1)=a1; end
if a2~=0, Ak(2,2,b2)=a2; end
if c1~=0, Ak(1,2,d1)=c1; end
if c2~=0, Ak(2,1,d2)=c2; end
Am=[];
for kk=1:p
    Am=[Am Ak(:,:,kk)];
end
% stability check
E=eye(M*p);AA=[Am;E(1:end-M,:)];lambda=eig(AA);lambdamax=max(abs(lambda));
if lambdamax>=1
    error('The simulated VAR process is not stable');
end


%% MSTE Analysis
% Information Dynamics matrix: storage on the diagonal, transfer off-diagonal
GCflt=nan*ones(M,M,nscales); GCdws=GCflt; %theoretical values of GC after filtering and downsampling
eGCflt=nan*ones(M,M,nscales,numsimu); eGCdws=eGCflt; %new estimate - SS-based
e0GCflt=nan*ones(M,M,nscales,numsimu); e0GCdws=e0GCflt; %"traditional" estimate - VAR-based
for s=1:nscales
    disp(['scale ' int2str(s)]);
    tau=tauv(s);
    
    % theoretical GC values at scale tau
    [GC2,GC1,b] = msgc(Am,Su,tau,ncoeff,whichfilter);
    GCflt(:,:,s)=GC1;
    GCdws(:,:,s)=GC2;
    
    
    %%% ESTIMATES
    for is=1:numsimu
        % uses functions of eMVAR toolbox to generate realizations of the VAR process
        U=eMVAR_InstModelfilter(N,Su,'StrictlyCausal'); % gaussian innovations with covariance Su
        Y=eMVAR_MVARfilter(Am,U);
        
        %%%STANDARD VAR ESTIMATION
        %%% Filtering
        Ytilda=filter(b,1,Y')';
        %model order selection
        if pcrit(1)=='a' || pcrit(1)=='b' 
            [pottaic,pottmdl,aic,mdl] = eMVAR_mos_idMVAR(Ytilda,pmax,0); %model order selection from eMVAR toolbox 
            if pcrit(1)=='a', p=pottaic; else p=pottmdl; end
        else
            p=pcrit;
        end
        % model identification
        GC1a=egc_gcMVAR(Ytilda,p); %Granger causality estimation from eGC toolbox
        e0GCflt(:,:,s,is)=GC1a;
        %%% Downsampling
        Ybar=downsample(Ytilda',tau)';
        if pcrit(1)=='a' || pcrit(1)=='b' 
            [pottaic,pottmdl,aic,mdl] = eMVAR_mos_idMVAR(Ybar,pmax,0); %model order selection from eMVAR toolbox 
            if pcrit(1)=='a', p=pottaic; else p=pottmdl; end
        else
            p=pcrit;
        end
        GC1d=egc_gcMVAR(Ybar,p); %Granger causality estimation from eGC toolbox
        e0GCdws(:,:,s,is)=GC1d;
        
        
        %%%NOVEL STATE SPACE ESTIMATION
        %model order selection
        if pcrit(1)=='a' || pcrit(1)=='b' 
            [pottaic,pottmdl,aic,mdl] = eMVAR_mos_idMVAR(Y,pmax,0); %model order selection from eMVAR toolbox 
            if pcrit(1)=='a', p=pottaic; else p=pottmdl; end
        else
            p=pcrit;
        end
        % model identification
        [eAm,eSu,Yp,Up]=eMVAR_idMVAR(Y,p,0); %model identification from eMVAR toolbox 
        % estimated GC values at scale tau
        [GC2,GC1,b] = msgc(eAm,eSu,tau,ncoeff,whichfilter);
        eGCflt(:,:,s,is)=GC1;
        eGCdws(:,:,s,is)=GC2;
  
    end % numsimu   

end % numscales

%% rearrange estimates - percentiles
eGCavg_d12=nan*ones(3,nscales); eGCavg_d21=eGCavg_d12;
eGCdws_d12=nan*ones(3,nscales); eGCdws_d21=eGCdws_d12;
eGC0avg_d12=nan*ones(3,nscales); eGC0avg_d21=eGC0avg_d12;
eGC0dws_d12=nan*ones(3,nscales); eGC0dws_d21=eGC0dws_d12;
for s=1:nscales
    eGCavg_d12(:,s)=prctile(squeeze(eGCflt(1,2,s,:)),[lo 50 hi])';
    eGCdws_d12(:,s)=prctile(squeeze(eGCdws(1,2,s,:)),[lo 50 hi])';
    eGCavg_d21(:,s)=prctile(squeeze(eGCflt(2,1,s,:)),[lo 50 hi])';
    eGCdws_d21(:,s)=prctile(squeeze(eGCdws(2,1,s,:)),[lo 50 hi])';
    
    eGC0avg_d12(:,s)=prctile(squeeze(e0GCflt(1,2,s,:)),[lo 50 hi])';
    eGC0dws_d12(:,s)=prctile(squeeze(e0GCdws(1,2,s,:)),[lo 50 hi])';
    eGC0avg_d21(:,s)=prctile(squeeze(e0GCflt(2,1,s,:)),[lo 50 hi])';
    eGC0dws_d21(:,s)=prctile(squeeze(e0GCdws(2,1,s,:)),[lo 50 hi])';
end


%% figs
switch whichfilter
case 'A'%%% AVERAGING FILTER
    labelfig='AVG';
case 'F' %%% FIR FILTER
    labelfig='FIR';
end
    
maxGCavg=1.5*max(GCflt(:));
maxGCdws=1.5*max(GCdws(:));
figure(1);
subplot(2,2,1); %AVG
plot(tauv,squeeze(GCflt(1,2,:)),'b.-'); hold on
plot(tauv,eGCavg_d12,'c.:');
plot(tauv,eGC0avg_d12,'y.:');
axis([1 max(tauv) 0 maxGCavg]);
title(['GC_{1\leftarrow2} , ' labelfig]); xlabel('scale');

subplot(2,2,3); %AVG
plot(tauv,squeeze(GCflt(2,1,:)),'b.-'); hold on
plot(tauv,eGCavg_d21,'c.:');
plot(tauv,eGC0avg_d21,'y.:');
axis([1 max(tauv) 0 maxGCavg]);
title(['GC_{2\leftarrow1} , ' labelfig]); xlabel('scale');

subplot(2,2,2); %DWS
plot(tauv,squeeze(GCdws(1,2,:)),'b.-'); hold on
plot(tauv,eGCdws_d12,'c.:');
plot(tauv,eGC0dws_d12,'y.:');
axis([1 max(tauv) 0 maxGCdws]);
title(['GC_{1\leftarrow2} , DWS']); xlabel('scale');

subplot(2,2,4); %DWS
plot(tauv,squeeze(GCdws(2,1,:)),'b.-'); hold on
plot(tauv,eGCdws_d21,'c.:');
plot(tauv,eGC0dws_d21,'y.:');
axis([1 max(tauv) 0 maxGCdws]);
title(['GC_{2\leftarrow1} , DWS']); xlabel('scale');

