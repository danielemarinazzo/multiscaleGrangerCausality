%% Multiscale GC analysis -test for multiscale VAR processes
% analysis at varying scale tau - simulation B of Faes et al. PRE 2017
% (numerical estimation of MSGC)
clear; close all; clc;

%% Parameters
tauv=(1:15)';
ncoeff=6; %  set the number of coeff of FIR filter
% ftauv=1./(2*tauv);

numsimu=3; %number of realizations
N=1000; %simulation length
pcrit='bic'; % 'aic', 'bic', or number for fixed model order
pmax=12; % max model order scanned by AIC/BIC

lo=25; hi=75; %percentiles (for visualization)

nscales=length(tauv);

% Simulation parameters
M=2; %generate two bivariate unidirectionally coupled AR processes,then merge
plag=2;
r1=0.95; f1=0; Su1=0.25; %Su1=0.1;
Su2=0.5;
r3=0.95; f3=0.075; Su3=1;
Su4=0.5;
c1=0.5; %c1=1; %x1->x2 at lag 1
c2=1; %x3->x4 at lag 1

%% theoretical coeffs
%%% theoretical coeffs 1->2
Ak=zeros(M,M,plag);
Ak(1,:,1)=[2*r1*cos(2*pi*f1) 0];
Ak(1,:,2)=[-r1^2 0];
Ak(2,:,1)=[c1 0];
Ak(2,:,2)=[0 0];
Am12=[];
for kk=1:plag
    Am12=[Am12 Ak(:,:,kk)];
end
% stability check
E=eye(M*plag);AA=[Am12;E(1:end-M,:)];lambda=eig(AA);lambdamax=max(abs(lambda));
if lambdamax>=1
    error('The simulated VAR process is not stable');
end

%%% theoretical coeffs 3->4
Ak=zeros(M,M,plag);
Ak(1,:,1)=[2*r3*cos(2*pi*f3) 0];
Ak(1,:,2)=[-r3^2 0];
Ak(2,:,1)=[c2 0];
Ak(2,:,2)=[0 0];
Am34=[];
for kk=1:plag
Am34=[Am34 Ak(:,:,kk)];
end
% stability check
E=eye(M*plag);AA=[Am34;E(1:end-M,:)];lambda=eig(AA);lambdamax=max(abs(lambda));
if lambdamax>=1
error('The simulated VAR process is not stable');
end


%% MSGC Analysis
eGCss=nan*ones(M,M,nscales,numsimu); %new estimate -SS-based
% eGCvar=nan*ones(M,M,nscales,numsimu); %"traditional" estimate - VAR-based
for is=1:numsimu
    clc; disp(['simu ' int2str(is) ' of ' int2str(numsimu)]);
    
    %%%%% Generate simulated data 1->2
    U=[sqrt(Su1)*randn(1,N); sqrt(Su2)*randn(1,N)];
    X12=eMVAR_MVARfilter(Am12,U);

    %%%%% Generate simulated data 3->4
    U=[sqrt(Su3)*randn(1,N); sqrt(Su4)*randn(1,N)];
    X34=eMVAR_MVARfilter(Am34,U);
    
    %%% instantaneous sum and identification
    Y=[X12(1,:)+X34(2,:); X12(2,:)+X34(1,:)];
    
    for s=1:nscales
        disp(['scale ' int2str(s)]);
        tau=tauv(s);
        
        %%% STATE SPACE ESTIMATION
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
        [GC2,GC1,b] = msgc(eAm,eSu,tau,ncoeff);
        eGCss(:,:,s,is)=GC2;        
    end

    
end


%% FIGURE
% rearrange estimates - percentiles
eGCss_d12=nan*ones(3,nscales); eGCss_d21=eGCss_d12;
% eGCvar_d12=nan*ones(3,nscales); eGCvar_d21=eGCvar_d12;
for s=1:nscales
    eGCss_d12(:,s)=prctile(squeeze(eGCss(1,2,s,:)),[lo 50 hi])';
    eGCss_d21(:,s)=prctile(squeeze(eGCss(2,1,s,:)),[lo 50 hi])';
end

figure(1);clf;
%one realization
subplot(3,4,[1 2 3]); plot(Y(1,:),'k'); ylabel('y_1')
subplot(3,4,[5 6 7]); plot(Y(2,:),'k'); ylabel('y_2')
subplot(3,4,[4 8]); 
plot(tauv,squeeze(eGCss(1,2,:,is)),'r.-'); hold on
plot(tauv,squeeze(eGCss(2,1,:,is)),'k.:');
xlim([1 max(tauv)]);
legend('GC_{1\leftarrow2}','GC_{2\leftarrow1}')
%distributions over numsimu realizations
subplot(3,4,[9 10]);
plot(tauv,eGCss_d12,'r.:'); hold on
xlim([1 max(tauv)]);
ylim([0 0.9]);
title('GC_{1\leftarrow2}'); xlabel('scale');
subplot(3,4,[11 12]);
plot(tauv,eGCss_d21,'k.:');
xlim([1 max(tauv)]);
ylim([0 0.9]);
title('GC_{2\leftarrow1}'); xlabel('scale');

 



