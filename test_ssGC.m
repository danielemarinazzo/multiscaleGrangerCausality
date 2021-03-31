%% computation of multiscale GC with state-space models for simple simulated VAR processes
clear; close all; clc;

%% Parameters
tau=1; %scale (NOTE: for tau=1 computes time-domain state-space GC, (ncoeff is irrelevant))
ncoeff=6; % number of coeff of FIR lowpass filter

N=300; % length of realization

%%% Simulation parameters
M=3;
p=2;
Ak=zeros(M,M,p);
Ak(1,1,1)=0.5;
Ak(2,1,2)=0.5;
Ak(3,1,1)=0.7;
Su=1*eye(M);
Am=[];
for kk=1:p
    Am=[Am Ak(:,:,kk)];
end
% stability check
E=eye(M*p);AA=[Am;E(1:end-M,:)];lambda=eig(AA);lambdamax=max(abs(lambda));
if lambdamax>=1
    error('The simulated VAR process is not stable');
end


%% GC Analysis
% theoretical GC values
GC = msgc(Am,Su,tau,ncoeff);

% estimated GC values
U=eMVAR_InstModelfilter(N,Su,'StrictlyCausal'); % gaussian innovations with covariance Su
Y=eMVAR_MVARfilter(Am,U); %simulated process

[eAm,eSu]=eMVAR_idMVAR(Y,p,0); %model identification from eMVAR toolbox 
eGCss = msgc(eAm,eSu,tau,ncoeff); %Granger causality estimated using ISS models

if tau==1
    
end

%% display
if tau>1
    clc;
    disp(['scale ' int2str(tau)])
    for ii=1:M
        for jj=ii+1:M
                disp(['GC ' int2str(ii) '->' int2str(jj) ', theor:' num2str(GC(jj,ii)) ', est.SS:' num2str(eGCss(jj,ii))]);
                disp(['GC ' int2str(jj) '->' int2str(ii) ', theor:' num2str(GC(ii,jj)) ', est.SS:' num2str(eGCss(ii,jj))]);
        end
    end
end

if tau==1
    eGCvar=egc_gcMVAR(Y,p); %Granger causality estimated from eGC toolbox (at scale 1)
    
    clc;
    for ii=1:M
        for jj=ii+1:M
                disp(['GC ' int2str(ii) '->' int2str(jj) ', theor:' num2str(GC(jj,ii)) ', est.SS:' num2str(eGCss(jj,ii)) ', est.VAR:' num2str(eGCvar(jj,ii))]);
                disp(['GC ' int2str(jj) '->' int2str(ii) ', theor:' num2str(GC(ii,jj)) ', est.SS:' num2str(eGCss(ii,jj)) ', est.VAR:' num2str(eGCvar(ii,jj))]);
        end
    end
end



   