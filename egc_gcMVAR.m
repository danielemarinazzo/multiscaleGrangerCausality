%% GRANGER CAUSALITY FROM STRICTLY CAUSAL MVAR MODEL: Y(n)=A(1)Y(n-1)+...+A(p)Y(n-p)+U(n)
% estimates Granger Causality in multiple time series from MVAR model fitted on data
% performs also row-by-row MVAR identification, equivalent to idMVAR (with vector least squares)

%%% input:
% Y, M*N matrix of time series (each time series is in a row)
% p, model order

%%% output:
% GC - M*M Granger Causality matrix - ij element is the GC from j to i
% p_val, matrix of p values of F test associated to GC causality
% Am=[A(1)...A(p)], M*pM matrix of the estimated MVAR model coefficients
% Su, estimated M*M input covariance matrix
% SigmaR, M*M matrix of restricted regression residual variances
% Ures, residuals of unrestricted model, dimension M*N 


function [GC,p_val,Am,Su,SigmaR,Ures]=egc_gcMVAR(Y,p)

M=size(Y,1);
tau=ones(1,M);
u=ones(1,M);
zerolag=zeros(1,M);
p_vett=p*ones(1,M);
[Vu]=egc_SetLag(p_vett,tau,u,zerolag);
Nu=size(Vu,1); %number of coeff for unrestricted regression

% Unrestricted regressions
Am=NaN*ones(M,p*M);
SigmaU=NaN*ones(M,M);
for jj=1:M
    [Sigma,Upu,coeff]=egc_LinReg(Y',jj,Vu); % execute linear regression 
    SigmaU(jj,:)=Sigma;
    Upred(:,jj)=Upu; % residuals of the unrestricted regression
    for k=1:size(Vu,1)
        Am(jj,M*(Vu(k,2)-1)+Vu(k,1))=coeff(k);
    end
end
Su=cov(Upred);
Ures=Upred';

% restricted regressions
SigmaR=NaN*ones(M,M); p_val=NaN*ones(M,M);
for jj=1:M
    for ii=1:M
        p_r=p_vett; p_r(ii)=0;
        Vr=egc_SetLag(p_r,tau,u,zerolag);
        [Sigma,Upr,~]=egc_LinReg(Y',jj,Vr);
        SigmaR(jj,ii)=Sigma;
        % Ftest for significance
        Nr=size(Vr,1);%number of coeff for restricted regression
        p_val(jj,ii) = egc_LinReg_Ftest(Upred(:,jj),Upr',Nu,Nr);
    end
end
GC=log(SigmaR./SigmaU);


end

