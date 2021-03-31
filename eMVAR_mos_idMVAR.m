%% Model Order Selection for identification strictly causal MVAR model

%%% input:
% Y, M*N matrix of time series (each time series is in a row)
% pmax, maximum tested model order
% idMode, determines estimation algorithm (0:builtin least squares, else other methods [see mvar.m from biosig package])

%%% output:
% pottaic: model order optimized with multichannel Akaike Information Criterion (AIC)
% pottaic: model order optimized with  Minimum Description Length (MDL) criterion
% aic: values of AIC index as a function of the model order
% mdl: values of MDL index as a function of the model order

function [pottaic,pottmdl,aic,mdl] = eMVAR_mos_idMVAR(Y,pmax,idMode)

N=size(Y,2);
M=size(Y,1); %dimensionalità serie

% figures of merit
aic=NaN*ones(pmax,1); mdl=aic;

for p=1:pmax
    
    [Am,S,Yp]=eMVAR_idMVAR(Y,p,idMode);

    %formula multivariate AIC 
    aic(p)=N*log(det(S))+2*M*M*p; % S matrice di covarianza
    
    %formula multivariate MDL
    mdl(p)=N*log(det(S))+log(N)*M*M*p; % S matrice di covarianza
    
        
end

pottaic=find(aic == min(aic));
pottmdl=find(mdl == min(mdl));