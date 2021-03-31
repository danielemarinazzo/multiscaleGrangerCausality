%% LINEAR REGRESSION

%%% INPUTS:
% data: N*M matrix of the M signals each having length N
% j: index (column) of the series considered as output, the one we want to describe
% V: two column vector of series (col 1) and lag (col 2) indexes

function [S,Up,Am]=egc_LinReg(data,j,V)

if isempty(V) %if no conditioning, ce will be the Entropy of B
    S=var(data(:,j));
    Up=data(:,j)-mean(data(:,j)); Am=[];
else % compute Conditional entropy
    B=egc_buildvectors(data,j,V); %% form the observation matrix
    
    % Linear Regression
    Yb=B(:,1)'; % inversion works with data organized in rows
    A=B; A(:,1)=[]; Z=A';
    Am=Yb/Z; % least squares!

    Yp=Am*Z; 
    Up=Yb-Yp;
    S=cov(Up');
    
end


