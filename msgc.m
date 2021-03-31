%% MULTISCALE GC computation

%%% inputs
% Am, Su: VAR parameters (Am: M x pM coeff matrix; Su: M x M innovation covariance matrix)
% tau: scale factor
% ncoeff: number of coefficients of FIR filter (if no averaging)
% whichfilter: 'F' for FIR (default), 'A' for averaging

%%% outputs
% GCdws: GC at scale tau - MxM matrix with GC i->j in position (j,i)
% GCflt: GC at scale tau after first step (filtering)
% b: filter coeffs

function [GCdws,GCflt,b] = msgc(Am,Su,tau,ncoeff,whichfilter)
    if nargin<5, whichfilter='F'; end;
    
    M=size(Am,1);
    GCflt=nan*ones(M,M);
    GCdws=nan*ones(M,M);
    
    % MA parameters resulting from the change of scale
    switch whichfilter
    case 'A'%%% AVERAGING FILTER
        epsi=0.001; ntau0=1+epsi*(tau-1); ntau=1-epsi; %have to put an epsilon to guarantee stability
        B0=ntau0/tau*eye(M);
        Bm=repmat(ntau/tau*eye(M),1,tau-1);
        b=(1/tau)*ones(1,tau); 
    case 'F' %%% FIR FILTER
        if tau==1
            q=0; b=1;
        else
            q=ncoeff; % number of filter coeffs
            ft=1/(2*tau); %cutoff frequency
            Wn=2*ft; %normalized cutoff frequency (fNyquist=1)
            b=fir1(q,Wn,'noscale'); %Hamming window, linear phase (symmetry of b coeffs)
        end
        Bk=zeros(M,M,q+1);
        for l=1:q+1
            Bk(:,:,l)=b(l)*eye(M);
        end
        Bm=[];
        for kk=1:q+1
            Bm=[Bm Bk(:,:,kk)];
        end
        B0=Bm(1:M,1:M);
        Bm=Bm(1:M,M+1:end);
    end
    
    % ISS parameters
    [A,C,K,V,Vy] = iss_varma2iss(Am,Bm,Su,B0); % max(abs(eig(A-K*C)))
    
    % FILTERING
    ret_flt = iss_PV(A,C,K,V);
    for jj=1:M
        for ii=1:M
            if ii~=jj
                GCflt(jj,ii)=log(ret_flt.Sigmaj_j(jj)/ret_flt.Sigmaj_ij(jj,ii));
            end
        end
    end
   
    %%% DOWNSAMPLING
    [Ad,Kd,Vd] = iss_ds(A,C,K,V,tau);
    Cd=C;
    
    ret_dws = iss_PV(Ad,Cd,Kd,Vd);
    for jj=1:M
        for ii=1:M
            if ii~=jj
                GCdws(jj,ii)=log(ret_dws.Sigmaj_jk(jj,ii)/ret_dws.Sigmaj_ijk(jj));
            end
        end
    end

end