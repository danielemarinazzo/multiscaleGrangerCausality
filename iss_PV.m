% Compute partial variances for a state space model from innovations form parameters 
% assumes one target and one source, and iterates
% A,C,K,V - innovations form state space parameters

% adapted from State-Space Granger Causality Matlab® Toolbox
% based on "Barnett, L. and Seth, A. K.,Granger causality for state-space models, 91(4),2015"

function ret = iss_PV(A,C,K,V,mostra)

if nargin<5, mostra='n'; end %display counters
M=size(V,1);

% determine the variance of the processes lambda0=E[Yn Yn'] for each target Y
O=dlyap(A,K*V*K');
lambda0=C*O*C'+V;
Sj=diag(lambda0);

% partial variance of target given past of all processes
Sj_ijk=diag(V); % j|all
Sj_j=nan*ones(M,1); % j|j
Sj_ij=nan*ones(M,M); % j|i,j
Sj_jk=nan*ones(M,M); % j|all but i
KVSQRT = K*chol(V,'lower');
Vrep=[]; %vector of reports
for jj=1:M  % set target
    % partial variance of target given its past only
    [~,tmp,rep] = iss_ss2iss(A,C(jj,:),KVSQRT*KVSQRT',V(jj,jj),K*V(:,jj));
    if isempty(tmp), tmp=NaN; end % per evitare che si incastri se la stima è "fallita"
    Sj_j(jj)=tmp; Vrep=[Vrep; rep];
    
    
    for ii=1:M
        if ii~=jj
            if ~strcmp(mostra,'n')
                clc; disp([mostra ', j=' int2str(jj) ', i=' int2str(ii)]);
            end
            % partial variance of target given its past and past of the source
            jjii=[jj ii];
            [~,tmp,rep] = iss_ss2iss(A,C(jjii,:),KVSQRT*KVSQRT',V(jjii,jjii),K*V(:,jjii));
            if isempty(tmp), tmp=NaN; end % per evitare che si incastri se la stima è "fallita"
            Sj_ij(jj,ii)=tmp(1,1); Vrep=[Vrep; rep];
            
            % partial variance of target given past of all processes except the source
            kk=1:M; kk([ii jj]) = []; % indices of remaining processes
            jjkk=[jj kk];
            [~,tmp,rep] = iss_ss2iss(A,C(jjkk,:),KVSQRT*KVSQRT',V(jjkk,jjkk),K*V(:,jjkk));
            Sj_jk(jj,ii)=tmp(1,1); Vrep=[Vrep; rep];
        end
    end
    
    for cnt=1:length(Vrep)
    rep=Vrep(cnt);
    if rep < 0
        if     rep == -1, warning('DARE: eigenvalues on/near unit circle');
        elseif rep == -2, warning('DARE: couldn''t find stablising solution');
        end
        ret.Sigmaj=nan*ones(M,1);
        ret.Sigmaj_j=nan*ones(M,1);
        ret.Sigmaj_ij=nan*ones(M,M);
        ret.Sigmaj_jk=nan*ones(M,M);
        ret.Sigmaj_ijk=nan*ones(M,1);
        return
    end
    if rep > sqrt(eps), warning('DARE: there were accuracy issues (relative residual = %e)',rep);
        ret.Sigmaj=nan*ones(M,1);
        ret.Sigmaj_j=nan*ones(M,1);
        ret.Sigmaj_ij=nan*ones(M,M);
        ret.Sigmaj_jk=nan*ones(M,M);
        ret.Sigmaj_ijk=nan*ones(M,1);
        return
    end
    end

end

ret.Sigmaj=Sj;
ret.Sigmaj_j=Sj_j;
ret.Sigmaj_ij=Sj_ij;
ret.Sigmaj_jk=Sj_jk;
ret.Sigmaj_ijk=Sj_ijk;

end

