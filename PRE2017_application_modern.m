%% Multiscale GC analysis of Temperature and CO2 data - modern climate data
% application of of Faes et al. PRE 2017
clc; clear; close all;

filtra='L1'; %'AR', 'diff', 'detrend', 'L1'
pfilter=0.92; % if 'AR'

% multiscale parameter
deltaT=12;
tauv=(deltaT:deltaT:1200)'; %vector of scales (in months)
whichfilter='F'; % 'A' for averaging, 'F' for FIR
ncoeff=6; % if FIR, set the number of coeff

% parameters for the analysis
pmax=20; %This is pmax akaike etc
p_crit='b'; % 'a'=Akaike, 'b'=BIC, 'c' p=pimp
pimp=10;

% IAAFT surrogates
numsurro=5; %no. of surrogates, 100 in the paper
lo=5; hi=95; %percentiles of surro distribution

% plottaGC='y'; % mettere 'y' se voglio calcolare GC, 'n' se voglio TE=GC/2

%% read, filter
load data_modern; data=[temp co2];
switch filtra
    case 'diff'
        data_f=diff(data);
    case 'AR'
        for m=1:size(data,2)
            [fia,fib]=AR_filter(data,m,pfilter);
            data_f(:,m)=fia;
            tendenza(:,m)=fib;
        end
    case 'detrend'
        tmp1=detrend(data(:,1));
        tmp2=detrend(data(:,2));
        data_f=[tmp1 tmp2];
    case 'L1'
        datatmp=data;
        load data_modern_L1
        data_f=data;
        data=datatmp;
    otherwise
        data_f=data;
end

temp=data(:,1); co2=data(:,2);
temp_f=data_f(:,1); co2_f=data_f(:,2);
temp_n=(temp-mean(temp))./std(temp);
co2_n=(co2-mean(co2))./std(co2);

temp_fn=(temp_f-mean(temp_f))./std(temp_f);
co2_fn=(co2_f-mean(co2_f))./std(co2_f);

Y=[temp_fn co2_fn]';
[M,N]=size(Y);

figure;
subplot(2,1,1);
plot(data(:,1),'k'); title('temp');zoom xon; xlim([1 N]);xlabel('months'); 
if strcmp(filtra,'AR'), hold on; plot(tendenza(:,1),'g'); end
subplot(2,1,2);
plot(data(:,2),'r'); title('co2');zoom xon; xlim([1 N]);xlabel('months'); 
if strcmp(filtra,'AR'), hold on; plot(tendenza(:,2),'g'); end

mesi=[1:1:708]';
outpergraf=[mesi temp_fn co2_fn];

figure;
subplot(2,1,1);
plot(temp_fn,'k','linewidth',1.5); title(['temp (filter ' filtra ')']); zoom xon; xlim([1 N]);xlabel('months')
subplot(2,1,2);
plot(co2_fn,'r','linewidth',1.5); title(['CO_2 (filter ' filtra ')']); zoom xon;xlim([1 N]);xlabel('months')

out.data=data;
out.datafilt=[temp_fn co2_fn];


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
        Ys(m,:,ns)=(surriaafft(Y(m,:)'))'; % function for IAAFT surrogates
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
    GC_s_d12(:,s)=prctile(squeeze(GC_s(1,2,s,:)),[lo 50 hi])';
    GC_s_d21(:,s)=prctile(squeeze(GC_s(2,1,s,:)),[lo 50 hi])';
end


%% plots
GC12=squeeze(GC(1,2,:));
GC21=squeeze(GC(2,1,:));

scale=tauv/12;

ymax=max([GC12; GC21]);
figure;
subplot(1,2,1);
plot(scale,GC12,'k.-');
hold on; plot(scale,GC_s_d12,':','color',[0.5 0.5 0.5 ],'linewidth',2);
title('GC_{temp \leftarrow CO2}');
legend('original', 'IAAFT surrogates');
ylim([0 1.1*ymax])

subplot(1,2,2);
plot(scale,GC21,'r.-');
hold on; plot(scale,GC_s_d21,':','color',[1 0.75 0.25],'linewidth',2);
title('GC_{CO2 \leftarrow Temp}');
legend('original', 'IAAFT surrogates');
ylim([0 1.1*ymax])

% legend('GC_{temp \leftarrow CO2}', 'GC_{CO2 \leftarrow Temp}');

disp(['model order p=' int2str(p)]);

out.GC=[tauv scale GC12 GC21];
out.surroGC=GC_s;
out.surroGC12_perctiles=GC_s_d12;
out.surroGC21_perctiles=GC_s_d21;




