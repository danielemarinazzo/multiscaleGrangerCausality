% genera iterative amplitude adjusted fourier tranform surrogates 
% algoritmo di Schreiber e Schmitz - Physical Review Letters 1996

% y: serie da surrogare
% nit: numero di iterazioni volute (default 7)
% stop: se metto 'spe' esce con lo spettro conservato, se metto 'dis' esce con la distribuzione conservata

function ys=surriaafft(y,nit,stop)

error(nargchk(1,3,nargin));%min e max di input arguments
if nargin < 3, stop='spe'; end %default matcha lo spettro
if nargin < 2, nit=7; end %default 7 iterazioni

% clear;close all;
% percorso='D:\johnny\lavoro\integrate_nlpred\elaborati_loo_si\';% percorso dei dati da analizzare
% nomefile='b-ca.prn';
% rs=load([percorso nomefile]);
% y=rs(:,1);
% y=(y-mean(y))/std(y);


%%
[ysorted,yindice]=sort(y);
my=abs(fft(y));
% ys=surrfft(y); %inizializzazione
ys=surrshuf(y); %inizializzazione


%% ciclo

for i=1:nit
    % step 1: impone lo spettro
    faseys=angle(fft(ys));
    fys=my.*(cos(faseys)+j*sin(faseys));
    ys=ifft(fys);ys=real(ys);
    ys=ys-mean(ys);

    % step 2: impone la distribuzione
    [yssorted,ysindice]=sort(ys);
    ypermuted=zeros(length(y),1);
    for i=1:length(y)
        ypermuted(ysindice(i))=ysorted(i);
    end
    ys=ypermuted;

end

%se volevo conservare lo spettro, faccio 1 altro mezzo giro dove impongo solo quello
if stop=='spe'
    faseys=angle(fft(ys));
    fys=my.*(cos(faseys)+j*sin(faseys));
    ys=ifft(fys);ys=real(ys);
    ys=ys-mean(ys);
end




