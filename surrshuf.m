%genera x surrogato con sample shuffling - distruggeanche gli autospettri

function xs=surrshuf(x)
sx=size(x);
p=randperm(sx(1));
xs=zeros(sx(1),1);
for k = 1:sx(1)
	xs(k)=x(p(k));
end
