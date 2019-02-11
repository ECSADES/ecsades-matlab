function F=rndsymlgs(alp,nD,n)
%generate random values from nD dimension symmetric logistics function with
%Frechet margins
%
%Algorithm 1.1 Alec Stephenson, Simulating Multivariate Extreme Value
%Distributions of Logistic Type, Extremes 6, 49-59, 2003
%INPUT
%alp is dependency parameter
%nD  is number of dimensions
%n is the number of realisations
%OUTPUT
%F random values

%step 1)  get T from unit exponentials
W=exprnd(1,[nD,n]);
T=bsxfun(@rdivide,W,sum(W));

%find mixture probabilities from recurrence relationship
p=NaN(nD,nD);
p(1,1)=1;

for iD=2:nD
    p(iD,1)=gamma(iD-alp)./(gamma(iD).*gamma(1-alp));
    for jD=2:(iD-1)
        t1=(iD-1-alp*jD).*p(iD-1,jD)+alp.*(jD-1).*p(iD-1,jD-1);
        p(iD,jD)=t1./(iD-1);
    end
    p(iD,iD)=alp.^(iD-1);
end
P=cumsum(p(end,:)');  %cumulative probability

%step2) Find k

U=rand(1,n);
k=sum(bsxfun(@gt,U,P))+1;
%step3) get Z from gamma(k,1);
Z=gamrnd(k,1,[1,n]);
%step4) find F on frechet margins
F=1./(bsxfun(@times,Z,(T.^alp)));