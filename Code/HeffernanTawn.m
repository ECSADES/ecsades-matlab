classdef HeffernanTawn
    
    %functions for fitting the conditional extremes model of Heffernan and Tawn (2004) "A conditional approach to modelling multivariate extreme values)" 

    properties        
        Prm   % 4 x nBoot, H&T model alpha, beta, mu, sig 
        Rsd   % nBoot x 1 (cell array), residuals from each bootstrap
        Thr   % nBoot x 1, H&T conditional threshold        
        NEP   % 1 x 1, non exceedence probability
        nBoot    % 1 x 1, number of bootstraps (from marginal)            
        X     % n x nBoot, Conditioned variable on Standard Margins
        Y     % n x nBoot, conditioning variable on Standard Margins
        RV % Conditional Return Value Y | X        
        n     % 1 x 1, number of observations
        nDmn  % 1 x 1, number of dimensions
        SmpLclRsdOn=true;   % 1 x 1, flag for sampling residuals locally from cvr bin (1) or globally (0)
    end
    
    properties(Hidden=true)
       FigureFolder='Figures';  %folder to save figures into (default Figures) 
      
    end
    
    properties(SetAccess=protected )
        nAlp     %number of alpha parmaeters
        nPrm      %number of model parameters
        nBin      %number of model nBins
        nRtr   %number of retrun periods
        NonStat   %flag for non stationary alpha
        A         %bin allocation        
        RsdInd   % nBoot x 1 (cell array), index of residuals from each bootstrap
        %% Cross Validation stuff TODO given options to change these defaults
        CVMth=0; %0 Only Cross Validate smoothness for original dataset (fast);
        %1 Cross Validate smoothness for every bootstrap resample (slow),
        nCV=10; %no. cross-validation groups
        nSmth=10;   %no. smoothnesses tried in CV
        SmthLB=-4;   %lower bound (log10)  for smmothness range
        SmthUB=4;   %upper bound (log10)  for smmothness range
        
        SmthSet     %Set of Candidate smoothness params in Cross Validation
        OptSmth     %Optimal smoothness
        CVLackOfFit %Lack of Fit of candidate smoothness params in CV
        MarginType = 'Laplace'; %Laplace or Gumbel
        TwoParamFlag = false;  %2P (alpha, beta) Hefferenan and Tawn model or 4P (alpha, beta , mu, sigma)       
    end
            
    methods
        function HT=HeffernanTawn(Mrg,HTNEP,NonStationary,CV,SmpLclRsdOn)
            %HT=HeffernanTawn(Mrg,CndInd,T...
            %INPUT
            % - Mrg 2 x 1, marignal model structure *output from stage 3
            % - HTNEP (scalar), Non-exceedance probability for conditional
            % - NonStationary flag for non stationary alpha
            %(Optional)
             %CV cross validation structure with control parameter
            
            %% Input checks
            if nargin == 0
                return
            end
            
            if ~isa(Mrg,'MarginalModel')
                error('Mrg should be a nDmn x 1 Marginal Model')
            end
            if numel(Mrg)<=1
                error('Mrg should be a nDmn x 1 Marginal Model')
            end
            HT.nDmn = numel(Mrg);
            
            %check have same number of observations for each margin
            for i = 2:HT.nDmn % loop over margins                
                if size(Mrg(1).Y,1)~=size(Mrg(i).Y,1)
                    error('Marginal 1 and %d should have the same number of observations. Mrg1 nObs = %d; Mrg%d nObs = %d',i,size(Mrg(1).Y,1),i,size(Mrg(i).Y,1));
                end                
            end
            %NEP must lead to postive threshold
            switch Mrg(1).MarginType
                case 'Gumbel' %Gumbel margin
                validateattributes(HTNEP, {'numeric'},{'<=',1,'>=',exp(-exp(0))},'HeffernanTawn','HTNEP',2);  %0<=Tau<=1
                case 'Laplace'%Laplace margin
                validateattributes(HTNEP, {'numeric'},{'<=',1,'>=',0.5},'HeffernanTawn','HTNEP',2);  %0<=Tau<=1
            end
            
            validateattributes(NonStationary,{'logical','numeric'},{'integer'},'HeffernanTawn','HTNEP',3)
            HT.NonStat=logical(NonStationary);  %flag for nonstationary HT model
            if nargin<5
                HT.SmpLclRsdOn=logical(SmpLclRsdOn); %flag for sampling residuals locally from bin or globally
            end
            %check have same number of bootstrap resamples for each margin
            %HT.NonEmptyBins = []; %preallocate
            for i = 2:HT.nDmn % loop over margins                
                %check have same number of bootstrap resamples for each margin
                if Mrg(1).nBoot~=Mrg(i).nBoot
                    error('Marginal 1 and %d should have the same number of bootstrap resamples. Mrg1 nBoot = %d; Mrg%d nBoot = %d',i,Mrg(1).nBoot,i,Mrg(i).nBoot);
                end
                %Check fitted each margin to same bootstrap resamples in                                
                if any(Mrg(1).BSInd(:)~=Mrg(i).BSInd(:))
                    error('Marginal 1 and %d should have the same bootstrap resamples',i);
                end
                %       HT.NonEmptyBins=unique([HT.NonEmptyBins;Mrg(i).NonEmptyBins;Mrg(j).NonEmptyBins]);                
            end
                      
            %% Smoothness parameters
            if nargin>=4 %% Cross Validation parameters (used for generalied Pareto fitting
                if isfield(CV,'CVMth')
                    validateattributes(CV.CVMth, {'numeric'},{'binary'},'MarginalModel','CV.CVMth',7);
                    HT.CVMth=CV.CVMth;
                end
                if isfield(CV,'nCV')
                    validateattributes(CV.nCV, {'numeric'},{'scalar','positive'},'MarginalModel','CV.nCV',7);
                    HT.nCV=CV.nCV;      %number cross-validation groups
                end
                if isfield(CV,'nSmth')
                    validateattributes(CV.nSmth, {'numeric'},{'scalar','positive'},'MarginalModel','CV.nSmth',7);
                    HT.nSmth=CV.nSmth;    %number smoothnesses tried in CV
                    
                    
                end
                if isfield(CV,'SmthLB')
                    validateattributes(CV.SmthLB, {'numeric'},{'scalar'},'MarginalModel','CV.SmthLB',7);
                    HT.SmthLB=CV.SmthLB;   %lower bound (log10)  for smmothness range
                end
                if isfield(CV,'SmthUB')
                    validateattributes(CV.SmthUB, {'numeric'},{'scalar'},'MarginalModel','CV.SmthUB',7);
                    HT.SmthUB=CV.SmthUB;   %upper bound (log10)  for smmothness range
                end
            end
            
            HT.FigureFolder=Mrg(1).FigureFolder;
            
            %% Pre-allocation
            HT.n=size(Mrg(1).BSInd,1);
            HT.nBoot=Mrg(1).nBoot;
            
            HT.NEP = rand(HT.nBoot,1).*(HTNEP(2)-HTNEP(1))+HTNEP(1);    
            HT.nBin=Mrg(1).Bn.nBin;
            if HT.NonStat %get number of model bins (number of alpha parameters).
                if HT.nBin == 1
                    HT.NonStat = logical(false);
                    warning('Cannot run non-stationary HT model with 1 bin, NonStationary set to false')
                end
                HT.nAlp= HT.nBin;  %nbin alpha                
            else  %stationary model case
                HT.nAlp=1;                
            end
            HT.nPrm= HT.nAlp+3;  %HT.nAlp + 3 other parametes
            
            HT.Prm=NaN(HT.nPrm,HT.nDmn-1,HT.nBoot);
            HT.Thr=NaN(HT.nBoot,HT.nDmn-1);
            HT.Rsd=cell(HT.nBoot,1); %different numbers of residuals for each bootstrap so store in cell
            HT.RsdInd=cell(HT.nBoot,1); %bin index of exceedences used in local sampling of residual;            
            
            HT.Y=NaN(HT.n,HT.nDmn-1,HT.nBoot);
            HT.A=NaN(HT.n,HT.nBoot);
            HT.OptSmth=NaN(HT.nBoot,1);
                        
            if HT.NonStat                
                HT.SmthSet=logspace(HT.SmthLB,HT.SmthUB,HT.nSmth); %try range smoothness penalties for sigma varying by bin
            else
                HT.nSmth=1;
                HT.SmthSet=0;    %Switch off smoothness              
            end
            HT.CVLackOfFit=NaN(HT.nSmth,HT.nBoot);
                        
            HT.X=Margins(Mrg(1));
           
            %% Fit H&T Model      
            for iBt=1:HT.nBoot           %loop over bootstrap resamples
                fprintf('Fitting for bootstrap sample %d of %d\n',iBt,HT.nBoot);
                %transform conditioned variable to Standard margins for iBt'th bootstrap
                HT.Thr(iBt)=Mrg(1).INV_Standard(HT.NEP(iBt)); 
                %HT.Thr(iBt)=quantile( HT.X(:,iBt),HT.NEP(iBt)); %find the H&T threshold
                IExc= HT.X(:,iBt)>HT.Thr(iBt);   %threshold exceedences
                %transform conditioning variable to Standard margins for iBt'th bootstrap                          
               
                J=Mrg(1).BSInd(:,iBt);   
                HT.A(:,iBt)=Mrg(1).Bn.A(J);
                                
                HT.RsdInd{iBt}=HT.A(IExc,iBt);
                for iDmn=2:HT.nDmn
                    HT.Y(:,iDmn-1,iBt)=Margins(Mrg(iDmn),iBt); 
                end
                %% Fit Model
                HT=Fit(HT,IExc,iBt);
            end
            
            %% Compute conditional Return Value
            HT=ConditionalReturnValue(HT,Mrg);                     
            
        end %HeffernanTawn constructor
                                       
        function Sml=Simulate(HT,Mrg,nRls)
            %Sml=Simulate(HT,Mrg,nRls)
            %simulate realisations under HT model
            %INPUTS:
            % - nDmn x 1 Mrg, MarginalModel class Mrg (used to transform back to orginal scale)
            % - 1 x 1 nRls, number of realisations under the model   
            %OUTPUT:
            % - Data structure Sml, data simulated from H&T model on Org,
            % Unif and Standard margins
                                              
            I=randi(HT.nBoot,nRls,1);%decide which bootstrap sample to use over all bootstraps HT.nBoot
            
            %% simulate covariate
   
            if Mrg(1).Bn.nBin>1
                Rat=Mrg(1).Rat(:,I);  %get rate of all observations                                
                RatCdf=cumsum(Rat)/sum(Rat); %get rate cdf
                
                %Simulate covariate with right rate
                Sml.A=sum(bsxfun(@gt,rand(1,nRls),RatCdf),1)'+1;  %bin allocation
                Sml.X=SampleCovariateFromBin(Mrg(1).Bn,Sml.A);  %Sample covariate value uniformly from  within bin
            else
                Sml.A=ones(nRls,1);
                Sml.X=rand(nRls,1)*360; %todo think about non periodic case here.
            end
            
            
            J=sub2ind([HT.nBin,HT.nBoot],Sml.A,I);  %common index across bin allocation and bootstrap
            
            %% simulate on standard scale
            
            U= rand(nRls,1); %sample uniform value
            G=INV_Standard(Mrg(1),U);            
            
            Sml.StnMrg=NaN(nRls,HT.nDmn);
            
            if HT.NonStat %non stationary alpha
                tAlp=NaN(nRls,HT.nDmn-1);
                for iDmn=1:HT.nDmn-1
                    t1=permute(HT.Prm(1:HT.nAlp,iDmn,:),[1,3,2]);  %nBin x nDmn - 1 x nBoot                   
                    tAlp(:,iDmn)=t1(J); %nRls x 1 
                end
                %%
            else
                tAlp=permute(HT.Prm(1,:,I),[3,2,1]); % nRls x nDmn-1
            end
            
            for iRls=1:nRls  %loop over realisations
                
                if mod(iRls,1e3)==0
                    fprintf('HT Simulation %d of %d\n',iRls,nRls) %progress printer
                end
                
                iBt=I(iRls); %bootstrap index for curent realisation
                if G(iRls)>HT.Thr(iBt)  %X value above the threshold                     
                    tA=Sml.A(iRls);
                    Z=HT.Rsd{iBt}(HT.RsdInd{iBt}==tA,:);  %nRsd x nDmn-1 joint residuals in each dimension
                    if size(Z,1)>0 && HT.nBin<10 %if binned residuals exist, and nBin not too large, sample resid from bin
                        ISmp= randi(size(Z,1),1);  %sample residual
                    else %too many bins - 
                        Z=HT.Rsd{iBt};  %nRsd x nDmn-1 joint residuals in each dimension
                        ISmp= randi(size(Z,1),1);  %sample residual
                    end
                    Z=Z(ISmp,:);
                    
                    Sml.StnMrg(iRls,1)=G(iRls); %use random gumbel/laplace
                    %get conditioning value from H&T model                   
                    
                    Sml.StnMrg(iRls,2:end)=tAlp(iRls,:).* G(iRls) + G(iRls).^HT.Prm(HT.nAlp+1,:,iBt).* ( HT.Prm(HT.nAlp+2,:,iBt)+ HT.Prm(HT.nAlp+3,:,iBt).*Z);
               
                else %resample below threshold
                    
                    IBlw=find(HT.X(:,iBt)<=HT.Thr(iBt)); %threshold on the gumbel scale
                    if any(IBlw)
                        if numel(IBlw)==1
                            ISmp=IBlw; %TODO check in nDmn
                        else
                            ISmp=randsample(IBlw,1); %TODO check in nDmn
                        end
                        Sml.StnMrg(iRls,1)=HT.X(ISmp,iBt);
                        Sml.StnMrg(iRls,2:end)=HT.Y(ISmp,:,iBt);
                    else
                        Sml.StnMrg(iRls,1)=HT.Thr(iBt);
                        Sml.StnMrg(iRls,2:end)=-Inf;
                    end
                end

            end

            %% Transform to uniform
            Sml.Unf=CDF_Standard(Mrg(1),Sml.StnMrg);          
            
            %% Transform back to original margin
            Sml.Org=NaN(size(Sml.Unf));            
            for iDmn=1:HT.nDmn %loop over dimension                
                Sml.Org(:,iDmn)=Mrg(iDmn).INV(Sml.Unf(:,iDmn),I,Sml.A);                
            end
                                                        
        end %simulate             
        
        function Sml=SimulateIS(HT,Mrg,nRls)
            %Sml=SimulateIS(HT,Mrg,nRls)
            %simulate uniformly space realisations under HT model (used in importance sampling
            %contours)
            %INPUTS:
            % - nDmn x 1 Mrg, MarginalModel class Mrg (used to transform back to orginal scale)
            % - 1 x 1 nRls, number of realisations under the model
            %OUTPUT:
            % - Data structure Sml, data simulated from H&T model on Org,
            % Unif and Standard margins
            
            Sml.I=randi(HT.nBoot,nRls,1);%decide which bootstrap sample to use over all bootstraps HT.nBoot       
            %Sml.I=ones(nRls,1);%decide which bootstrap sample to use over all bootstraps HT.nBoot       
            %% simulate covariate                       
            RatNrm=Mrg(1).Rat./sum(Mrg(1).Rat,1); %nBin x nBoot
            %Simulate covariate with right rate
            Sml.A=randi(HT.nBin,nRls,1);  %bin allocation
            Sml.X=SampleCovariateFromBin(Mrg(1).Bn,Sml.A);  %Sample covariate value uniformly from  within bin
            
            if HT.nBin==1
                J=Sml.I;
                W=RatNrm(J)'; %nRls x 1 %probability weight for chosen bin;
            else                
                J=sub2ind([HT.nBin,HT.nBoot],Sml.A,Sml.I);  %common index across bin allocation and bootstrap
                W=RatNrm(J); %nRls x 1 %probability weight for chosen bin;
            end
            
          
            %% simulate on original scale
            %gp upper endpoint in each dimension   
            
            Sml.Org= rand(nRls,HT.nDmn); %simulation resluts on original scale
            Sml.Unf= NaN(nRls,HT.nDmn);   %simulation resluts on uniform scale
            Sml.StnMrg= NaN(nRls,HT.nDmn); %simulation resluts on standard margins
            Sml.fog=NaN(nRls,HT.nDmn);  %importance weights
            
            %% HT alp parameter
            if HT.NonStat %non stationary alpha
                tAlp=NaN(nRls,HT.nDmn-1);
                for iDmn=1:HT.nDmn-1
                    t1=permute(HT.Prm(1:HT.nAlp,iDmn,:),[1,3,2]);  %nBin x nDmn - 1 x nBoot
                    tAlp(:,iDmn)=t1(J); %nRls x 1 %TODO Check is right in nDmn
                end
                %%
            else
                tAlp=permute(HT.Prm(1,:,Sml.I),[3,2,1]); % nRls x nDmn-1
            end
            
            %% simulated data
            for iDmn=1:HT.nDmn
                %% Find sensible range to draw values over
                [UL,LL,Rng] = Mrg(iDmn).makeRange(Sml.I,Sml.A);
                Sml.Org(:,iDmn)= Sml.Org(:,iDmn).*Rng+LL; %sample uniform value                
                             
                %transform to uniform scale;
                Sml.Unf(:,iDmn)=Mrg(iDmn).CDF(Sml.Org(:,iDmn),Sml.A,Sml.I);   
                %transform to standard margins;
                Sml.StnMrg(:,iDmn)=INV_Standard(Mrg(iDmn),Sml.Unf(:,iDmn));
                            
                %% Density of computed points
                if iDmn==1  %first dimension standard marginal denisty
                    %compute density of chosen points (other dimensions use the HT density)
                    Sml.fog(:,1)=Mrg(1).PDF(Sml.Org(:,1),Sml.A,Sml.I).*Rng.*W; %importance weights
                    fx = Mrg(1).PDF_Standard(Sml.StnMrg(:,1)); %density of values on standard margins
                    IExc=HT.Thr(Sml.I)<Sml.StnMrg(:,1); %index of exceedences
                    nE=sum(IExc); %number of exceedences
                else %conditional density of each subsequent dimension.
                    iAsc=iDmn-1;
                    DimSet=[1,iDmn];                           
                    
                    fyx=NaN(nRls,1);
                    %% below threshold 2D kernel density estimation 
                    Q=[HT.X(:,1),HT.Y(:,iAsc,1)];
                    Z=Sml.StnMrg(:,DimSet);
                    fyx(~IExc) = ksdensity(Q,Z(~IExc,:)); %joint empirical density of data.
                    fyx(~IExc)= fyx(~IExc)./fx(~IExc);
                    fyx(any(isnan(Z),2))=0;
                    %% above threshold  HT density                                                                                
                    %HT parameters;
                    IBtExc=Sml.I(IExc); %bootstrap index of exceedences only                    
                    
                    tBet=squeeze(HT.Prm(HT.nAlp+1,iAsc,IBtExc));
                    tMu=squeeze(HT.Prm(HT.nAlp+2,iAsc,IBtExc));
                    tSgm=squeeze(HT.Prm(HT.nAlp+3,iAsc,IBtExc));
                    
                    % Y= aX+X^b (m+s Z);
                    %  ((Y- aX) -m X ^b)./ s X^b =Z
                    
                    Xpb=Z(IExc,1).^tBet;                                        
                    mu=bsxfun(@minus,bsxfun(@minus,Z(IExc,2),tAlp(IExc,iAsc).*Z(IExc,1)),Xpb.*tMu);                    
                    sgm=Xpb.*tSgm;
                    Z=bsxfun(@rdivide,mu,sgm);
                    %residual desnity comes from kernel deinsty estimation
                    fygx_HT=NaN(nE,1);
                   
                    for iBt=1:HT.nBoot %loop over bootstraps
                        I=IBtExc==iBt; %
                        if any(I)                                                                                    
                            %TODO add localisation of residuals!!                           
                            %TODO jointly sample residuals!!
                            %conditional density
                            fygx_HT(I)=ksdensity(HT.Rsd{iBt}(:,iAsc),Z(I)); %f(y|x) on standard scale
                        end
                    end
                                                                               
                    %f(x,y) = f(y|x)f(x)
                    fyx(IExc)=fygx_HT; %conditional density f(y|x) on laplace scale
                    
                    Sml.fog(:,iDmn)=fyx.*Rng; %importance weights
                end                                
            end %loop over dimensions                                                                     
            
            if any(isnan(Sml.fog(:)))
               error('NaNs detected in importance weights') 
            end
            
        end %SimulateIS
        
        function HT=ConditionalReturnValue(HT,Mrg,nRls)
            %HT=ConditionalReturnValue(HT)
            %compute conditional return value Y | X using monte carlo simuation use return period defined in marginal model.
            if nargin<=2
                nRls=1000; %default number of MC draws
            end
                      
            HT.nRtr=numel(Mrg(1).RtrPrd);
            %preallocate
            nAsc=HT.nDmn-1;
            RV.X_Stn=NaN(HT.nBin,nRls,HT.nRtr);
            RV.X=NaN(HT.nBin,nRls,HT.nRtr);
            RV.Y_Stn=NaN(HT.nBin,nAsc,nRls,HT.nRtr);
            RV.Y=NaN(HT.nBin,nAsc,nRls,HT.nRtr);
               
            I=randi(HT.nBoot,nRls,1);  %bootstrap samples to use
            %draw random resdiuals
           
            if ~HT.SmpLclRsdOn   %if too many bins (few obs per bin), sample residuals globally
                Z=cell2mat(cellfun(@(x)x(randi(numel(x),HT.nBin,1))',HT.Rsd(I,:),'uniformoutput',false))'; %TODO check multivariate cases
                Z=permute(Z,[1,3,2]);
                
            else   %when your bins are big enough (decent no. of obs per bin), sample residuals locally from bin
                Z=NaN(nRls,nAsc,HT.nBin);
                for iBin=1:HT.nBin
                    tRsd=cellfun(@(x,y)y(x==iBin,:),HT.RsdInd(I,:),HT.Rsd(I,:),'uniformoutput',false);
                    J=cellfun(@length,tRsd)>0;
                    Z(J,:,iBin)=cell2mat(cellfun(@(x)x(randi(size(x,1)),:),tRsd(J,:),'uniformoutput',false));
                end
                Z=permute(Z,[3,2,1]);
            end
                           
            for iRtr=1:HT.nRtr %loop over return periods
                            
                %% Sample from return value distribution within each bin     
                rho=Mrg(1).Rat(:,I); %annual rate of occurence
                LT=rho*Mrg(1).RtrPrd(iRtr); %poisson Rate                
                                
                UX=rand(HT.nBin,nRls); %U should be in the range [ P0, 1] where P0 is the non occurence rate.                
                P=1+log(UX)./(LT);    %adjust for return period  (inv of exp(-L*T*(1-C));
                P(bsxfun(@lt,P,HT.NEP(I)'))=NaN;  %this is the non-exceedence rate on standard scale
                %P(P<0)=NaN;  %this is the non-exceedence rate on standard scale
                RV.X(:,:,iRtr)=Mrg(1).INV(P,I); %RVX values on original scale
                
                %transform form uniform to standard margins using CDF                                                                
                RV.X_Stn(:,:,iRtr)=INV_Standard(Mrg(1),P);
                               
                %compute Y using conditional model given X 
                tX=permute(RV.X_Stn(:,:,iRtr),[1,3,2]);
                RV.Y_Stn(:,:,:,iRtr)=bsxfun(@times,HT.Prm(1:HT.nAlp,:,I),tX) + bsxfun(@power,tX,HT.Prm(HT.nAlp+1,:,I)).*bsxfun(@plus,HT.Prm(HT.nAlp+2,:,I),bsxfun(@times,HT.Prm(HT.nAlp+3,:,I),Z));
                                                                                                                 
                for iAsc=1:nAsc
                    %transform Y from standard to uniform margins using CDF   
                    UY=permute(Mrg(iAsc+1).CDF_Standard(RV.Y_Stn(:,iAsc,:,iRtr)),[1,3,2]); %nBin x nAsc x nRls x nRtr
                    %Transform Y to back original margins
                    RV.Y(:,iAsc,:,iRtr)=permute(Mrg(iAsc+1).INV(UY,I),[1,3,2]);
                end                                  
            end
            
            %% Get Omni Value (covariate free) return value and its associated conditions
            if HT.nBin > 1
                
                XOmni_Stn=max(RV.X_Stn,[],1);
                [XOmni,J]=max(RV.X,[],1);    %J is index of location of max
                
                YOmni=NaN(1,nAsc,nRls,HT.nRtr);
                YOmni_Stn=NaN(1,nAsc,nRls,HT.nRtr);
                for iRtr=1:HT.nRtr
                    I=sub2ind([HT.nBin,nRls],J(:,:,iRtr),(1:nRls)); %need composite index of realisation and max location
                    %original margins
                    tYOmni=reshape(permute(RV.Y(:,:,:,iRtr),[1,3,2]),HT.nBin*nRls,nAsc);
                    YOmni(1,:,:,iRtr)=permute(tYOmni(I,:),[3,2,1]);
                    %standard margins
                    tYOmni_Stn=reshape(permute(RV.Y_Stn(:,:,:,iRtr),[1,3,2]),HT.nBin*nRls,nAsc);
                    YOmni_Stn(1,:,:,iRtr)=permute(tYOmni_Stn(I,:),[3,2,1]);
                end
                
                RV.Y=cat(1,RV.Y,YOmni);
                RV.Y_Stn=cat(1,RV.Y_Stn,YOmni_Stn);
                RV.X=cat(1,RV.X,XOmni);
                RV.X_Stn=cat(1,RV.X_Stn,XOmni_Stn);
                
            end
            
            %% Store Return Value simulation
            HT.RV=RV;
            HT.RV.nRls=nRls;
                                    
        end %ConditionalReturnValue
        
        function Plot(HT,Mrg)
            %% Store bin start and end points for plot labels
%             BinSt=Mrg(1).DrcEdg;
%             BinEnd=circshift(Mrg(1).DrcEdg,-1);
%             if Mrg(1).DrcEdg(1) > 0 %if there is a bin which straddles 0, make it bin 1
%                 BinSt=circshift(BinSt,1);
%                 BinEnd=circshift(BinEnd,1);
%             end 
            %% simulate under model
            nRls=length(Mrg(1).Y)*10;
            Sml=Simulate(HT,Mrg,nRls); %simulate 10 times length of original data
            
            %% Plot  data and simulation
            figure;
            clf;  
            for iDmn = 2:HT.nDmn
                subplot(2,HT.nDmn-1,iDmn-1)
                
                plot(Sml.StnMrg(:,1),Sml.StnMrg(:,iDmn),'r.')
                hold on
                plot(HT.X(:,1),HT.Y(:,iDmn-1,1),'k.','markersize',10)
                xlabel(sprintf('%s: Conditioning variable',Mrg(1).RspLbl))
                ylabel(sprintf('%s: Conditioned variable',Mrg(iDmn).RspLbl));
                title(sprintf('%s|%s: H&T simulation on %s margins',Mrg(iDmn).RspLbl,Mrg(1).RspLbl,Mrg(1).MarginType))
                legend('Simulation','Data','location','NorthWest')
                axis tight
                box on
                grid on
            end
            
            for iDmn = 2:HT.nDmn
                subplot(2,HT.nDmn-1,HT.nDmn-1+iDmn-1)
                plot(Sml.Org(:,1),Sml.Org(:,iDmn),'r.')
                hold on
                plot(Mrg(1).Y,Mrg(iDmn).Y,'k.','markersize',10)
                
                xlabel(sprintf('%s: Conditioning variable',Mrg(1).RspLbl))
                ylabel(sprintf('%s: Conditioned variable',Mrg(iDmn).RspLbl));
                title(sprintf('%s|%s: H&T simulation on original scale',Mrg(iDmn).RspLbl,Mrg(1).RspLbl))
                axis tight
                grid on
            end
            
            
            savePics(fullfile(HT.FigureFolder,'Stg4_HT_1_SmlvsData'))
                        
            %% Residual Diagnostics orignal sample
            figure;
            clf;
            tRsd=cell2mat(HT.Rsd);    
                     
            %histogram of redisuals in ith variable
            for iDmn = 1:(HT.nDmn-1)
                subplot((HT.nDmn-1),2,2*iDmn-1)
                if verLessThan('Matlab','8.5')
                    hist(tRsd(:,iDmn));
                    h = findobj(gca,'Type','patch');
                    set(h,'FaceColor',[1 1 1]*0.5);
                    set(h,'linestyle','none')
                else
                    histogram(tRsd(:,iDmn),'edgecolor','none','facecolor','k')
                end
                axis tight
                title(sprintf('%s|%s: resid hist',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl))
                xlabel('Residual')                
            end
            %qq plot of redisuals in ith variable
            for iDmn = 1:(HT.nDmn-1)

                subplot(HT.nDmn-1,2,2*iDmn);
                h=qqplot(tRsd(:,iDmn));
                if verLessThan('Matlab','8.5')                    
                    h = findobj(gca,'Type','patch');
                    set(h,'Marker','.');
                    set(h,'MarkerEdgeColor',[0 0 0]);
                    set(h,'linestyle','--')
                    set(h,'LineWidth',1)
                else
                    h(1).Marker='.';
                    h(1).MarkerEdgeColor=[0,0,0];
                end
                ylabel('Empirical quantiles of residuals')
                xlabel('Standard normal quantiles')
                axis tight
                box on
                title(sprintf('%s|%s: resid QQ-plot',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl)) 
            end
            savePics(fullfile(HT.FigureFolder,'Stg4_HT_2_ResidualDiagnostics1'))
   
            
            %% Residual diagnostics: residuals against covariates
            figure;clf;
            iP = 0; %counter on subplot number
            for iDmn = 1:(HT.nDmn-1)
                for iC=1:Mrg(1).nCvr
                    iP = iP +1;
                    subplot(HT.nDmn-1,Mrg(1).nCvr,iP)
                    IExc=HT.X(:,1)>HT.Thr(1);   %threshold exceedences
                    IExc=IExc(Mrg(1).BSInd(:,1)>0);
                    plot(Mrg(1).X(IExc,iC),HT.Rsd{1}(:,iDmn),'k.')
                    axis tight
                    box on
                    grid on
                    xlabel(Mrg(1).CvrLbl(iC))
                    title(sprintf('%s|%s: Residuals on %s',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl,Mrg(1).CvrLbl{iC}))
                end
            end
        
            
            savePics(fullfile(HT.FigureFolder,'Stg4_HT_2_ResidualDiagnostics2'))
            
            %% Parameter Plot (bootstrap uncertainty)
            LgnLbl=cell(HT.nDmn-1,1);
            for iDmn = 1:(HT.nDmn-1)
                LgnLbl{iDmn} = sprintf('%s|%s',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl);
            end
            Lbl={'\alpha','\beta','m','s'};
            
            for iDmn = 1:(HT.nDmn-1) %different parameters figure for each dimension
                figure;
                clf;
                
                for i=1:4          %loop over parameters
                    j=i+HT.nAlp-1;
                    if HT.nAlp>1 && i==1
                        for iC=1:Mrg(1).nCvr
                            subplot(Mrg(1).nCvr,2,(iC-1)*2+1)
                            hold on
                            
                            tAlp=squeeze(HT.Prm(1:HT.nAlp,iDmn,:));
                            PlotParameter(Mrg(1).Bn,tAlp,iC,'color',[0,0,0],'linewidth',2)
                            
                            xlabel(Mrg(1).CvrLbl(iC))
                            ylabel(Lbl{i})
                            box on
                            grid on
                        end
                    else
                        if i==1
                            subplot(1,2,1)
                        else
                            subplot(3,2,(i-1)*2)
                        end
                        
                        if verLessThan('Matlab','8.5')
                            hist(squeeze(HT.Prm(j,iDmn,:)));
                            h = findobj(gca,'Type','patch');
                            set(h,'FaceColor',[1 1 1]*0.5);
                            set(h,'linestyle','none')
                        else
                            histogram(squeeze(HT.Prm(j,iDmn,:)),'edgecolor','none','facecolor','k')
                        end
                        hold on
                        xlabel(Lbl{i})
                        if i==1 %\alpha
                            switch HT.MarginType
                                case 'Gumbel'
                                    xlim([0,1]);
                                case 'Laplace'
                                    xlim([-1,1]);
                            end
                        end
                        if i==3 %\mu
                            
                            t1=max(abs(HT.Prm(j,iDmn,:)));
                            if t1 < eps
                                t1=0.1;
                            end
                            xlim([-t1,t1]);  %make mu plot centered around 0
                        end
                        if i==4
                            legend(LgnLbl{iDmn},'location','best')
                        end
                    end
                end
                axes('position',[0.1300    0.1100    0.7750    0.8150]);
                axis off
                title(sprintf('%s|%s: Histograms of H&T parameter uncertainty',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl))
                
                savePics(fullfile(HT.FigureFolder,sprintf('Stg4_HT_3_Parameters_%s',Mrg(iDmn+1).RspLbl)))
            end
            
            %% Threshold plot  (alpha as function of NEP)
            %Alpha varies by covariate bin, so to avoid over-crowded plots; have a different version of this plot for each dimension.    
            nPlt1=ceil(sqrt(Mrg(1).Bn.nBin)); %max size nPlt x nPlt
            nPlt2=ceil(Mrg(1).Bn.nBin./nPlt1);
            
 
            
            for iDmn = 1:(HT.nDmn-1)
                figure;
                clf;
                for iAlp = 1:HT.nAlp
                    subplot(nPlt2,nPlt1,iAlp)
                    plot(HT.NEP,squeeze(HT.Prm(iAlp,iDmn,:)),'k.','markersize',20)
                    tAlp=squeeze(HT.Prm(iAlp,iDmn,:));
                    
                    if Mrg(1).nBoot>=20
                        hold on
                        nNEP=10; %number of bins for threshold diagnostic plot
                        NEPEdg=linspace(min(HT.NEP),max(HT.NEP),nNEP+1);
                        NEPBin=(NEPEdg(1:end-1)+NEPEdg(2:end))/2;
                        [~,tA] = histc(HT.NEP,NEPEdg);
                        tA(tA>nNEP)=nNEP;
                        
                        
                        S=accumarray(tA,tAlp,[nNEP,1],@median,NaN);
                        Slb=accumarray(tA,tAlp,[nNEP,1],@(x)quantile(x,0.025),NaN);
                        Sub=accumarray(tA,tAlp,[nNEP,1],@(x)quantile(x,0.975),NaN);
                        
                        NanI=isnan(S);
                        plot(NEPBin(~NanI),S(~NanI),'r-','linewidth',2)
                        plot(NEPBin(~NanI),Slb(~NanI),'r--','linewidth',2)
                        plot(NEPBin(~NanI),Sub(~NanI),'r--','linewidth',2)
                    end
                    grid on

                    if iAlp == 1
                        xlabel('HT NEP')
                        ylabel('\alpha')
                        if HT.nAlp==1 %stationary case
                            title(sprintf('%s|%s: H&T parameter stability by threshold',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl))
                        else
                            title(sprintf('%s|%s: H&T parameter stability by threshold \n Bin %s',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl,Mrg(1).Bn.BinLbl{iAlp}))
                        end
                    else
                        title(sprintf('Bin %s',Mrg(1).Bn.BinLbl{iAlp}))
                    end
                end
                hold off;
                savePics(fullfile(HT.FigureFolder,sprintf('Stg4_HT_4_AlphaThresholdStability_%s',Mrg(iDmn+1).RspLbl)))
            end
           
               
            %% Threshold plot (beta as function of NEP)
            figure; clf;            %All dimensions on one plot: beta is stationary so don't have subplots for covariate bins
            nPlt1=ceil(sqrt(HT.nDmn-1)); 
            nPlt2=ceil((HT.nDmn-1)./nPlt1);
            for iDmn = 1:(HT.nDmn-1)
                subplot(nPlt2,nPlt1,iDmn);
                plot(HT.NEP,squeeze(HT.Prm(HT.nAlp+1,iDmn,:)),'k.','markersize',20)
                tBet=squeeze(HT.Prm(HT.nAlp+1,iDmn,:));
                if Mrg(1).nBoot>20
                    hold on
                    nNEP=10;
                    NEPEdg=linspace(min(HT.NEP),max(HT.NEP),nNEP+1);
                    NEPBin=(NEPEdg(1:end-1)+NEPEdg(2:end))/2;
                    [~,tA] = histc(HT.NEP,NEPEdg);
                    tA(tA>nNEP)=nNEP;
                    n=accumarray(tA,tA,[nNEP,1],@numel,NaN);%#ok<*PROPLC>
                    S=accumarray(tA,tBet,[nNEP,1],@median,NaN);
                    Slb=accumarray(tA,tBet,[nNEP,1],@(x)quantile(x,0.025),NaN);
                    Sub=accumarray(tA,tBet,[nNEP,1],@(x)quantile(x,0.975),NaN);
                    plot(NEPBin(n>0),S(n>0),'r-','linewidth',2)
                    plot(NEPBin(n>0),Slb(n>0),'r--','linewidth',2)
                    plot(NEPBin(n>0),Sub(n>0),'r--','linewidth',2)
                end
                grid on
                xlabel('HT NEP')
                ylabel('\beta')
                title(sprintf('%s|%s: H&T parameter stability by threshold',Mrg(iDmn+1).RspLbl,Mrg(1).RspLbl))
            end
            
            savePics(fullfile(HT.FigureFolder,'Stg4_HT_4_BetaThresholdStability'))
                                   
            %% Lack of fit plot
            if HT.nAlp>1
                figure;
                clf;
                if HT.nBoot>1 && HT.CVMth==1
                   plot(log10(HT.SmthSet),nanmedian(HT.CVLackOfFit,2),'k-','linewidth',2) 
                   hold on
                   plot(log10(HT.SmthSet),quantile(HT.CVLackOfFit,0.025,2),'k--','linewidth',2) 
                   plot(log10(HT.SmthSet),quantile(HT.CVLackOfFit,0.975,2),'k--','linewidth',2) 
                else
                   plot(log10(HT.SmthSet),HT.CVLackOfFit(:,1),'k-','linewidth',2)
                end                
                axis tight
                hold on
                grid on
                plot(log10(median(HT.OptSmth))*[1,1],ylim,'r--','linewidth',2)
                ylabel('Lack of fit')
                xlabel('$\log_{10}(\tilde{\lambda)}$','Interpreter','latex')              
                title(sprintf('HT cross-validation lack of fit'))  
                savePics(fullfile(HT.FigureFolder,'Stg4_HT_5_SectorGoodnessOfFit'))
            end
            
            %% Cond. RV CDFs
            ColMat=hsv(HT.nBin);

            figure;
            clf;
            c=0;
            for iDmn=2:HT.nDmn %loop over associated variables
                for iRtr=1:HT.nRtr
                    c=c+1;
                    
                    subplot(HT.nDmn-1,HT.nRtr,c)
                    grid on
                    % Per bin: -----------
                    if HT.nBin<16  & HT.nBin > 1  %plot directional CDFs when there are no more than 15 bins
                        try %sort() different in matlab versions 
                            tX=sort(permute(HT.RV.Y(1:HT.nBin,iDmn-1,:,iRtr),[1,3,2]),2,'MissingPlacement','first');
                        catch
                            tX=sort(permute(HT.RV.Y(1:HT.nBin,iDmn-1,:,iRtr),[1,3,2]),2);
                        end
                        for iBin = 1:HT.nBin
                            plot(tX(iBin,:),linspace(0,1,HT.RV.nRls),'linewidth',2,'color',ColMat(iBin,:))
                            hold on
                        end
                    end
                    
                    %Omni: --------------
                    try
                        tX=sort(permute(HT.RV.Y(end,iDmn-1,:,iRtr),[1,3,2]),2,'MissingPlacement','first');
                    catch
                        tX=sort(permute(HT.RV.Y(end,iDmn-1,:,iRtr),[1,3,2]),2);
                    end
                    if HT.nBin == 1
                        plot(tX,linspace(0,1,HT.RV.nRls),'k-','linewidth',2)
                    else
                        plot(tX,linspace(0,1,HT.RV.nRls),'--k','linewidth',2)
                    end
                    
                    ylabel('Cumulative Probability')
                    xlabel(sprintf('%s | %s',Mrg(iDmn).RspLbl, Mrg(1).RspLbl))
                    grid on
                    
                    if (HT.nBin<16  & HT.nBin >1 ) & (c == (HT.nDmn-1)*HT.nRtr)
                        legend([Mrg(1).Bn.BinLbl(:);'Omni'],'location','best');
                    end
                    
                    if iDmn==2
                        title(sprintf('Cnd RV CDF %g Years',Mrg(iDmn).RtrPrd(iRtr)))                       
                    end
                end                               
            end
            savePics(fullfile(HT.FigureFolder,'Stg4_HT_6_ConditionalReturnValueCDF'))
        end %plot
        
    end %methods
    methods (Access = private)
        function HT=Fit(HT,IExc,iBt)
            % HT=Fit(HT,IExc,iBt)
            % INPUTS:
            % - n x 1 IExc logical array of exceedances
            % - 1 x 1 iBt bootstrap number
            % OUTPUTS:
            % - HT.Prm(:,:,iBt) fitted H&T parameters for iBt'th bootstrap
            % - HT.Rsd{:,iBt} residuals for iBt'th bootstrap
            
            %% Starting Solution
            p0=repmat([0.5*ones(HT.nAlp,1);0.5;0;1],1,HT.nDmn-1); %starting value (alpha, beta, mu, sigma)
            opts=optimset('display','off');
            
            %% Get exceedances
            nExc=sum(IExc);
            tX=HT.X(IExc,iBt); %exceedances for iBt'th bootstrap
            tY=HT.Y(IExc,:,iBt);  %conditioning value given exceedance in conditioned value
            if HT.NonStat %bin allocation
                tA=HT.A(IExc,iBt);
                if sum(isnan(tA))>0
                   warning('NaNs in allocation') 
                end
            else
                tA=ones(nExc,1);
            end             
            
            %% Flag for do Cross Validation
            if (HT.CVMth==0 &&  iBt>1) || HT.nSmth==1 
                %% Cross validation off
                if HT.nSmth==1  %only one smoothness
                   HT.OptSmth(iBt)=HT.SmthSet;
                else  %use first bootstrap instead of cross validating every time
                    HT.OptSmth(iBt)=HT.OptSmth(1);
                end                                                                                
            else
                fprintf('Starting Cross Validation to estimate HT alpha smoothness:\n');
                LackOfFit=NaN(HT.nCV,HT.nSmth); %lack of fit
                ICV=randi(HT.nCV,nExc,1); %split nExc observations into nCV cross validation groups
                for iCV=1:HT.nCV  %cross validation loop
                    fprintf('.')
                    %Pull out fit set
                    tXFit=tX(ICV~=iCV); %main response
                    tYFit=tY(ICV~=iCV,:); %associated variable
                    tAFit=tA(ICV~=iCV,:);
                    %Pull out prediction set
                    tXPrd=tX(ICV==iCV);
                    tYPrd=tY(ICV==iCV,:);
                    tAPrd=tA(ICV==iCV,:);
                    for iLam=1:HT.nSmth %Smoothness penalties
                        if HT.TwoParamFlag %2 parameter HT
                            P12=fminsearch(@(p)HeffernanTawn.likelihood(tXFit,tYFit,tAFit,p,HT.MarginType,HT.SmthSet(iLam),HT.TwoParamFlag),p0(1:2),opts);
                            PrmHat=[P12;zeros(HT.nDmn-1,1)';ones(nDmin-1,1)'];
                            LackOfFit(iCV,iLam)=HeffernanTawn.likelihood(tXPrd,tYPrd,tAPrd,PrmHat,HT.MarginType);
                        else %4 parameter HT
                            PrmHat=fminsearch(@(p)HeffernanTawn.likelihood(tXFit,tYFit,tAFit,p,HT.MarginType,HT.SmthSet(iLam),HT.TwoParamFlag),p0,opts);
                            LackOfFit(iCV,iLam)=HeffernanTawn.likelihood(tXPrd,tYPrd,tAPrd,PrmHat,HT.MarginType);
                        end
                    end %SigPen
                end %CV
                fprintf('\n')
                HT.CVLackOfFit(:,iBt)=sum(LackOfFit,1)';
                [~,tI]=min(HT.CVLackOfFit(:,iBt));
                HT.OptSmth(iBt)=HT.SmthSet(tI)';                 
            end
            
            %% fit model: find optimal parameter values for heffernan and tawn model         
            if HT.TwoParamFlag %2 parameter HT
                tHTPrm=fminsearch(@(p)HeffernanTawn.likelihood(tX,tY,tA,p,HT.MarginType,HT.OptSmth(iBt),HT.TwoParamFlag),p0(1:2),opts);
                tHTPrm=[tHTPrm;zeros(HT.nDmn-1,1)';ones(nDmin-1,1)'];
            else
                tHTPrm=fminsearch(@(p)HeffernanTawn.likelihood(tX,tY,tA,p,HT.MarginType,HT.OptSmth(iBt),HT.TwoParamFlag),p0,opts);
            end            
            HT.Prm(:,:,iBt)=tHTPrm;
            
            %% keep residuals
            tAlp=HT.Prm(1:HT.nAlp,:,iBt);
            tXpB=tX.^HT.Prm(HT.nAlp+1,:,iBt);
            mn=tY-tAlp(tA,:).*tX-HT.Prm(HT.nAlp+2,:,iBt).*tXpB;
            std=HT.Prm(HT.nAlp+3,:,iBt).*tXpB;
            HT.Rsd{iBt}=mn./std;
            if any(~isreal(HT.Rsd{iBt}(:)))
                warning('Complex Residuals found!!\n')
            end
            
        end %fit
        
    end %Methods private
    methods(Static)

       
        
        function PLOGL=likelihood(X,Y,A,p,MrgTp,L,TwoPrmFlg)
            %function PLOGL=likelihood(X,Y,A,p,MrgTp,L,TwoPrmFlg)
            %compute penalised Heffernan and Tawn likelihood for data on Gumbel scale.
            %INPUT
            %-X         n x 1 conditioned value
            %-Y         n x 1 x (nDmn-1) conditioning value
            %-A         n x 1 x (nDmn-1) bin allocation
            %-p         nBin+3 x (nDmn-1) parameter values
            %-L         smoothness parameter (alpha)
            %-MrgTp     Gumbel or Laplace
            %-TwoPrmFlag TwoParameters HT Model if true, else 4 parameter     
            if nargin<=5
                L=0; %no penalty
            end
                        
            if nargin<=6
                TwoPrmFlg=false;
            end
            
            n = size(p,1);
            nDmn = size(p,2)+1;
            
            if TwoPrmFlg
                nBin=n-1;
                Alp=p(1:nBin,:);
                Bet=p(nBin+1,:);
                M=zeros(nDmn-1,1);
                S=ones(nDmn-1,1);
            else
                nBin=n-3;
                Alp=p(1:nBin,:); %nBin x nD-1
                Bet=p(nBin+1,:); %1 xnD-1
                M=p(nBin+2,:);  %1 x nD -1
                S=p(nBin+3,:); % 1 x nD-1
            end
            
            %nD checks need unwrap Alp(:)
            if any(Alp(:)>1) || any(Alp(:)<-1) || any(S<0) || any(Bet>1)  %constraints
                PLOGL=Inf;
                return
            end
            
            if strcmp(MrgTp,'Gumbel')
                if any(Alp(:)<0)
                    PLOGL=Inf;
                    return
                end
            end
            
            Xb=bsxfun(@power,X,Bet);  %X^b
            Std=bsxfun(@times,S,Xb); %standard deviation
            
            Y = reshape(Y,[],nDmn-1);
            Al = Alp(A,:);
                       
            NLOGL=sum(sum(0.5*((Y-bsxfun(@times,Al,X)-bsxfun(@times,M,Xb)).^2./(Std).^2)+log(Std)));
                         
            PLOGL=NLOGL+L.*sum(sum((Alp-mean(Alp,1)).^2,1),2);     %Penalised Likelhiood
        end %likelihood
        
         
    end
end