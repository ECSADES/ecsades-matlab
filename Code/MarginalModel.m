classdef MarginalModel
    %Fit Piecewise constant EVA Model
    %
    % o Estimates piecewise constant threshold using local bin quantile given non-exceedance
    % probability (NEP).
    % o Annual rate calculated using local Possion MLE for each bin = (no. observations)/(no. years of data
    % o Generalised Pareto distribution fitted above the threshold:
    %     - Shape assumed constant across all bins
    %     - Scale can varies across covariate bins (constant within bins)
    %     - Variation in scale param across bins controlled by smoothness penalty
    %     estimated using cross validation
    
    properties
        %% Data
        X        %nDat x nCvr, covariate direction data (X) [0,360]
        Y        %nDat x 1, response data (Y)
        Yrs=1;   %1 x 1, number of years of data
        RspLbl      %1 x 1, (string), reponse label
        RspSavLbl   %1 x 1 (string), reponse label for saving plots (Aviod special characters)
        CvrLbl      %nCvr x 1, (string), covariate label vecotr
        nB=100;  %1 x 1, number of bootstrap resesamples
        RtrPrd=100; %nRtr x 1 Return Period  in years
        
        %% Parameters
        Bn   % Bin class defining binning properties
        Scl      %nBin x 1, Fitted Generalised Paraeto Scale nBin x nB
        Shp      %1 x nBoot, Fitted Generalised Paraeto Shape (stationary shape from PPC marginal fitting)
        %nBin x nBoot (nonstationary shape)
        Alp      %nBin x nBoot   gamma parameter
        Bet      %nBin x nBoot   gamma parameter
        GmmLct
        
        NEP      %nBoot x 1, Non Exceedence Probability Quantile threshold
        Thr      %nBin x 1, Extreme value threshold nBin x nB
        RatExc   %nBin x nBoot, count no. exceedence observations in each bin nBin x nB
        RatBlw   %nBin x nBoot, count no. non-exceedence observations in each bin nBin x nB
        BSInd    %nDat x nBoot, bootstrap index;
        
        nCvr      % 1 x 1, number of covariates
        nDat      % 1 x 1, number of data obs
        nRtr=1     % 1 x 1 number of return periods
        %% Return Value
        RVPrb    %(nBin+1) x nRVX return value cumulative probabilities for specified return period
        RVX      %nRVX x 1 return value response value
        RVMed    %(nBin+1) x 1 median return value (interpolate CDF)
       
        nRVX=400;  %number of points at which the return value is computed.
        MarginType = 'Laplace' %Laplace or Gumbel
    end
    
    properties(Hidden=true)
        FigureFolder='Figures';  %folder to save figures into (default Figures)
        BnMax         %bin max (used in check for upper endpoint)
    end
    
    properties (Access=protected)
        IsPrd   %nCvr x 1  flag for if it is a periodic covariate
        %% Cross Validation (CV) options
        CVMth=0;      %1 x 1 (boolean), CV method
        %=0 Cross Validate smoothness parameter on original dataset only (fast);
        %=1 Cross Validate smoothness for every bootstrap resample (slow),
        nCV=10;       %1 x 1, no. cross-validation groups
        nSmth=10;     %1 x 1, no. smoothnesses tried in CV
        SmthLB=-4;    %1 x 1, lower bound (log10)  for smmothness range
        SmthUB=4;     %1 x 1, upper bound (log10)  for smmothness range
        SmthSet       %nSmth x 1, Set of Candidate smoothness params in Cross Validation
        OptSmth       %1 x 1, Optimal smoothness from SmthSet
        CVLackOfFit   %nSmth x nBoot, Lack of Fit of candidate smoothness params in CV
        
    end
    
    methods
        function MM=MarginalModel(Dat,iDmn,NEP,Bn,nB,Yrs,RtrPrd,CV,MarginType)
            %MM=MarginalModel(X,NEP,Y,BinEdg,IsPrd,nB,Yrs)
            %INPUTS:
            % -Dat structure  (from stage 1)
            %     - Dat.X     nDat x nCvr  covariate values
            %     - Dat.Y     nDat x nDmn  response data (Y)
            %     - Dat.IsPrd   nCvr x 1 vector if covairte is periodic or not
            %     - Dat.CvrLbl    char string for reponse label
            %     - Dat.RspLbl    char string for reponse label
            % - iDmn  Dimension of Y to use as reponse
            % - NEP   1x1 Non Exceedence Probability Quantile threshold
            % -Bn CovariateBinning class containing binning information
            %(Optional)
            % - nB number of bootstrap resamples (assumed 1 if not specified)
            % - Yrs  number of years of data (assumed 1 if not specified)
            %- RtrPrd  return period (years) (assumed to be 100 if non speicifed)
            % - CV cross validation structure with control parameters
            % - MarginType, 'Gumbel', 'Laplace'
            %OUTPUT:
            % - MM, Marginal Model class containing details of data and
            % fitted piecewise constant marginal model
            
            if nargin == 0 %empty constructor
                return;
            end
            
            %% Check inputs
            [MM.nDat,MM.nCvr]=size(Dat.X);
            validateattributes(Dat.IsPrd,{'numeric','logical'},{'numel',MM.nCvr,'integer'},'BinAllocation','IsPrd',1);
            MM.IsPrd=Dat.IsPrd;
            
            for iC=1:MM.nCvr
                if MM.IsPrd(iC)  %periodic check
                    validateattributes(Dat.X(:,iC), {'numeric'},{'vector','<=',360,'>=',0},'MarginalModel','X',1);
                end
            end
            MM.X=Dat.X;
            validateattributes(iDmn,{'numeric'},{'numel',1,'integer','positive'},'BinAllocation','Edg',2);
            validateattributes(Dat.Y(:,iDmn), {'numeric'},{'vector','numel', MM.nDat},'MarginalModel','Y',1);
            MM.Y=Dat.Y(:,iDmn);
            MM.RspLbl=Dat.RspLbl{iDmn};
            MM.RspSavLbl=Dat.RspLbl{iDmn};
            MM.CvrLbl=Dat.CvrLbl;
            
            
            validateattributes(NEP, {'numeric'},{'<=',1,'>=',0},'SimulateMarginal','nDat',3);  %0<=Tau<=1 nep range
            if ~isa(Bn,'CovariateBinning')
                error('Input Bn, should be of class: CovariateBinning');
            end
            MM.Bn=Bn;
            
            %% Optional inputs
            if nargin>=5
                validateattributes(nB, {'numeric'},{'scalar','nonnegative'},'MarginalModel','nB',5);
                MM.nB=nB;  %number of bootstraps
                if MM.nB==0  %case where bootstrapping is off!!
                    MM.nB=1;
                end
            end
            MM.NEP = [range(NEP)/2+min(NEP) ;sort(rand(MM.nB-1,1)).*range(NEP)+min(NEP)];  %sample NEPs over range with middle of range first
            if nargin>=6
                validateattributes(Yrs, {'numeric'},{'scalar','positive'},'MarginalModel','Yrs',6);
                MM.Yrs=Yrs;
            end
            if nargin>=7
                validateattributes(RtrPrd, {'numeric'},{'vector','positive'},'MarginalModel','RV',7);
                MM.RtrPrd=sort(RtrPrd);
                MM.nRtr=numel(RtrPrd);
            end
            
            if nargin>=9
                validateattributes(MarginType,{'string','char'},{},'MarginalModel','MarginType',9)
                MM.MarginType = validatestring(MarginType,{'Gumbel','Laplace'});
            end
            
            %% Check cross-validation inputs
            if nargin>=8 %% Cross Validation parameters (used for generalied Pareto fitting
                if isfield(CV,'CVMth')
                    validateattributes(CV.CVMth, {'numeric'},{'binary'},'MarginalModel','CV.CVMth',8);
                    MM.CVMth=CV.CVMth;
                end
                if isfield(CV,'nCV')
                    validateattributes(CV.nCV, {'numeric'},{'scalar','positive'},'MarginalModel','CV.nCV',8);
                    MM.nCV=CV.nCV;      %number cross-validation groups
                end
                if isfield(CV,'nSmth')
                    validateattributes(CV.nSmth, {'numeric'},{'scalar','positive'},'MarginalModel','CV.nSmth',8);
                    MM.nSmth=CV.nSmth;    %number smoothnesses tried in CV
                end
                if isfield(CV,'SmthLB')
                    validateattributes(CV.SmthLB, {'numeric'},{'scalar'},'MarginalModel','CV.SmthLB',8);
                    MM.SmthLB=CV.SmthLB;   %lower bound (log10)  for smmothness range
                end
                if isfield(CV,'SmthUB')
                    validateattributes(CV.SmthUB, {'numeric'},{'scalar'},'MarginalModel','CV.SmthUB',8);
                    MM.SmthUB=CV.SmthUB;   %upper bound (log10)  for smmothness range
                end
            end
            
            
            %% Preallocate Output Arrays
            MM.Scl=NaN(MM.Bn.nBin,MM.nB);  %Fitted Generalised Pareto Scale nBin x nB
            MM.Shp=NaN(MM.nB,1);    %Fitted Generalised Paraeto Shape  1 x nB
            MM.Alp=NaN(MM.Bn.nBin,MM.nB); %Gamma Parameter
            MM.Bet=NaN(MM.Bn.nBin,MM.nB); %Gamma Parameter nBin x nB
            MM.GmmLct=NaN(MM.Bn.nBin,1); %Gamma Location Parameter nBin x 1
            MM.Thr=NaN(MM.Bn.nBin,MM.nB); %Extreme value threshold nBin x nB
            MM.RatExc=NaN(MM.Bn.nBin,MM.nB); %Annual Exceedence Rate no. exceedence observations in each bin nBin x nB
            MM.RatBlw=NaN(MM.Bn.nBin,MM.nB);  %Annual non-exceedence rate observations in each bin nBin x nB
            MM.OptSmth=NaN(MM.nB,1);  %Fitted smoothness parameter (smoothness in scale param across bins)
            MM.CVLackOfFit=NaN(MM.nSmth,MM .nB); %Lack of fit associated with different smoothness parameters
            MM.BSInd=NaN(numel(MM.Y),MM.nB);  %Bootstrap sample indices
            MM.SmthSet=logspace(MM.SmthLB,MM.SmthUB,MM.nSmth); %try range smoothness penalties for sigma varying by bin
            
            %% Fit model
            MM = Fit(MM);
            
            %% Return Value
            MM=ReturnValueCDF(MM);
            %% Plots
            if ~exist('Figures','dir')
                mkdir('Figures');
            end
            Plot(MM)
            
        end %MarginalModel constructor
        
        function MM=Fit(MM)
            %MM=Fit(MM)
            %Bootstrap fit the threshold, rate and generalised Pareto model.
            %Generalised Pareto estimates smoothness of scale using k-fold cross
            %validation
            
            rng(1); %reset random seed to get same bootstrap sample for all marginals
            MM.BSInd=[(1:length(MM.Y))',randi(length(MM.Y),length(MM.Y),MM.nB-1)];
            rng('shuffle');
            MM.BnMax=accumarray(MM.Bn.A,MM.Y,[MM.Bn.nBin,1],@max,NaN);
                        
            MM.GmmLct=accumarray(MM.Bn.A,MM.Y,[MM.Bn.nBin,1],@(x)min(x)-range(x)*0.01,0); %threshold (simple quantile in each bin)
                            
            %% Bootstrap loop
            for iBt = 1:MM.nB %will end up with MM.nB+1 sample parameters (1 for orig data, MM.nB boostrap resampled data)
                if MM.nB>1
                    fprintf('Fitting for bootstrap sample %d of %d\n',iBt,MM.nB);
                else
                    fprintf('Fitting sample\n');
                end
                %% Bootstrap Allocation
                [tY,A]=GetBootstrapSample(MM,MM.BSInd(:,iBt));
                                
                %% Fit gamma distribution
                for iBn=1:MM.Bn.nBin
                    I=A==iBn;
                    if any(I)
                        p=gamfit(tY(I)-MM.GmmLct(iBn));
                        
                        MM.Alp(iBn,iBt)=p(1);
                        MM.Bet(iBn,iBt)=p(1).*p(2); %orthogonal parameterisation.
                    end
                end

                %% Threshold                
                MM.Thr(:,iBt)=MM.gaminv(MM.NEP(iBt),MM.Alp(:,iBt),MM.Bet(:,iBt),MM.GmmLct); %threshold (simple quantile in each bin)
                
                %% Get Exceedences
                IExc=(tY>MM.Thr(A,iBt));   %Index of exceedenses
                AExc= A(IExc);  %bins allocation of excdeences
                
                Rsd=tY(IExc)-MM.Thr(AExc,iBt);   %Y-u above threshold used for gpfit
                ObsMax=MM.BnMax(AExc)-MM.Thr(AExc,iBt); %observation max used in upper end point of gp
                
                %% Rate of occurence above and below threshold
                MM.RatExc(:,iBt)=accumarray(AExc,AExc,[MM.Bn.nBin,1],@numel)./MM.Yrs;
                MM.RatBlw(:,iBt)=accumarray(A(~IExc),A(~IExc),[MM.Bn.nBin,1],@numel)./MM.Yrs;
                
                %% Generalised Pareto Fit
                MM=GPCrossValidation(MM,Rsd,ObsMax,AExc,iBt);
                
            end
            
            
        end %BootMargModel
        
        function [YMrg,YUnif]=Margins(MM,iBt)
            %[YMrg,YUnif]=Margins(MM,iBt)
            %% Transform response Y to Gumbel Margins (GP -> Unif -> Gumbel)
            %INPUT
            % - (Optional) 1 x 1 iBt, index on bootstrap resample (default 1)            
            %OUTPUT
            % - nDat x 1 YGmbl,  response data (Y) on Gumbel margins
            % - nDat x 1 YUnif,  response data (Y) on Uniform margins
            
            if nargin==1
                iBt=1:MM.nB;  %default to orginal data
            end        
            
            YUnif=NaN(numel(MM.Y),numel(iBt));
            
            for i=1:numel(iBt)     
                jBt=iBt(i);
                tY=MM.Y(MM.BSInd(:,jBt));
                tA=MM.Bn.A(MM.BSInd(:,jBt));
                 
                if MM.Shp(jBt)<0
                    if max(tY-MM.Thr(tA,jBt)+MM.Scl(tA,jBt)./MM.Shp(jBt))>0
                        error('Upper end point invalid');
                    end
                end                      
                YUnif(:,i)=CDF(MM,tY,tA,jBt);
            end
            %% transform from uniform to standard margins            
            YMrg=INV_Standard(MM,YUnif);                        
            
            if any(isinf(YMrg(:)))
                error('Bad transformation to %s',MM.MarginType)
            end
        end %Margins
        
        function MM=ReturnValueCDF(MM,RtrPrd)
            %MM=ReturnValueCDF(MM);
            %Marignal Return Value CDF for PPC model
            %OUTPUTS:
            % - MM.RVPrb nBin+1 x 100 CDF values for each covariate bin;
            % final row is omni-directional
            % - MM.RVX 100 x 1 Response values at which CDF calculated
            if nargin < 2
                RtrPrd = MM.RtrPrd;
            end
            
            %% CDF range
            CDF_X=permute(linspace(min(min(MM.Y(:)),0),max(MM.Y(:))*3,MM.nRVX),[3,1,2]); %choose range over which to compute return value;                   
         
            %% Empirical_GP Cdf
            P=CDF(MM,CDF_X);
            
            MM.RVX=CDF_X(:);
            MM.RVPrb=NaN(MM.Bn.nBin+1,MM.nRVX,MM.nRtr);
            MM.RVMed=NaN(MM.Bn.nBin+1,MM.nRtr);%return value median 
            MM.RVMedStn=NaN(MM.Bn.nBin+1,MM.nRtr); %return value median on standard scale
            %% Return Value CDF
            for iRtr=1:MM.nRtr
                R=exp(-bsxfun(@times,(MM.RatExc+MM.RatBlw).*RtrPrd(iRtr),(1-P))); %Directional RV CDF
                R(end+1,:,:)=prod(R,1); %#ok %get omni return value
               
                MM.RVPrb(:,:,iRtr)=squeeze(mean(R,2));  %posterior predictive return value                
               
                for iBin=1:MM.Bn.nBin+1
                    %original margins
                    [tP,I]=unique(MM.RVPrb(iBin,:,iRtr));
                    tX=MM.RVX(I);
                    MM.RVMed(iBin,iRtr)=interp1(tP,tX,0.5);                  
                end                
            end
            
                                                          
        end %ReturnValue
               
        function P=CDF(MM,X,A,I)            
            %CDF function for marginal model  (empirical below threshold - GP above)
            %CDF(MM,X) compute CDF for all bins and bootstraps using original data X is the locations to compute cdf at
            %CDF(MM,X,A) compute CDF for all bins and bootstraps using original data at specific
            %bins A
            %CDF(MM,X,A,I) compute CDF for all bins and bootstraps using original data at specific
            %bins and bootstraps indexes I                         
            if nargin==3  %Y not specified but A specifed                
                I=1:MM.nB; 
            end                     
           
            if nargin<=2
                P=MarginalModel.gamgpcdf(X,MM.Shp',MM.Scl,MM.Thr,MM.Alp,MM.Bet,MM.GmmLct,MM.NEP');             
            else %specified specfic bins and bootstraps           
                if numel(A)==numel(I) && numel(A)>1
                    J=sub2ind([MM.Bn.nBin,MM.nB],A,I);
                    P=MarginalModel.gamgpcdf(X,MM.Shp(I),MM.Scl(J),MM.Alp(J),MM.Bet(J),MM.GmmLct(A),MM.NEP(I));
                else
                    P=MarginalModel.gamgpcdf(X,MM.Shp(I)',MM.Scl(A,I),MM.Thr(A,I),MM.Alp(A,I),MM.Bet(A,I),MM.GmmLct(A),MM.NEP(I)');
                end
            end
            
        end %CDF
        
        function P=PDF(MM,X,A,I)
            %PDF function for marginal model  (empirical below threshold - PDF above)
            %PDF(MM,X) compute PDF for all bins and bootstraps using original data X is the locations to compute cdf at
            %PDF(MM,X,A) compute PDF for all bins and bootstraps using original data at specific
            %bins A
            %PDF(MM,X,A,I) compute PDF for all bins and bootstraps using original data at specific
            %bins and bootstrap indexes I                         
            if nargin==3  %Y not specified but A specifed                
                I=1:MM.nB; 
            end                     
           
            if nargin<=2
                P=MarginalModel.gamgppdf(X,MM.Shp',MM.Scl,MM.Thr,MM.Alp,MM.Bet,MM.GmmLct,MM.NEP');
            else %specified specifc bins and bootstraps
                if numel(A)==numel(I) && numel(A)>1
                    J=sub2ind([MM.Bn.nBin,MM.nB],A,I);
                    P=MarginalModel.gamgppdf(X,MM.Shp(I),MM.Scl(J),MM.Alp(J),MM.Bet(J),MM.GmmLct(A),MM.NEP(I));
                else
                    P=MarginalModel.gamgppdf(X,MM.Shp(I)',MM.Scl(A,I),MM.Thr(A,I),MM.Alp(A,I),MM.Bet(A,I),MM.GmmLct(A),MM.NEP(I)');
                end
            end
            
        end %PDF
        
        function X=INV(MM,P,I,A)
            %Inverse CDF for marginal model  (empirical below threshold - GP above)
            %P probability           
            %I index of bootstraps to use
            %A index of bins to use
            %if I scalar --> case where finding inverse CDF in single bin
            %if I vector --> case where inverse CDF in across sampled bins and bootstraps
            %if P matrix --> case where finding inverse CDF for all bootstraps and bins
            
            X=NaN(size(P));
            p=size(P);
            if numel(I)==1
                Cs=1; %scalar
            elseif prod(p(2:end))>1
                Cs=3; %matrix                
            else
                Cs=2; %vector
            end
                                 
            switch Cs
                case 1 %I scalar --> case where finding inverse CDF in single bin
                    tNEP=MM.NEP(I);  %non exceedence probability
                    IExc=P>tNEP;  %index of exceedences          
                    tP=(P(IExc)-tNEP)./(1-tNEP); %rescale NEP for conditional distribution
                    X(IExc)=MarginalModel.gpinv(tP,MM.Shp(I),MM.Scl(A,I),MM.Thr(A,I));
                    
                    X(~IExc)=MarginalModel.gaminv(P(~IExc),MM.Alp(A,I),MM.Bet(A,I),MM.GmmLct(A,I));
                case 2  %I vector --> case where inverse CDF in across sampled bins and bootstraps
                    tNEP=MM.NEP(I);  %non exceedence probability
                    IExc=P>tNEP;  %index of exceedences
                    tP=(P(IExc)-tNEP(IExc))./(1-tNEP(IExc));
                    
                    if MM.Bn.nBin==1
                        X(IExc)=MarginalModel.gpinv(tP,MM.Shp(I(IExc)),MM.Scl(I(IExc))',MM.Thr(I(IExc))');
                        X(~IExc)=MarginalModel.gaminv(P(~IExc),MM.Alp(I(~IExc))',MM.Bet(I(~IExc))',MM.GmmLct);
                    else
                        J=sub2ind([MM.Bn.nBin,MM.nB],A,I);
                        JExc=J(IExc); %pull out indices from common bin allc x bootstrap corresponding to exceedences
                        JBlw=J(~IExc); %pull out indices from common bin allc x bootstrap corresponding to exceedences
                        
                        
                        X(IExc)=MarginalModel.gpinv(tP,MM.Shp(I(IExc)),MM.Scl(JExc),MM.Thr(JExc));
                        X(~IExc)=MarginalModel.gaminv(P(~IExc),MM.Alp(JBlw),MM.Bet(JBlw),MM.GmmLct(A(~IExc)));
                    end
                                      
                case 3 %I matrix --> case where finding inverse CDF for all bootstraps and bins
                    tNEP=MM.NEP(I)';  %non exceedence probability
                    IExc=bsxfun(@gt,P,tNEP);  %index of exceedences
                    tP=bsxfun(@rdivide,bsxfun(@minus,P,tNEP),(1-tNEP));    
                    tP(isnan(tP) | tP<0 )=0;
                    X=MarginalModel.gpinv(tP,MM.Shp(I)',MM.Scl(:,I),MM.Thr(:,I));                         
                    %% Resample below threshold
                    X2=MarginalModel.gaminv(tP,MM.Alp(:,I),MM.Bet(:,I),MM.GmmLct);                                             
                    for iBin=1:MM.Bn.nBin
                        if any(~IExc(iBin,:))
                            X(iBin,~IExc(iBin,:))=X2(iBin,~IExc(iBin,:));
                        end
                    end                    
            end
        end %INV
        
        function X=INV_Standard(MM,P)
            %transform from uniform to standard margins
            %using inverse CDF
            switch MM.MarginType
                case 'Gumbel'
                    X = -log(-log(P));
                case 'Laplace'
                    X = sign(0.5-P).*log(2*min(1-P,P)); 
                otherwise
                    error('Margin Type not recognised')
            end
                    
        end %INV_Standard
            
        function P=CDF_Standard(MM,X)
            %transform from standard to uniform margins using CDF
            switch MM.MarginType
                case 'Gumbel'
                    P = exp(-exp(-X));
                case 'Laplace'
                    P = (X>0)-.5*sign(X).*exp(-abs(X));
                otherwise
                    error('Margin Type not recognised')
            end
        end %CDF_Standard
        
        function F=PDF_Standard(MM,X)
            %get probability density on standard margins
            switch MM.MarginType
                case 'Gumbel'
                    F = exp(-(X +exp(-X)));
                case 'Laplace'
                    F =0.5.*exp(-abs(X));
                otherwise
                    error('Margin Type not recognised')
            end
            F(isnan(F))=0;
          end %PDF_Standard
        
        function Plot(MM)
            
            %Plot raw data on all margins you have calculated so far
            mThr=median(MM.Thr,2);
            IExc=(MM.Y>mThr(MM.Bn.A));   %Index of exceedenses
            
            %% Transform Orignal Sample to Standard Margins
            [YMrg,YUnif]=Margins(MM,1);
            nSubPlt = 1+2*(~isempty(YUnif)); %number of subplots; if have transformed raw data onto gumbel margins,
            
            
            figure(1) %Raw response by (Drc) bins
            clf;
            for iC=1:MM.Bn.nCvr
                subplot(MM.Bn.nCvr,nSubPlt,1+3*(iC-1))
                plot(MM.X(IExc,iC),MM.Y(IExc),'k.')
                hold on
                plot(MM.X(~IExc,iC),MM.Y(~IExc),'.','color',[1,1,1]*0.7)
                if iC==1
                    title(sprintf('%s: Raw data', MM.RspLbl))
                end
                
                PlotBinEdge(MM.Bn,iC);
                PlotParameter(MM.Bn,MM.GmmLct,iC,'color','r','linewidth',2);
                PlotParameter(MM.Bn,MM.Thr,iC,'color','b','linewidth',2);
                xlabel(MM.CvrLbl(iC))
                ylabel(MM.RspLbl)
                axis tight
                
                if nSubPlt > 1
                    %Data on uniform margins
                    subplot(MM.Bn.nCvr,nSubPlt,2+3*(iC-1))
                    plot(MM.X(:,iC),YUnif,'k.')
                    hold on
                    if iC==1
                        title(sprintf('On uniform margins \n NEP = %.2f', MM.NEP(1)))
                    end
                    xlabel(MM.CvrLbl(iC))
                    ylabel(sprintf('%s-uniform scale',MM.RspLbl))
                    
                    PlotBinEdge(MM.Bn,iC);
                    
                    set(gca,'xtick',0:45:360)
                    ylim([0,1])
                    xlim([0,360])
                    
                    %Data on standard margins
                    subplot(MM.Bn.nCvr,nSubPlt,3*iC)
                    plot(MM.X(:,iC),YMrg,'k.')
                    hold on
                    if iC==1
                        title(sprintf('On %s margins \n NEP = %.2f',MM.MarginType, MM.NEP(1)))
                    end
                    xlabel(MM.CvrLbl(iC))
                    ylabel(sprintf('%s-%s scale',MM.RspLbl,MM.MarginType))
                    
                    PlotBinEdge(MM.Bn,iC)
                    
                    xlim([0,360])
                    set(gca,'xtick',0:45:360)
                end
            end
            savePics(fullfile(MM.FigureFolder,sprintf('Stg3_%s_1_DataTransform',MM.RspSavLbl)))
            
            
            %% Diagnostic plots for GP fit
            %Piecewise Sigma: fitted GP scale (black)
            figure(2)
            clf;
            
            for iC=1:MM.Bn.nCvr
                %GP Scale
                subplot(MM.Bn.nCvr,3,(iC-1)*3+1)
                PlotParameter(MM.Bn,MM.Scl,iC,'color','k','linewidth',2);
                PlotBinEdge(MM.Bn,iC);
                xlabel(MM.CvrLbl(iC))
                ylabel('\sigma')
                title(sprintf('%s: GP scale',MM.RspLbl))
            
                %Gam Alpha
             
                subplot(MM.Bn.nCvr,3,(iC-1)*3+2)
                PlotParameter(MM.Bn,MM.Alp,iC,'color','k','linewidth',2);
                PlotBinEdge(MM.Bn,iC);
                xlabel(MM.CvrLbl(iC))
                ylabel('\alpha')
                title(sprintf('%s: Gam alpha',MM.RspLbl))
                
                %Gam Beta
                subplot(MM.Bn.nCvr,3,(iC-1)*3+3)
                PlotParameter(MM.Bn,MM.Bet,iC,'color','k','linewidth',2);
                PlotBinEdge(MM.Bn,iC);
                xlabel(MM.CvrLbl(iC))
                ylabel('\beta')
                title(sprintf('%s: Gam beta ',MM.RspLbl))

            end
            
            savePics(fullfile(MM.FigureFolder,sprintf('Stg3_%s_2_Parameters',MM.RspSavLbl)))
            
            
            %Constant Xi: Fitted GP shape (black)  and true shape (green)
            figure(3);
            clf
            if verLessThan('Matlab','8.5')
                hist(MM.Shp)
                h = findobj(gca,'Type','patch');
                set(h,'FaceColor',[1 1 1]*0.5);
                set(h,'linestyle','none')
            else
                histogram(MM.Shp,'edgecolor','none','facecolor','k')
            end
            xlabel('\xi')
            
            
            title(sprintf('%s: GP shape',MM.RspLbl))
            
            savePics(fullfile(MM.FigureFolder,sprintf('Stg3_%s_3_ParametersShape',MM.RspSavLbl)))
            
            %% Lack of fit against different smoothness Lambda
            if ~all(isnan(MM.CVLackOfFit(:)))
                figure(4);
                clf;
                if MM.nB>1 && MM.CVMth==1
                    plot(MM.SmthSet,nanmedian(MM.CVLackOfFit,2),'k-','linewidth',2)
                    hold on
                    plot(MM.SmthSet,quantile(MM.CVLackOfFit,0.025,2),'k--','linewidth',2)
                    plot(MM.SmthSet,quantile(MM.CVLackOfFit,0.975,2),'k--','linewidth',2)
                else
                    plot(MM.SmthSet,MM.CVLackOfFit(:,1),'k-','linewidth',2)
                end
                axis tight
                hold on
                plot(median(MM.OptSmth)*[1,1],ylim,'r--','linewidth',2)
                ylabel('Lack of fit')
                xlabel('\lambda')
                set(gca,'xscale','log');
                title(sprintf('%s: GP cross-validation lack of fit',MM.RspLbl))
                savePics(fullfile(MM.FigureFolder,sprintf('Stg3_%s_4_CV',MM.RspSavLbl)))
            end
            
            %% Q-Q plots for each bin
            figure(5);
            clf           
            nPlt1=ceil(sqrt(MM.Bn.nBin)); %max size nPlt x nPlt
            nPlt2=ceil(MM.Bn.nBin./nPlt1);
            nQ=100;
            Q=permute(linspace(0,max(MM.Y)*1.2,nQ),[1,3,2]);
            C=MarginalModel.gpcdf(Q,MM.Shp',MM.Scl,MM.Thr);
            Q=squeeze(Q);
            C=permute(C,[3,1,2]);
            
            if MM.nB>1
                qC=quantile(C,[0.025,0.5,0.975],3);
            end
            
            
            for iB=1:MM.Bn.nBin
                subplot(nPlt2,nPlt1,iB)
                I=IExc & MM.Bn.A==iB;
                if any(I)
                    P=linspace(0,1,sum(I));
                    plot(sort(MM.Y(I)),log10(1-P),'r.')  %data
                    axis tight
                    hold on
                    if MM.nB>1
                        plot(Q,log10(1-qC(:,iB,2)),'k-') %fitted CDF
                        plot(Q,log10(1-qC(:,iB,1)),'k--') %fitted CDF
                        plot(Q,log10(1-qC(:,iB,3)),'k--') %fitted CDF
                    else
                        plot(Q,log10(1-C(:,iB)),'k-') %fitted CDF
                    end
                    xlim([0,max(MM.Y(I))])
                    
                    title(sprintf('%s: Bin %s, nExc %g',MM.RspLbl,MM.Bn.BinLbl{iB},sum(I)))
                    
                    ylabel('log(1-p)')
                    xlabel(MM.RspLbl)
                else %if no data in bin, leave plot empty
                    title(sprintf('%s: Bin %s, nExc %g',MM.RspLbl,MM.Bn.BinLbl{iB},sum(I)))
                    box on
                end
            end
            hold off;
            savePics(fullfile(MM.FigureFolder,sprintf('Stg3_%s_5_SectorGoodnessOfFit',MM.RspSavLbl)))
            
            % Overall QQ plot
            figure(6)
            clf;
            P=linspace(0,1,sum(IExc));
            plot(sort(MM.Y(IExc)),log10(1-P),'r.')
            axis tight
            hold on
            %DR This was removed on 4/11/2017 don't think its needed anymore
            %but cannot confirm!!
            %w=bsxfun(@rdivide,MM.RatExc(MM.NonEmptyBins,:),sum(MM.RatExc(MM.NonEmptyBins,:)));  %propn of exceedence data falling in each bin
            w=bsxfun(@rdivide,MM.RatExc,sum(MM.RatExc));  %propn of exceedence data falling in each bin
            COmni=sum(bsxfun(@times,shiftdim(w,-1),C),2); %sum ( prob in bin .* CDF(:,iBin) )
            
            COmni(COmni>=1)=1;
            if MM.nB>1
                qCOmni=squeeze(quantile(COmni,[0.025,0.5,0.975],3));
                plot(Q,log10(1-qCOmni(:,2)),'k-')
                hold on
                plot(Q,log10(1-qCOmni(:,[1,3])),'k--')
            else
                plot(Q,log10(1-COmni),'k-')
            end
            xlim([0,max(MM.Y)])
            title(sprintf('%s: Overall goodness of fit',MM.RspLbl))
            ylabel('log(1-p)')
            xlabel(MM.RspLbl)
            savePics(fullfile(MM.FigureFolder,sprintf('Stg3_%s_6_OverallGoodnessOfFit',MM.RspSavLbl)))
            
            %% Threshold uncertianty
            if MM.nB>1 && numel(unique(MM.NEP))>1
                figure(7);
                clf;
                if MM.nB>50
                    nNEP = 10; %number of bins for threshold diagnostic plot
                    NEPEdg=linspace(min(MM.NEP),max(MM.NEP),nNEP+1);
                    NEPBin=(NEPEdg(1:end-1)+NEPEdg(2:end))/2;
                    
                    [~,A] = histc(MM.NEP,NEPEdg);
                    A(A>nNEP)=nNEP;
                    
                    S=accumarray(A,MM.Shp,[nNEP,1],@median,NaN);
                    Slb=accumarray(A,MM.Shp,[nNEP,1],@(x)quantile(x,0.025),NaN);
                    Sub=accumarray(A,MM.Shp,[nNEP,1],@(x)quantile(x,0.975),NaN);
                    NanI=isnan(S);
                    plot(NEPBin(~NanI),S(~NanI),'r-','linewidth',2)
                    hold on
                    plot(NEPBin(~NanI),Slb(~NanI),'r--','linewidth',2)
                    plot(NEPBin(~NanI),Sub(~NanI),'r--','linewidth',2)               
                end
                plot(MM.NEP,MM.Shp,'k.','markersize',20)
                xlabel('NEP')
                ylabel('GP shape \xi')
                title(sprintf('%s: GP shape stability by threshold',MM.RspLbl))
                savePics(fullfile(MM.FigureFolder,sprintf('Stg3_%s_7_ThresholdStability',MM.RspSavLbl)))
            end
            
            %% Return Value CDF plot
            
            figure(8);
            clf;
            for iRtr=1:MM.nRtr
                subplot(MM.nRtr,1,iRtr)
                plot(MM.RVX,MM.RVPrb(1:MM.Bn.nBin,:,iRtr),'linewidth',2)
                hold on
                plot(MM.RVX,MM.RVPrb(end,:,iRtr),'k','linewidth',2)
                
                I0 = max(MM.RVPrb(:,:,1))>1e-3;
                I1 = min(MM.RVPrb(:,:,end))<1-1e-3;
                xlim([min(MM.RVX(I0)),max(MM.RVX(I1))])
                hold on
                plot(xlim,[0.5,0.5],'--k')
                plot(xlim,[exp(-1),exp(-1)],'--k')
                
                ylabel('Cumulative Probability')
                xlabel(MM.RspLbl)
                
                legend([MM.Bn.BinLbl(:);'Omni'],'location','best');
                
                title(sprintf('Return-Value CDF %g Years',MM.RtrPrd(iRtr)))
            end
            savePics(fullfile(MM.FigureFolder,sprintf('Stg3_%s_8_ReturnValueCDF',MM.RspSavLbl)))
            
        end %PlotDiagnGP
        
    end %methods
    
    methods (Access = private)
        
        function [Y,A,I]=GetBootstrapSample(MM,I)
            % [X,Y,I]=GetBootstrapSample(MM,iBt)
            % if iBt==1 use original sample
            Y = MM.Y(I);
            A = MM.Bn.A(I);
        end %GetBootstrapSample
        
        function MM=GPCrossValidation(MM,Rsd,ObsMax,AExc,iBt)
            % MM=GPCrossValidation(MM,tRsd,AExc,iBt)
            % fit piecewise constant GP to threshold exceedances
            % INPUTS:
            % - Rsd
            % - AExc
            % - iBt
            % OUTPUT:
            
            %% Constant Starting Solution (method of moments)
            xbar=mean(Rsd);
            s2=var(Rsd);
            xi0 =0;
            sigma0 = .5 .* xbar .* (xbar.^2 ./ s2 + 1)*2;
            p0=[xi0;log(ones(MM.Bn.nBin,1)*sigma0)];
            
            nExc=length(Rsd);
            opts=optimset('display','off');
            
            %% Flag for do Cross Validation
            if (MM.CVMth==0 &&  iBt>1) || MM.nSmth==1
                %% Cross validation off
                if MM.nSmth==1  %only one smoothness
                    MM.OptSmth(iBt)=MM.SmthSet;
                else  %use first bootstrap instead of cross validating every time
                    MM.OptSmth(iBt)=MM.OptSmth(1);
                end
            else
                fprintf('Starting Cross Validation to estimate GP Scale smoothness:\n');
                LackOfFit=NaN(MM.nCV,MM.nSmth); %lack of fit
                ICV=randi(MM.nCV,nExc,1); %split nExc observations into nCV cross validation groups
                
                for iCV=1:MM.nCV  %cross validation loop
                    fprintf('.')
                    %Pull out fit set
                    ObsMaxFit=  ObsMax(ICV~=iCV);
                    RsdFit=Rsd(ICV~=iCV);
                    AExcFit=AExc(ICV~=iCV);
                    %Pull out prediction set
                    ObsMaxPrd=  ObsMax(ICV==iCV);
                    RsdPrd=Rsd(ICV==iCV);
                    AExcPrd=AExc(ICV==iCV);
                    for iLam=1:MM.nSmth %Smoothness penalties
                        InitialNLL=MarginalModel.GPLikeNonStationary(p0,RsdFit,AExcFit,ObsMaxFit,MM.SmthSet(iLam));
                        if isinf(InitialNLL)
                            error('bad starting value')
                        end
                        PrmHat=fminsearch(@(p)MarginalModel.GPLikeNonStationary(p,RsdFit,AExcFit,ObsMaxFit,MM.SmthSet(iLam)),p0,opts);
                        LackOfFit(iCV,iLam)=MarginalModel.GPLikeNonStationary(PrmHat,RsdPrd,AExcPrd,ObsMaxPrd);
                    end %SigPen
                end %CV
                fprintf('\n')
                MM.CVLackOfFit(:,iBt)=sum(LackOfFit,1)';
                [~,tI]=min(MM.CVLackOfFit(:,iBt));
                MM.OptSmth(iBt)=MM.SmthSet(tI)';
            end
            
            %% Fit Model using optimal smoothness            
            InitialNLL=MarginalModel.GPLikeNonStationary(p0,Rsd,AExc,ObsMax,MM.OptSmth(iBt));
            if isinf(InitialNLL)
                error('bad starting value')
            end
            %fit model using opt smoothness
            PrmHat = fminsearch(@(p)MarginalModel.GPLikeNonStationary(p,Rsd,AExc,ObsMax,MM.OptSmth(iBt)),p0,opts);
            MM.Shp(iBt)=PrmHat(1);
            MM.Scl(:,iBt)=exp(PrmHat(2:end));
                          
        end %GPCrossValidation
                                      
    end %methods private
    
    methods (Static)
        function PLOGL=GPLikeNonStationary(PARAMS,Dat,BinAlc,ObsMax,SigPen)
            %PLOGL=gplikenonstationary(PARAMS,Dat,BinAllocation,L)
            %Piecewise constant gplike
            %INPUTS:
            % - (nBins+1) x 1 PARAMS, piecewise constant marginal model
            % parameters(= [xi, sig_1 ,..., sig_nBins])
            % - nDat x 1 Dat, GP data (= threshold exceedences - threshold)
            % - nDat x 1 BinAlc, bin allocation index on data
            % - 1 x 1 SigPen, penalty imposing smoothness in GP scale across covariate
            % - max seen in  each bin used for global upper end point constraint
            %OUTPUT:
            % - 1 x 1 PLOGL, Penalised negative log likelihood
            
            if size(Dat,1) ~= size(BinAlc,1)
                error('Input dimension mismatch: dmn Dat =?= dmn BinAllocation')
            end
            
            if nargin<=4
                SigPen=0;  %unpenalised
            end          
            Xi=PARAMS(1);
            Sig=exp(PARAMS(2:end));
            nBins =numel(Sig);
            
            NLOGL=0;
            for iBin=1:nBins %Additive likelihood over covariate bins
                I=BinAlc==iBin;
                if any(I)
                    tGPLike = gplike([Xi,Sig(iBin)],Dat(I));
                    if isinf(tGPLike) || any(Xi<-0.5) || (Xi<0 && any(ObsMax(I)>-Sig(iBin)./Xi))
                        NLOGL=Inf;
                    else
                        NLOGL=NLOGL+tGPLike;
                    end
                end
            end
            
            PLOGL=NLOGL+SigPen.*sum((Sig-mean(Sig)).^2);
            
        end %likelihood
                                    
        function X=gpinv(P,Xi,Sgm,Thr)
            %    X=gpinv(P,Xi,Sgm,Thr) returns the inverse of generalized Pareto (GP)
            %     cdf with tail index (shape) parameter Xi, scale parameter Sgm,
            %     and threshold (location) parameter Thr, evaluated at the values in X.
            %     The size of P is the common size of the input arguments.
            
            t1=bsxfun(@rdivide,bsxfun(@power,1-P,-Xi)-1,Xi);
            
            %%Gumbel case
            I=abs(Xi)<1e-5;
            if any(I(:))
%                 t1(:,I)=-log(1-P(:,I));
                t1(I)=-log(1-P(I));
            end
            
            X=bsxfun(@plus,bsxfun(@times,Sgm,t1), Thr);
        end %gpinv
        
        function F=gppdf(X,Xi,Sgm,Thr)
            %     F = gppdf(X,Xi,Sgm,Thr) returns the pdf of the generalized Pareto (GP)
            %     distribution with tail index (shape) parameter Xi, scale parameter Sgm,
            %     and threshold (location) parameter THETA, evaluated at the values in X.
            %     The size of P is the common size of the input arguments.
            
            Z=bsxfun(@rdivide,bsxfun(@minus,X,Thr),Sgm);
            t1=1+bsxfun(@times,Xi,Z);
            t1(t1<=0)=0;  %deal with upper end point of distribution
            
            F=bsxfun(@times,1./Sgm,bsxfun(@power,t1,-1./Xi-1));
            
            %% Gumbel case
            I=abs(Xi)<1e-5;
            if any(I(:))
                Xi=bsxfun(@times,Xi,ones(size(Z)));
                I=abs(Xi)<1e-5;
                Sgm=bsxfun(@times,Sgm,ones(size(Z)));
                
                F(I)=bsxfun(@times,1./Sgm(I),exp(-Z(I)));
            end
            
            F(bsxfun(@lt,X,Thr))=0; %density of zero below the threshold
            F(Z<0)=0;
            
        end %gppdf
        function P=gpcdf(X,Xi,Sgm,Thr)
            %     P = CumulativeDensityFunction(obj,X) returns the cdf of the generalized Pareto (GP)
            %     distribution with tail index (shape) parameter Xi, scale parameter Sgm,
            %     and threshold (location) parameter THEThrTA, evaluated at the values in X.
            %     The size of P is the common size of the input arguments.
            Z=bsxfun(@rdivide,bsxfun(@minus,X,Thr),Sgm);
            t1=1+bsxfun(@times,Xi,Z);
            t1(t1<=0)=0;  %deal with upper end point of distribtion
            P=1-bsxfun(@power,t1,-1./Xi);
            %% Gumbel case
            I=abs(Xi)<1e-5;
            if any(I(:))
                Xi=bsxfun(@times,Xi,ones(size(Z)));
                I=abs(Xi)<1e-5;
                
                P(I)=1-exp(-Z(I));
            end
            P(bsxfun(@le,X,Thr))=NaN; %density of zero below the threshold
            P(Z<=0)=NaN;
            
        end %gpcdf
        
        function F=gampdf(X,Alp,Bet,GmmLct)
            %     F = gampdf(X,Alp,Bet,GmmLct) returns the pdf of the gamma distribution using orthog0nal
            %     parameterisation
            %density
            Z=bsxfun(@minus,X,GmmLct);
            Z(Z<0)=0;
            
            F = bsxfun(@times,(Bet./Alp).^(-Alp)./gamma(Alp),bsxfun(@power,Z,(Alp-1))).*...
                exp(-bsxfun(@times,Alp./Bet,Z));
        end %gampdf
        
        function P=gamcdf(X,Alp,Bet,GmmLct)
            %     P = gpcdf(X,Alp,Bet,GmmLct) returns the cdf of the gamma distribution using orthognal
            %     parameterisation
            Z=bsxfun(@minus,X,GmmLct);
            Z(Z<0)=0;
            Z = bsxfun(@times,Alp./Bet,Z);
            P = bsxfun(@gammainc,Z, Alp);
        end %gamcdf
        
        function X=gaminv(P,Alp,Bet,GmmLct)
            %     X = gaminv(P,Alp,Bet,GmmLct)returns the inverse the gamma distribution using orthognal
            %     parameterisation
            
            %icdf
            
            P=bsxfun(@times,P,ones(size(Alp)));
            Alp=bsxfun(@times,Alp,ones(size(P)));
            
            I=~isnan(P);
            q=NaN(size(P));
            
            q(I) = gammaincinv(P(I),Alp(I));
            X = bsxfun(@plus,bsxfun(@times,q,Bet./Alp),GmmLct);
        end %gaminv
                                                     
        function P=gamgpcdf(X,Xi,Sgm,Thr,Alp,Bet,GmmLct,Tau)
            %     P=gamgpcdf(X,Xi,Sgm,Thr,Alp,Bet,GmmLct,Tau) returns the cdf of the gamma-gp distribution
            
            %threshold
%             Thr=MarginalModel.gaminv(Tau,Alp,Bet,GmmLct);            
            IBlw=bsxfun(@le,X,Thr);
            %gamma part
            P1=MarginalModel.gamcdf(X,Alp,Bet,GmmLct);
            % mjj 08/11/18: I think tau scaling for Gamma was missing:
            % added it.
            %gp part
            P2=bsxfun(@plus,bsxfun(@times,MarginalModel.gpcdf(X,Xi,Sgm,Thr),(1-Tau)),Tau);
            %combine                        
            P=NaN(size(P1));
            P(IBlw)=P1(IBlw);
            P(~IBlw)=P2(~IBlw);
            
        end %gamgpcdf
        
        
        
        function f=gamgppdf(X,Xi,Sgm,Thr,Alp,Bet,GmmLct,Tau)
            % f=gamgpcdf(X,Xi,Sgm,Alp,Bet,GmmLct,Tau) returns the pdf of the gamma-gp distribution
           
            %threshold
%             Thr=MarginalModel.gaminv(Tau,Alp,Bet,GmmLct);
            IBlw=bsxfun(@le,X,Thr);
            %gamma part
            f1=MarginalModel.gampdf(X,Alp,Bet,GmmLct);            
            %gp part
            f2=bsxfun(@times,MarginalModel.gppdf(X,Xi,Sgm,Thr),(1-Tau));
            %combine                     
            f=NaN(size(f1));
            f(IBlw)=f1(IBlw);
            f(~IBlw)=f2(~IBlw);
            
        end %gamgpcdf
        
        
        function X=gamgpinv(P,Xi,Sgm,Thr,Alp,Bet,GmmLct,Tau)
            %X=gamgpinv(P,Xi,Sgm,Thr,Alp,Bet,GmmLct,Tau) returns the inv of the gamma-gp distribution
                              
            %gamma part
            X1=MarginalModel.gaminv(P,Alp,Bet,GmmLct);
            
            %gp part (adjust probabilities)
            PGP=bsxfun(@rdivide,(P-Tau),(1-Tau));
            PGP(PGP<0)=0;
            X2=MarginalModel.gpinv(PGP,Xi,Sgm,Thr);                                    
                        
            X=NaN(size(X1));
            X(IBlw)=X1(IBlw);
            X(~IBlw)=X2(~IBlw);
        end %gamgpcdf
                
        
    end %methods
end %class
