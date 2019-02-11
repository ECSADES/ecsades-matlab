classdef Contour
    %class of functions for fitting extreme value contours
    %Contour methods currently implemented:
    % - 'Exc' constant exceedence contour, requires long simulation
    %under model to get reliable results;
    % - 'EmpDns' uses an empirical density on Original Scale to define contour requires long simulation
    %under model to get reliable results;
    % - 'HTDns' c   onstant density contour of StandardMargins Margins, density form of H & T to get contour
    % - 'RadQnt' radial quantile estimate

    properties
        nPon=100;  %1 x 1, how many points from which to draw contour
        XRng %nPon x (nBin+1) x Cnt.nLvln, conditioned values for contour
        YVal; %nMth x 1, (cell array) contour value in conditioning variable
        XY;
        Mth;  %nMth x 1, (cell array) of contour methods used
        nMth; %1 x 1, number of contouring methods
        nBin; %1 x 1, number of covariate bins
        nLvl; %1 x 1, number of contour levels chosen
        nAsc; %1 x 1 number of associated variables
        Sml %structure importance sampled simulation under the model
    end
    
    properties
        DnsLvl  %(nBin+1) x nLvl, level used by contour plot
        HTDnsLvl  %(nBin+1) x nLvl, level used by contour plot
        DnsEdg   %(nBin+1) x 2 (cell array), grid edges used by contour
        HTDnsEdg  %(nBin+1) x 2 (cell array), grid edges used by contour
        RadQntX  %nAng x (nBin+1) x nLvl, X values used in radial quantile
<<<<<<< HEAD
        HusAng  
        PltOn=true; %flag to switch on explanatory diagrams
        LvlOrg; %(nBin+1) x nLvl x 2, contour level on orginal scale on conditioned variable
        LvlStnMrg; %(nBin+1) x nLvl x 2, contour level on StandardMargins scale ion conditioned variable
=======
        HusAng
        PltOn=false; %flag to switch on explanatory diagrams
        LvlOrg; %(nBin+1) x nLvl x 2, contour level on orginal scale on conditioned variable       
>>>>>>> CEVAContours
    end
    
    methods
        function Cnt=Contour(HT,Mrg,Mth)
            %Cnt=CEVA_Contour(HT,Mrg,Lvl,Mth)
            %INPUTS:
            % - HT Hefferenan and Tawn class output
            % - Mrg Marginal model class output
            % - RtrPrd nLvl x 1 probability level (1-prb)
            % - Mth contour method to be used can be more than one (cell array)
            % - Typ  (either Exc  YLvl=median(Y | X>=h)
            %                 Eq  YLvl=median(Y | X=h)
            %OUTPUT:
            % - Cnt Contours class containing settings/options for contours
            % to be calculated
            % - Sml data simulated from fitted H&T model
            
            %% Check inputs
            if nargin==0
                return
            end
            if ~isa(Mrg,'MarginalModel')
                error('Mrg should be a nDmn x 1 Marginal Model')
            end
            
            if ~isa(HT,'HeffernanTawn')
                error('HT should be a HeffernanTawn class')
            end
            Mth=unique(Mth); %check if method in more than once.
            
            %% Pre-allocation
            Cnt.nLvl=numel(Mrg(1).RtrPrd);
            Cnt.nMth=numel(Mth);
            Cnt.YVal=cell(Cnt.nMth,1);
            Cnt.Mth=Mth;
            Cnt.nBin=Mrg(1).Bn.nBin;
            Cnt.nAsc=HT.nDmn-1;
            Cnt.nPon=100;  %number of points at which contour is defined
            
            %% simulate
            Cnt.Sml=SimulateIS(HT,Mrg,1e6);
            
            
            
        end %CEVA_Contour constructor
        
        function Cnt = makeContours(Cnt,Mrg,Mth,A)
            % makeContour: construct the contours (of required)
            
            if nargin<=3 %PPC case
                A=Cnt.Sml.A;
                Cnt=GetLockPoint(Cnt,HT,Mrg,A);
            end
            % unique bins
            uqA = unique(A(A>0));
            % # unique bins
            n_uqA = numel(uqA);
            % # bins + omni
            if n_uqA>1
                nBinsplsOmni = n_uqA+1;
            else
                nBinsplsOmni = 1;
            end
            
            %% Get Contour levels from Marginal and Conditional analysis
            Cnt.XRng=NaN(Cnt.nPon,nBinsplsOmni, Cnt.nLvl);  %X range for contours
            
            %% Compute lock ponit for each bin
            
            %Initialise XVal
            for iBin=1:nBinsplsOmni
                for iLvl=1:Cnt.nLvl
                    Cnt.XRng(:,iBin,iLvl)=linspace(min(Mrg(1).Y),Cnt.LvlOrg(iBin,iLvl,1),Cnt.nPon);
                end
            end
            %Initialise YVal
            Cnt.XY = cell(Cnt.nMth,1);
            
            
            %% Constant Exceedence contour (invariant to margins);
            if any(strcmp(Mth,'Exc'))
                Cnt=Exceedence(Cnt,A);
            end
            %% Huseby Contour Method  (Original Scale);
            if any(strcmp(Mth,'Hus'))
                Cnt=HusebyContour(Cnt,A);
            end
            %% Constant Hefferenan and Tawn Density contour  (must be on StandardMargins margins);
            if any(strcmp(Mth,'HTDns'))
                Cnt=HTDensity(Cnt,A);
            end
        end
        
        function Cnt=GetLockPoint(Cnt,HT,Mrg,~)
            % unique bins
            uqA = unique(A(A>0));
            % # unique bins
            n_uqA = numel(uqA);
            % # bins + omni
            if n_uqA>1
                nBinsplsOmni = n_uqA+1;
            else
                nBinsplsOmni = 1;
            end
            
            %In PPC code contours just computed at the bin so nothing special needed
            Cnt.LvlOrg=NaN(nBinsplsOmni, Cnt.nLvl,HT.nDmn); %contour level on standard scale
            Cnt.LvlOrg(:,:,1)= Mrg(1).RVMed;  %Mrg(1).RVMed; %range in X for contour on oroginal scale
            Cnt.LvlOrg(:,:,2:end)=permute(nanmedian(HT.RV.Y,3),[1,4,3,2]); %contour level on standard scale                        
        end
        
        function Cnt=Exceedence(Cnt,A)
            %Cnt=Exceedence(Cnt,Sml)
            %% Constant Exceedence
            % INPUTS:
            % - Cnt Contours class
            % - Sml data simulated from fitted H&T model
            % OUTPUT:
            % - Cnt.XVal,Cnt.YVal for iMth = 'Exc' data to draw contours
                        
            if nargin<=1
                A=Cnt.Sml.A;
            end
            
            IMth=strcmp(Cnt.Mth,'Exc');
            
            % # points to use for Y CDF grid
            nG = 1000;
            
             % unique bins
            uqA = unique(A(A>0));
            % # unique bins
            n_uqA = numel(uqA);
            % # bins + omni
            if n_uqA>1
                nBinsplsOmni = n_uqA+1;
            else
                nBinsplsOmni = 1;
            end
            
            % initialise
            Cnt.XY{IMth}=NaN(2.*Cnt.nPon,2,nBinsplsOmni,Cnt.nLvl,Cnt.nAsc);
            
                        
            for iQ=1:Cnt.nLvl %contour levels
                for iB=1:nBinsplsOmni
                    % pre-allocate
                    YDw = nan(Cnt.nPon,Cnt.nAsc);
                    YUp = nan(Cnt.nPon,Cnt.nAsc);
                    % loop over associated variables
                    for iAsc = 1:Cnt.nAsc
                        %lock point not defined for this bin/return level combo
                        if any(isnan(Cnt.LvlOrg(iB,iQ,:)))
                           continue 
                        end
                        if iB==nBinsplsOmni %omni case
                            SmlX = Cnt.Sml.Org(:,1);
                            fogX = Cnt.Sml.fog(:,1);
                            SmlY = Cnt.Sml.Org(:,iAsc+1);
                            fogY = Cnt.Sml.fog(:,iAsc+1);
                        else     %pick out simulations from each bin in turn
                            SmlX = Cnt.Sml.Org(A==iB,1);
                            fogX = Cnt.Sml.fog(A==iB,1);
                            SmlY = Cnt.Sml.Org(A==iB,iAsc+1);
                            fogY = Cnt.Sml.fog(A==iB,iAsc+1);
                        end
                        % compute overall f(.)/g(.)
                        fog = fogX.*fogY;
                        
                        % compute P(X>x_{loc},Y>y_{loc})
                        I_gtLoc_Up = (SmlX>Cnt.LvlOrg(iB,iQ,1))&...
                            (SmlY>Cnt.LvlOrg(iB,iQ,iAsc+1));
                        I_gtLoc_Dwn = (SmlX>Cnt.LvlOrg(iB,iQ,1))&...
                            (SmlY<=Cnt.LvlOrg(iB,iQ,iAsc+1));
                        p_Loc_Up = sum(fog(I_gtLoc_Up))./sum(fog);
                        p_Loc_Dwn = sum(fog(I_gtLoc_Dwn))./sum(fog);
                        
                        % grid in main and associated
                        x_grd = Cnt.XRng(:,iB,iQ);
                        y_grd = linspace(min(SmlY),max((SmlY)),nG)';
                        
                        % P(X>x)
                        P_x = sum((x_grd>SmlX').*fog',2)./sum(fog);
                        % P(Y>y|X>x)
                        P_ygx = nan(nG,Cnt.nPon);
                        for iG = 1:Cnt.nPon
                            I_x = SmlX(:,1)>x_grd(iG);
                            [srtY,srtI] = sort(SmlY(I_x));
                            tfog = fog(I_x);
                            srtfog = tfog(srtI);
                            P_ygx(:,iG) = ReturnLevel.Dissipation.WeightedCDF(y_grd,srtY,srtfog);
%                             tP_ygx = sum((y_grd>SmlY(I_x)').*fog(I_x)',2)./sum(fog(I_x));
                        end
                        
                        % find y_Dw closest to this probability level
                        [~,I_min_dw] = min(abs((P_ygx).*(1-P_x')-p_Loc_Dwn),[],1);
                        YDw(:,iAsc) = y_grd(I_min_dw);
                        [~,I_min_up] = min(abs((1-P_ygx).*(1-P_x')-p_Loc_Up),[],1);
                        YUp(:,iAsc) = y_grd(I_min_up);
                    end
                    
                    %store the 2 parts of the contour (up and down) together
                    for iAsc=1:Cnt.nAsc
                        Cnt.XY{IMth}(:,2,iB,iQ,iAsc) = [permute(YDw(:,iAsc),[1,3,4,5,2]);flipud(permute(YUp(:,iAsc),[1,3,4,5,2]))]; %TO FINISH: last dimention is Associated
                        Cnt.XY{IMth}(:,1,iB,iQ,iAsc) = [Cnt.XRng(:,iB,iQ);flipud(Cnt.XRng(:,iB,iQ))];  %what should X be??
                    end
                end
            end
            
            
            %% diagram
            if Cnt.PltOn
                clf;
                
                rectangle('position',[Cnt.LvlOrg(iB,iQ,1),Cnt.LvlOrg(iB,iQ,2),20,20],'curvature',0,'edgecolor','none','facecolor','b');
                hold on;
                rectangle('position',[Cnt.LvlOrg(iB,iQ,1),0,20,Cnt.LvlOrg(iB,iQ,2)],'curvature',0,'edgecolor','none','facecolor','r')
                plot(Cnt.Sml.Org(:,1),Cnt.Sml.Org(:,2),'k.','markersize',5)
                hold on
                plot(Cnt.XY{1}(1:Cnt.nPon,1,end,iQ,1),Cnt.XY{1}(1:Cnt.nPon,2,end,iQ,1),'r-','linewidth',2)
                plot(Cnt.XY{1}((Cnt.nPon+1):end,1,end,iQ,1),Cnt.XY{1}((Cnt.nPon+1):end,2,end,iQ,1),'b-','linewidth',2)
                
                plot(Cnt.LvlOrg(iB,iQ,1),Cnt.LvlOrg(iB,iQ,2),'go','markersize',10,'linewidth',2)
                axis tight
                plot([Cnt.LvlOrg(iB,iQ,1),Cnt.LvlOrg(iB,iQ,1)*2],Cnt.LvlOrg(iB,iQ,2)*[1,1],'k--','linewidth',2)
                plot(Cnt.LvlOrg(iB,iQ,1)*[1,1],ylim,'k--','linewidth',2)
                
                xlim([min(Cnt.Sml.Org(:,1)),max(Cnt.Sml.Org(:,1))]);
                ylim([min(Cnt.Sml.Org(:,2)),max(Cnt.Sml.Org(:,2))]);
                
                savePics('ExcDiagram')
            end
            
        end %Exceedence
        
        function Cnt=HTDensity(Cnt,A)
            % Cnt=HTDensity(Cnt,Mrg,HT)
            %% Constant density
            % INPUTS:
            % - Cnt Contours class
            % - Mrg Marginal model class output
            % - HT Hefferenan and Tawn class output
            % OUTPUT:
            % - Cnt.XVal,Cnt.YVal for iMth = 'HTDns' data to draw contours                        
            if nargin<=1
                A=Cnt.Sml.A;
            end
            %% Parse inputs
            % method
            IMth=strcmp(Cnt.Mth,'HTDns');          
             % unique bins
            uqA = unique(A(A>0));
            % # unique bins
            n_uqA = numel(uqA);
            % # bins + omni
            if n_uqA>1
                nBinsplsOmni = n_uqA+1;
            else
                nBinsplsOmni = 1;
            end
            
            nG=[200,200]; %number of grid points for binned density                    
            
            % initialise            
            Cnt.XY{IMth}=cell(nBinsplsOmni,Cnt.nLvl,Cnt.nAsc);                                    
           
            % loop over associated variables
            for iAsc = 1:Cnt.nAsc
                fog=prod(Cnt.Sml.fog(:,[1,iAsc+1]),2);
                % loop over bins
                for iBin = 1:nBinsplsOmni  
                   
                    %% use contour function to get iso-line.
                    %setup gridded to bin density into
                    Lmt=squeeze(Cnt.LvlOrg(iBin,end,:));
                    if iBin>n_uqA %omni case!!
                        Edgx=linspace(min(Cnt.Sml.Org(:,1)),Lmt(1)*1.2,nG(1)+1);
                        Edgy=linspace(min(Cnt.Sml.Org(:,2)),max(Cnt.Sml.Org(:,2)),nG(2)+1);
                    else
                        Edgx=linspace(min(Cnt.Sml.Org(A==iBin,1)),Lmt(1)*1.2,nG(1)+1);
                        Edgy=linspace(min(Cnt.Sml.Org(A==iBin,2)),max(Cnt.Sml.Org(A==iBin,2)),nG(2)+1);                        
                    end
                    Grdx=(Edgx(1:end-1)+Edgx(2:end))/2;
                    Grdy=(Edgy(1:end-1)+Edgy(2:end))/2;
                    
                    [GX,GY]=ndgrid(Grdx,Grdy);
                    G=[GX(:),GY(:)];
                    
                    Ax=discretize(Cnt.Sml.Org(:,1),Grdx);
                    Ay=discretize(Cnt.Sml.Org(:,2),Grdy);
                    GrdInd=sub2ind(nG,Ax,Ay);
                    
                     %% find gridded density estimate for current bin
                    if iBin>n_uqA %omni case!!
                        I=~isnan(GrdInd); %in current bin and not a nan
                    else
                        I=~isnan(GrdInd) & A==iBin; %in current bin and not a nan
                    end
                    %% which bin is lock point in?
                    Lck=shiftdim(Cnt.LvlOrg(iBin,:,:),1);
                    Lx=discretize(Lck(:,1),Grdx);
                    Ly=discretize(Lck(:,2),Grdy);
                    L=sub2ind(nG,Lx,Ly);
                                                          
                    %f=reshape(accumarray(GrdInd(I),fog(I),[prod(nG),1],@mean,0),nG);
                    %%smooth points onto grid
                    tf=fog(I);
                    [IDX, D]= knnsearch(Cnt.Sml.Org(I,:),G,'K',2000,'Distance','seuclidean');
                    Rho=exp(-D);
                    Rho=Rho./sum(Rho,2);
                    f=reshape(sum(Rho.*tf(IDX),2),nG);
                                
                    %% find lock poing density value for each return level
                    IL=find(~isnan(L));
                    Lvl=f(L(IL));
                    if numel(Lvl)==1
                        Lvl=Lvl.*[1,1];
                    end
                   % contour(Grdx,Grdy,f',Lvl);
                    % compute contour using low level function
                    C =contourc(Grdx,Grdy,f',Lvl);
                    
                    [I,J]=ismember(C(1,:),Lvl);
                    Ind=J(I); %which contour does each segment belong to
                    Count=cumsum(C(2,I)+1);
                                                
                    c=0;
                    for i=1:numel(Ind) %loop over line segments
                        iLvl=IL(Ind(i)); %contour level;
                        tI=c+2:Count(i);
                        c=Count(i);
                        %assign line segement to right level;
                        if isempty(Cnt.XY{IMth}{iBin,iLvl,iAsc})
                            Cnt.XY{IMth}{iBin,iLvl,iAsc}=C(:,tI);
                        else
                            Cnt.XY{IMth}{iBin,iLvl,iAsc}=[Cnt.XY{IMth}{iBin,iLvl,iAsc},[NaN;NaN],C(:,tI)];
                        end                        
                    end                    
                end
            end%loop over bins                                                     
            
        end %HTDensity
        
        function Cnt=HusebyContour(Cnt,A)
            %O=HusebyContour(x, varargin) Drawing Huseby contours from a bivariate sample
            %
            %INPUT
            %Cnt   %contour structure
            %Mrg  %marginal model (only used for return period)
            %Sml  %set of simulations to draw contour through
            %A op[tional redefined bins!!
            
            if nargin<=1
                A=Cnt.Sml.A;
            end
            
            %% Parse inputs
            % method
            IMth=strcmp(Cnt.Mth,'Hus');                        
            
            % unique bins
            uqA = unique(A(A>0));
            % # unique bins
            n_uqA = numel(uqA);
            % # bins + omni
            if n_uqA>1
                nBinsplsOmni = n_uqA+1;
            else
                nBinsplsOmni = 1;
            end
            % initialise
            Cnt.XY{IMth}=NaN(Cnt.nPon,2,nBinsplsOmni,Cnt.nLvl,Cnt.nAsc);
            
            % set of angles for contour eval.
            angles = linspace(-pi,pi,Cnt.nPon); %set of angles
            
            % loop over associated variables
            for iAsc = 1:Cnt.nAsc
                % loop over bins
                for iBin = 1:nBinsplsOmni
                    
                    J=[1,iAsc+1];
                    if (iBin==nBinsplsOmni)&&(nBinsplsOmni>1)
                        x=Cnt.Sml.Org(:,J);
                        fog = prod(Cnt.Sml.fog(:,J),2);
                    else
                        x=Cnt.Sml.Org(A==iBin,J);
                        fog = prod(Cnt.Sml.fog(A==uqA(iBin),J),2);
                    end
                    
                    % E[X]
                    E_X = (sum(x.*fog,1)./sum(fog))';
                    % Cov[X]
                    Cov_X = (x-E_X')'*(fog.*(x-E_X'))./sum(fog);
                    % chol(Cov[X])
                    cholCov_X = cholcov(Cov_X);
                    
                    % Z = \Sig^{-1/2}*(x-\mu)                     
                    z = (x-E_X')/cholCov_X;
                    
                    % compute prob. associated with lock point
                    IGd=find(~isnan(Cnt.LvlOrg(iBin,:,1)) & Cnt.LvlOrg(iBin,:,1)>E_X(1));
                    
                    p_c = sum(bsxfun(@gt,x(:,1),Cnt.LvlOrg(iBin,IGd,1)).*fog,1)./sum(fog,1);
                    
                    %k = floor(size(x,1)*(1-pf)); %number of points above probiblity level?
                    
                    % set up knots for the IS density
                    Z_knt = linspace(min(z(:))-1e-2,max(z(:))+1e-2,20000)';

                    
                    %% Calculate the function C_theta
                    C_theta=NaN(Cnt.nPon, numel(p_c)); %preallocate
                    for iAng=1:length(angles)
                        % X_p = c*X
                        Z_p = z(:,1)*cos(angles(iAng)) + z(:,2)*sin(angles(iAng));
                        % Compute importance sampled density
                        %P_Zp = sum((Z_knt<=Z_p').*fog',2)./sum(fog);
                        [srtZ,srtI] = sort(Z_p);
                        srtfog = fog(srtI);
                        P_Zp = 1-ReturnLevel.Dissipation.WeightedCDF(Z_knt,srtZ,srtfog);
%                         old_P_Zp = 1-sum((Z_knt>Z_p').*fog',2)./sum(fog);
                        % find closest point to req. probability
                        [~,I_min] = min(abs(P_Zp-p_c),[],1);
                        % get C(\theta)
                        C_theta(iAng,:) = Z_knt(I_min);
                    end
                    C_theta=movmean(C_theta,3); %deal with fase moving corners better!!
                    
                    %% Calculate the intersection points
                    
                    for iRtr=1:numel(p_c) %loop over return periods
                        for iAng =1:Cnt.nPon-1  %loop over poitns
                            % compute contour points for standardised
                            % variables
                            tXY = nan(1,2);
                            tXY(1) = (sin(angles(iAng+1))*C_theta(iAng, iRtr)-sin(angles(iAng))*C_theta(iAng+1,iRtr))./...
                                    (sin(angles(iAng+1))*cos(angles(iAng))-sin(angles(iAng))*cos(angles(iAng+1)));
                            tXY(2) = (-cos(angles(iAng+1))*C_theta(iAng, iRtr)+cos(angles(iAng))*C_theta(iAng+1, iRtr))./...
                                    (sin(angles(iAng+1))*cos(angles(iAng))-sin(angles(iAng))*cos(angles(iAng+1)));
                            % transform back to original scale
                            Cnt.XY{IMth}(iAng,:,iBin,IGd(iRtr),iAsc) = tXY*cholCov_X + E_X';
                            %Cnt.XY{IMth}(iAng,:,iBin,iRtr,iAsc) = tXY + E_X';
                        end
                    end
%                     Cnt.XY{IMth}(end,:,iBin,:,iAsc) = Cnt.XY{IMth}(1,:,iBin,:,iAsc); %final point is start of contour again
                    
                end %loop over bins
            end
            
        end % HusebyContour
        
        function Plot(Cnt,Mrg)
            
            %% plot omni contour
            figure(1);
            clf;
            for iAsc=1:Cnt.nAsc
                subplot(1,Cnt.nAsc,iAsc);
                
                hold on
                plot(Mrg(1).Y,Mrg(iAsc+1).Y,'k.','markersize',10,'handlevisibility','off')
                hold on
                
                xlabel(sprintf('%s: Conditioned variable',Mrg(1).RspLbl))
                ylabel(sprintf('%s: Conditioning variable',Mrg(iAsc+1).RspLbl))
                C=lines(Cnt.nMth);
                
                for iQ=1:Cnt.nLvl %return period
                    for iCnt=1:Cnt.nMth  %method                        
                        %General plotting code:
                        if strcmp(Cnt.Mth(iCnt),'HTDns')
                            if any(Cnt.XY{iCnt}{end,iQ,iAsc}(:))
                                plot(Cnt.XY{iCnt}{end,iQ,iAsc}(1,:),Cnt.XY{iCnt}{end,iQ,iAsc}(2,:),'-','color',C(iCnt,:),'linewidth',2);
                                   hold on
                            end
                        else
                            plot(Cnt.XY{iCnt}(:,1,end,iQ,iAsc),Cnt.XY{iCnt}(:,2,end,iQ,iAsc),'-','color',C(iCnt,:),'linewidth',2);
                            hold on
                        end
                    end
                    plot(Cnt.LvlOrg(end,:,1),Cnt.LvlOrg(end,:,iAsc+1),'go','markersize',8, 'LineWidth',2);
                end
                
                title(sprintf('%s|%s: Omni directional contour',Mrg(iAsc+1).RspLbl,Mrg(1).RspLbl))
                
                legend(Cnt.Mth,'location','best')
            end
            
            savePics('Figures/Stg5_Contour_1_Omni');
            %% plot contour by bin
            for iAsc=1:Cnt.nAsc
                figure(iAsc+1);
                clf;
                nPlt1=ceil(sqrt(Cnt.nBin)); %max size nPlt x nPlt
                nPlt2=ceil(Cnt.nBin./nPlt1);
                for iBin=1:Cnt.nBin
                    
                    subplot(nPlt2,nPlt1,iBin)

                    hold on
                    J=Mrg(1).Bn.A==iBin;
                    plot(Mrg(1).Y(J),Mrg(iAsc+1).Y(J),'k.','markersize',10,'handlevisibility','off')
                    xlabel(sprintf('%s: Conditioned',Mrg(1).RspLbl))
                    ylabel(sprintf('%s: Conditioning',Mrg(iAsc+1).RspLbl))
                    C=lines(Cnt.nMth);
                    
                    for iQ=1:Cnt.nLvl %return period
                        for iCnt=1:Cnt.nMth %method
                            if strcmp(Cnt.Mth(iCnt),'HTDns')
                                if any(Cnt.XY{iCnt}{iBin,iQ,iAsc}(:))
                                    plot(Cnt.XY{iCnt}{iBin,iQ,iAsc}(1,:),Cnt.XY{iCnt}{iBin,iQ,iAsc}(2,:),'-','color',C(iCnt,:),'linewidth',2);
                                       hold on
                                end
                            else
                                plot(Cnt.XY{iCnt}(:,1,iBin,iQ,iAsc),Cnt.XY{iCnt}(1:end,2,iBin,iQ,iAsc),'color',C(iCnt,:),'linewidth',2);
                                hold on
                            end
                        end
                        
                    end
                    plot(Cnt.LvlOrg(iBin,:,1),Cnt.LvlOrg(iBin,:,iAsc+1),'go','markersize',8, 'LineWidth',2);  %BUG: NaN in Cnt.LvlOrg(iBin,:,iAsc+1)
                    
                    if iBin==1
                        title(sprintf('Bin %s',Mrg(1).Bn.BinLbl{iBin}))
                        legend(Cnt.Mth,'location','best')
                    else
                        title(sprintf('Bin %s',Mrg(1).Bn.BinLbl{iBin}))
                    end
                    
                end
                savePics(sprintf('Figures/Stg5_Contour_2_Binnned_%s',Mrg(iAsc+1).RspLbl));
            end
            
        end %plot
    end%methods
    
end %end class


