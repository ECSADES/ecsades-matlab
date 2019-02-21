clear; clc; close all;
addpath('..\Code');
%% Stage4_FitH&T

FM=cellstr(ls('Output/MM*.mat'));
for iM=1:numel(FM) %load all marginal models
    load(sprintf('Output/MM%g',iM),'MM')  %Load margin 2
    if iM==1
        Mrg=MM;
    else
        Mrg=cat(1,Mrg,MM);
    end
    clear MM;
end

%% Fit heffernan and Tawn model
HTNEP=[0.7,0.85];  %Conditional Extreme NEP

NonStationary= true;  %if true: estimate the H&T alpha parameter in a piecewise constant non-stationary way using the same bins as the Marginal analysis;
                      %if false: estimate stationary H&T alpha
SampleLocalResid = true; %if true: when simulating under H&T model, resample residuals locally from same covariate sector; if false: resample residuals from any sector
                         %Note: if you have any bins with very few observations; set this to false.  

%% Cross Validation defaults (Optional can specify some or all of these)
CV.CVMth=0;     %0 Only Cross Validate smoothness for original dataset (fast);
                %1 Cross Validate smoothness for every bootstrap resample (slow),
CV.nCV=2;      %number cross-validation groups
CV.nSmth=10;    %number smoothnesses tried in CV
CV.SmthLB=-8;   %lower bound (log10)  for smmothness range
CV.SmthUB=12;   %upper bound (log10)  for smmothness range

%% Fit model
HT=HeffernanTawn(Mrg,HTNEP,NonStationary,CV,SampleLocalResid); %Fit Heffernan & Tawn joint exceedance model

%% Save result
if ~exist('Output','dir')
   mkdir('Output') 
end
save('Output/HT','HT')

Plot(HT,Mrg)







