clear; clc; close all;

addpath('..\..\Code');
%% Stage5_Contour
%%%%%%%%%%%%%%%%%%%%%%%%%
%load marginal model
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
%load Heffernan and Tawn model
load('Output/HT','HT');

%% Contour
Mth={'Hus','Exc','HTDns'};   %cell array of contour methods to be used 

    %'Exc' constant exceedence contour, 
  %'Hus'  Huseby contour method
  %'HTDns' constant density contour on standard Margins, uses density form of H&T to get contour         
   
nSml = 1e6; %number of simulations under H&T model (may need to increase this when you have lots of bins, or see non-smooth Huseby contours)  

%% Estimate Contour
Cnt=Contour(HT,Mrg,Mth,nSml);
Cnt = makeContours(Cnt,Mrg,Mth,Cnt.Sml.A,HT);
Plot(Cnt,Mrg)


%% Save result
if ~exist('Output','dir')
   mkdir('Output') 
end
save('Output/Cnt','Cnt')


