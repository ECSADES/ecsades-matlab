clear; clc; close all;
addpath('..\..\Code');
%% Stage2_BinData
% Choose covariate bins

%% Inputs 
load('Output/Data','Dat')

%Choose bins for each covariate dimension
BinEdg={[25,230,275,315]',[145,270]'}; %In single covariate case: input in format {[]'}. In multiple covariate case: input in format {[]',[]'}.
           %Note: covariate bins must be suitable for all dimensions (main and associated)
                                
%Note: if using non-periodic covariate (you set 'IsPrdCvr' = 0 in Stage1),
%you must provide the endpoint/maximum value of the covariate as last bin-edge

%Note: covariate bins must be suitable for all dimensions (main and associated)
                                
%Note: if using non-periodic covariate (you set 'IsPrdCvr' = 0 in Stage1),
%you must provide the endpoint/maximum value of the covariate as last bin-edge
%% Allocate Data to Bins and Plot
Bn=CovariateBinning(Dat.X,BinEdg,Dat.IsPrd,Dat.Y,Dat.RspLbl,Dat.CvrLbl);

%% Save bin Edges
save('Output/Bin','Bn')
