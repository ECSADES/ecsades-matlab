clear; clc; close all;
addpath('PPC\Code');

%% Stage 1: Peak-picking
% Extract peaks data 

load('PPC\CNS_mo_response')

%% Parameters
RspLbl={'Hs','Tp'}; %main and associated variable labels
CvrLbl={'Direction'}; %covariate labels

Rsp=wave.hs;  %main response data

Cvr=[wave.dm];  %covariate data (e.g. direction, season, location)
IsPrdCrv=1;  %flag dictating periodicity of covariate(s). If 1, covariate data loops on 360. When have >1 covariate, this is a vector input e.g. [1,0].   

Asc=[wave.tp_sea];  %associated variable(s)   

NEP=0.7; %non-exceedence quantile level used to set peak picking threshold 

%% Peak Picking

if ~exist(fullfile(cd,'Figures'),'dir')
    mkdir('Figures')
end
    
Dat=PeakPicking(Rsp,Cvr,Asc,IsPrdCrv,NEP,RspLbl,CvrLbl);

%% Save
if ~exist('Output','dir')
    mkdir('Output')
else
    F=cellstr(ls('Output/*.mat'));
    if numel(F)>0
       warning('Existing files in output directory') 
    end
end
save('Output\Data','Dat');
