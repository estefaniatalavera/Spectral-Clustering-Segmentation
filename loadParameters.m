%%%%%
%
% File for loading the main parameters for Spectral-Clustering
%
%%%%%%%%%%%%%%%%%%%%%%

% clear all, close all

%% Load paths
%addpath('Adwin;Data_Loading;Evaluation;Features_Preprocessing');
%addpath('GCMex;GraphCuts;PCA;Tests;Utils')
addpath('Data_Loading;Evaluation;Features_Preprocessing;Utils');
addpath('SpectralFunctions;PCA;SimilaritiesMatrix;ResultsFM&JI');


%% Data loading
% directorio_im = 'D:/LIFELOG_DATASETS'; % SHARED PC
%directorio_im = '/Volumes/SHARED HD/Video Summarization Project Data Sets/R-Clustering'; % MARC PC
directorio_im='/Users/estefaniatalaveramartinez/Desktop/LifeLogging/IbPRIA/Sets'; % EST PC

% directorio_im = ''; % put your own datasets location

%% -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%
% To plot
%
Dataset_ops={'All','Narrative','SenseCam'};
%Dataset_ops={'All'}

% directorio_results = 'D:/R-Clustering_Results'; % SHARED PC
%directorio_results = '/Volumes/SHARED HD/R-Clustering Results'; % MARC PC
directorio_results = '/Users/estefaniatalaveramartinez/Desktop/LifeLogging/IbPRIA/Results'; % EST PC
% directorio_results = ''; % put your own results location


%% R-Clustering parameters
clus_type = 'Clustering'; % Clustering type used before the GraphCuts. 
                        % It can take the following values:
                        %           'Clustering'
                        %           'Both'
paramsPCA.minVarPCA=0.95;
paramsPCA.standarizePCA=false;
paramsPCA.usePCA_Spect = true;
                        


%% Spectral parameters

paramsSpec.NN = false;
NN=5;%2:2:20;
paramsSpec.Sig = true;
Sig=1;%0.5:0.5:10;
paramsSpec.Eps = false;
Eps=1;%0.1:0.1:1.5;

sim_matrix={'Sigma','NN','Epsilon'};
%sim_matrix={'Epsilon'};
k_values=6:4:36;

%% Evaluation parameters
tol=5; % tolerance for the final evaluation
minImCl=0; % (deprecated)