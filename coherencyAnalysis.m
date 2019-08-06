%% Coherency analysis
% This code estimates Welch Spectrum and global coherence of raw data.
clc; close all; clear;

%% to add Fieldtrip and Chronux toolbox to MATLAB path
addpath(genpath('C:\Users\CBDL\Documents\bikash\toolbox\fieldtrip\'));
ft_defaults;
addpath(genpath('C:\Users\CBDL\Documents\bikash\toolbox\chronux\'));

subIDs = dlmread('CamCAN_list.txt'); % CamCAN_list.txt contains the subject IDs
subIDs = subIDs(:);
subIDs(subIDs==0) = [];

subList = [];
for sI = 1:650
    % directories to save the results must be created before running PARFOR loop
    saveDir = ['C:\Users\CBDL\Documents\bikash\coherencyAnalysis\CC' num2str(subIDs(sI)) '\']; 
    if ~exist(saveDir); mkdir(saveDir); end
    if ~exist([saveDir 'cohrGrad.mat']); subList = [subList sI]; end
end

%% main analysis (parallelized using PARFOR)
subList = 1;
for (sI = 1:length(subList))
    tic
    disp(subIDs(sI));
    dataDir = ['C:\Users\CBDL\Documents\bikash\CamCAN\CC' num2str(subIDs(subList(sI))) '\'];
    saveDir = ['C:\Users\CBDL\Documents\bikash\coherencyAnalysis\CC' num2str(subIDs(subList(sI))) '\'];
    fileName = 'transdef_mf2pt2_rest_raw.fif'; % filename of data that has been head position aligned across subjects
    
    hdr = ft_read_header([dataDir fileName]); % reads the header present in the '.fif' file
    data = ft_read_data([dataDir fileName]); % read data from the file
    data = downsample(data', 4); % downsample data from original sampling freq 1000 Hz to 250 Hz
    eog = data(:, 321:322); % channel number 321 and 322 recorded EOG
    ecg = data(:, 323); % channel number 323 recorded ECG
    data = data(:,1:306); % gradiometer and magnetometer data
    numSamples = size(data,1); 
    samplingRate = 250;
    
    magnetometer = find(strcmp(hdr.chantype, 'megmag')); % channel indices of the magnetometers
    gradiometer = find(strcmp(hdr.chantype, 'megplanar')); % channel indices of the gradiometers
    % separate 'save' function need to be defined to save data inside a parallel loop
    saveData([saveDir 'datasetInfo.mat'], eog, ecg, numSamples, samplingRate, magnetometer, gradiometer);
    
    %% Welch's spectrum of the data
    % calculates Welch's periodogram on nonoverlapping windows of length 20s using 20*Fs=5000 point FFT
    [spectrumData, frequency] = pwelch(data, 20*samplingRate, 0, 20*samplingRate, samplingRate); 
    saveWelchSpectrum([saveDir 'welchSpectrum.mat'], spectrumData, frequency);
    
    %% Global coherency
    % global coherency has been estimated using Chronux
    params = struct();
    params.tapers = [2 3]; % time-bandwidth prodct = 2, number of tapers = 3
    params.Fs = samplingRate;
    params.pad = -1; % -1: no padding, 0: padding to the immediate 2^x
    params.fpass = [0 40]; % frequency range of interest
    winSize = 5; % window size for estimation
    numWin = floor(size(data,1)/winSize/samplingRate); 
    dataMag = data(1:numWin*winSize*samplingRate,magnetometer); % reshaping the resting state data into segments of 5s
    dataGrad = data(1:numWin*winSize*samplingRate,gradiometer);
    % variable 'totalCohrM' contains the global coherence at different frequencies
    [crossSpecM, cohrM, totalCohrM, eigVecM, altCohrM, frequency] = CrossSpecMatc(dataMag, winSize, params);
    [crossSpecG, cohrG, totalCohrG, eigVecG, altCohrG, frequency] = CrossSpecMatc(dataGrad, winSize, params);
    saveGlobalCohr([saveDir 'cohrMag.mat'], crossSpecM, cohrM, totalCohrM, eigVecM, altCohrM, frequency);
    saveGlobalCohr([saveDir 'cohrGrad.mat'], crossSpecG, cohrG, totalCohrG, eigVecG, altCohrG, frequency);
    toc
end
% Defining save function to save results from Welch spectrum analysis
function saveWelchSpectrum(fName, spectrumData, frequency)
    save(fName, 'spectrumData', 'frequency');
end

% Defining save function to save ECG, EOG data and other information
function saveData(fName, eog, ecg, numSamples, samplingRate, magnetometer, gradiometer)
    save(fName, 'eog', 'ecg', 'numSamples', 'samplingRate', 'magnetometer', 'gradiometer');
end

% Defining save function to save results from global coherence analysis
function saveGlobalCohr(fName, crossSpecC, cohrC, totalCohrC, eigVecC, altCohrC, frequency)
    save(fName, 'crossSpecC', 'cohrC', 'totalCohrC', 'eigVecC', 'altCohrC', 'frequency');
end