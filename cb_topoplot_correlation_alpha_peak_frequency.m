%% Sensor topography of correlation between peak alpha frequency and age in each sensor

clc; close all; clear;

%load('../key_analyses/roiLayout.mat', 'layout'); % load the sensor layout i.e. coordinates for Neuromag magnetometer
load('C:/Bikash scripts/key_analyses/welchSpectrumTopomap.mat', 'peakAlphaFreq', 'megAge'); % load the peak alpha frequency values and age values

% Spearman's rank correlation between age and peak alpha frequency
[C, P] = corr(peakAlphaFreq', megAge, 'type', 'spearman'); 

fS = cbl_figure([0.03 0.05 0.3 0.43]); % create a figure with the mentioned size
cbl_topoplot(layout, [C'; nan(2,1)], 'both'); axis square off;
colorbar; colormap jet;
saveas(fS, 'cb_topotmap_peak_alpha_frequency_correlation', 'png'); % save the plot in PNG format
[~, fdrP] = fdr_bh(P);