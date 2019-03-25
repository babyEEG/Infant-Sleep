function paramData = fit_model_to_data()
% Anton Tokariev (University of Helsinki, BABA Center, Finland)
% James Roberts  (QIMR Berghofer, Brain Modelling Group, Australia)
%
% This script fits synthetic connectivity data (model) to empirical data
%
% Input:
% Define the path to empirical connectivity data
% Example: load('Fit model to data\Connectivity_data_alpha.mat');
% AND to the file with model parameters
% Example: load('Model parameters\model_params.mat');
%
% Output: paramData - mode weights (a1, a3, a4) for the best model fit 
%
% Copyright (C) 2019
% 
% This is a free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation
% 
% This program is distributed WITHOUT ANY WARRANTY
% (without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE)
%
% See the GNU General Public License: http://www.gnu.org/licenses

  %close all
  %clearvars

% Load broadband AAC connectivity for alpha range (empirical data)
  load('Fit model to data\Connectivity_data_alpha.mat');

% Load mean/std of simulated conn. matr. (for cost function)
% FCmodel_mean: mean of funct.conn.(FC) matrix
% FCmodel_std:  std  of funct.conn.(FC) matrix
% params: matrix with coeffs. grid (for grid search)
  load('Model parameters\model_params_short.mat');

% Restore original combinations of coeffs (aj)
% paramList_full - original full grid (over 8M combinations)
% paramList_full - full grid can be computed using 'compute_Model.m'
  paramList_full = combvec(params(:, 1)', params(:, 2)', params(:, 3)');   %#ok<NODEF>
  paramList_full = fliplr(paramList_full');
% Fragment of the full grid (containing optimal aj for test data) for demo 
  paramList = paramList_full(1250000:1260000, :);
  
% Load FidelityOperator (for conn.matr. correction)  
  load ('Head Model\FidelityOperator.mat'); 
   
% Load eigenmodes for infant cx
  load ('Eigenmodes\eigModes.mat');
    
  eigModes = eigModes';                                                    %#ok<NODEF>

  
% Empirical data   
   Data_mean = mean(nonzeros(triu(Connectivity_data)));
    Data_std = std(nonzeros(triu(Connectivity_data)));

   
   
% FITTING  
% ======================================================================= %  
% Find min of cost f = (meanModel-meanData)^2 + (stdModel-stdData)^2        
  [~, ind] = min((FCmodel_mean - Data_mean).^2 + (FCmodel_std - Data_std).^2);
% Find corresponding mode coeffs (aj) in the list     
  paramData = paramList(ind, :);
      
% Empirical data weights
  FCdat = Connectivity_data;
 
% Model best fit weights  
  FCmod = get_modelFC(paramData, eigModes, FidelityOperator); 

% Show Model-to-Data fitting  
  plot_hist_fit(FCdat, FCmod);

end 


function FCfit = get_modelFC(x, eigModes, FidelityOperator)

  Nch = size(eigModes, 2); % # nodes  
   
% Signal
  L = 30000;     % length
 Fs = 100;       % sampling rate
  t = (1:L)'/Fs; % time
  w = 2*pi*10;   % carrier frequency (10 Hz)
 nu = 2*pi*0.1;  % frequency of amplitude modulation (0.1 Hz)
 
  signal = hilbert(cos(nu*t).*cos(w*t));
 
% alpha Noise
  sigma = 1;
 
  [b1, a1] = butter(7, 13/(Fs/2), 'low');
  [b2, a2] = butter(7, 8/(Fs/2), 'high');

  noise = hilbert(filtfilt(b1, a1, filtfilt(b2, a2, sigma*randn(L, Nch)))); % time-by-ch
  
  z = (x(1)*eigModes(1, :) + x(2)*eigModes(3, :) + x(3)*eigModes(4, :)) .* signal + noise;
  
  FCfit = corr(abs(z)).*FidelityOperator;
 
end


function plot_hist_fit(FCdat, FCmod)
% Anton Tokariev (University of Helsinki, BABA Center, Finland)

    FCdat = nonzeros(triu(FCdat));
    FCmod = nonzeros(triu(FCmod));
    
    xhist = -0.5:0.02:0.5;
       
    figure;
    set(gcf, 'Color', 'w');
    hold on 
  
    h1 = histogram(FCdat, xhist, 'Normalization', 'probability'); 
    h2 = histogram(FCmod, xhist, 'Normalization', 'probability'); 
   
    set(h1, 'Facecolor', [0.50 0.50 0.50]); % gray - data
    set(h2, 'Facecolor', [1.00 0.35 0.00]); % orange - model
 
    set(h1, 'Facealpha', 0.6);
    set(h2, 'Facealpha', 0.6);
  
  % For this example  
    xlim([-0.1 0.4]);
    ylim([0 0.2]);
    
   %[rho, ~] = corr(FCdat, FCmod, 'Type', 'Pearson');
   % title(['edge-to-edge correlation R = ' num2str(roundn(rho, -2))]);
   
    legend('data', 'model')
   
    title('Fitting Model to Data');   
    xlabel('weights', 'Fontweight', 'bold');
    ylabel('probability', 'Fontweight', 'bold');
    
    hold off
 
end