function compute_Model()
% Anton Tokariev (University of Helsinki, BABA Center, Finland)
% James Roberts  (QIMR Berghofer, Brain Modelling Group, Australia)
%
% Eigenmodes for infant cortical surface were computed using original
% Matlab script developed by Moo K. Chung and colleagues.
% The original code can be found here:
% http://brainimaging.waisman.wisc.edu/~chung/lb/
%
% 'eigenModes.mat' provided with current code contains eigenvalues of the
% first four strongest modes (columns = modes; rows = parcel centroids) 
%
% Model uses 3 modes: a1(mode1) + a3(mode3) + a4(mode4); a2 = 0 (no mode2)
% aj - coefficients (weights) of modes 
% This Model simulates brain activity in 'alpha' frequency range
% For more details see original paper
%
% Input: define the range for eigenmode weights in 'Eigenmode coefficients'
% Example:
%           a1 = linspace(0, 0.5, 201); % range [0...0.5]; step = 0.0025 
%           a3 = linspace(0, 0.5, 201); % range [0...0.5]; step = 0.0025
%           a4 = linspace(0, 0.5, 201); % range [0...0.5]; step = 0.0025
%
% Output: mean (FCmodel_mean) and std (FCmodel_std) values for synthetic
%         connectivity matrices obtained with different combinations of
%         mode weights. These data will be used for the model fitting to
%         data in 'fit_model_to_data'
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
  
% load eigenmodes for infant cx surface
  load('Eigenmodes\eigModes.mat');
  
% Visualize first 4 strongest eigenmodes for infant cx  
  load('Head Model\Atlas.mat'); % coordinates of parcel centroids 
  
  for k = 1:size(eigModes, 2) %#ok<NODEF>
      figure;
      set(gcf, 'Color', 'w');
      hold on
      scatter3(Atlas.Centroids(:,1), Atlas.Centroids(:,2), Atlas.Centroids(:,3), 80, eigModes(:, k), 'filled');
     %colormap('cool');
      colorbar
      title(['Mode ' num2str(k)]);
      hold off
  end
  
   
  eigModes = eigModes';
    
  Nch = size(eigModes, 2); % Number of parcels
   
        
% Generate MODEL: a1 = var; a2 = 0; a3 = var; a4 = var
% ----------------------------------------------------------------------- %
% Eigenmode coefficients (final range) 
  a1 = linspace(0, 0.5, 201);
  a3 = linspace(0, 0.5, 201);
  a4 = linspace(0, 0.5, 201);

% Grid of coefficients (parameter space)  
  params = [a1; a3; a4]'; %#ok<NASGU>
    
  paramList = combvec(a1, a3, a4); % combinations of coeffs (aj) as a list
  paramList = fliplr(paramList');
   
% Signal
  L = 30000;     % length
 Fs = 100;       % sampling rate        
  t = (1:L)'/Fs; % time   
  w = 2*pi*10;   % carrier frequency (10 Hz)
 nu = 2*pi*0.1;  % frequency of amplitude modulation (0.1 Hz)
 
  signal = hilbert(cos(nu*t).*cos(w*t)); % complex values

% alpha 'noise'
  sigma = 1;
 
% 'fast' bandpass filter for 'alpha' range  
 [b1, a1] = butter(7, 13/(Fs/2), 'low');
 [b2, a2] = butter(7, 8/(Fs/2), 'high');

  noise = hilbert(filtfilt(b1, a1, filtfilt(b2, a2, sigma*randn(L, Nch)))); % time-by-ch

% Plot example of simulated activity  
  figure;
  set(gcf, 'Color', 'w');
  subplot(311);
  plot(t, real(signal), 'r');
  ylabel('signal');
      subplot(312);
      plot(t, real(noise(:,1)), 'Color', [.5 .5 .5]);
      ylabel('noise');
          subplot(313);
          plot(t, real(signal)+real(noise(:,1)), 'b');
          ylabel('combination');
          xlabel('time, sec');
  
  
% Binary mask for correction of Connectivity matrix (see Tokariev 2018)  
  load('Head Model\FidelityOperator.mat'); 
% if no FidelityOperator available, generate it as:
% FidelityOperator = ones(Nch);   
  
% init output 
  FCmodel_mean = zeros(size(paramList, 1), 1);
  FCmodel_std  = zeros(size(paramList, 1), 1);

  tic
  for k = 1:10000 %size(paramList, 1) % check all coeffs (aj) combinations
         
    % Apply combination of coefficients (aj) from parameter list: a(j)mode(j) x signal + noise  
      z = (paramList(k, 1)*eigModes(1, :) + paramList(k, 2)*eigModes(3, :) + paramList(k, 3)*eigModes(4, :)) .* signal + noise;
      
    % Remove non-reliable edges from simulated connectivity matrix
      FCbuf = corr(abs(z)).*FidelityOperator;
    
    % Mean of simulated connectivity matrix   
      FCmodel_mean(k, 1) = mean(nonzeros(triu(FCbuf, 1)));
      
    % STD of simulated connectivity matrix
      FCmodel_std(k, 1) = std(nonzeros(triu(FCbuf, 1)));
      
  end
  toc 
 
  save('Model parameters\model_params.mat', 'FCmodel_mean', 'FCmodel_std', 'params');
  
end


