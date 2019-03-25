% Pipeline compilation from 25/03/2019
% by Anton Tokariev (Baby Brain Activity Center/University of Helsinki)
%
% Original paper:
%'Large-scale brain modes reorganize between infant sleep states and carry
% prognostic information for preterms'
%
% This pipeline computes amplitude-amplitude connectivity for cortical
% signals from multi-channel EEG data
%
% Input: define the input EEG file (which is in .mat format)
%        see 'INPUT EEG DATA' section below
%        Example: load('Data example\filename.mat');]
%
% Output: broadband (AAC_bb) and narrowband (AAC_nb) connectivity
%
% For detailed instructions how to use the script, see README file
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


function compute_AAC()

 tic
      
 %close all 
 %clearvars

% INPUT EEG DATA
% ======================================================================= %
% eeg_data.eeg    - input EEG [channels x samples]
% eeg_data.Fs     - sampling rate
% eeg_data.Labels - channel names
  load('Data example\infant_eeg.mat'); %[channels x samples] 

% Visualization:  
  plot_signals(eeg_data.eeg', eeg_data.Labels, eeg_data.Fs, 'EEG example');               
     
     
% LOAD FILTERS
% ======================================================================= %
% N = 4  EEG filters (delta, theta, alpha, beta)
 %load('Filters\filters_eeg_4_bands.mat');   % top: HPF; bottom: LPF
   
% N = 21 EEG filters (extended set)
  load('Filters\filters_eeg_21_bands.mat');  % top: HPF; bottom: LPF     
   
 % N = 15 Amplitude filters 
   load('Filters\filters_amp_15_bands.mat'); % top: HPF; bottom: LPF 
   
   N_flt_eeg = size(flt_eeg, 2); % N of filters for EEG
   N_flt_amp = size(flt_amp, 2); % N of filters for amplitudes
   
   disp(['N = ' num2str(N_flt_eeg) ' EEG filters vs. N = ' num2str(N_flt_amp) ' amplitude filters']); 
 
% LOAD HEAD MODEL
% ======================================================================= %         
  load('Head Model\Atlas.mat');
  load('Head Model\InverseOperator.mat');   
  load('Head Model\CollapseOperator.mat');
  

% ANALYSIS
% ======================================================================= %
% 1. Filter EEG
% ---------------------------------------------------
  disp('Filter EEG...');
  eeg_filtered = filter_eeg(eeg_data.eeg', flt_eeg);
% Output > eeg_filtered: cell{1 x N filters}(samples x channels)
 
% Visualization:
  %plot_signals(eeg_filtered{1, 3}, eeg_data.Labels, eeg_data.Fs, 'Filtered EEG');

  
% 2. Compute parcel signals
% ---------------------------------------------------
  disp('Compute cortical signals...');
  parcels = get_parcel_signals(eeg_filtered, InverseOperator, CollapseOperator, Atlas);
% Output > parcels: cell{1 x N filters}(samples x parcels)

% Visualization
  %plot_signals(parcels{1, 3}(:, 1:20), Atlas.Areas(1:20), eeg_data.Fs, 'First 20 (of 58) parcel signals');
  
  
% 3. Hilbert transform
% ---------------------------------------------------
  disp('Hilbert transform...');
  edge = 500;
  for n = 1:N_flt_eeg
      L = size(parcels{1, n}, 1);                                                                        % signal length
      buf = [flipud(parcels{1, n}(1:edge, :)); parcels{1, n}; flipud(parcels{1, n}(end-edge+1:end, :))]; % control edge effects in Hilbert
      buf = hilbert(buf);                                                                                % get analytic signal  (complex)
      parcels{1, n} = buf(edge+1:edge+L, :);                                                             % cut edges
  end
% Output > parcels: cell{1 x N filters}(samples x parcels): complex
      


  AAC_bb{1, N_flt_eeg} = [];         % init conn.matr. for fullband AAC 
     
  AAC_nb{N_flt_amp, N_flt_eeg} = []; % init conn.matr. for narrowband AAC
  
  disp('Computing functional connectivity... Wait...');
  disp(repmat('=', 1, 50));
  
  for n = 1:N_flt_eeg % loop for fr.bands
      
    disp(['Frequency band i = ' num2str(n) ' out of ' num2str(N_flt_eeg)]);
   
  % 4. Signal orthogonalization
  % ---------------------------------------------------  
    A_orth = orthogonalize_signals(parcels{1, n}, 200);
  % Output > A_orth [ch x ch x samples] envelopes of orthog. signals
  
  
  % 5. Compute broadband (bb) amp-amp connectivity  (AAC)
  % ---------------------------------------------------
    AAC_bb{1, n} = compute_AAC_bb(A_orth, parcels{1, n});
  % Output > AAC_bb cell{1 x N filters} AAC data [ch x ch]
  
 
  % 6. Compute narrowband (nb) amp-amp connectivity (AAC)
  % ---------------------------------------------------
    AAC_nb(:, n) = compute_AAC_nb(A_orth, parcels{1, n}, flt_amp);
  % Output > AAC_nb cell{N amp.filters x N filters} AAC data [ch x ch]
  
    disp('Done');

  end
   
  disp(repmat('=', 1, 50));
   
  % Visualization of step (5) (fr = 3: alpha band)
    %fr = 3; imagesc(AAC_bb{1, fr});colorbar; hold on; set(gcf, 'Color', 'w'); xlabel('parcels'); ylabel('parcels'); title('Broadband AAC (alpha)', 'Fontweight', 'bold');
    
% 7. Connectivity matrix correction
%   (see Tokariev et al., 2019, Cer.Cor. for details)  
% ---------------------------------------------------
  disp('Correction of connectivity matrices...');
  
  load('Head Model\FidelityOperator.mat'); % mask to remove 'non-reliable' edges
  
% correct broadband  
  for n = 1:N_flt_eeg
      AAC_bb{1, n} = AAC_bb{1, n}.*FidelityOperator;
  end
  
% correct narrowband  
  for n1 = 1:N_flt_amp
      for n2 = 1:N_flt_eeg
          AAC_nb{n1, n2} = AAC_nb{n1, n2}.*FidelityOperator;
      end
  end 
  
  
  
  disp('Saving connectivity matrices to "AAC output" folder...');
  
% Save AAC data     
  save('AAC_output\AAC_bb.mat', 'AAC_bb'); % broadband AAC
  save('AAC_output\AAC_nb.mat', 'AAC_nb'); % narrowband AAC

 toc 

end % end of main function


function plot_signals(s, lables, Fs, ttl)
% Anton Tokariev (University of Helsinki, BABA Center, Finland)
%
% Input:
%      s: signal        [channels x samples]
% labels: channel names [channels x 1]
%     Fs: sampling rate
%    ttl: title (char, example ttl = 'Title') 

    L = size(s, 1); % length
   ch = size(s, 2); % number of channels
        
   shift = 5*(max(std(s))); % b/w different channels for visualization
   
   baseline = (1:1:ch)*shift;
   baseline = fliplr(repmat(baseline, L, 1));
  
 % time   
   t =(0:1:L-1)/Fs;
   t = repmat(t, ch, 1);
   t = t';
     
% add baseline to each channel for visualization
  sig_plot = s + baseline;
  
% plot signals
  figure;
  set(gcf, 'Color', 'w');
  hold on
  
  plot(t, sig_plot, '-k');
  
  xlabel('Time, sec', 'Fontweight', 'bold');
  ylabel('Channels',  'Fontweight', 'bold');
  
  set(gca, 'YTick', fliplr(baseline(1,:)), 'YTickLabel', fliplr(lables));
  
  title(ttl); % from input
    
  hold off

end



function signal_flt = filter_eeg(signal, flt)
% Anton Tokariev (University of Helsinki, BABA Center, Finland)
%
% Input:
% signal: [samples x ch]
%    flt: cell array with filter objects (top - HPF, bottom - LPF)
%    Fs: sampling rate
%
% Output:
% signal_flt: cell array{1, N of fr.bands} of filtered data [samples x ch]
 
 % Signal Length
   L = size(signal, 1);

 % Flipped signal     
   signal_ud = flipud(signal);
   
 % Add pieces to signal (to fix edge effects)
   signal_big = [signal_ud; signal; signal_ud];                                                 
  
 % Number of bandpass filters (= fr.bands); columns = bandpass filters  
   N_flt = size(flt, 2); 
 
 % Init cell array for bandpass filtered signals   
   signal_flt{1, N_flt} = []; 
   
% Filter signals (bandpass filter = HPF + LPF)  
  for fr = 1:N_flt
      
      buf = [];                                                            %#ok<NASGU>
      buf = filtfilt(flt{1, fr}, signal_big); % HPF/cutoff = 0.85xFc 
      buf = filtfilt(flt{2, fr}, buf);        % LPF/cutoff = 1.15xFc
      
      signal_flt{1, fr} = buf(L+1:L+L, :);    % cut signal_big >> orig.sig.
      
  end
    
end



function signal_src = get_parcel_signals(signal, InverseOperator, CollapseOperator, MyAtlas)
 % Anton Tokariev (University of Helsinki, BABA Center, Finland)
 %
 % see also Tokariev et al., 2018, Cerebral Cortex
 %
 % Input:
 % signal: cell array {1 x N freq.} of filtered EEG data [samples x ch]
 % InverseOperator: Inverse solution for infant head with N = 19 trodes
 % CollapseOperator: weights for src signals within 'host' parcels
 % MyAtlas: assignment of src/verticies to parcels (in MyAtlas.Parcels)
 %
 % Output:
 % signal_src: cell array {1 x N freq.} of filtered parcel signals [samples x ch] 
 
   N_fr = size(signal, 2);       % number of frequencies
   
   Np = length(MyAtlas.Parcels); % number of parcels
   
   L = size(signal{1, 1}, 1);    % EEG length 
   
   signal_src{1, N_fr} = [];     % init output array
   
   CollapseOperator = repmat(CollapseOperator, 1, L);
   
   for k = 1:N_fr
      
     % source signals  
       src_buf = (InverseOperator * signal{1, k}').* CollapseOperator;     % [sources x samples]
       
     % parcel/cortical signals  
       parcel_buf = zeros(Np, L);
  
       for j = 1:Np
           parcel_buf(j, :) = mean(src_buf(MyAtlas.Parcels{j, 1}, :));     % [parcels x samples] 
       end
        
       signal_src{1, k} = parcel_buf';                                     % [samples x parcels] 
       
   end

end



function A_orth = orthogonalize_signals(A, WindowLength)
% Anton Tokariev     (University of Helsinki, BABA Center, Finland)
% Alexander Zhigalov (University of Birmingham, School of Psychology, UK)
%
% Implemented based on Brookes et al., 2012, Neuroimage 
%
% Input:
% A: signals [samples x channels] (complex values)
% WindowLength: non-overlapping windows for orthogonalization [samples]
%
% Output:
% A_orth: array with envelopes of mutually orthogonalized signals [channels x channels x samples]
%
% Example:
% A_orth(1,2,:) contains envelope of signal A(:,2) orth. re signal A(:,1) 
% A_orth(2,1,:) contains envelope of signal A(:,1) orth. re signal A(:,2)
  
    N_Windows = fix(size(A, 1)/WindowLength);
  
    A(WindowLength*N_Windows+1:end, :) = []; % cut extra samples at the end
    
      L = size(A, 1);
    Nch = size(A, 2); 
      
    A_orth = zeros(Nch, Nch, L); % init 3D array 
   
    for ChA = 1:Nch % dim 1
     
         X = A(:, ChA);     % [samples x 1]
                
         for ChB = 1:Nch % dim 2
             
             Y = A(:, ChB); % [samples x 1]
             
             if ChA ~=ChB % ignore diagonal
                     
              % Reshape arrays: 1 Column = 1 window  
                X_blocks = reshape(X, WindowLength, N_Windows);
                Y_blocks = reshape(Y, WindowLength, N_Windows);
                
              % Orthogonalize Y re X >> Y_orth
                B = sum((real(X_blocks)./repmat(sum(real(X_blocks).^2, 1), WindowLength, 1) ).*real(Y_blocks), 1);
       
                Y_orth_blocks = Y_blocks - repmat(B, WindowLength, 1) .* X_blocks;
       
              % Put Y_orth to A_orth(ChA, ChB, :) 
                A_orth(ChA, ChB, :) = reshape(Y_orth_blocks, L, 1);
             
             end % end of if
                       
         end % ChB loop
         
    end % ChA loop
    
    A_orth = abs(A_orth); % envelopes from analytic signal
      
end



function oCC = compute_AAC_bb(A, B)
% Anton Tokariev (University of Helsinki, BABA Center, Finland)
%
% Computes correlations between broadband amplitude envelopes of signal B
% and signal A (A was preliminarily orthogonalized relative to signal B)
%
% Input:
% A: envelopes of mutually orthogonalized signals [ch x ch x samples]
% B: original signals                             [samples x ch]: complex
%
% Output:
% oCC: connectivity matrix containing orthogonalized corr. coeffs. (oCC)
% oCC [channels x channels]
 
  d1 = size(A, 1); % N of ch 
  d2 = size(A, 2); % N of ch
   
  B = abs(B); % envelopes from complex bandpass filtered signals
   
  oCC = zeros(d1, d2); % init conn.matrix
   
  for ChA = 1:d1
       
      X = B(:, ChA); % orig. signal
        
      for ChB = 1:d2
          
          if ChA ~=ChB % ignore diagonal
           
             Y = squeeze(A(ChA, ChB, :)); % orth. signal
           
             oCC(ChA, ChB) = corr(X, Y);
           
          end

       end % ChB
           
  end % ChA
  
% Make oCC symmetric  
  oCC = (oCC + oCC')./2; 
  
end



function oCC = compute_AAC_nb(A, B, flt_amp)
% Anton Tokariev (University of Helsinki, BABA Center, Finland)
%
% see also Hipp et al., 2012, Nature Neuroscience
%
% Computes correlations between narrowband (= bandpass filtered) amplitude
% envelopes of signal B and signal A
% Note, that signal A was preliminarily orthogonalized relative to signal B
%
% Input:
% A: envelopes of mutually orthogonalized signals [ch x ch x samples]
% B: original signals                             [samples x ch]: complex
% flt_amp: cell array of filter objects for amp.env. (top:HPF, bottom: LPF)
%
% Note!
% Amplitude envelopes in this script are downsampled from 100 Hz to 20 Hz
% So, filters for amplitudes were designed for sampling rate Fs = 20 Hz
%
% Output:
% oCC: cell array {1 x N of amplitude filters}, where each cell contains
% connectivity matrix [channels x channels] of orthogonalized corr. coeffs.
% corresponding to bandpass filtered amplitude envelopes (within fr.range
% defined by filters in flt_amp: columns = bandpass filters, upper cell is
% high-pass filter, lower cell is low-pass filter)
 
  oCC{size(flt_amp, 2), 1} = [];
   
   d1 = size(A, 1);
   d2 = size(A, 2);
   d3 = size(A, 3); % length
   
   B = abs(B); % envelope from complex banpass filtered signals
   
 % Reshape to 2D
   A_2D = reshape(A, d1*d2, d3, 1);
   A_2D = A_2D';                          % [samples x ChA-ChB]
   
 % Pick non-zero columns (to speed up computation)
   A_2D_good = A_2D(:, all(A_2D, 1));
   
 % Resampling (with anti-aliasing FIR LPF)  
   A_2D_good = resample(A_2D_good, 1, 5); % 100 Hz >> 20 Hz
   
   d3 = size(A_2D_good, 1); % new length
   
   B = resample(B, 1, 5);                 % 100 Hz >> 20 Hz   
     
 % Mirrored signal to edges (fix edge effects)
   edge = 1800; % 90 sec (Fs = 20 Hz)  

   A_2D_good = [flipud(A_2D_good(1:edge, :)); A_2D_good; flipud(A_2D_good(end-edge+1:end, :))];
     
   B = [flipud(B(1:edge, :)); B; flipud(B(end-edge+1:end, :))];
   
   for k = 1:size(flt_amp, 2) % amp.freq. loop (flt_amp for 20Hz!!!)
       
       buf_A = [];                                                         %#ok<NASGU>
       buf_A = filtfilt(flt_amp{1, k}, A_2D_good); % Hd1/HPF
       buf_A = filtfilt(flt_amp{2, k}, buf_A);     % Hd2/LPF
       
       buf_B = [];                                                         %#ok<NASGU>
       buf_B = filtfilt(flt_amp{1, k}, B);         % Hd1/HPF
       buf_B = filtfilt(flt_amp{2, k}, buf_B);     % Hd2/LPF
   
     % Reconstruct to 3D array for orth signals
       A_flt = [];                                                         %#ok<NASGU>
       
       A_flt = zeros(d3, d1*d2);                         % [length x ch*ch]
   
       A_flt(:, all(A_2D, 1)) = buf_A(edge+1:edge+d3, :);% filtered 
  
       A_flt = A_flt';                                   % [ch*ch x length]                   
   
       A_flt = reshape(A_flt, d1, d2, d3);               % [ch x ch x length]
   
     % cut filtered original envelopes
       B_flt = [];                                                         %#ok<NASGU>
       
       B_flt = buf_B(edge+1:edge+d3, :);
  
  
     % Compute oCC for amp.frequency(k) 
       buf_oCC = zeros(d1, d2);
   
       for ChA = 1:d1
       
           X = B_flt(:, ChA);
        
           for ChB = 1:d2
           
                Y = squeeze(A_flt(ChA, ChB, :));
           
                if ChA ~= ChB % ignore diagonal
                    
                   buf_oCC(ChA, ChB) = corr(X, Y);
                   
                end
           
           end % ChB
           
       end % ChA
       
       oCC{k, 1} = (buf_oCC + buf_oCC')./2; % make oCC symmetric and store 
       
   end % fr.loop 
           
end









 
