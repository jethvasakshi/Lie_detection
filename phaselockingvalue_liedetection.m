
clc; 
clear all;

%% frequency band of interest 
freq_band = [29 46];  % (e.g., alpha: [8, 12])
fs = 128;  

%% band-pass filter
[b, a] = butter(2, freq_band / (fs / 2), 'bandpass');  % 2nd order Butterworth filter


%% Truth condition

folder = "C:\Users\saksh\Desktop\LIE WAVES-PAUSE\Truth cond_ATAR";
files = dir(fullfile(folder, '*.csv'));  

for k = 1:numel(files)
    % read each file and store the data
    F = fullfile(folder, files(k).name);   
    files(k).data = readmatrix(F);
    disp(['Processing Subject ', num2str(k), '...']);
    sub = (files(k).data)';  % transposing
    
    count = 0;
    clear alltrails;
    
    % data segmentation
    for i = 1:384:size(sub, 2)
        count = count + 1;
        temp = sub(:, i:i+383);  
        baseline_segment = temp(:, 1:128);  
        baseline_value = mean(baseline_segment, 2);  
        temp_baseline_corrected = temp - baseline_value;  
        
        % band-pass filtering to the segment
        temp_filtered = filtfilt(b, a, temp_baseline_corrected')'; 
        
        alltrails(:, :, count) = temp_filtered;  
    end
    
    % number of trails, channels, and time points
    num_trails = size(alltrails, 3);  
    num_channels = size(alltrails, 1);  
    time_points = size(alltrails, 2);   
    
    % PLV computation
    for trail_idx = 1:num_trails
        trail_data = alltrails(:, :, trail_idx);
        n_channels = size(trail_data, 1);  
        channels_phase = zeros(size(trail_data));  
        
        % Hilbert transform
        for p = 1:n_channels
            channels_phase(p, :) = angle(hilbert(trail_data(p, :)));
        end 
        
        plv_val = zeros(n_channels, n_channels);
        
        % PLV between all pairs of channels
        for c = 1:n_channels
            for d = c+1:n_channels
                phase_diff = channels_phase(c, :) - channels_phase(d, :);
                plv_val(c, d) = abs(mean(exp(1i * phase_diff)));
            end
        end
        
        % PLV values for the current trail
        trails_plv(:, :, trail_idx) = plv_val;
    end
    
    % Average PLV across all trails for the current subject
    mean_trail_plv = mean(trails_plv, 3);  % Average over trails
    
    all_sub_plv(:, :, k) = mean_trail_plv;
end
    
% averaged PLV for all subjects 
truth_plv = mean(all_sub_plv, 3);

% Plotting PLV matrix for truth condition
figure;
heatmap(truth_plv);
colorbar;
title('PLV Matrix for Truth Condition');
xlabel('Channels');
ylabel('Channels');


%% Lie condition


folder = "C:\Users\saksh\Desktop\LIE WAVES-PAUSE\Lie cond_ATAR";
files = dir(fullfile(folder, '*.csv')); 

for k = 1:numel(files)
    F = fullfile(folder, files(k).name);   
    files(k).data = readmatrix(F);
    disp(['Processing Subject ', num2str(k), '...']);
    sub = (files(k).data)';  
    
    count = 0;
    clear alltrails;
    
    for i = 1:384:size(sub, 2)
        count = count + 1;
        temp = sub(:, i:i+383);  
        baseline_segment = temp(:, 1:128);  
        baseline_value = mean(baseline_segment, 2);  
        temp_baseline_corrected = temp - baseline_value;  
        temp_filtered = filtfilt(b, a, temp_baseline_corrected')';  
        alltrails(:, :, count) = temp_filtered;  
    end
    
    num_trails = size(alltrails, 3);  
    num_channels = size(alltrails, 1);  
    time_points = size(alltrails, 2);   
    
    for trail_idx = 1:num_trails
        trail_data = alltrails(:, :, trail_idx);
        n_channels = size(trail_data, 1); 
        channels_phase = zeros(size(trail_data));  
        for p = 1:n_channels
            channels_phase(p, :) = angle(hilbert(trail_data(p, :)));
        end 
        
        plv_val = zeros(n_channels, n_channels);
        
        for c = 1:n_channels
            for d = c+1:n_channels
                phase_diff = channels_phase(c, :) - channels_phase(d, :);
                plv_val(c, d) = abs(mean(exp(1i * phase_diff)));
            end
        end
       
        trails_plv(:, :, trail_idx) = plv_val;
    end
   
    mean_trail_plv = mean(trails_plv, 3); 
    all_sub_plv(:, :, k) = mean_trail_plv;
end

lie_plv = mean(all_sub_plv, 3);

% Plotting PLV matrix for lie condition
figure;
heatmap(lie_plv);
colorbar;
title('PLV Matrix for Lie Condition');
xlabel('Channels');
ylabel('Channels');


%% Paired sample t-test

% matrices to vectors
truth_vector = truth_plv(:);
lie_vector = lie_plv(:);

% filtering non-zero
nonZeroIdx1 = truth_vector ~= 0;
nonZeroIdx2 = lie_vector ~= 0;
validIdx = nonZeroIdx1 & nonZeroIdx2;
truth_pairs = truth_vector(validIdx);
lie_pairs = lie_vector(validIdx);

% paired t-test
[h, p, ci, stats] = ttest(truth_pairs, lie_pairs);

% results display
disp(['Paired t-test p-value: ', num2str(p)]);
disp('T-test statistics:');
disp(stats);

%% Difference matrix

difference_PLV_matrix= abs(truth_plv-lie_plv)

figure;
imagesc(difference_PLV_matrix)
title('Imaginary Coherence (Truth-Lie)');
xlabel('Channels');
ylabel('Channels');

figure;
heatmap(difference_PLV_matrix)
title('Imaginary Coherence (Truth-Lie)');
xlabel('Channels');
ylabel('Channels');
