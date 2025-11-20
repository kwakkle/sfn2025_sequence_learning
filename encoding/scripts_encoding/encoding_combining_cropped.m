%% WT Mice - Combined Crop and Euclidean Distances

close all
clearvars
clc

% WT mouse 
wt_mice = {'OE12', 'OE15', 'OE24', 'OE35', 'OE39', 'OE40', 'OE45', 'OE46', 'OE47'};

% Initialize arrays
C_dff = NaN(5, length(wt_mice), 50, 50);
C_dprime = NaN(5, length(wt_mice), 50, 50);
all_dprime_dists = NaN(length(wt_mice), 5);  % for d′
all_dff_dists = NaN(length(wt_mice), 5);     % for df/F

for i = 1:length(wt_mice)
    mouse_id = wt_mice{i};

    % Load cropped maps
    maps_file = sprintf('cropped_maps_mouse_%s_naive.mat', mouse_id);  % change if expert or naive
    load(maps_file, 'cropped_dff', 'cropped_dprime');

    % Load Euclidean distances
    dist_file = sprintf('euclidean_dists_%s_naive.mat', mouse_id);  % change if expert or naive
    load(dist_file, 'euclidean_dists', 'euclidean_dff_dists');

    for stim = 1:5
        C_dff(stim, i, :, :) = squeeze(cropped_dff(stim, :, :));
        C_dprime(stim, i, :, :) = squeeze(cropped_dprime(stim, :, :));
        all_dprime_dists(i, stim) = euclidean_dists(stim);        % d′ distance
        all_dff_dists(i, stim) = euclidean_dff_dists(stim);       % df/F distance
    end
end

% Compute WT averages
avg_dff = squeeze(nanmean(C_dff, 2));
avg_dprime = squeeze(nanmean(C_dprime, 2));
avg_dprime_dists = mean(all_dprime_dists, 1);
avg_dff_dists = mean(all_dff_dists, 1);

% Plot 5 rows × 2 columns with distance text
stim_labels = {'Pre', 'A', 'B', 'C', 'D'};
figure;
for stim = 1:5
    % df/F plot
    subplot(5, 2, (stim-1)*2 + 1);
    imagesc(squeeze(avg_dff(stim, :, :)), [-0.01 0.02]);
    axis image; colormap('jet'); colorbar;
    title(['Cropped dF/F ' stim_labels{stim}]);
    text(80, 45, sprintf('Avg Euclidean\ndF/F: %.2f px', avg_dff_dists(stim)), ...
         'Color', 'black', 'FontSize', 9, 'FontWeight', 'bold', ...
         'BackgroundColor', 'white', 'VerticalAlignment', 'bottom');

    % d′ plot
    subplot(5, 2, (stim-1)*2 + 2);
    imagesc(squeeze(avg_dprime(stim, :, :)), [0 2]);
    axis image; colormap('jet'); colorbar;
    hold on;
    text(80, 45, sprintf('Avg Euclidean\nd′: %.2f px', avg_dprime_dists(stim)), ...
         'Color', 'white', 'FontSize', 9, 'FontWeight', 'bold', ...
         'BackgroundColor', 'black', 'VerticalAlignment', 'bottom');
    title(['Cropped d′ ' stim_labels{stim}]);
end
set(gcf, 'Position', [50 50 900 1000]);
sgtitle('WT Naive Average Cropped Maps – df/F and d′ with Euclidean Distances');

% === Separate figure for stimulus A only ===
stim_to_plot = find(strcmp(stim_labels, 'A')); % find index for A

figure;

% dF/F map
subplot(1, 2, 1);
imagesc(squeeze(avg_dff(stim_to_plot, :, :)), [-0.01 0.02]);
axis image;
set(gca, 'XTick', [], 'YTick', []);   % hide ticks but keep axes
box on;                               % now MATLAB draws the border
colormap('jet'); 
colorbar('eastoutside');
text(size(avg_dff,3)+25, size(avg_dff,2)/2, 'df/F', 'FontSize', 12, ...
     'Rotation', 90, 'HorizontalAlignment', 'center');

% d′ map
subplot(1, 2, 2);
imagesc(squeeze(avg_dprime(stim_to_plot, :, :)), [0 2]);
axis image;
set(gca, 'XTick', [], 'YTick', []);   % hide ticks but keep axes
box on;                               % now MATLAB draws the border
colormap('jet'); 
colorbar('eastoutside');
text(size(avg_dprime,3)+20, size(avg_dprime,2)/2, 'd-prime', 'FontSize', 12, ...
      'Rotation', 90, 'HorizontalAlignment', 'center');

set(gcf, 'Position', [100 400 700 300]);
sgtitle('WT Naive – Stimulus A Only');

%% FXS Mice - Combined Crop and Euclidean Distances

close all
clearvars
clc

% FXS mouse 
fxs_mice = {'OE48', 'OE49', 'OE50'};

% Initialize arrays
C_dff = NaN(5, length(fxs_mice), 50, 50);
C_dprime = NaN(5, length(fxs_mice), 50, 50);
all_dprime_dists = NaN(length(fxs_mice), 5);
all_dff_dists = NaN(length(fxs_mice), 5);

for i = 1:length(fxs_mice)
    mouse_id = fxs_mice{i};

    % Load cropped maps
    maps_file = sprintf('cropped_maps_mouse_%s_naive.mat', mouse_id);  % change if expert or naive
    load(maps_file, 'cropped_dff', 'cropped_dprime');

    % Load Euclidean distances
    dist_file = sprintf('euclidean_dists_%s_naive.mat', mouse_id); % change if expert or naive
    load(dist_file, 'euclidean_dists', 'euclidean_dff_dists');

    for stim = 1:5
        C_dff(stim, i, :, :) = squeeze(cropped_dff(stim, :, :));
        C_dprime(stim, i, :, :) = squeeze(cropped_dprime(stim, :, :));
        all_dprime_dists(i, stim) = euclidean_dists(stim);
        all_dff_dists(i, stim) = euclidean_dff_dists(stim);
    end
end

% Compute FXS averages
avg_dff = squeeze(nanmean(C_dff, 2));
avg_dprime = squeeze(nanmean(C_dprime, 2));
avg_dprime_dists = mean(all_dprime_dists, 1);
avg_dff_dists = mean(all_dff_dists, 1);

% Plot
stim_labels = {'Pre', 'A', 'B', 'C', 'D'};
figure;
for stim = 1:5

    % df/F plot
    subplot(5, 2, (stim-1)*2 + 1);
    imagesc(squeeze(avg_dff(stim, :, :)), [-0.01 0.02]);
    axis image; colormap('jet'); colorbar;
    title(['Cropped dF/F ' stim_labels{stim}]);
    text(80, 45, sprintf('Avg Euclidean\ndF/F: %.2f px', avg_dff_dists(stim)), ...
         'Color', 'black', 'FontSize', 9, 'FontWeight', 'bold', ...
         'BackgroundColor', 'white', 'VerticalAlignment', 'bottom');

    % d′ plot
    subplot(5, 2, (stim-1)*2 + 2);
    imagesc(squeeze(avg_dprime(stim, :, :)), [0 2]);
    axis image; colormap('jet'); colorbar;
    hold on;
    text(80, 45, sprintf('Avg Euclidean\nd′: %.2f px', avg_dprime_dists(stim)), ...
         'Color', 'white', 'FontSize', 9, 'FontWeight', 'bold', ...
         'BackgroundColor', 'black', 'VerticalAlignment', 'bottom');
    title(['Cropped d′ ' stim_labels{stim}]);
end
set(gcf, 'Position', [50 50 900 1000]);
sgtitle('FXS Naive Average Cropped Maps – df/F and d′ with Euclidean Distances');

% === Separate figure for stimulus A only ===
stim_to_plot = find(strcmp(stim_labels, 'A')); % find index for A

figure;

% dF/F map
subplot(1, 2, 1);
imagesc(squeeze(avg_dff(stim_to_plot, :, :)), [-0.01 0.02]);
axis image;
set(gca, 'XTick', [], 'YTick', []);   % hide ticks but keep axes
box on;                               % now MATLAB draws the border
colormap('jet'); 
colorbar('eastoutside');
text(size(avg_dff,3)+25, size(avg_dff,2)/2, 'df/F', 'FontSize', 12, ...
     'Rotation', 90, 'HorizontalAlignment', 'center');

% d′ map
subplot(1, 2, 2);
imagesc(squeeze(avg_dprime(stim_to_plot, :, :)), [0 2]);
axis image;
set(gca, 'XTick', [], 'YTick', []);   % hide ticks but keep axes
box on;                               % now MATLAB draws the border
colormap('jet'); 
colorbar('eastoutside');
text(size(avg_dprime,3)+20, size(avg_dprime,2)/2, 'd-prime', 'FontSize', 12, ...
      'Rotation', 90, 'HorizontalAlignment', 'center');

set(gcf, 'Position', [100 400 700 300]);
sgtitle('Fmr1 KO Naive – Stimulus A Only');
