%% pre_stimulus_aligned_response_combined

% Naive

clear all 
clc
close all


WT_naive_resp_avg_combined = {
    'OE12_240216_naive_respAlign_ABCD_10back.mat'
    'OE15_240312_naive_respAlign_ABCD_10back.mat'
    'OE24_240913_naive_respAlign_ABCD_10back.mat'
    'OE35_250208_naive_respAlign_ABCD_10back.mat'
    'OE39_250209_naive_respAlign_ABCD_10back.mat'
    'OE40_250210_naive_respAlign_ABCD_10back.mat'
    'OE45_250319_naive_respAlign_ABCD_10back.mat'
    'OE47_250319_naive_respAlign_ABCD_10back.mat'
};


all_maps = [];
for i = 1:length(WT_naive_resp_avg_combined)         
    tmp = load(WT_naive_resp_avg_combined{i});        
    if isempty(all_maps)
        [r, c, f] = size(tmp.ABCD_resp_avg);
        all_maps = nan(r, c, f, length(WT_naive_resp_avg_combined));  
    end
    all_maps(:,:,:,i) = tmp.ABCD_resp_avg;
end

group_avg = mean(all_maps, 4, 'omitnan');

% plotting
n_back = size(group_avg, 3);         
image_time_res = 0.1;                 

figure
for k = 1:n_back
    subplot(3,4,k)
    imagesc(group_avg(:,:,k));        % <-- plot the combined average
    colormap('jet'); 
    colorbar
    axis image
    box on
    caxis([-0.01 0.02])
    rel_t = -(n_back - k + 0) * image_time_res;   % e.g., -1.0, -0.9, ..., -0.1
    title(sprintf('t = %.1fs', rel_t))
end
set(gcf,'Position',[10 10 900 700])

% === Visualization: only last pre-response frame ===
last_frame = n_back;  % this is the frame immediately before the lick (t = -0.1 s)

figure
imagesc(group_avg(:,:,last_frame));
axis image;
colormap('jet');
colorbar;
caxis([-0.01 0.02]);

rel_t = -(n_back - last_frame + 0) * image_time_res;  % compute relative time (e.g., -0.1 s)
title(sprintf('Wild Type Naive Frame before lick (t = %.1fs)', rel_t));

set(gcf, 'Position', [100 100 500 400]);  % optional: make window smaller



%% pre_stimulus_aligned_response_combined

% Expert

clear all 
clc
close all


WT_expert_resp_avg_combined = {
    'OE12_240913_expert_respAlign_ABCD_10back.mat'
    'OE15_240316_expert_respAlign_ABCD_10back.mat'
    'OE24_240921_expert_respAlign_ABCD_10back.mat'
    'OE35_250213_expert_respAlign_ABCD_10back.mat'
    'OE39_250213_expert_respAlign_ABCD_10back.mat'
    'OE40_250214_expert_respAlign_ABCD_10back.mat'
    'OE45_250323_expert_respAlign_ABCD_10back.mat'
    'OE46_250323_expert_respAlign_ABCD_10back.mat'
    'OE47_250325_expert_respAlign_ABCD_10back.mat'
};

all_maps = [];
for i = 1:length(WT_expert_resp_avg_combined)          % <-- fix var name
    tmp = load(WT_expert_resp_avg_combined{i});        % <-- fix var name
    if isempty(all_maps)
        [r, c, f] = size(tmp.ABCD_resp_avg);
        all_maps = nan(r, c, f, length(WT_expert_resp_avg_combined));  % <-- fix
    end
    all_maps(:,:,:,i) = tmp.ABCD_resp_avg;
end

group_avg = mean(all_maps, 4, 'omitnan');

% plotting
n_back = size(group_avg, 3);         
image_time_res = 0.1;                 

figure
for k = 1:n_back
    subplot(3,4,k)
    imagesc(group_avg(:,:,k));        % <-- plot the combined average
    colormap('jet'); 
    colorbar
    axis image
    box on
    caxis([-0.01 0.02])
    rel_t = -(n_back - k + 0) * image_time_res;   % e.g., -1.0, -0.9, ..., -0.1
    title(sprintf('t = %.1fs', rel_t))
end
set(gcf,'Position',[10 10 900 700])


% === Visualization: only last pre-response frame ===
last_frame = n_back;  % this is the frame immediately before the lick (t = -0.1 s)

figure
imagesc(group_avg(:,:,last_frame));
axis image;
colormap('jet');
colorbar;
caxis([-0.01 0.02]);

rel_t = -(n_back - last_frame + 0) * image_time_res;  % compute relative time (e.g., -0.1 s)
title(sprintf('Wild Type Expert Frame before lick (t = %.1fs)', rel_t));

set(gcf, 'Position', [100 100 500 400]);  % optional: make window smaller


%% pre_stimulus_aligned_response_combined

% fmr1ko Naive

clear all 
clc
close all


FMR1KO_naive_resp_avg_combined = {
    'OE48_250407_naive_respAlign_ABCD_10back.mat'
    'OE49_250407_naive_respAlign_ABCD_10back.mat'
    'OE50_250407_naive_respAlign_ABCD_10back.mat'
};


all_maps = [];
for i = 1:length(FMR1KO_naive_resp_avg_combined)         
    tmp = load(FMR1KO_naive_resp_avg_combined{i});        
    if isempty(all_maps)
        [r, c, f] = size(tmp.ABCD_resp_avg);
        all_maps = nan(r, c, f, length(FMR1KO_naive_resp_avg_combined));  
    end
    all_maps(:,:,:,i) = tmp.ABCD_resp_avg;
end

group_avg = mean(all_maps, 4, 'omitnan');

% plotting
n_back = size(group_avg, 3);         
image_time_res = 0.1;                 

figure
for k = 1:n_back
    subplot(3,4,k)
    imagesc(group_avg(:,:,k));        % <-- plot the combined average
    colormap('jet'); 
    colorbar
    axis image
    box on
    caxis([-0.01 0.02])
    rel_t = -(n_back - k + 0) * image_time_res;   % e.g., -1.0, -0.9, ..., -0.1
    title(sprintf('t = %.1fs', rel_t))
end
set(gcf,'Position',[10 10 900 700])

% === Visualization: only last pre-response frame ===
last_frame = n_back;  % this is the frame immediately before the lick (t = -0.1 s)

figure
imagesc(group_avg(:,:,last_frame));
axis image;
colormap('jet');
colorbar;
caxis([-0.01 0.02]);

rel_t = -(n_back - last_frame + 0) * image_time_res;  % compute relative time (e.g., -0.1 s)
title(sprintf('FMR1KO Naive Frame before lick (t = %.1fs)', rel_t));

set(gcf, 'Position', [100 100 500 400]);  % optional: make window smaller

%% pre_stimulus_aligned_response_combined

% fmr1ko Expert

clear all 
clc
close all


FMR1KO_expert_resp_avg_combined = {
    'OE48_250414_expert_respAlign_ABCD_10back.mat'
    'OE49_250412_expert_respAlign_ABCD_10back.mat'
    'OE50_250412_expert_respAlign_ABCD_10back.mat'
};

all_maps = [];
for i = 1:length(FMR1KO_expert_resp_avg_combined)          % <-- fix var name
    tmp = load(FMR1KO_expert_resp_avg_combined{i});        % <-- fix var name
    if isempty(all_maps)
        [r, c, f] = size(tmp.ABCD_resp_avg);
        all_maps = nan(r, c, f, length(FMR1KO_expert_resp_avg_combined));  % <-- fix
    end
    all_maps(:,:,:,i) = tmp.ABCD_resp_avg;
end

group_avg = mean(all_maps, 4, 'omitnan');

% plotting
n_back = size(group_avg, 3);         
image_time_res = 0.1;                 

figure
for k = 1:n_back
    subplot(3,4,k)
    imagesc(group_avg(:,:,k));        % <-- plot the combined average
    colormap('jet'); 
    colorbar
    axis image
    box on
    caxis([-0.01 0.02])
    rel_t = -(n_back - k + 0) * image_time_res;   % e.g., -1.0, -0.9, ..., -0.1
    title(sprintf('t = %.1fs', rel_t))
end
set(gcf,'Position',[10 10 900 700])

% === Visualization: only last pre-response frame ===
last_frame = n_back;  % this is the frame immediately before the lick (t = -0.1 s)

figure
imagesc(group_avg(:,:,last_frame));
axis image;
colormap('jet');
colorbar;
caxis([-0.01 0.02]);

rel_t = -(n_back - last_frame + 0) * image_time_res;  % compute relative time (e.g., -0.1 s)
title(sprintf('FMR1KO Expert Frame before lick (t = %.1fs)', rel_t));

set(gcf, 'Position', [100 100 500 400]);  % optional: make window smaller