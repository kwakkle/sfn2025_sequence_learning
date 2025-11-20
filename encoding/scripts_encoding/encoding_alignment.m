%% Encoding Alignment Code  

clc
clear all
close all

% Wild type mice ------


load("OE12_240219_Stim_movie.mat") % 16 and 19 and 21
load("OE12_240219_trial_history.mat")

% load("OE15_240316_Stim_movie.mat") % 12 and 16 and 17
% load("OE15_240316_trial_history.mat")
 
% load("OE24_240921_Stim_movie.mat") % 13 and 21 and 24
% load("OE24_240921_trial_history.mat")

% load("OE35_250213_Stim_movie.mat") % 08 and 13 and 14
% load("OE35_250213_trial_history.mat")

% load("OE39_250213_Stim_movie.mat") % 09 and 13 and 14 
% load("OE39_250213_trial_history.mat")
 
% load("OE40_250214_Stim_movie.mat") % 10 and 14 and 15
% load("OE40_250214_trial_history.mat")

% load("OE45_250323_Stim_movie.mat") % 19 and 23 and 24
% load("OE45_250323_trial_history.mat")

% load("OE46_250323_Stim_movie.mat") % 19 and 23 and 24
% load("OE46_250323_trial_history.mat")

% load("OE47_250325_Stim_movie.mat") % 19 and 25 and 26
% load("OE47_250325_trial_history.mat")

% FXS Mice -------

% load("OE48_250414_Stim_movie.mat") % 07 and 14 and 15
% load("OE48_250414_trial_history.mat")

% load("OE49_250412_Stim_movie.mat") % 07 and 12 and 14
% load("OE49_250412_trial_history.mat")

% load("OE50_250412_Stim_movie.mat") % 07 and 12 and 14 
% load("OE50_250412_trial_history.mat")


% load (
%if reversal session set to 1

rev = 0;

ABCD_trial=1;
A_CD_trial=2;
AB_D_trial=3;
AAAA_trial=4;
DCBA_trial=5;

num_frames = 6000;
num_per_stim=30;
num_trials=num_frames/num_per_stim
image_time_res=0.1;

% low pass filter movie
hsize=[7, 7];
sigma = 3;
t=fspecial('gaussian',hsize,sigma);
movieLP=zeros(size(movie));

for i=1:num_frames
    movieLP(:,:,i)=imfilter(movie(:,:,i),t);
end

prestim_F=zeros(size(movieLP(:,:,num_trials),1),size(movieLP(:,:,num_trials),2),num_trials);

for i=1:num_trials
    start=num_per_stim*(i-1)+1; %%% TO UNDO DIFF, CHANGE 2 to 1
    range=start:start+9; %%% TO UNDO DIFF, CHANGE 8 to 9
    prestim_F(:,:,i)=mean(movieLP(:,:,range),3); %%% TO UNDO DIFF, CHANGE HERE
end
    
prestim_dFF=zeros(size(movieLP(:,:,num_trials),1),size(movieLP(:,:,num_trials),2),num_trials.*num_per_stim);
index_trials=ceil((1:num_trials*num_per_stim)./num_per_stim);

for i=1:size(prestim_dFF,3)
    prestim_dFF(:,:,i)=(movieLP(:,:,i)-prestim_F(:,:,index_trials(i)))./prestim_F(:,:,index_trials(i)); %%% TO UNDO DIFF, CHANGE HERE
end

%% Define trial types

All_trials = 1:num_trials;
rows=size(movie,1);
columns=size(movie,2);

All_movie=zeros(rows,columns,num_per_stim*length(All_trials));
All_indiv_movie=zeros(rows,columns,num_per_stim,length(All_trials));

for i=1:length(All_trials)
    for j=1:num_per_stim        
        All_movie(:,:,num_per_stim*(i-1)+j)=prestim_dFF(:,:,num_per_stim*(All_trials(i)-1)+j);
    end
    All_indiv_movie(:,:,:,i)=prestim_dFF(:,:,num_per_stim*(All_trials(i)-1)+1:num_per_stim*(All_trials(i)-1)+num_per_stim);
end

% implay(R_movie,5)

All_avg_movie=zeros(rows,columns,num_per_stim);
for i=1:num_per_stim
    frames=(0:length(All_trials)-1)*num_per_stim+i;
    All_avg_movie(:,:,i)=mean(All_movie(:,:,frames),3);
end

% ABCD trials only

ABCD_trials=find(trial_history(1,1:200)==1);
ABCD_avg_movie=zeros(rows,columns, num_per_stim); 
for i=1:num_per_stim
    temp_frames=(0:length(All_trials)-1)*num_per_stim+i;
    frames=temp_frames(ABCD_trials);
    ABCD_avg_movie(:,:,i)=mean(All_movie(:,:,frames),3);
end

figure
for i=1:num_per_stim
    subplot(6,5,i)%change 5 to 4
    imagesc(ABCD_avg_movie(:,:,i))
    colormap('jet')
    colorbar
    % caxis([-0.05 0.05])
    caxis([-0.01 0.02])
    title(['frame', num2str(i)])
end

% % AAAA_trials only
% 
% %AAAA_trials=find(event_history(1,:)==4);
% AAAA_trials=find(trial_history(1,:)==4);
% AAAA_avg_movie=zeros(rows,columns, num_per_stim); %This is the issue
% for i=1:num_per_stim
%     temp_frames=(0:length(All_trials)-1)*num_per_stim+i;
%     frames=temp_frames(AAAA_trials);
%     AAAA_avg_movie(:,:,i)=mean(All_movie(:,:,frames),3);
% end
% 
% figure
% title('AAAA');
% for i=1:num_per_stim
%     subplot(6,5,i)%change 5 to 4
%     imagesc(AAAA_avg_movie(:,:,i))
%     colormap('jet')
%     colorbar
%   %  caxis([-0.05 0.05])
%     caxis([-0.01 0.02])
%     title(['frame', num2str(i)])
% end
% set(gcf,'Position',[10 10 1104 847])
% 
% % savefig([ 'Stim_movie_AAAA.fig']);
% % saveas(gcf,[ 'Stim_movie_AAAA.png']);
% 
% % DCBA_trials only
% 
% %DCBA_trials=find(event_history(1,:)==5);
% DCBA_trials=find(trial_history(1,:)==5);
% DCBA_avg_movie=zeros(rows,columns, num_per_stim); %This is the issue
% for i=1:num_per_stim
%     temp_frames=(0:length(All_trials)-1)*num_per_stim+i;
%     frames=temp_frames(DCBA_trials);
%     DCBA_avg_movie(:,:,i)=mean(All_movie(:,:,frames),3);
% end
% 
% figure
% title('DCBA');
% for i=1:num_per_stim
%     subplot(6,5,i)%change 5 to 4
%     imagesc(DCBA_avg_movie(:,:,i))
%     colormap('jet')
%     colorbar
% %    caxis([-0.05 0.05])
%      caxis([-0.01 0.02])
%     title(['frame', num2str(i)])
% end
% set(gcf,'Position',[10 10 1104 847])
% 
% % savefig([ 'Stim_movie_DCBA.fig']);
% % saveas(gcf,[ 'Stim_movie_DCBA.png']);

mask=0;

%% Encoding analyses

stim_frame=[7 12 17 22 27]

figure
for i=1:length(stim_frame)
    subplot(length(stim_frame),2,(i-1)*2+1)
    imagesc(ABCD_avg_movie(:,:,stim_frame(i))-ABCD_avg_movie(:,:,stim_frame(i)-1))
    colormap('jet')
    colorbar
    caxis([-0.01 0.02])
    title(['frame', num2str(stim_frame(i))])
end
set(gcf,'Position',[10 10 1104 847])

encoding_frames=zeros(length(stim_frame),rows,columns,length(ABCD_trials));
for i=1:length(stim_frame)
    temp=(ABCD_trials-1)*num_per_stim+stim_frame(i);
    encoding_frames(i,:,:,:)=All_movie(:,:,temp)-All_movie(:,:,temp-1);
end

hist_bins=(min(min(min(min(encoding_frames)))):0.001:max(max(max(max(encoding_frames)))));

all_encoding=zeros(length(stim_frame),rows,columns);

for m=1:length(stim_frame)
    for i=1:rows
        for j=1:columns
            stim_hist=histcounts(reshape(encoding_frames(m,i,j,:),1,[]),hist_bins);
            not_stim_hist=histcounts(reshape(encoding_frames(setdiff(1:length(stim_frame),m),i,j,:),1,[]),hist_bins);
            stim_sum=4*cumsum(stim_hist(end:-1:1));
            not_stim_sum=cumsum(not_stim_hist(end:-1:1));
            AUC=trapz(not_stim_sum,stim_sum)./(not_stim_sum(end).*stim_sum(end));
            all_encoding(m,i,j)=sqrt(2).*norminv(AUC);
        end
    end
end


for i=1:length(stim_frame)
    subplot(length(stim_frame),2,i*2)
    imagesc(squeeze(all_encoding(i,:,:)))
    colormap('jet')
    colorbar
    caxis([0 2])
    title(['frame', num2str(stim_frame(i))])
end
set(gcf,'Position',[10 10 600 900])

%% Creating Cropped images 

stim_frame = [7 12 17 22 27];
stim_labels = {'Pre', 'A', 'B', 'C', 'D'};

[rows, cols, ~] = size(ABCD_avg_movie);

cropped_dff = NaN(5, 50, 50);
cropped_dprime = NaN(5, 50, 50);

% mask for sensory area 
mask_rect = false(rows, cols);
mask_rect(:, 1:50) = true;

figure;
for i = 1:5
    % full maps
    dff_diff = ABCD_avg_movie(:,:,stim_frame(i)) - ABCD_avg_movie(:,:,stim_frame(i)-1);
    dprime_map = squeeze(all_encoding(i,:,:));

    % detect brightest pixel detection within mask
    dff_masked = dff_diff; dff_masked(~mask_rect) = -Inf;
    dprime_masked = dprime_map; dprime_masked(~mask_rect) = -Inf;
% notes: -Inf used to mask out regions not wanted in peak search

    [y_dff, x_dff] = find(dff_masked == max(dff_masked(:)), 1);
    [y_dprime, x_dprime] = find(dprime_masked == max(dprime_masked(:)), 1);

% find -> brightest pixel (MAX) within masked area 

    % df/F cropping

    half_size = 25;
    x_min = x_dff - half_size + 1; x_max = x_dff + half_size;
    y_min = y_dff - half_size + 1; y_max = y_dff + half_size;
    
    x_in = max(1, x_min):min(cols, x_max);
    y_in = max(1, y_min):min(rows, y_max);
    x_out = x_in - x_min + 1;
    y_out = y_in - y_min + 1;
    
    temp_crop_dff = NaN(50, 50);
    temp_crop_dff(y_out, x_out) = dff_diff(y_in, x_in);
    cropped_dff(i,:,:) = temp_crop_dff;

    % d' cropping
    x_min = x_dprime - half_size + 1; x_max = x_dprime + half_size;
    y_min = y_dprime - half_size + 1; y_max = y_dprime + half_size;
    
    x_in = max(1, x_min):min(cols, x_max);
    y_in = max(1, y_min):min(rows, y_max);
    x_out = x_in - x_min + 1;
    y_out = y_in - y_min + 1;
    
    temp_crop_dprime = NaN(50, 50);
    temp_crop_dprime(y_out, x_out) = dprime_map(y_in, x_in);
    cropped_dprime(i,:,:) = temp_crop_dprime;

% plotting -------------

    % plot full df/F
    subplot(5, 4, (i-1)*4 + 1);
    imagesc(dff_diff, [-0.01 0.02]); axis image; colormap('jet'); colorbar;
    hold on; plot(x_dff, y_dff, 'wx', 'MarkerSize', 10, 'LineWidth', 2);
    title(['Full dF/F ' stim_labels{i}]);

    % plot cropped df/F
    subplot(5, 4, (i-1)*4 + 2);
    imagesc(squeeze(cropped_dff(i,:,:)), [-0.01 0.02]);
    axis image; colormap('jet'); colorbar;
    title(['Cropped dF/F ' stim_labels{i}]);

    % plot full d'
    subplot(5, 4, (i-1)*4 + 3);
    imagesc(dprime_map, [0 2]); axis image; colormap('jet'); colorbar;
    hold on; plot(x_dprime, y_dprime, 'wx', 'MarkerSize', 10, 'LineWidth', 2);
    title(['Full d′ ' stim_labels{i}]);

    % crop full d'
    subplot(5, 4, (i-1)*4 + 4);
    imagesc(squeeze(cropped_dprime(i,:,:)), [0 2]);
    axis image; colormap('jet'); colorbar;
    title(['Cropped d′ ' stim_labels{i}]);
end

set(gcf, 'Position', [50 50 1600 900]);


%% Saving for averaging later

% save('cropped_maps_mouse_OE24_expert.mat', 'cropped_dprime', 'cropped_dff');
% saveas(gcf, 'cropped_maps_mouse_OE24_expert.png', 'cropped_dprime', 'cropped_dff');

