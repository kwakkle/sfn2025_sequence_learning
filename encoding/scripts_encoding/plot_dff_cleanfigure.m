%% 
clc
clear all
close all

% Wild type mice ------
% 
% load("OE12_240219_Stim_movie.mat") % 16 and 19 and 21
% load("OE12_240219_trial_history.mat")

load("OE15_240316_Stim_movie.mat") % 12 and 16 and 17
load("OE15_240316_trial_history.mat")
% 
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

%if reversal session set to 1

rev = 0;

ABCD_trial=1;
A_CD_trial=2;
AB_D_trial=3;
AAAA_trial=4;
DCBA_trial=5;

% num_frames = size(movie,3);
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
ABCD_avg_movie=zeros(rows,columns, num_per_stim); %This is the issue(?)
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

%% Encoding Analyses - dF/F only and labeled A, B, C, D

stim_frame = [12 17 22 27];  % Assuming these correspond to A, B, C, D

% Preallocate
encoding_frames = zeros(length(stim_frame), rows, columns, length(ABCD_trials));
for i = 1:length(stim_frame)
    temp = (ABCD_trials - 1) * num_per_stim + stim_frame(i);
    encoding_frames(i,:,:,:) = All_movie(:,:,temp) - All_movie(:,:,temp - 1);  % dF/F difference
end

% Define bins for histogram-based AUC estimate
hist_bins = min(encoding_frames(:)) : 0.001 : max(encoding_frames(:));

% Initialize encoding result
all_encoding = zeros(length(stim_frame), rows, columns);

% Compute encoding strength per pixel
for m = 1:length(stim_frame)
    for i = 1:rows
        for j = 1:columns
            stim_vals = reshape(encoding_frames(m, i, j, :), 1, []);
            not_stim_vals = reshape(encoding_frames(setdiff(1:length(stim_frame), m), i, j, :), 1, []);
            stim_hist = histcounts(stim_vals, hist_bins);
            not_stim_hist = histcounts(not_stim_vals, hist_bins);
            
            stim_sum = 4 * cumsum(stim_hist(end:-1:1));
            not_stim_sum = cumsum(not_stim_hist(end:-1:1));
            
            AUC = trapz(not_stim_sum, stim_sum) / (not_stim_sum(end) * stim_sum(end));
            all_encoding(m,i,j) = sqrt(2) * norminv(AUC);
        end
    end
end

% Plot encoding maps with labels A, B, C, D
figure;
colormap jet;
labels = {'A', 'B', 'C', 'D'};
for i = 1:length(stim_frame)
    subplot(length(stim_frame), 1, i);
    imagesc(squeeze(all_encoding(i, :, :)));
    axis image off;
    colorbar;
    caxis([0 2]);  % Adjust if needed
    text(-5, 20, labels{i}, 'FontSize', 16, 'FontWeight', 'bold');
end
set(gcf, 'Position', [100, 100, 500, 900]);
%%

figure;
labels = {'A', 'B', 'C', 'D'};

for i = 1:length(stim_frame)
    ax = subplot(length(stim_frame), 1, i);
    
    imagesc(ABCD_avg_movie(:, :, stim_frame(i)) - ABCD_avg_movie(:, :, stim_frame(i) - 1));
    colormap('jet');
    caxis([-0.01 0.02]);

    % BIG bold A–D labels
    text(-45, 70, labels{i}, 'FontSize', 30, 'FontWeight', 'bold');

    % Keep square plot, remove ticks but keep black border
    axis square;
    box on;
    xticks([1 size(ABCD_avg_movie,2)]);
    yticks([1 size(ABCD_avg_movie,1)]);
    set(gca, 'XTickLabel', {}, 'YTickLabel', {}, 'LineWidth', 1);

    % Customize colorbar
    c = colorbar;
    c.Ticks = [-0.01, 0, 0.01, 0.02];
    c.Label.String = 'df/F';
    c.Label.FontSize = 16;
    c.FontSize = 14;

    % Position subplot to center
    set(ax, 'Units', 'normalized');
    set(ax, 'Position', [0.25, 1 - i * 0.22, 0.5, 0.18]);
end

set(gcf, 'Position', [100, 100, 1000, 800]);
