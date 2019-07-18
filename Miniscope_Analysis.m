%% Load file
root_path = 'C:\Users\Virginia Casasnovas\Documents\Virginia\Data\MINISCOPES\HeadFreeTreadmillMiniscope\';
% root_path = 'C:\Users\Virginia Casasnovas\Documents\Virginia\Data\MINISCOPES\tmp\';
addpath(genpath(root_path));
tracking_folder = 'TM TRACKING FILES';
processed_folder = 'TM PROCESSED FILES';
registered_folder = 'TM PROCESSED IMAGING FILES';
segmented_folder = 'TM SEGMENTED IMAGING FILES';
raw_folder = 'TM RAW FILES';

experiment_folder = '20190507 - miniscopes catch slow';
mouse_folder = 'MC2370';

trial = 18;
trial_folder = ['miniscope\T',num2str(trial)];

trial_video_data_folder = [raw_folder,'\',experiment_folder,'\',mouse_folder,'\',trial_folder,'\']; % Folder with imaging timestamp file (.dat)
trial_video_folder = [registered_folder,'\',experiment_folder,'\',mouse_folder,'\']; % Folder with registered imaging files (.tiff)
trial_neuron_folder = [segmented_folder,'\',experiment_folder,'\',mouse_folder,'\']; % Folder with edited segmented cell signal files (.mat)
trial_tracking_folder = [tracking_folder,'\',experiment_folder,'\',mouse_folder,'\']; % Folder with tracking files (.mat)
trial_processed_folder = [processed_folder,'\',experiment_folder,'\',mouse_folder,'\']; % Folder with processed tracking files (.mat)

video_data = importdata([trial_video_data_folder,'timestamp.dat']); % Loading imaging timestamp file
sysClock = video_data.data(:,3); % sysClock contains timestamps of miniscope video
missing_frames_idx = find(diff(sysClock)>=40); % Check where timestamps (frames) are missing

%% Loading and processing
% Load tiff
video_file = dir([root_path,trial_video_folder,mouse_folder,'_*_',num2str(trial),'_m_Ready.tiff']);
video_name = video_file.name;
info = imfinfo(video_name);

image_stack = [];
num_images = length(info);
for k = 1:num_images
    current_image = imread(video_name, k, 'Info', info);
    image_stack(:,:,k) = current_image;
end

if length(sysClock)~=num_images
    fprintf('\nNumber of frames in video file does not match the number of timestamps.\n');
else
    % Segment cells
    if exist([root_path,trial_neuron_folder,'neurons_',num2str(trial),'.mat'],'file') % Segmented file does not exist
        options.nPCs = 1000;
        options.used_PCs = 1:200; % If this value is empty the algorithm prompts the user for manual selection
        options.mu = 0.5;
        rng(1);
        
        cells_sort_file = [root_path,trial_video_folder,video_name];
        muk = get_mukamel_rois(cells_sort_file, [], options.nPCs, options.used_PCs, options.mu, []);
        
        % convert mukamel rois to carey neuron
        cn = convert_mukamel_to_carey_neuron(muk.ica_segments, muk.seg_centroid);
        
        [cn_data, ~, avg] = get_carey_neurons_mean_intensity(cells_sort_file, cn);
        cn.intensity = cn_data;
        
        % Deconvolute to get spikes
        ops.fs = 30;
        ops.recomputeKernel = 0;
        ops.sensorTau = 0.7;
        ops.estimateNeuropil = 1;
        ops.deconvType = 'L0';
        
        [sp, ca] = wrapperDECONV(ops, cn_data);
        cn.sp = sp;
        cn.sp_thres = sp>max(sp(:))*0.1;
        cn.ca = ca;
        
        % Saving mukamel and carey neurons
        if ~exist([root_path,trial_neuron_folder])
            mkdir([root_path,trial_neuron_folder]);
        end
        
        save([root_path,trial_neuron_folder,'neurons_',num2str(trial)],'cn','muk','options');
        
        % Plot segmentation
        figure('Name','Cell contours');
        imshow(image_stack(:,:,200),[]); hold on;
        for i = 1:cn.n_cells
            imcontour(cn.mask{i})
            hold on;
        end
    else
        % Load segmented cells
        load([root_path,trial_neuron_folder,'neurons_',num2str(trial)]);
        
        % Plot segmentation
        figure('Name','Cell contours');
        imshow(image_stack(:,:,200),[]); hold on;
        for i = 1:cn.n_cells_final
            imcontour(cn.mask_final{i})
            hold on;
        end
    end
end
   
%% Edit segmentation
[cn] = edit_segmentation_mukamel(cn,image_stack);

% fs = 30;
% cn.time = (0:size(cn.intensity_final,1)-1)/fs;
cn.time = sysClock/1000;
close all;

if exist([root_path,trial_neuron_folder,'neurons_',num2str(trial),'.mat'],'file')
    save_neurons = input('\nOverwrite neuron data files? [y/n] ','s');
    
    if strcmp(save_neurons,'y')
        save([root_path,trial_neuron_folder,'neurons_',num2str(trial)],'cn','muk','options');
    else
        fprintf('\nNot saving neuron data.\n')
    end
end

%% Plot final segmentation
figure('Name','Cell contours final');
imshow(image_stack(:,:,200),[]); hold on;
for i = 1:cn.n_cells_final
    imcontour(cn.mask_final{i})
    hold on;
end

figure('Name','Cell centroids final');
imshow(image_stack(:,:,200),[]); hold on;
for i = 1:cn.n_cells
    plot(cn.centroid{i}(1),cn.centroid{i}(2),'r*');
end
for i = 1:cn.n_cells_final
    plot(cn.centroid_final{i}(1),cn.centroid_final{i}(2),'b*');
    hold on;
end

figure('Name','Cell intensity final');
for i = 1:cn.n_cells_final
    plot(cn.time,cn.intensity_final(:,i)+20*i);
    hold on;
end
set(gca,'YTickLabel',[]);
xlabel('Time [s]');

%% Plot tiff
% figure;
% colormap(gray);
% 
% key = '';
% i = 0;
% 
% while (~strcmp(key, 'space'))&&(i<size(image_stack,3))
%     key = get(gcf,'currentkey');
%     i = i+1;
%     
%     imagesc(image_stack(:,:,i)); hold on;
% %   imcontour(I)
%     title(['Frame: ',num2str(i)]);
%     pause(0.0001);
% end

%% Tracking
tracking_file = dir([root_path,trial_tracking_folder,mouse_folder,'_*_',num2str(trial),'.mat']);
tracking_name = tracking_file.name;

tracking_data = load([root_path,trial_tracking_folder,tracking_name]);

processed_file = dir([root_path,trial_processed_folder,mouse_folder,'_*_',num2str(trial),'.mat']);
if ~isempty(processed_file) && ~isfield(tracking_data,'ptsraw')
    processed_name = processed_file.name;
    processed_data = load([root_path,trial_processed_folder,'\',processed_name]);
    tracking_data.time = (0:size(tracking_data.final_tracks,3)-1)/processed_data.trial_data.fps;
    tracking_data.strides = processed_data.strides;
end

tracking_data.st_sw_frames = get_st_sw_frames(tracking_data);

[loco_frame,loco_time,loco_step] = loco_epochs(tracking_data); % Extracting locomotion epochs
ipaw = 1;

%% Plot locomotion epochs with tracking
figure;
yl = [100 max(max(tracking_data.final_tracks(1,1:4,:)))+10];
for j = 1:size(loco_frame{ipaw},1)
    area([loco_frame{ipaw}(j,1) loco_frame{ipaw}(j,2)]/tracking_data.fps,[yl(2) yl(2)],'FaceColor',[0.9,0.8,0.9],'FaceAlpha',0.4,'EdgeColor','none');
    hold on;
end
plot(tracking_data.time,squeeze(tracking_data.final_tracks(1,1:4,:))')
ylim(yl);
xlabel('Time [s]'); ylabel('x [mm]');

%% Plot locomotion epochs with spikes
figure;
yl = [100 max(max(tracking_data.final_tracks(1,1:4,:)))+10];
for j = 1:size(loco_frame{ipaw},1)
    area([loco_frame{ipaw}(j,1) loco_frame{ipaw}(j,2)]/tracking_data.fps,[yl(2) yl(2)],'FaceColor',[0.9,0.8,0.9],'FaceAlpha',0.4,'EdgeColor','none');
    hold on;
end

plot(tracking_data.time,squeeze(tracking_data.final_tracks(1,1:4,:))')
offset = yl(1)+100;
for i = 1:cn.n_cells_final
    plot([cn.time(cn.sp_thres_final(:,i)) cn.time(cn.sp_thres_final(:,i))]',...
        repmat([0 5],size(cn.time(cn.sp_thres_final(:,i))))'+offset,'k')
    hold on;
    offset = offset+10;
end
ylim(yl);

%% Plot locomotion epochs with raw intensity, convoluted intensity and spikes
figure;
ylim([0 100*cn.n_cells_final]);
yl = ylim;
for j = 1:size(loco_frame{ipaw},1)
    area([loco_frame{ipaw}(j,1) loco_frame{ipaw}(j,2)]/330,[yl(2) yl(2)],'FaceColor',[0.9,0.8,0.9],'FaceAlpha',0.4,'EdgeColor','none');
    hold on;
end

offset = 20;
for i = 1:cn.n_cells_final
    plot(cn.time,cn.intensity_final(:,i)+offset,'b') 
    hold on;
    plot(cn.time,cn.ca_final(:,i)+offset,'k') 
    plot([cn.time(cn.sp_thres_final(:,i)) cn.time(cn.sp_thres_final(:,i))]',repmat([0 2],size(cn.time(cn.sp_thres_final(:,i))))'+offset-15,'k')
    offset = offset+max(cn.ca_final(:,i))+20;
    % plot(cn.time,cn.intensity_final(:,i)+offset);
    % offset = offset+max(cn.intensity_final(:,i))+10;
end
ylim([0 offset-10]);

% for i = 1:cn.n_cells_final
%     plot(cn.time,cn.intensity_final(:,i)+20*i);
%     hold on;
% end
set(gca,'YTickLabel',[]);
xlabel('Time [s]');

%% ISI
intersp_times = {};
for i = 1:cn.n_cells_final
    intersp_times{i} = diff(cn.time(cn.sp_thres_final(:,i)));
    figure;
    histogram(intersp_times{i},'BinWidth',0.2);
end

%% Figure with locomotion epochs and intensity, segmentation masks and ISI
for i = 1:cn.n_cells_final
    figure;
    subplot(2,2,3:4)
    
    ylim([0 100]); yl = ylim;
    for j = 1:size(loco_frame{ipaw},1)
        area([loco_frame{ipaw}(j,1) loco_frame{ipaw}(j,2)]/330,[yl(2) yl(2)],'FaceColor',[0.9,0.8,0.9],'FaceAlpha',0.4,'EdgeColor','none');
        hold on;
    end
    
    h1 = plot(cn.time,cn.intensity_final(:,i)+10,'b');
    hold on;
    h2 = plot(cn.time,cn.ca_final(:,i)+10,'k');
    plot([cn.time(cn.sp_thres_final(:,i)) cn.time(cn.sp_thres_final(:,i))]',repmat([0 2],size(cn.time(cn.sp_thres_final(:,i))))','k')
    ylim([0 max(cn.ca_final(:,i))+15]);
    legend([h1,h2],{'Raw','Deconv'})
    xlabel('Time [s]'); %title(['Imaging signal']);
    
    subplot(2,2,1)
    imshow(image_stack(:,:,200),[]);
    hold on;
    imcontour(cn.mask_final{i});
    %title('Contour')
       
    subplot(2,2,2)
    hst = histogram(intersp_times{i},'BinWidth',0.2);
    set(gca,'XTick',0:hst.BinLimits(2),'xlim',hst.BinLimits);
    xlabel('Time [s]'); %title('Histogram')
end

%% Trial data struct with: tracking and imaging fields
% Trials to be analysed
% trials = [2, 3, 5, 6, 8, 9, 11, 14, 17, 18];
% trials = [4, 5, 6, 9, 12, 15, 17, 18];
% trials = 12;
trials = [2, 3, 5, 6];

trial_data = [];

for itrial = 1:length(trials)
    
    trial = trials(itrial);
    trial_data(itrial).trial = trial;
    
    tracking_file = dir([root_path,trial_tracking_folder,mouse_folder,'_*_',num2str(trial),'.mat']);
    tracking_name = tracking_file.name;
    
    tracking_data = load([root_path,trial_tracking_folder,tracking_name]);
    
    processed_file = dir([root_path,trial_processed_folder,mouse_folder,'_*_',num2str(trial),'.mat']);
    if ~isempty(processed_file) && ~isfield(tracking_data,'ptsraw')
        processed_name = processed_file.name;
        processed_data = load([root_path,trial_processed_folder,processed_name]);
        trial_data(itrial).tracking.strides = processed_data.strides;
    end
    
    trial_data(itrial).speed_R = tracking_data.speed_R;
    trial_data(itrial).speed_L = tracking_data.speed_L;
    
    if tracking_data.speed_R == tracking_data.speed_L
        trial_data(itrial).trial_type = 'tied';
    else
        trial_data(itrial).trial_type = 'split';
    end
    
    trial_data(itrial).tracking.time = (0:size(tracking_data.final_tracks,3)-1)/tracking_data.fps;
    trial_data(itrial).tracking.st_sw_frames = get_st_sw_frames(tracking_data);
    
    trial_data(itrial).tracking.ptsraw = tracking_data.ptsraw;
    trial_data(itrial).tracking.final_tracks = tracking_data.final_tracks;
    [trial_data(itrial).tracking.loco_frame,trial_data(itrial).tracking.loco_time,trial_data(itrial).tracking.loco_step] = loco_epochs(trial_data(itrial).tracking);
    
    % Final cell segmentation signals saved in neurons_trialnumber.mat file
    trial_data(itrial).imaging = load([root_path,trial_neuron_folder,'neurons_',num2str(trial)]);
    trial_data(itrial).imaging.time = trial_data(itrial).imaging.cn.time;  
end

[trial_data] = match_loco_frames(trial_data);
[trial_data] = loco_spike_count(trial_data);
[trial_data] = match_img_tacking_frames(trial_data);

nbins = 10; % Number of bins for stride modulation
[trial_data] = stride_modulation(trial_data,nbins); % Stride modulation analysis

N = 30; % Set higher than maximum number of segmented cells
red = [0, 0, 1];
pink = [255, 192, 203]/255;
colors_p = [linspace(red(1),pink(1),N)', linspace(red(2),pink(2),N)', linspace(red(3),pink(3),N)'];

%% Tracking stride modulation
paw_lab = {'FR','HR','FL','HL'};
figure;
k = 0;
for itrial = 1:length(trial_data)
    for ipaw = 1:4
        k = k+1;
        subplot(length(trial_data),4,k);
        plot(trial_data(itrial).stride.stride_bins_mean{ipaw}(:,1),'*-'); hold on;
        %         plot(trial_data(itrial).stride_bins{ipaw}(:,1,1),'*-');
        errorbar(1:10,trial_data(itrial).stride.stride_bins_mean{ipaw}(:,1),trial_data(itrial).stride.stride_bins_sem{ipaw}(:,1),'.')
        if itrial == 1
            title(paw_lab{ipaw})
        end
        if ipaw == 1
            ylabel(trial_data(itrial).trial_type)
        end
    end
end

%% Imaging stride modulation 
paw_lab = {'FR','HR','FL','HL'};
figure;
k = 0;
for itrial = 1:length(trial_data)
    for ipaw = 1:4
        k = k+1;
        subplot(length(trial_data),4,k);
        for icell = 1:size(trial_data(itrial).stride.stride_bins_image_mean{ipaw},2)
            plot(trial_data(itrial).stride.stride_bins_image_mean{ipaw}(:,icell),'*-','Color',colors_p(icell,:));
            hold on;
            errorbar(1:10,trial_data(itrial).stride.stride_bins_image_mean{ipaw}(:,icell),trial_data(itrial).stride.stride_bins_image_sem{ipaw}(:,icell),'Color',colors_p(icell,:))
        end
        if itrial == 1
            title(paw_lab{ipaw})
        end
        if ipaw == 1
            ylabel(trial_data(itrial).trial_type)
        end
    end
end

%% Stride modulation - tracking and intensity for individual segmented cells
paw_lab = {'FR','HR','FL','HL'};
figure;
k = 0;
for itrial = 1:length(trial_data)
    for ipaw = 1:4
        k = k+1;
        subplot(length(trial_data),4,k);
        
        [hax,li1,li2] = plotyy(1:10,trial_data(itrial).stride.stride_bins_mean{ipaw}(:,1),1:10,trial_data(itrial).stride.stride_bins_image_mean{ipaw});
             
        for i = 1:2
            hax(i).XLim = ([0 11]);    
        end
        hold(hax(1), 'on');
        errorbar(hax(1), 1:10, trial_data(itrial).stride.stride_bins_mean{ipaw}(:,1),...
            trial_data(itrial).stride.stride_bins_sem{ipaw}(:,1),'k.');
        
        li1.LineWidth = 1.5;        
        li1.Color = 'k';
        li1.Marker = '.';
        
        hold(hax(2), 'on');
        for i = 1:length(li2)
            li2(i).Marker = '.';
            li2(i).Color = colors_p(i,:);
            errorbar(hax(2), 1:10, trial_data(itrial).stride.stride_bins_image_mean{ipaw}(:,i),...
                trial_data(itrial).stride.stride_bins_image_sem{ipaw}(:,i),'.','Color',colors_p(i,:));       
        end
           
        if itrial == 1
            title(paw_lab{ipaw})
        elseif itrial == length(trial_data)
            xlabel('% stride')
        end
        if ipaw == 1
            ylabel(hax(1),{['\bf\color{black}',trial_data(itrial).trial_type,' (',num2str(trial_data(itrial).trial),')']; ['\rmR = ',num2str(trial_data(itrial).speed_R)] ;[' L = ', num2str(trial_data(itrial).speed_L)];...
                ['\rm\color[rgb]{',num2str(hax(1).YColor),'}x [mm]']},'Interpreter','tex')
        elseif ipaw == 4
            ylabel(hax(2),{'area';'\DeltaF/F [-]'})
        end
    end
end

%% Stride modulation - tracking and mean intensity of segmented cells
paw_lab = {'FR','HR','FL','HL'};
figure;
k = 0;
for itrial = 1:length(trial_data)
    for ipaw = 1:4
        k = k+1;
        subplot(length(trial_data),4,k);
        
        [hax,li1,li2] = plotyy(1:10,trial_data(itrial).stride.stride_bins_mean{ipaw}(:,1),1:10,trial_data(itrial).stride.stride_bins_image_mean_cell{ipaw});
             
        for i = 1:2
            hax(i).XLim = ([0 11]);    
        end
        hold(hax(1), 'on');
        errorbar(hax(1), 1:10, trial_data(itrial).stride.stride_bins_mean{ipaw}(:,1),...
            trial_data(itrial).stride.stride_bins_sem{ipaw}(:,1),'k.');
        
        li1.LineWidth = 1.5;        
        li1.Color = 'k';
        li1.Marker = '.';
        
        hold(hax(2), 'on');
        for i = 1:length(li2)
            li2(i).Marker = '.';
            li2(i).Color = colors_p(i,:);
            errorbar(hax(2), 1:10, trial_data(itrial).stride.stride_bins_image_mean_cell{ipaw}(:,i),...
                trial_data(itrial).stride.stride_bins_image_sem_cell{ipaw}(:,i),'.','Color',colors_p(i,:));       
        end
           
        if itrial == 1
            title(paw_lab{ipaw})
        elseif itrial == length(trial_data)
            xlabel('% stride')
        end
        if ipaw == 1
            ylabel(hax(1),{['\bf\color{black}',trial_data(itrial).trial_type,' (',num2str(trial_data(itrial).trial),')']; ['\rmR = ',num2str(trial_data(itrial).speed_R)] ;[' L = ', num2str(trial_data(itrial).speed_L)];...
                ['\rm\color[rgb]{',num2str(hax(1).YColor),'}x [mm]']},'Interpreter','tex')
        elseif ipaw == 4
            ylabel(hax(2),{'area';'\DeltaF/F [-]'})
        end
    end
end

%% Plot all strides 
ipaw = 1;
figure;
for itrial = 1:length(trial_data)
    subplot(length(trial_data),1,itrial);
    for i = 1:length(trial_data(itrial).stride.stance_vals{ipaw})
        plot(trial_data(itrial).stride.stance_vals{ipaw}{i}(1,:)-(trial_data(itrial).stride.stance_vals{ipaw}{i}(1,1)))
        hold on;
    end
end

%% Plot tracking and spikes
figure;
for itrial = 1:length(trial_data)
    subplot(length(trial_data),1,itrial)
    ylim([0 100*trial_data(itrial).imaging.cn.n_cells_final]);
    yl = ylim;
    for j = 1:size(trial_data(itrial).tracking.loco_frame{ipaw},1)
        area([trial_data(itrial).tracking.loco_frame{ipaw}(j,1) trial_data(itrial).tracking.loco_frame{ipaw}(j,2)]/330,[yl(2) yl(2)],'FaceColor',[0.9,0.8,0.9],'FaceAlpha',0.4,'EdgeColor','none');
        hold on;
    end
    
    offset = 20;
    for i = 1:trial_data(itrial).imaging.cn.n_cells_final
        plot(trial_data(itrial).imaging.cn.time,trial_data(itrial).imaging.cn.intensity_final(:,i)+offset,'b')
        hold on;
        plot(trial_data(itrial).imaging.time,trial_data(itrial).imaging.cn.ca_final(:,i)+offset,'k')
        plot([trial_data(itrial).imaging.time(trial_data(itrial).imaging.cn.sp_thres_final(:,i)) ...
            trial_data(itrial).imaging.time(trial_data(itrial).imaging.cn.sp_thres_final(:,i))]',...
            repmat([0 2],size(trial_data(itrial).imaging.cn.time(trial_data(itrial).imaging.cn.sp_thres_final(:,i))))'+offset-15,'k')
        offset = offset+max(trial_data(itrial).imaging.cn.ca_final(:,i))+20;
    end
    ylim([0 offset-10]);
    set(gca,'YTickLabel',[]);
    title(['Trial ',num2str(trial_data(itrial).trial),': ',trial_data(itrial).trial_type, ' [R = ', ...
        num2str(trial_data(itrial).speed_R), '; L = ', num2str(trial_data(itrial).speed_L),']'])
    if itrial == length(trial_data)
        xlabel('Time [s]');
    end
end

%% Plot spike count
xlab = {};

figure;
for itrial = 1:length(trial_data)
    plot(itrial,mean(trial_data(itrial).imaging.loco_spike_count{ipaw}),'b*');
    hold on;
%     errorbar(itrial,mean(trial_data(itrial).loco_spike_count),...
%         std(trial_data(itrial).loco_spike_count)/sqrt(length(trial_data(itrial).loco_spike_count)),'b.')
    errorbar(itrial,mean(trial_data(itrial).imaging.loco_spike_count{ipaw}),...
        std(trial_data(itrial).imaging.loco_spike_count{ipaw}),'b.')
    xlab{itrial} = [trial_data(itrial).trial_type, '\newline[R = ', ...
        num2str(trial_data(itrial).speed_R), '; L = ', num2str(trial_data(itrial).speed_L),']'];
end
xlim([0.5 3.5]);
ylabel('Mean Spikes');
set(gca,'XTick',1:3);
% set(gca,'HorizontalAlignment', 'center');
set(gca,'XTickLabels',xlab,'TickLabelInterpreter', 'tex');

%% Plot intensity on top of tracks
figure;
for itrial = 1:length(trial_data)
    subplot(3,1,itrial)
    yl = [100 max(max(trial_data(itrial).tracking.final_tracks(1,1:4,:)))+10];
    for j = 1:size(trial_data(itrial).tracking.loco_frame{ipaw},1)
        area([trial_data(itrial).tracking.loco_frame{ipaw}(j,1) trial_data(itrial).tracking.loco_frame{ipaw}(j,2)]/330,[yl(2) yl(2)],'FaceColor',[0.9,0.8,0.9],'FaceAlpha',0.4,'EdgeColor','none');
        hold on;
    end
    plot(trial_data(itrial).tracking.time,squeeze(trial_data(itrial).tracking.final_tracks(1,1:4,:))');
    offset = yl(1)+350;
    
    for i = 1:trial_data(itrial).imaging.cn.n_cells_final
        %     plot([trial_data(itrial).imaging.cn.time(trial_data(itrial).imaging.cn.sp_thres_final(:,i)) trial_data(itrial).imaging.cn.time(trial_data(itrial).imaging.cn.sp_thres_final(:,i))]',...
        %         repmat([0 5],size(trial_data(itrial).imaging.cn.time(trial_data(itrial).imaging.cn.sp_thres_final(:,i))))'+offset,'k');
        
        plot(trial_data(itrial).imaging.cn.time,trial_data(itrial).imaging.cn.intensity_final(:,i)+offset,'b')
        
        hold on;
        offset = offset+10;
    end
    
    ylim(yl);
end

%% Trying 
% %% Select image ROI
% figure;
% imshow(image_stack(:,:,200),[]);
%     
% roi = imfreehand();
% roi_mask = createMask(roi);
% 
% %% Intensity signal
% roi_mask_mat = repmat(roi_mask,[1 1 size(image_stack,3)]);
% imageStack_mask = image_stack.*roi_mask_mat;
% imageStack_mask_mean = mean(imageStack_mask(:));
% roi_mean = (squeeze(mean(mean(imageStack_mask)))-imageStack_mask_mean)/imageStack_mask_mean;
% 
% % fc = 0.1; fs = 30;
% % [b,a] = butter(2,fc/(fs/2),'high');
% % roi_mean_filt = filtfilt(b,a,roi_mean);
% 
% ops.fs = 30;
% ops.recomputeKernel = 0;
% ops.sensorTau = 0.7;
% ops.estimateNeuropil = 1;
% ops.deconvType = 'L0';
% 
% [sp, ca] = wrapperDECONV(ops, roi_mean);
% 
% 
% imageStack_bkg_mask = image_stack.*(~roi_mask_mat);
% imageStack_proc = image_stack-imageStack_bkg_mask*mean(image_stack(:));
% imageStack_mask = imageStack_proc;%.*roi_mask_mat;
% imageStack_mask_mean = mean(imageStack_mask(:));
% roi_mean_proc = (squeeze(mean(mean(imageStack_mask)))-imageStack_mask_mean)/imageStack_mask_mean;
% 
% ops.sensorTau = 0.7;
% [sp, ca_proc] = wrapperDECONV(ops, roi_mean_proc);
% 
% % fc = 8; fs = 30;
% % [b,a] = butter(8,fc/(fs/2),'high');
% % roi_mean_filt = filtfilt(b,a,ca);
% 
% %% Plot
% figure;
% plot(roi_mean(2:end),'b')
% hold on;
% plot(roi_mean_proc(2:end),'r')
% % plot(ca(2:end),'r--')
% % plot(ca_proc(2:end),'b--')
% legend({'Raw','Bkg corr'})

