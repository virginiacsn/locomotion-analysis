%% Load file
root_path = 'C:\Users\Virginia Casasnovas\Documents\Virginia\Data\MINISCOPES\HeadFreeTreadmillMiniscope\';
addpath(genpath(root_path));
raw_folder = 'TM RAW FILES';
tracking_folder = 'TM TRACKING FILES';
processed_folder = 'TM PROCESSED FILES';
registered_folder = 'TM PROCESSED IMAGING FILES';
segmented_folder = 'TM SEGMENTED IMAGING FILES';

experiment_folder = '20190507 - miniscopes catch slow';

%% Loading and processing
mouse_folder = dir([root_path,registered_folder,'\',experiment_folder,'\']);
mouse_folder = mouse_folder(3:end);

for imouse = 1:length(mouse_folder)
    seg_files =  rdir([root_path,segmented_folder,'\',experiment_folder,'\',mouse_folder(imouse).name,'\','*_seg.mat']);
    
    trial_data = [];
    k = 0;
    for itrial = 1:length(seg_files)
        
        cn = load(seg_files(itrial).name,'cn');
        if ~isempty(cn.cn.intensity)
            k = k+1;
            trial_str = strsplit(seg_files(itrial).name,{'_','.'});
            trial = str2num(trial_str{end-3});
            
            trial_data(k).trial = trial;
            
            tracking_file = dir([root_path,tracking_folder,'\',experiment_folder,'\',mouse_folder(imouse).name,'\*_',num2str(trial),'.mat']);
            tracking_name = tracking_file.name;
            
            tracking_data = load([root_path,tracking_folder,'\',experiment_folder,'\',mouse_folder(imouse).name,'\',tracking_name]);
            
            processed_file = dir([root_path,processed_folder,'\',experiment_folder,'\',mouse_folder(imouse).name,'_*_',num2str(trial),'.mat']);
            if ~isempty(processed_file) && ~isfield(tracking_data,'ptsraw')
                processed_name = processed_file.name;
                processed_data = load([root_path,processed_folder,'\',experiment_folder,'\',processed_name]);
                trial_data(k).tracking.strides = processed_data.strides;
            end
            
            trial_data(k).speed_R = processed_data.trial_data.speed_R;
            trial_data(k).speed_L = processed_data.trial_data.speed_L;
            
            if processed_data.trial_data.speed_R == processed_data.trial_data.speed_L
                trial_data(k).trial_type = 'tied';
            else
                trial_data(k).trial_type = 'split';
            end
            
            trial_data(k).tracking.time = (0:size(tracking_data.final_tracks,3)-1)/processed_data.trial_data.fps;
            trial_data(k).tracking.st_sw_frames = get_st_sw_frames(processed_data);
            
            %     trial_data(itrial).tracking.ptsraw = tracking_data.ptsraw;
            trial_data(k).tracking.final_tracks = tracking_data.final_tracks;
            [trial_data(k).tracking.loco_frame,trial_data(k).tracking.loco_time,trial_data(k).tracking.loco_step] = loco_epochs(trial_data(k).tracking);
            
            % Final cell segmentation signals saved in neurons_trialnumber.mat file
            trial_data(k).imaging = load(seg_files(itrial).name);
            if ~isfield(trial_data(k).imaging.cn,'time')
                timestamp = importdata([root_path,raw_folder,'\',experiment_folder,'\',mouse_folder(imouse).name,'\miniscope\T',num2str(trial),'\timestamp.dat']); % Loading imaging timestamp file
                trial_data(k).imaging.cn.time = timestamp.data(:,3)/1000; % sysClock contains timestamps of miniscope video
            end
            trial_data(k).imaging.time = trial_data(k).imaging.cn.time;
        end
    end
    
    [trial_data] = match_loco_frames(trial_data);
    [trial_data] = loco_spike_count(trial_data);
    [trial_data] = match_img_tacking_frames(trial_data);
    
    nbins = 10; % Number of bins for stride modulation
    [trial_data] = stride_modulation(trial_data,nbins); % Stride modulation analysis
    
    [~,sort_idx] = sort([trial_data.trial]);
    trial_data = trial_data(sort_idx);
    
    plotStrideModulation(mouse_folder(imouse).name,trial_data);
end

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