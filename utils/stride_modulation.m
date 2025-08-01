function[trial_data] = stride_modulation(trial_data,nbins,signal)
% Function for stride modulation analysis.
%
% Inputs:
%   - trial_data: struct containing tracking_data and cn imaging structs.
%   - nbins: number of bins to group data.
% Outputs:
%   Collected in trial_data.stride substruct. Cell for each paw.
%   - stride_vals: cells with tracking data for each stride.
%   - stride_frames: cells with tracking data frames for each stride.
%
%   - img_vals: cells with imaging raw intensity data for each stride.
%   - img_frames: cells with imaging data frames for each stride.
%   - img_vals_interp: cells with interpolated imaging data for each
%   stride (in order to bin imaging data, more data points are 
%   needed, artificially creating them).
%   - img_frames_interp: cells with extended vector of imaging data frames 
%   for each stride.
%
%   - stride_pts_bins: mean tracking data points of each bin for each stride.
%   - img_pts_bins: mean imaging data points of each bin for each stride.
%
%   - stride_bins_mean: mean tracking data binned points across strides.
%   - stride_bins_sem: sem of mean tracking data binned points across strides. 
%
%   - img_bins_mean: mean binned imaging data across strides.
%   - img_bins_sem: sem of mean binned imaging data across strides.
%   - img_bins_mean_cell: mean binned imaging data across cells.
%   - img_bins_sem_cell: sem of mean binned imaging data across cells.

if length(trial_data) == 1
    imaging_data = trial_data.imaging;
    tracking_data = trial_data.tracking;
else
    error('Change trial_data format.');
end
    
for itrial_img =  1:length(imaging_data)
    for ipaw = 1:4
        k = 0;
        
        stride_idx = [];
        stride_vals = {};
        stride_frames = {};
        
        img_frames = {};
        img_vals = {};
        img_frames_interp = {};
        img_vals_interp = {};
        
%         median_phase_times = [];
                
        itrial_trk = find([tracking_data.trial_num] == imaging_data(itrial_img).trial_num);

        for istance = 1:size(tracking_data(itrial_trk).st_sw_frames{ipaw},1)-1
            % Frame of stance onset, stride start
            stride_ini = tracking_data(itrial_trk).st_sw_frames{ipaw}(istance,1);
            % Frame of next stance onset - 1, stride end
            stride_end = tracking_data(itrial_trk).st_sw_frames{ipaw}(istance+1,1)-1;
            
            % Include frames within locomotion epochs only and up to imaging frame end
            if any(stride_ini >= tracking_data(itrial_trk).loco_frame{ipaw}(:,1) & ...
                    stride_end <= tracking_data(itrial_trk).loco_frame{ipaw}(:,2) & ...
                    istance < length(imaging_data(itrial_img).stance_frame{ipaw}))
                
                k = k+1;
                
                stride_idx(k) = istance;
                
                stride_vals{k} = squeeze(tracking_data(itrial_trk).final_tracks(:,ipaw,stride_ini:stride_end));
                stride_frames{k} = stride_ini:stride_end;
                
                img_frames{k} = imaging_data(itrial_img).stance_frame{ipaw}(istance):imaging_data(itrial_img).stance_frame{ipaw}(istance+1);
                img_vals{k} = imaging_data(itrial_img).(signal)(img_frames{k},:);
                
%                 fprintf('\n trial: %d, paw: %d, stance: %d',itrial,ipaw,istance);
            
                % Make sure image frames are not a single number, else
                % create artificial vectors
                if img_frames{k}(1) ~= img_frames{k}(end)
                    img_frames_interp{k} = img_frames{k}(1):0.1:img_frames{k}(end);
                    img_vals_interp{k} = interp1(img_frames{k},img_vals{k},img_frames_interp{k});
                    if size(img_vals_interp{k},1) == 1
                       img_vals_interp{k} = img_vals_interp{k}';
                    end
                else
                    img_frames_interp{k} = zeros(1,20);
                    img_vals_interp{k} = zeros(20,size(imaging_data(itrial_img).(signal),2));
                end
                
                % Find corresponding swing onset frame
                iswing = tracking_data(itrial_trk).st_sw_frames{ipaw}((tracking_data(itrial_trk).st_sw_frames{ipaw}(:,2)>stride_ini & tracking_data(itrial_trk).st_sw_frames{ipaw}(:,2)<stride_end),2);
                % Store number of frames for stance and swing
%                 median_phase_times(k,:) = [iswing-stride_ini stride_end-iswing];
            end
        end
        
        trial_data.tracking(itrial_trk).stride_idx{ipaw} = stride_idx;
        trial_data.tracking(itrial_trk).stride_vals{ipaw} = stride_vals;
        trial_data.tracking(itrial_trk).stride_frames{ipaw} = stride_frames;
        
        trial_data.imaging(itrial_img).img_vals{ipaw} = img_vals;
        trial_data.imaging(itrial_img).img_frames{ipaw} = img_frames;
        
        trial_data.imaging(itrial_img).img_vals_interp{ipaw} = img_vals_interp;
        trial_data.imaging(itrial_img).img_frames_interp{ipaw} = img_frames_interp;
        
        % Median number of frames for stance and swing across strides
%         trial_data.tracking(itrial_trk).median_phase{ipaw} = median(median_phase_times,1);
        % trial_data(itrial).stride_ratio{ipaw} = trial_data(itrial).median_phase{ipaw}(1)/trial_data(itrial).median_phase{ipaw}(2);
        
        % Number of bins that correspond to stance
%         trial_data.tracking(itrial_trk).nstance_bins{ipaw} = round(trial_data.imaging(itrial_img).median_phase{ipaw}(1),-1);
        
        % Binning tracking data
        for istride = 1:length(stride_vals)
            
            stride_length = size(stride_vals{istride},2);
            pts_bin = floor(stride_length/nbins);
            
            % For each bin mean the tracking data points for that stride
            if pts_bin>0
                for ibin = 0:nbins-2
                    trial_data.tracking(itrial_trk).stride_pts_bins{ipaw}(ibin+1,:,istride) = mean(stride_vals{istride}(:,1+ibin*pts_bin:(ibin+1)*pts_bin),2)';
                end
                trial_data.tracking(itrial_trk).stride_pts_bins{ipaw}(ibin+2,:,istride) = mean(stride_vals{istride}(:,1+(ibin+1)*pts_bin:end),2)';
            else
                trial_data.tracking(itrial_trk).stride_pts_bins{ipaw}(:,:,istride) = zeros(nbins,4);
            end
        end
        
        for irow = 1:size(trial_data.tracking(itrial_trk).stride_pts_bins{ipaw},1)
            for icol = 1:size(trial_data.tracking(itrial_trk).stride_pts_bins{ipaw},2)
                
                trial_data.tracking(itrial_trk).stride_bins_mean{ipaw}(irow,icol) = ...
                    mean(trial_data.tracking(itrial_trk).stride_pts_bins{ipaw}(irow,icol,~isnan(trial_data.tracking(itrial_trk).stride_pts_bins{ipaw}(irow,icol,:))));
                
                N = length(trial_data.tracking(itrial_trk).stride_pts_bins{ipaw}(irow,icol,~isnan(trial_data.tracking(itrial_trk).stride_pts_bins{ipaw}(irow,icol,:))));
                
                trial_data.tracking(itrial_trk).stride_bins_sem{ipaw}(irow,icol) = ...
                    std(trial_data.tracking(itrial_trk).stride_pts_bins{ipaw}(irow,icol,~isnan(trial_data.tracking(itrial_trk).stride_pts_bins{ipaw}(irow,icol,:))))/sqrt(N);
            end
        end
        
        % Binning imaging data
        for istride = 1:length(img_frames_interp)
            stride_length = size(img_frames_interp{istride},2);
            pts_bin = floor(stride_length/nbins);
            stride_time = (0:pts_bin-1)';
            
            % For each bin caculate the area over the curve of the raw
            % intensity signal using trapz for that stride
            for ibin = 0:nbins-2
                stride_bin_vals = trial_data.imaging(itrial_img).img_vals_interp{ipaw}{istride}(1+ibin*pts_bin:(ibin+1)*pts_bin,:);
                trial_data.imaging(itrial_img).img_pts_bins{ipaw}(ibin+1,:,istride) = trapz(stride_time,stride_bin_vals,1);
                % trial_data(itrial).img_pts_bins{ipaw}(ibin+1,:,istride) = mean(stride_bin_vals);
            end
            stride_time = (1:stride_length-(ibin+1)*pts_bin)';
            stride_bin_vals = trial_data.imaging(itrial_img).img_vals_interp{ipaw}{istride}(1+(ibin+1)*pts_bin:end,:);
            trial_data.imaging(itrial_img).img_pts_bins{ipaw}(ibin+2,:,istride) = trapz(stride_time,stride_bin_vals,1);
            % trial_data(itrial).img_pts_bins{ipaw}(ibin+2,:,istride) = mean(stride_bin_vals);
        end
        
        trial_data.imaging(itrial_img).img_bins_mean{ipaw} = mean(trial_data.imaging(itrial_img).img_pts_bins{ipaw},3);
        trial_data.imaging(itrial_img).img_bins_sem{ipaw} = std(trial_data.imaging(itrial_img).img_pts_bins{ipaw},[],3)/sqrt(size(trial_data.imaging(itrial_img).img_pts_bins{ipaw},3));
        
        trial_data.imaging(itrial_img).img_bins_mean_cell{ipaw} = mean(trial_data.imaging(itrial_img).img_bins_mean{ipaw},2);
        trial_data.imaging(itrial_img).img_bins_sem_cell{ipaw} = std(trial_data.imaging(itrial_img).img_bins_mean{ipaw},[],2)/sqrt(size(trial_data.imaging(itrial_img).img_bins_mean{ipaw},2));       
    end
end
end
