function[trial_data] = stride_binning(trial_data,param,nbins)
% Function for stride modulation analysis.
%
% Inputs:
%   - trial_data: struct containing tracking_data and cn imaging structs.
%   - nbins: number of bins to group data.
% Outputs:
%   Collected in trial_data.stride substruct. Cell for each paw.
%   - stance_vals: cells with tracking data for each stride.
%   - stance_frames: cells with tracking data frames for each stride.
%
%   - stride_pts: mean tracking data points of each bin for each stride.
%
%   - stride_bins_mean: mean tracking data binned points across strides.
%   - stride_bins_sem: sem of mean tracking data binned points across strides. 

if isempty(param)
    param = 'final_tracks';
end

for itrial =  1:length(trial_data)
    for ipaw = 1:4
        k = 0;
        
        stance_vals = {};
        stance_frames = {};
    
        for istance = 1:size(trial_data(itrial).st_sw_frames{ipaw},1)-1
            % Frame of stance onset, stride start
            stride_ini = trial_data(itrial).st_sw_frames{ipaw}(istance,1);
            % Frame of next stance onset - 1, stride end
            stride_end = trial_data(itrial).st_sw_frames{ipaw}(istance+1,1)-1;
            
            % Include frames within locomotion epochs only and up to imaging frame end
            if any(stride_ini >= trial_data(itrial).loco_frame{ipaw}(:,1) & ...
                    stride_end <= trial_data(itrial).loco_frame{ipaw}(:,2))
                
                k = k+1;
                
                stance_vals{k} = squeeze(trial_data(itrial).(param)(:,ipaw,stride_ini:stride_end));
                stance_frames{k} = stride_ini:stride_end;

                % Find corresponding swing onset frame
                iswing = trial_data(itrial).st_sw_frames{ipaw}((trial_data(itrial).st_sw_frames{ipaw}(:,2)>stride_ini & trial_data(itrial).st_sw_frames{ipaw}(:,2)<stride_end),2);
            end
        end
        
        trial_data(itrial).stance_vals{ipaw} = stance_vals;
        trial_data(itrial).stance_frames{ipaw} = stance_frames;
        
        % Binning tracking data
        for istride = 1:length(stance_vals)
            
            stride_length = size(stance_vals{istride},2);
            pts_bin = floor(stride_length/nbins);
            
            % For each bin mean the tracking data points for that stride
            if pts_bin>0
                for ibin = 0:nbins-2
                    trial_data(itrial).stride_pts{ipaw}(ibin+1,:,istride) = mean(stance_vals{istride}(:,1+ibin*pts_bin:(ibin+1)*pts_bin),2)';
                end
                trial_data(itrial).stride_pts{ipaw}(ibin+2,:,istride) = mean(stance_vals{istride}(:,(ibin+1)*pts_bin:end),2)';
            else
                trial_data(itrial).stride_pts{ipaw}(:,:,istride) = zeros(nbins,4);
            end
        end
        
        for irow = 1:size(trial_data(itrial).stride_pts{ipaw},1)
            for icol = 1:size(trial_data(itrial).stride_pts{ipaw},2)
                
                trial_data(itrial).stride_bins_mean{ipaw}(irow,icol) = ...
                    mean(trial_data(itrial).stride_pts{ipaw}(irow,icol,~isnan(trial_data(itrial).stride_pts{ipaw}(irow,icol,:))));
                
                N = length(trial_data(itrial).stride_pts{ipaw}(irow,icol,~isnan(trial_data(itrial).stride_pts{ipaw}(irow,icol,:))));
                
                trial_data(itrial).stride_bins_sem{ipaw}(irow,icol) = ...
                    std(trial_data(itrial).stride_pts{ipaw}(irow,icol,~isnan(trial_data(itrial).stride_pts{ipaw}(irow,icol,:))))/sqrt(N);
            end
        end     
    end
end
end
