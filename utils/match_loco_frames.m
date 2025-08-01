function[trial_data] = match_loco_frames(trial_data)
% Function to match imaging frames to tracking frames in order to get
% frames of locomotion epoch start and end in imaging frame reference.
%
% Inputs: 
%   - trial_data: struct containing tracking_data and cn imaging structs.
% Outputs: 
%   Collected in trial_data.imaging substruct. Cell for each paw.
%   - loco_frame: Nx2 matrix with frame number of locomotion epoch start
%   and end ([start end]).
%   - loco_time: loco_frame converted to timestamps.
%   - loco_step: Step signal with 1 when within locomotion epoch and 0
%   otherwise.
if length(trial_data) == 1
    imaging_data = trial_data.imaging;
    tracking_data = trial_data.tracking;
else
    error('Change trial_data format.');
end
    
for ipaw = 1:4
    for itrial_img = 1:length(imaging_data)
        loco_frame_aux = [];
        loco_time_aux = [];
        loco_step_aux = zeros(size(imaging_data(itrial_img).time,1),1);
                    
        img_time = imaging_data.time;
        
        itrial_trk = find([tracking_data.trial_num] == imaging_data(itrial_img).trial_num);
        
        for itime = 1:size(tracking_data(itrial_trk).loco_time{ipaw},1)
            track_time = tracking_data(itrial_trk).loco_time{ipaw}(itime,:);
            
            if track_time(2) <= img_time(end)
                idx_min = [];
                for i = 1:length(track_time)
                    [~,idx_min(i)] = min(abs(img_time-track_time(i)));
                end
                
                loco_frame_aux(itime,:) = idx_min;
                loco_time_aux(itime,:) =  img_time(idx_min);
                
                loco_step_aux(idx_min(1):idx_min(2)) = 1;
            end
        end
        trial_data.imaging(itrial_img).loco_frame{ipaw} = loco_frame_aux;
        trial_data.imaging(itrial_img).loco_time{ipaw} = loco_time_aux;
        trial_data.imaging(itrial_img).loco_step{ipaw} = loco_step_aux;
        
    end
end
end