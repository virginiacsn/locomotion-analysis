function[trial_data] = match_img_tacking_frames(trial_data)
% Function to match imaging frames to tracking frames in order to get
% frames of stance and swing onset in imaging frame reference.
%
% Inputs: 
%   - trial_data: struct containing tracking_data and cn imaging structs.
% Outputs: 
%   Collected in trial_data.imaging substruct. Cell for each paw.
%   - stance_frame: Best match of imaging frame to tracking data stance onset frame.
%   - stance_time: stance_frame converted to timestamp.
%   - swing_frame: Best match of imaging frame to tracking data swing onset frame.
%   - swing_time: swing_frame converted to timestamp.

if length(trial_data) == 1
    imaging_data = trial_data.imaging;
    tracking_data = trial_data.tracking;
else
    error('Change trial_data format.');
end
    
for ipaw = 1:4
    for itrial_img = 1:length(imaging_data)
        stance_frame_aux = [];
        stance_time_aux = [];
        
        img_time = imaging_data(itrial_img).time;
        
        itrial_trk = find([tracking_data.trial_num] == imaging_data(itrial_img).trial_num);
        
        for istance = 1:size(tracking_data(itrial_trk).st_sw_frames{ipaw},1)
            % Find tracking timestamp for given stance onset frame
            track_time = tracking_data(itrial_trk).time(tracking_data(itrial_trk).st_sw_frames{ipaw}(istance,1));
            
            % Only if given tracking timestamp is smaller than last imaging
            % timestamp, as imaging recording could stop before the end of
            % the trial
            if track_time <= img_time(end)
                % Find smallest difference time index in imaging timestamp
                % vector between this and tracking stance onset timestamp
                [~,idx_min] = min(abs(img_time-track_time));
                
                stance_frame_aux(istance) = idx_min;
                stance_time_aux(istance) = img_time(idx_min);
            end
        end
        
        trial_data.imaging(itrial_img).stance_frame{ipaw} = stance_frame_aux;
        trial_data.imaging(itrial_img).stance_time{ipaw} = stance_time_aux;
        
        
        swing_frame_aux = [];
        swing_time_aux = [];
        
        for iswing = 1:size(tracking_data(itrial_trk).st_sw_frames{ipaw},1)
            track_time = tracking_data(itrial_trk).time(tracking_data(itrial_trk).st_sw_frames{ipaw}(iswing,2));
            
            if track_time <= img_time(end)
                [~,idx_min] = min(abs(img_time-track_time));
                
                swing_frame_aux(iswing) = idx_min;
                swing_time_aux(iswing) = img_time(idx_min);
            end
        end
        
        trial_data.imaging(itrial_img).swing_frame{ipaw} = swing_frame_aux;
        trial_data.imaging(itrial_img).swing_time{ipaw} = swing_time_aux;       
    end
end
end