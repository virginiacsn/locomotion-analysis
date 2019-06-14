function[loco_frame,loco_time,loco_step] = loco_epochs(tracking_data)
% Function to identify clean locomotion epochs. 
%
% Inputs: 
%   - tracking_data: struct loaded from tracking data file. 
% Outputs: 
%   Cell for each paw.
%   - loco_frame: Nx2 matrix with frame number of locomotion epoch start
%   and end ([start end]).
%   - loco_time: loco_frame converted to timestamps.
%   - loco_step: Step signal with 1 when within locomotion epoch and 0
%   otherwise.

for ipaw = 1:4
    % Column 5 corresponding to frame number of stance onset point
    stance_pts_frame = tracking_data.st_sw_frames{ipaw}(:,1);
    
    loco_frame_aux = diff(stance_pts_frame);
    % Identifing break in locomotion if time between two consecutive
    % stance onset frames was longer than 40% more than the mean time
    % between two consecutive stance onset frames.
    loco_frame_idx = find(loco_frame_aux>1.4*mean(loco_frame_aux));
    
    % To get clean and long enough locomotion epochs, only consider these 
    % if containing more than 5 strides
    max_stride_inbet = 5;
    strides_inbet = diff(loco_frame_idx);
    stride_gap_idx_end = loco_frame_idx(find(strides_inbet>=max_stride_inbet)+1);
    stride_gap_idx_ini = loco_frame_idx(strides_inbet>=max_stride_inbet);
    
    loco_frame_ini_final = [stance_pts_frame(1); stance_pts_frame(stride_gap_idx_ini+1); stance_pts_frame(loco_frame_idx(end)+1)];
    loco_frame_end_final = [stance_pts_frame(loco_frame_idx(1)); stance_pts_frame(stride_gap_idx_end); stance_pts_frame(end)];
    if loco_frame_ini_final(end)~=loco_frame_end_final(end)
        loco_frame{ipaw} = [loco_frame_ini_final loco_frame_end_final];
    else
        loco_frame{ipaw} = [loco_frame_ini_final(1:end-1) loco_frame_end_final(1:end-1)];
    end
    
    loco_time{ipaw} = loco_frame{ipaw}*60/19800;

    loco_step_aux = zeros(size(tracking_data.final_tracks,3),1);
    for i = 1:size(loco_frame{ipaw},1)
        loco_step_aux(loco_frame{ipaw}(i,1):loco_frame{ipaw}(i,2)) = 1;
    end
    loco_step{ipaw} = loco_step_aux;
end
end