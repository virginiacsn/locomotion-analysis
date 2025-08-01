function[trial_data_tracking] = get_body_kin(trial_data_tracking)

fs = 5;
[b,a] = butter(3, fs/330*2,'low');

for itrial = 1:length(trial_data_tracking)
    
    body_pos = squeeze(mean(trial_data_tracking(itrial).final_tracks(1,1:4,:)));
    time_pos = trial_data_tracking(itrial).time;
    
    time_pos(isnan(body_pos)) = [];
    body_pos(isnan(body_pos)) = [];
    
    body_pos_interp = interp1(time_pos,body_pos,trial_data_tracking(itrial).time,'spline');
    trial_data_tracking(itrial).body_pos = filtfilt(b,a,body_pos_interp);   
    
    dt = trial_data_tracking(itrial).time(2)-trial_data_tracking(itrial).time(1);
    
    trial_data_tracking(itrial).body_vel = filtfilt(b,a,diff(body_pos_interp)/dt); 
    trial_data_tracking(itrial).body_acc = filtfilt(b,a,diff(trial_data_tracking(itrial).body_vel)/dt);
end

end