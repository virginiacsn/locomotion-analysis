function[trial_data_output] = assign_col_speed(group_speed,trial_data_output)

blue = [0, 0, 1];
pink = [255, 192, 203]/255;
magenta = [1, 0, 1];
red = [1,0,0];
orange = [255, 165, 0]/255;
green = [0,1,0];

nspd = length(group_speed);
cols = [linspace(blue(1),red(1),nspd)', linspace(blue(2),red(2),nspd)', linspace(blue(3),red(3),nspd)'];
cols = [red; green; orange; magenta; blue; pink];

for ipaw = 1:4
    col_spd_stride{ipaw} = zeros(size(trial_data_output.stride.trial_num{ipaw},1),3);
end
col_spd_trial = zeros(size(trial_data_output.trial.trial_num,1),3);

for ispd = 1:nspd
    
    trials_spd = group_speed(ispd).trial_num;
    
    for ipaw = 1:4
        
        for itrial = 1:length(trials_spd)
            idx_trials = trial_data_output.stride.trial_num{ipaw}==trials_spd(itrial);
            col_spd_stride{ipaw}(idx_trials,:) = repmat(cols(ispd,:),length(find(idx_trials)),1);
        end
    end
    
    for itrial = 1:length(trials_spd)
        idx_trials = trial_data_output.trial.trial_num==trials_spd(itrial);
        col_spd_trial(idx_trials,:) = repmat(cols(ispd,:),length(find(idx_trials)),1);
    end
    
end

trial_data_output.stride.col_spd = col_spd_stride;
trial_data_output.trial.col_spd = col_spd_trial;


end