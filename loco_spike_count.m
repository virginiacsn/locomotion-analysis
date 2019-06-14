function[trial_data] = loco_spike_count(trial_data)

for ipaw = 1:4
    for itrial = 1:length(trial_data)
        trial_data(itrial).imaging.loco_spike_count{ipaw} = sum(trial_data(itrial).imaging.cn.sp_thres.*...
            repmat(trial_data(itrial).imaging.loco_step{ipaw},1,size(trial_data(itrial).imaging.cn.sp_thres,2)));
    end
end
end