function[TD] = group_by_speed(TD,trk_params,img_params,out_params)

speeds = cat(1,[TD.tracking.speed_L],[TD.tracking.speed_R])';
uni_speeds = unique(speeds,'rows');

for ispd = 1:size(uni_speeds,1)
    
    TD.group_speed(ispd).speed_L = uni_speeds(ispd,1);
    TD.group_speed(ispd).speed_R = uni_speeds(ispd,2);
    
    % Tracking
    if ~isempty(trk_params)
        trk_trials = [TD.tracking.trial_num];
        trk_trial_spd = trk_trials([TD.tracking.speed_L]==uni_speeds(ispd,1)&[TD.tracking.speed_R]==uni_speeds(ispd,2));
        TD.group_speed(ispd).tracking.trial_num = trk_trial_spd;
        
        for iparam = 1:length(trk_params)
            
            k = 0;
            for itrial = trk_trial_spd
                if ~isempty(TD.tracking(itrial).(trk_params{iparam}{1}))
                    k = k+1;
                    for ipaw = 1:4
                        if length(trk_params{iparam}) < 2
                            TD.group_speed(ispd).tracking.(trk_params{iparam}){ipaw}(:,k) = TD.tracking(itrial).(trk_params{iparam}){ipaw}(:,1);
                        else
                            TD.group_speed(ispd).tracking.(trk_params{iparam}{2}){ipaw}(:,k) = TD.tracking(itrial).(trk_params{iparam}{1}).(trk_params{iparam}{2}){ipaw}(:,1);
                        end
                    end
                end
                
            end
            for ipaw = 1:4
                TD.group_speed(ispd).tracking.([trk_params{iparam}{end},'_mean']){ipaw} = mean(TD.group_speed(ispd).tracking.(trk_params{iparam}{end}){ipaw},2);
                TD.group_speed(ispd).tracking.([trk_params{iparam}{end},'_sem']){ipaw} = std(TD.group_speed(ispd).tracking.(trk_params{iparam}{end}){ipaw},[],2)/sqrt(length(trk_trial_spd));
            end
        end
        
        % Imaging
        if ~isempty(img_params)
            img_trials = [TD.imaging.trial_num];
            img_trial_spd = find(ismember(img_trials,trk_trial_spd));
            
            for iparam = 1:length(img_params)
                
                k = 0;
                for itrial = img_trial_spd
                    k = k+1;
                    if length(img_params{iparam}) < 2
                        TD.group_speed(ispd).imaging.(img_params{iparam})(:,k) = TD.imaging(itrial).(img_params{iparam});
                    else
                        if iscell(TD.imaging(itrial).(img_params{iparam}{1}).(img_params{iparam}{2}))
                            for ipaw = 1:4
                                TD.group_speed(ispd).imaging.(img_params{iparam}{2}){ipaw}(:,k) = TD.imaging(itrial).(img_params{iparam}{1}).(img_params{iparam}{2}){ipaw}(:,1);
                            end
                        else
                            TD.group_speed(ispd).imaging.(img_params{iparam}{2})(:,k) = TD.imaging(itrial).(img_params{iparam}{1}).(img_params{iparam}{2})(:,1);
                        end
                    end
                end
                if iscell(TD.imaging(itrial).(img_params{iparam}{1}).(img_params{iparam}{2}))
                    for ipaw = 1:4
                        TD.group_speed(ispd).imaging.([img_params{iparam}{end},'_mean']){ipaw} = mean(TD.group_speed(ispd).imaging.(img_params{iparam}{end}){ipaw},2);
                        TD.group_speed(ispd).imaging.([img_params{iparam}{end},'_sem']){ipaw} = std(TD.group_speed(ispd).imaging.(img_params{iparam}{end}){ipaw},[],2)/sqrt(length(img_trial_spd));
                    end
                else
                    TD.group_speed(ispd).imaging.([img_params{iparam}{end},'_mean']) = mean(TD.group_speed(ispd).imaging.(img_params{iparam}{end}),2);
                    TD.group_speed(ispd).imaging.([img_params{iparam}{end},'_sem']) = std(TD.group_speed(ispd).imaging.(img_params{iparam}{end}),[],2)/sqrt(length(img_trial_spd));
                end
            end
        end
        
        % Output.trial
        if ~isempty(out_params)
            out_trials = [TD.output.trial.trial_num];
            out_trial_spd = out_trials([TD.output.trial.speed_L]==uni_speeds(ispd,1)&[TD.output.trial.speed_R]==uni_speeds(ispd,2));
            for iparam = 1:length(out_params)
                
                k = 0;
                for itrial = out_trial_spd
                    k = k+1;
                    for ipaw = 1:4
                        if length(out_params{iparam}) < 2
                            TD.group_speed(ispd).output.trial.(out_params{iparam}){ipaw}(:,k) = TD.output.trial(itrial).(out_params{iparam})(:,ipaw);
                        else
                            TD.group_speed(ispd).output.trial.(out_params{iparam}{2})(:,k) = TD.output.trial(itrial).(out_params{iparam}{1}).(out_params{iparam}{2})(:,1);
                            
                        end
                    end
                end
                for ipaw = 1:4
                    TD.group_speed(ispd).output.trial.([out_params{iparam}{end},'_mean']){ipaw} = mean(TD.group_speed(ispd).output.trial.(out_params{iparam}{end}){ipaw},2);
                    TD.group_speed(ispd).output.trial.([out_params{iparam}{end},'_sem']){ipaw} = std(TD.group_speed(ispd).output.trial.(out_params{iparam}{end}){ipaw},[],2)/sqrt(length(out_trial_spd));
                end
            end
        end
    end
end