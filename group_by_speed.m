function[group_speed] = group_by_speed(trial_data,group_speed,unique_spd,trial_spd,params)


for iparam = 1:length(params)
    for ispd = 1:size(unique_spd,1)
        group_speed(ispd,:).spd = unique_spd(ispd,:);
        group_speed(ispd).trial_num = trial_spd{ispd};
        
        if iscell(trial_data(1).trial_num)
            trial_num = [trial_data.trial_num{1}];
        else
            trial_num = [trial_data.trial_num];
        end
        trials = find(ismember(trial_num,trial_spd{ispd}));
        
        for itrial = trials
            if ~isempty(trial_data(itrial).(params{iparam}))
                itrial_field = itrial;
            end
        end
        
        if iscell(trial_data(itrial_field).(params{iparam}))
            for ipaw = 1:4
                if iscell(trial_data(1).trial_num)
                    trial_num = [trial_data.trial_num{ipaw}];
                else
                    trial_num = [trial_data.trial_num];
                end
                
                trials = find(ismember(trial_num,trial_spd{ispd}));
                
                for ivar = 1:size(trial_data(itrial_field).(params{iparam}){ipaw},2)
                    varb = [];
                    for itrial = trials
                        if ~isempty(trial_data(itrial).(params{iparam}))
                            varb = [varb, trial_data(itrial).(params{iparam}){ipaw}(:,ivar)];
                        end
                    end
                    group_speed(ispd).(params{iparam}){ipaw}{ivar} = varb;
                    group_speed(ispd).([params{iparam},'_spd_mean']){ipaw}(:,ivar) = mean(varb,2);
                    group_speed(ispd).([params{iparam},'_spd_sem']){ipaw}(:,ivar) = std(varb,[],2)/sqrt(size(varb,2));
                end
            end
        else
            
            trial_num = [trial_data.trial_num];
            trials = find(ismember(trial_num,trial_spd{ispd}));
            
            for ivar = 1:size(trial_data(itrial_field).(params{iparam}),2)
                varb = [];
                for itrial = trials
                    if ~isempty(trial_data(itrial).(params{iparam}))
                        varb = [varb, trial_data(itrial).(params{iparam})(:,ivar)];
                    end
                end
                group_speed(ispd).(params{iparam}){ivar} = varb;
                group_speed(ispd).([params{iparam},'_spd_mean'])(:,ivar) = mean(varb,2);
                group_speed(ispd).([params{iparam},'_spd_sem'])(:,ivar) = std(varb,[],2)/sqrt(size(varb,2));
            end
        end
    end
end