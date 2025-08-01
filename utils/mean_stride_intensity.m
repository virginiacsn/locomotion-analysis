function[trial_data] = mean_stride_intensity(trial_data)
% trial_data is trial_data.imaging struct

if ~isfield(trial_data,'img_vals')
    error('Necessary fields missing.')
else
    %     [~,tr_idx] = sort([trial_data.trial_num]);
    for ipaw = 1:4
        for itrial = 1:length(trial_data)
            for istride = 1:length(trial_data(itrial).img_vals{ipaw})
                trial_data(itrial).img_vals_stride{ipaw}(istride,:) = mean(trial_data(itrial).img_vals{ipaw}{istride},1);
            end
        end
    end
end
end
