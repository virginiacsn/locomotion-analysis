function[trial_data] = stance_phase_intensity(trial_data)

tracking_data = trial_data.tracking;
imaging_data = trial_data.imaging;
output_data = trial_data.output.stride;

for itrial_img = 1:length(imaging_data)
    %     itrial_trk = find([tracking_data.trial_num] == imaging_data(itrial_img).trial_num);
    nstrides = min(length(imaging_data(itrial_img).img_vals{1}),length(imaging_data(itrial_img).img_vals{3}));
    
    for istride = 1:nstrides-1
        idx_paw = find(imaging_data(itrial_img).img_frames{1}{istride}==imaging_data(itrial_img).img_frames{3}{istride}(1));
        if isempty(idx_paw)
            idx_paw = find(imaging_data(itrial_img).img_frames{1}{istride}==imaging_data(itrial_img).img_frames{3}{istride+1}(1));
        end
        if ~isempty(idx_paw)
            imaging_data(itrial_img).img_vals_stance_mean(istride,:) = mean(imaging_data(itrial_img).img_vals{1}{istride}(1:idx_paw,:),1);
            imaging_data(itrial_img).img_vals_stance_area(istride,:) = trapz(imaging_data(itrial_img).img_frames{1}{istride}(1:idx_paw),imaging_data(itrial_img).img_vals{1}{istride}(1:idx_paw,:),1);
        else
            imaging_data(itrial_img).img_vals_stance_mean(istride,:) = nan(1,size(imaging_data(itrial_img).img_vals{1}{istride},2));
            imaging_data(itrial_img).img_vals_stance_area(istride,:) = nan(1,size(imaging_data(itrial_img).img_vals{1}{istride},2));
        end
    end
end

trial_data.imaging = imaging_data;
end