function[trial_data] = param_imaging(trial_data)

if ~isfield(trial_data,'intensity_bkg')
    error('Necessary fields missing.')
else
    for itrial = 1:length(trial_data)
        
        trial_data(itrial).int_det = detrend(trial_data(itrial).intensity,'constant');
        trial_data(itrial).mean_int_det = mean(trial_data(itrial).int_det,1);
        trial_data(itrial).area_int_det = trapz(trial_data(itrial).time,trial_data(itrial).int_det,1);
        trial_data(itrial).interp_int_det =  interp1(trial_data(itrial).time,trial_data(itrial).int_det,[1:19800]'/330);
        
        trial_data(itrial).int_bkg_det = detrend(trial_data(itrial).intensity_bkg,'constant');
        trial_data(itrial).mean_int_bkg_det = mean(trial_data(itrial).int_bkg_det,1);
        trial_data(itrial).area_int_bkg_det = trapz(trial_data(itrial).time,trial_data(itrial).int_bkg_det,1);
        trial_data(itrial).interp_int_bkg_det =  interp1(trial_data(itrial).time,trial_data(itrial).int_bkg_det,[1:19800]'/330);
    end
end