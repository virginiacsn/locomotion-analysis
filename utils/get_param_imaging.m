function[trial_data] = get_param_imaging(trial_data)

if ~isfield(trial_data,'intensity_bkg')
    error('Necessary fields missing.')
else
    for itrial = 1:length(trial_data)
        
        int_aux = trial_data(itrial).intensity;
        int_aux(1,:) = int_aux(2,:);
        trial_data(itrial).intensity = int_aux;
        
        int_bkg_aux = trial_data(itrial).intensity_bkg;
        int_bkg_aux(1,:) = int_bkg_aux(2,:);
        trial_data(itrial).intensity_bkg = int_bkg_aux;
        
        trial_data(itrial).int_det = detrend(trial_data(itrial).intensity,'constant');
        trial_data(itrial).mean_int_det = mean(trial_data(itrial).int_det,1);
        trial_data(itrial).area_int_det = trapz(trial_data(itrial).time,trial_data(itrial).int_det,1);
        
        trial_data(itrial).int_bkg_det = detrend(trial_data(itrial).intensity_bkg,'constant');
        trial_data(itrial).mean_int_bkg_det = mean(trial_data(itrial).int_bkg_det,1);
        trial_data(itrial).area_int_bkg_det = trapz(trial_data(itrial).time,trial_data(itrial).int_bkg_det,1);
        
        data_min = trial_data(itrial).intensity-mean(trial_data(itrial).intensity(:));
        trial_data(itrial).int_norm = data_min./std(trial_data(itrial).intensity(:));
        
        %         data_min = trial_data(itrial).intensity_bkg-mean(trial_data(itrial).intensity_bkg(:));
        %         trial_data(itrial).int_bkg_norm = data_min./std(trial_data(itrial).intensity_bkg(:));
        
%         ncells = size(trial_data(itrial).intensity_bkg,2);
        npts = size(trial_data(itrial).intensity_bkg,1);
        data_min = trial_data(itrial).intensity_bkg-repmat(mean(trial_data(itrial).intensity_bkg,1),npts,1);
        trial_data(itrial).int_bkg_norm = data_min./repmat(std(trial_data(itrial).intensity_bkg,[],1),npts,1);
        
        trial_data(itrial).mean_int_norm = mean(trial_data(itrial).int_norm,1);
        trial_data(itrial).area_int_norm = trapz(trial_data(itrial).time,trial_data(itrial).int_norm,1);
        
        trial_data(itrial).mean_int_bkg_norm = mean(trial_data(itrial).int_bkg_norm,1);
        trial_data(itrial).area_int_bkg_norm = trapz(trial_data(itrial).time,trial_data(itrial).int_bkg_norm,1);
    end
end