function[trial_data] = norm_intensity(trial_data)

npts = size(trial_data.int_det,1);
data_min = trial_data.int_det-repmat(min(trial_data.int_det,[],1),npts,1);
trial_data.int_norm = data_min./repmat(max(data_min,[],1),npts,1);

end