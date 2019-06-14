function[TD] = assign_col_trial(TD)

if ~isfield(TD,'trial_num')
    error('Trial numbers not available in input struct.')
end

blue = [0, 0, 1];
pink = [255, 192, 203]/255;

red = [1,0,0];
orange = [255, 165, 0]/255;
green = [0,1,0];

if isempty(strfind([TD.trial_type{:,1}],'split'))
    trial_speeds = [TD.speed_L];
    unique_speeds = unique(trial_speeds);
    Ns = length(unique_speeds);
    cols_speed = [linspace(blue(1),red(1),Ns)', linspace(blue(2),red(2),Ns)', linspace(blue(3),red(3),Ns)'];
end

for ipaw = 1:size(TD.trial_num,2)
    
    if iscell(TD.trial_num)
        trial_num = TD.trial_num{ipaw};
    else
        trial_num = TD.trial_num;
    end
    
    Nt = max(trial_num);
    cols = [linspace(blue(1),pink(1),Nt)', linspace(blue(2),pink(2),Nt)', linspace(blue(3),pink(3),Nt)'];
    
    %TD.col_trial{ipaw} = zeros(length(trial_num),1);
    for i = 1:Nt
        ind_col = find(trial_num == i);
        
        if iscell(TD.trial_num) % Stride data
            if strcmp(TD.trial_type{ipaw}(ind_col(1)),'split')
                type_col = red;
            else
                if ind_col(1) < Nt/2
                    type_col = green;
                else
                    type_col = orange;
                end
            end
            TD.col_trial{ipaw}(ind_col,:) = repmat(cols(i,:),length(ind_col),1);
            TD.col_type{ipaw}(ind_col,:) = repmat(type_col,length(ind_col),1);
                        
            if isempty(strfind([TD.trial_type{:,1}],'split'))
                stride_speeds = TD.speed_L{ipaw};
                TD.col_speed{ipaw}(ind_col,:) = repmat(cols_speed(unique_speeds == stride_speeds(ind_col(1)),:),length(ind_col),1);
            end
        else
            if strcmp(TD.trial_type(ind_col(1)),'split')
                type_col = red;
            else
                if ind_col(1) < Nt/2
                    type_col = green;
                else
                    type_col = orange;
                end
            end
            TD.col_trial(ind_col,:) = repmat(cols(i,:),length(ind_col),1);
            TD.col_type(ind_col,:) = repmat(type_col,length(ind_col),1);
            
            if isempty(strfind([TD.trial_type{:,1}],'split'))
                TD.col_speed(ind_col,:) = cols_speed(unique_speeds == trial_speeds(ind_col),:);
            end
        end
    end
end
end