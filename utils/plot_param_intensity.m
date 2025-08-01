function[] = plot_param_intensity(trial_data,param,plt_stride,plt_trial,col_spd)

paw_lab = {'FR','HR','FL','HL'};
paw_col = {'r','m','b','c'};
spd_mrk = {'+','o','*','v','x','<','>'};

% nspd = length(group_speed);
% lilliac = [76, 0, 153]/255;
% purple = [204, 153, 253]/255;
% spd_col = [linspace(lilliac(1),purple(1),nspd)', linspace(lilliac(2),purple(2),nspd)', linspace(lilliac(3),purple(3),nspd)'];

paw_order = [3,1,4,2];

if plt_stride
    figure;
    k = 0;
    for ipaw = paw_order
        k = k+1;
        for itrial = 1:length(trial_data.imaging)
            
            img_trial = trial_data.imaging(itrial).trial_num;
            
            stride_intensity = trial_data.imaging(itrial).img_vals_stride{ipaw};
            
            stride_idx = trial_data.tracking([trial_data.tracking.trial_num] == img_trial).stride_idx{ipaw};
            
            out_trial = find(trial_data.output.stride.trial_num{ipaw} == img_trial);
            stride_param = trial_data.output.stride.(param){ipaw}(out_trial(stride_idx));
            stride_col =  trial_data.output.stride.col_trial{ipaw}(out_trial(stride_idx),:);
            
            subplot(2,2,k)
            scatter(stride_param,mean(stride_intensity,2),15,stride_col)
            hold on;
        end
        if k == 4 || k == 3
        xlabel(strrep(param,'_',' '));
        end
        if k == 1 || k == 3
            ylabel('mean \DeltaF [-]');
        end
        title(paw_lab{ipaw});
    end
end
if plt_trial
    
    img_trials = sort([trial_data.imaging.trial_num]);
    [uni,ia,~] = unique(trial_data.output.trial.col_spd,'rows','stable');
    
    figure;
    k = 0;
    for ipaw = paw_order
        k = k+1; m = 0; hl = []; lg = {}; memcol = zeros(length(ia),1);
        for itrial = 1:length(trial_data.imaging)
            
            img_trial = trial_data.imaging(itrial).trial_num;
            
            % mean of trial for all neurons
            trial_intensity = mean(trial_data.imaging(itrial).img_vals_stride{ipaw}(:),1);
            
            % out param computed for non-clean locomotion
            out_trial = find(trial_data.output.trial.trial_num == img_trial);
            stride_param = trial_data.output.trial.(param)(out_trial,ipaw);
            if col_spd
                trial_col =  trial_data.output.trial.col_spd(out_trial,:);
            else
                trial_col =  trial_data.output.trial.col_type(out_trial,:);
            end
            
            subplot(2,2,k)            
            if any(ismember(uni,trial_col,'rows')) && ipaw == 1 && memcol(ismember(uni,trial_col,'rows'))==0 && col_spd
                memcol = memcol + ismember(uni,trial_col,'rows');
                m = m+1;
                hs = plot(stride_param,trial_intensity,'o');
                %hs.MarkerSize = 20;
                hs.MarkerFaceColor = trial_col; hs.MarkerEdgeColor = trial_col;
                hl = [hl, hs];
                lg{m} = ['L: ',num2str(trial_data.output.trial.speed_L(out_trial)),'; R: ',num2str(trial_data.output.trial.speed_R(out_trial)),' m/s'];
            else
                scatter(stride_param,trial_intensity,20,trial_col,'o','Filled')
            end
            hold on;
            text(stride_param+0.001,trial_intensity-0.001,num2str(img_trial),'FontSize',8);
        end
        if ipaw == 1 && col_spd
            legend(hl,lg);
        end
        if k == 4 || k == 3
        xlabel(strrep(param,'_',' '));
        end
        if k == 1 || k == 3
            ylabel('mean \DeltaF [-]');
        end
        title(paw_lab{ipaw});
    end
end
end

