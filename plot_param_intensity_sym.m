function[] = plot_param_intensity_sym(trial_data,param,col_spd)

% nspd = length(group_speed);
% lilliac = [76, 0, 153]/255;
% purple = [204, 153, 253]/255;
% spd_col = [linspace(lilliac(1),purple(1),nspd)', linspace(lilliac(2),purple(2),nspd)', linspace(lilliac(3),purple(3),nspd)'];

paw_leg = {'Front','Hind'};

[uni,ia,~] = unique(trial_data.output.trial.col_spd,'rows','stable');

figure;
m = 0; hl = []; lg = {}; memcol = zeros(length(ia),1);
for ipaw = 1:2
    for itrial = 1:length(trial_data.imaging)
        
        img_trial = trial_data.imaging(itrial).trial_num;
        
        % mean of trial for all neurons
        trial_intensity = mean(trial_data.imaging(itrial).img_vals_stride{ipaw}(:),1);
        
        % out param computed for non-clean locomotion
        out_trial = find(trial_data.output.trial.trial_num == img_trial);
        stride_param = trial_data.output.trial.sym.(param)(out_trial,ipaw);
        if col_spd
            trial_col =  trial_data.output.trial.col_spd(out_trial,:);
        else
            trial_col =  trial_data.output.trial.col_type(out_trial,:);
        end
        
        subplot(1,2,ipaw);
        
        if any(ismember(uni,trial_col,'rows')) && ipaw == 2 && memcol(ismember(uni,trial_col,'rows'))==0
            memcol = memcol + ismember(uni,trial_col,'rows');
            m = m+1;
            hs = plot(stride_param,trial_intensity,'o');
            hs.MarkerFaceColor = trial_col; hs.MarkerEdgeColor = trial_col;
            hl = [hl, hs];
            lg{m} = ['L: ',num2str(trial_data.output.trial.speed_L(out_trial)),'; R: ',num2str(trial_data.output.trial.speed_R(out_trial)),' m/s'];
        else
            scatter(stride_param,trial_intensity,20,trial_col,'o','Filled')
        end
        hold on;
        text(stride_param+0.001,trial_intensity-0.001,num2str(img_trial),'FontSize',8);
    end
    if ipaw == 2
        legend(hl,lg);
    end
    xlabel(strrep(param,'_',' '));
    ylabel('mean \DeltaF [-]');
    title(paw_leg{ipaw});
    axis square
end
end
