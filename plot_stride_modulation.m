function[] = plotStrideModulation(mouse_name,trial_data)
%% Stride modulation figs
N = 300; % Set higher than maximum number of segmented cells
red = [0, 0, 1];
pink = [255, 192, 203]/255;
colors_p = [linspace(red(1),pink(1),N)', linspace(red(2),pink(2),N)', linspace(red(3),pink(3),N)'];

nfig = floor(length(trial_data)/5);

paw_lab = {'FR','HR','FL','HL'};

%% Stride modulation - tracking and intensity for individual segmented cells
for ifig = 1:nfig
    figure('Name',mouse_name);
    k = 0;
    
    for itrial = (ifig-1)*5+1:(ifig*5)%length(trial_data)
        for ipaw = 1:4
            k = k+1;
            
            % subplot(length(trial_data),4,k);
            subplot(5,4,k);
            
            [hax,li1,li2] = plotyy(1:10,trial_data(itrial).tracking.stride_bins_mean{ipaw}(:,1),1:10,trial_data(itrial).imaging.stride.stride_bins_image_mean{ipaw});
            
            for i = 1:2
                hax(i).XLim = ([0 11]);
            end
            hold(hax(1), 'on');
            errorbar(hax(1), 1:10, trial_data(itrial).stride.stride_bins_mean{ipaw}(:,1),...
                trial_data(itrial).stride.stride_bins_sem{ipaw}(:,1),'k.');
            
            li1.LineWidth = 1.5;
            li1.Color = 'k';
            li1.Marker = '.';
            
            hold(hax(2), 'on');
            for i = 1:length(li2)
                li2(i).Marker = '.';
                li2(i).Color = colors_p(i,:);
                errorbar(hax(2), 1:10, trial_data(itrial).stride.stride_bins_image_mean{ipaw}(:,i),...
                    trial_data(itrial).stride.stride_bins_image_sem{ipaw}(:,i),'.','Color',colors_p(i,:));
            end
            
            if itrial == 1
                title(paw_lab{ipaw})
            elseif itrial == length(trial_data)
                xlabel('% stride')
            end
            if ipaw == 1
                ylabel(hax(1),{['\bf\color{black}',trial_data(itrial).trial_type,' (',num2str(trial_data(itrial).trial),')']; ['\rmR = ',num2str(trial_data(itrial).speed_R)] ;[' L = ', num2str(trial_data(itrial).speed_L)];...
                    ['\rm\color[rgb]{',num2str(hax(1).YColor),'}x [mm]']},'Interpreter','tex')
            elseif ipaw == 4
                ylabel(hax(2),{'area';'\DeltaF/F [-]'})
            end
        end
    end
end
%% Stride modulation - tracking and mean intensity of segmented cells

for ifig = 1:nfig
    figure('Name',mouse_name);
    k = 0;
    for itrial = (ifig-1)*5+1:(ifig*5)%length(trial_data)
        for ipaw = 1:4
            k = k+1;
            % subplot(length(trial_data),4,k);
            subplot(5,4,k);
            
            [hax,li1,li2] = plotyy(1:10,trial_data(itrial).stride.stride_bins_mean{ipaw}(:,1),1:10,trial_data(itrial).stride.stride_bins_image_mean_cell{ipaw});
            
            for i = 1:2
                hax(i).XLim = ([0 11]);
            end
            hold(hax(1), 'on');
            errorbar(hax(1), 1:10, trial_data(itrial).stride.stride_bins_mean{ipaw}(:,1),...
                trial_data(itrial).stride.stride_bins_sem{ipaw}(:,1),'k.');
            
            li1.LineWidth = 1.5;
            li1.Color = 'k';
            li1.Marker = '.';
            
            hold(hax(2), 'on');
            for i = 1:length(li2)
                li2(i).Marker = '.';
                li2(i).Color = colors_p(i,:);
                errorbar(hax(2), 1:10, trial_data(itrial).stride.stride_bins_image_mean_cell{ipaw}(:,i),...
                    trial_data(itrial).stride.stride_bins_image_sem_cell{ipaw}(:,i),'.','Color',colors_p(i,:));
            end
            
            if itrial == 1
                title(paw_lab{ipaw})
            elseif itrial == length(trial_data)
                xlabel('% stride')
            end
            if ipaw == 1
                ylabel(hax(1),{['\bf\color{black}',trial_data(itrial).trial_type,' (',num2str(trial_data(itrial).trial),')']; ['\rmR = ',num2str(trial_data(itrial).speed_R)] ;[' L = ', num2str(trial_data(itrial).speed_L)];...
                    ['\rm\color[rgb]{',num2str(hax(1).YColor),'}x [mm]']},'Interpreter','tex')
            elseif ipaw == 4
                ylabel(hax(2),{'area';'\DeltaF/F [-]'})
            end
        end
    end
end
end