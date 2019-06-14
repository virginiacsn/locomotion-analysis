function[] = plot_stride_md(trial_data,itrial,plt_cell,plt_mean_cell)

paw_lab = {'FR','HR','FL','HL'};
paw_col = {'r','m','b','c'};

ncells = size(trial_data.imaging([trial_data.imaging.trial_num] == itrial).intensity,2);
% green = [102, 0, 54]/255; pink = [255, 192, 203]/255;
lilliac = [76, 0, 153]/255;
purple = [204, 153, 253]/255;
cell_col = [linspace(lilliac(1),purple(1),ncells)', linspace(lilliac(2),purple(2),ncells)', linspace(lilliac(3),purple(3),ncells)'];

paw_order = [3,1,4,2];

tracking_data = trial_data.tracking([trial_data.tracking.trial_num] == itrial);
imaging_data = trial_data.imaging([trial_data.imaging.trial_num] == itrial);

speed_L = trial_data.tracking([trial_data.tracking.trial_num] == itrial).speed_L;
speed_R = trial_data.tracking([trial_data.tracking.trial_num] == itrial).speed_R;

if plt_cell
    figure('Name',['Trial: ',num2str(itrial);'; L: ',num2str(speed_L),', R: ',num2str(speed_R)]);
    k = 0;
    for ipaw = paw_order
        k = k+1;
        subplot(2,2,k)
        
        [hax,li1,li2] = plotyy(1:10,tracking_data.stride_bins_mean{ipaw}(:,1),1:10,imaging_data.img_bins_mean{ipaw});
        
        for i = 1:2
            hax(i).XLim = ([1 10]);
            hax(i).YColor = 'k';
            hax(i).YLim = [hax(i).YLim(1)-0.2 hax(i).YLim(2)+0.2];
        end
        
        hold(hax(1), 'on');
        %     errorbar(hax(1), 1:10, tracking_data.stride_bins_mean{ipaw}(:,1),...
        %         tracking_data.stride_bins_sem{ipaw}(:,1),[paw_col{ipaw},'.']);
        plot_error_shade(hax(1), 1:10, [], [], tracking_data.stride_bins_mean{ipaw}(:,1),...
            tracking_data.stride_bins_sem{ipaw}(:,1),paw_col{ipaw});
        
        li1.LineWidth = 2;
        li1.Color = paw_col{ipaw};
        li1.Marker = '.';
        
        axes(hax(2));
        hold(hax(2), 'on');
        for i = 1:length(li2)
            li2(i).LineWidth = 2;
            li2(i).Marker = '.';
            li2(i).Color = cell_col(i,:);
            %         errorbar(hax(2), 1:10,imaging_data.img_bins_mean{ipaw}(:,i),...
            %             imaging_data.img_bins_sem{ipaw}(:,i),'.','Color',cell_col(i,:));
            plot_error_shade(hax(2), 1:10, [], [], imaging_data.img_bins_mean{ipaw}(:,i),...
                imaging_data.img_bins_sem{ipaw}(:,i),cell_col(i,:));
        end
        xlabel('stride bin [-]');
        ylabel(hax(1),'x [mm]');
        ylabel(hax(2),'area \DeltaF/F [-]');
        title(paw_lab{ipaw});
    end
end
if plt_mean_cell
    figure('Name',['Trial: ',num2str(itrial),'; L: ',num2str(speed_L),', R: ',num2str(speed_R)]);
    k = 0;
    for ipaw = paw_order
        k = k+1;
        subplot(2,2,k);
        
        [hax,li1,li2] = plotyy(1:10,tracking_data.stride_bins_mean{ipaw}(:,1),1:10,imaging_data.img_bins_mean_cell{ipaw});
        
        for i = 1:2
            hax(i).XLim = [1 10];
            hax(i).YColor = 'k';
            hax(i).YLim = [hax(i).YLim(1)-0.2 hax(i).YLim(2)+0.2];
        end
        
        hold(hax(1), 'on');
        %     errorbar(hax(1), 1:10, tracking_data.stride_bins_mean{ipaw}(:,1),...
        %         tracking_data.stride_bins_sem{ipaw}(:,1),[paw_col{ipaw},'.']);
        plot_error_shade(hax(1), 1:10, [], [], tracking_data.stride_bins_mean{ipaw}(:,1),...
            tracking_data.stride_bins_sem{ipaw}(:,1),paw_col{ipaw});
        
        li1.LineWidth = 2;
        li1.Color = paw_col{ipaw};
        li1.Marker = '.';
        
        axes(hax(2));
        hold(hax(2), 'on');
        li2.LineWidth = 2;
        li2.Color = purple;
        li2.Marker = '.';
        
        %     errorbar(hax(2), 1:10,imaging_data.img_bins_mean_cell{ipaw},...
        %         imaging_data.img_bins_sem_cell{ipaw},'.','Color',purple);
        plot_error_shade(hax(2), 1:10, [], [], imaging_data.img_bins_mean_cell{ipaw},...
            imaging_data.img_bins_sem_cell{ipaw},purple);
        
        xlabel(hax(1),'stride bin [-]');
        ylabel(hax(1),'x [mm]');
        ylabel(hax(2),'area \DeltaF/F [-]');
        title(paw_lab{ipaw});
    end
end
end