function[] = plot_stride_md_spd_sep(group_speed)

paw_lab = {'FR','HR','FL','HL'};
paw_col = {'r','m','b','c'};
spd_mrk = {'+','o','*','v','x','<','>'};
spd_line = {'-','--'};
type_leg = {'split','tied'};

nspd = length(group_speed);
lilliac = [76, 0, 153]/255;
purple = [204, 153, 253]/255;
spd_col = [linspace(lilliac(1),purple(1),nspd)', linspace(lilliac(2),purple(2),nspd)', linspace(lilliac(3),purple(3),nspd)'];

paw_order = [3,1,4,2];
nbins = length(group_speed(1).stride_bins_mean_spd_mean{1}(:,1));

figure;
k = 0; li2 = [];
for ipaw = paw_order
    k = k+1; hl = []; li = 1; ylimit = [];
    subplot(2,2,k)
    
    for ispd = 1:length(group_speed)
        
        hl(ispd) = plot(1:nbins,group_speed(ispd).stride_bins_mean_spd_mean{ipaw}(:,1),spd_line{ispd},'Color',paw_col{ipaw},'LineWidth',2);
        hold on;
        ax = gca;
        plot_error_shade(ax, 1:nbins, [], [], group_speed(ispd).stride_bins_mean_spd_mean{ipaw}(:,1),...
            group_speed(ispd).stride_bins_mean_spd_sem{ipaw}(:,1),paw_col{ipaw});
        
        ax.XLim = ([1 nbins]);
        
        ylimit(ispd,:) = [min(group_speed(ispd).stride_bins_mean_spd_mean{ipaw}(:,1)) max(group_speed(ispd).stride_bins_mean_spd_mean{ipaw}(:,1))];        
    end
    ax.YLim = [round(min(ylimit(:,1)))-5 5+round(max(ylimit(:,2)))];
    ax.YTick = round(linspace(ax.YLim(1),ax.YLim(2), 5),1);
    
    if k == 2
        legend(hl,type_leg)
    end
    if k == 1 || k == 3
        ylabel('x [mm]');
    end
    if k == 3 || k == 4
        xlabel('stride bin [-]');
    end
    box off;
    title(paw_lab{ipaw});
    set(ax,'FontSize',13);
end

figure;
k = 0; li2 = [];
for ipaw = paw_order
    k = k+1; hl = []; li = 1; ylimit = [];
    subplot(2,2,k)
    for ispd = 1:length(group_speed)
        
        hl(ispd) = plot(1:nbins,group_speed(ispd).img_bins_mean_cell_spd_mean{ipaw}(:,1),spd_line{ispd},'Color',spd_col(ispd,:),'LineWidth',2);
        hold on;
        ax = gca;
        plot_error_shade(ax, 1:nbins, [], [], group_speed(ispd).img_bins_mean_cell_spd_mean{ipaw},...
            group_speed(ispd).img_bins_mean_cell_spd_sem{ipaw},spd_col(ispd,:));
        
        ax.XLim = ([1 nbins]);
        ylimit(ispd,:) = [min(group_speed(ispd).img_bins_mean_cell_spd_mean{ipaw}(:)) max(group_speed(ispd).img_bins_mean_cell_spd_mean{ipaw}(:))];
        
    end
    ax.YLim = [round(min(ylimit(:,1)))-0.5 0.5+round(max(ylimit(:,2)))];
    ax.YTick = round(linspace(ax.YLim(1),ax.YLim(2), 5),1);
    if k == 2
        legend(hl,type_leg)
    end
    if k == 1 || k == 3
        ylabel(ax,'area \DeltaF [-]');
    end
    if k == 3 || k == 4
        xlabel('stride bin [-]');
    end
    title(paw_lab{ipaw});
    box off;
    set(ax,'FontSize',13);
end
end
