function[] = plot_intensity_spd(group_speed,signal)

if ~isfield(group_speed,['mean_',signal]) && ~isfield(group_speed,['area_',signal])
    error('Necessary fields missing.')
else
    nspd = length(group_speed);
    
    ncells = size(group_speed(1).(['mean_',signal]),2);
    lilliac = [76, 0, 153]/255;
    purple = [204, 153, 253]/255;
    cell_col = [linspace(lilliac(1),purple(1),ncells)', linspace(lilliac(2),purple(2),ncells)', linspace(lilliac(3),purple(3),ncells)'];
    
    figure('Name','Cell mean intensity'); 
    for ispd = 1:nspd
        for icell = 1:ncells
            plot(ispd,group_speed(ispd).(['mean_',signal,'_spd_mean'])(icell),'o', 'MarkerEdgeColor',cell_col(icell,:));
            hold on;
%             errorbar(ispd,group_speed(ispd).mean_intensity_spd_mean(icell),group_speed(ispd).mean_intensity_spd_sem(icell),'Color',cell_col(icell,:));
        end
        
        plot(ispd,mean(group_speed(ispd).(['mean_',signal,'_spd_mean'])),'ko', 'MarkerFaceColor','k');
        errorbar(ispd,mean(group_speed(ispd).(['mean_',signal,'_spd_mean'])),std(group_speed(ispd).(['mean_',signal,'_spd_mean']),[],2)/sqrt(ncells),'Color','k');
        
        xlabel('belt speed [m/s]');
        ylabel('mean \DeltaF/F [-]');
        
    end
    ax = gca;
    ax.XLim = [0 nspd+1];
    ax.XTick = [1:nspd];
    for ispd = 1:nspd
        ax.XTickLabel{ispd} = ['L = ',num2str(group_speed(ispd).spd(1)),'; R = ',num2str(group_speed(ispd).spd(2))];
    end
    
    figure('Name','Cell area intensity');
    for ispd = 1:nspd
        for icell = 1:ncells
            plot(ispd,group_speed(ispd).(['area_',signal,'_spd_mean'])(icell),'o', 'MarkerEdgeColor',cell_col(icell,:));
            hold on;
        end
        
        plot(ispd,mean(group_speed(ispd).(['area_',signal,'_spd_mean'])),'ko', 'MarkerFaceColor','k');
        errorbar(ispd,mean(group_speed(ispd).(['area_',signal,'_spd_mean'])),std(group_speed(ispd).(['area_',signal,'_spd_mean']),[],2)/sqrt(ncells),'Color','k');
        
        xlabel('belt speed [m/s]');
        ylabel('area \DeltaF/F [-]');
        
    end
    ax = gca;
    ax.XLim = [0 nspd+1];
    ax.XTick = [1:nspd];
    for ispd = 1:nspd
        ax.XTickLabel{ispd} = ['L = ',num2str(group_speed(ispd).spd(1)),'; R = ',num2str(group_speed(ispd).spd(2))];
    end
    
%     figure('Name','Mean cell mean intensity');
%     for ispd = 1:nspd
%         
%         plot(ispd,mean(group_speed(ispd).mean_intensity_spd_mean),'o', 'MarkerEdgeColor',purple,'MarkerFaceColor',purple);
%         hold on;
%         errorbar(ispd,mean(group_speed(ispd).mean_intensity_spd_mean),std(group_speed(ispd).mean_intensity_spd_mean,[],2)/sqrt(ncells),'Color',purple);
%         
%         xlabel('belt speed [m/s]');
%         ylabel('mean \DeltaF/F [-]');
%         
%     end
%     ax = gca;
%     ax.XLim = [0 nspd+1];
%     ax.XTick = [1:nspd];
%     for ispd = 1:nspd
%         ax.XTickLabel{ispd} = ['L = ',num2str(group_speed(ispd).spd(1)),'; R = ',num2str(group_speed(ispd).spd(2))];
%     end
end
