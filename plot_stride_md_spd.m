function[] = plot_stride_md_spd(group_speed,plt_cell,plt_mean_cell)

if ~isfield(group_speed,'stride_bins_mean') && ((~isfield(group_speed,'img_bins_mean')&& plt_cell)|| (~isfield(group_speed,'img_bins_mean_cell')&& plt_mean_cell))
    error('Necessary fields missing.')
else
    
    paw_lab = {'FR','HR','FL','HL'};
    paw_col = {'r','m','b','c'};
    spd_mrk = {'+','o','*','v','x','<','>'};
    
    nspd = length(group_speed);
    lilliac = [76, 0, 153]/255;
    purple = [204, 153, 253]/255;
    spd_col = [linspace(lilliac(1),purple(1),nspd)', linspace(lilliac(2),purple(2),nspd)', linspace(lilliac(3),purple(3),nspd)'];
    
    paw_order = [3,1,4,2];
    nbins = length(group_speed(1).stride_bins_mean_spd_mean{1}(:,1));
    
    if plt_cell
        ncells = size(group_speed(1).img_bins_mean_spd_mean{1},2);
        
        figure;
        k = 0;
        for ipaw = paw_order
            k = k+1; hl = []; li = 1; ylimit = [];
            subplot(2,2,k)
            
            for ispd = 1:length(group_speed)
                
                if ispd == 1
                    [hax,li1,li2] = plotyy(1:nbins,group_speed(ispd).stride_bins_mean_spd_mean{ipaw}(:,1),1:nbins,group_speed(ispd).img_bins_mean_spd_mean{ipaw});
                    li1.Marker = spd_mrk{ispd};
                    for i = 1:length(li2)
                        li2(i).Marker = spd_mrk{ispd};
                        li2(i).Color = spd_col(ispd,:);
                    end
                else
                    li1(end+1) = plot(hax(1),1:nbins,group_speed(ispd).stride_bins_mean_spd_mean{ipaw}(:,1));
                    li2 = [li2; plot(hax(2),1:nbins,group_speed(ispd).img_bins_mean_spd_mean{ipaw})];
                    
                    li1(end).Marker = spd_mrk{ispd};
                    for i = (ispd-1)*ncells:ispd*ncells
                        li2(i).Marker = spd_mrk{ispd};
                        li2(i).Color = spd_col(ispd,:);
                    end
                end
                hold(hax(1),'on');
                plot_error_shade(hax(1), 1:nbins, [], [], group_speed(ispd).stride_bins_mean_spd_mean{ipaw}(:,1),...
                    group_speed(ispd).stride_bins_mean_spd_sem{ipaw}(:,1),paw_col{ipaw});
                hold(hax(2),'on');
                
                ylimit(li,:) = [min(group_speed(ispd).stride_bins_mean_spd_mean{ipaw}(:,1)) max(group_speed(ispd).stride_bins_mean_spd_mean{ipaw}(:,1))];
                ylimit(li+1,:) = [min(group_speed(ispd).img_bins_mean_spd_mean{ipaw}(:)) max(group_speed(ispd).img_bins_mean_spd_mean{ipaw}(:))];
                li = li+2;
                
                for i = 1:2
                    hax(i).XLim = ([1 nbins]);
                    hax(i).YColor = 'k';
                end
                for i = 1:length(li1)
                    li1(i).LineWidth = 1.2;
                    li1(i).Color = paw_col{ipaw};
                end
                for i = 1:length(li2)
                    li2(i).LineWidth = 1.2;
                end
                if k == 2
                    strl{ispd} = ['L: ',num2str(group_speed(ispd).spd(1)),'; R: ',num2str(group_speed(ispd).spd(2)),' m/s'];
                    hl{ispd} = scatter(NaN,NaN,'k',spd_mrk{ispd},'LineWidth',1.2);
                end
            end
            for i = 1:2
                hax(i).YLim = [round(min(ylimit(i:2:end,1)))-0.2 0.2+round(max(ylimit(i:2:end,2)))];
                hax(i).YTick = round(linspace(hax(i).YLim(1), hax(i).YLim(2), 5), 1);
            end
            
            if k == 2
                legend([hl{:}],strl)
            end
            xlabel('stride bin [-]');
            ylabel(hax(1),'x [mm]');
            ylabel(hax(2),'area F [-]');
            title(paw_lab{ipaw});
            
        end
    end
    if plt_mean_cell
        figure;
        k = 0; li2 = [];
        for ipaw = paw_order
            k = k+1; hl = []; li = 1; ylimit = [];
            subplot(2,2,k)
            
            for ispd = 1:length(group_speed)
                
                if ispd == 1
                    [hax,li1,li2] = plotyy(1:nbins,group_speed(ispd).stride_bins_mean_spd_mean{ipaw}(:,1),1:nbins,group_speed(ispd).img_bins_mean_cell_spd_mean{ipaw});
                    li1.Marker = spd_mrk{ispd};
                    li2.Marker = spd_mrk{ispd};
                    li2.Color = spd_col(ispd,:);
                else
                    li1(end+1) = plot(hax(1),1:nbins,group_speed(ispd).stride_bins_mean_spd_mean{ipaw}(:,1));
                    li2(end+1) = plot(hax(2),1:nbins,group_speed(ispd).img_bins_mean_cell_spd_mean{ipaw});
                    
                    li1(end).Marker = spd_mrk{ispd};
                    li2(end).Marker = spd_mrk{ispd};
                    li2(end).Color = spd_col(ispd,:);
                end
                
                hold(hax(1),'on');
                plot_error_shade(hax(1), 1:nbins, [], [], group_speed(ispd).stride_bins_mean_spd_mean{ipaw}(:,1),...
                    group_speed(ispd).stride_bins_mean_spd_sem{ipaw}(:,1),paw_col{ipaw});
                
                hold(hax(2),'on');
                plot_error_shade(hax(2), 1:nbins, [], [], group_speed(ispd).img_bins_mean_cell_spd_mean{ipaw},...
                    group_speed(ispd).img_bins_mean_cell_spd_sem{ipaw},spd_col(ispd,:));
                
                ylimit(li,:) = [min(group_speed(ispd).stride_bins_mean_spd_mean{ipaw}(:,1)) max(group_speed(ispd).stride_bins_mean_spd_mean{ipaw}(:,1))];
                ylimit(li+1,:) = [min(group_speed(ispd).img_bins_mean_cell_spd_mean{ipaw}(:)) max(group_speed(ispd).img_bins_mean_cell_spd_mean{ipaw}(:))];
                li = li+2;
                
                for i = 1:2
                    hax(i).XLim = ([1 nbins]);
                    hax(i).YColor = 'k';
                end
                for i = 1:length(li1)
                    li1(i).LineWidth = 1.2;
                    li1(i).Color = paw_col{ipaw};
                end
                for i = 1:length(li2)
                    li2(i).LineWidth = 1.2;
                end
                if k == 2
                    strl{ispd} = ['L: ',num2str(group_speed(ispd).spd(1)),'; R: ',num2str(group_speed(ispd).spd(2)),' m/s'];
                    hl{ispd} = scatter(NaN,NaN,'k',spd_mrk{ispd},'LineWidth',1.2);
                end
            end
            for i = 1:2
                hax(i).YLim = [round(min(ylimit(i:2:end,1)))-0.2 0.2+round(max(ylimit(i:2:end,2)))];
                hax(i).YTick = round(linspace(hax(i).YLim(1), hax(i).YLim(2), 5),1);
            end
            if k == 2
                legend([hl{:}],strl)
            end
            xlabel('stride bin [-]');
            ylabel(hax(1),'x [mm]');
            ylabel(hax(2),'area \DeltaF [-]');
            title(paw_lab{ipaw});
            set(hax(1),'FontSize',12);
            set(hax(2),'FontSize',12);            
        end
    end
end
end