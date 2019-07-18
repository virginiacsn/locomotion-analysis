function[] = plot_mean_int_belt_speed(trial_data)

N = length(trial_data);
blue = [0, 0, 1];
pink = [255, 192, 203]/255;
colors_p = [linspace(blue(1),pink(1),N)', linspace(blue(2),pink(2),N)', linspace(blue(3),pink(3),N)'];

figure;
for itrial = 1:length(trial_data)
    if trial_data(itrial).speed_R == 0.15
        x_R = 1;
    else
        x_R = 2;
    end
    if trial_data(itrial).speed_L == 0.15
        x_L = 3;
    else
        x_L = 4;
    end
    
    area_int = mean(trapz(trial_data(itrial).imaging.time,trial_data(itrial).imaging.cn.intensity_final,1));
    %     area_int_sem(itrial) = mean(trapz(trial_data(itrial).imaging.time,trial_data(itrial).imaging.cn.intensity_final,1));
    
    plot(x_R,area_int,'Color',colors_p(itrial,:),'MarkerFaceColor',colors_p(itrial,:),'Marker','o','MarkerSize',7);
    hold on,
    plot(x_L,area_int,'Color',colors_p(itrial,:),'MarkerFaceColor',colors_p(itrial,:),'Marker','o','MarkerSize',7);
end

ax = gca;
ax.XTick = [1:4];
ax.XTickLabel = {'R = 0.15','R = 0.3','L = 0.15','L = 0.3'};
xlim([0 5]);
end