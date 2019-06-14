%% Figures
%% Tracking stride modulation
%% Subplot per paw
paw_lab = {'FR','HR','FL','HL'};
figure;
kproc = 0;
plot_trials = [13,14,15,16,17,18,19,20];

N = length(plot_trials); % Set higher than maximum number of segmented cells
red = [0, 0, 1];
pink = [255, 192, 203]/255;
colors_p = [linspace(red(1),pink(1),N)', linspace(red(2),pink(2),N)', linspace(red(3),pink(3),N)'];

for imouse = 2%:length(TD)
    for itrial = 1:length(plot_trials)
        for ipaw = 1:4
            kproc = kproc+1;
            subplot(3,4,ipaw);
            plot(TD(imouse).tracks(plot_trials(itrial)).stride_bins_mean{ipaw}(:,1),'o-','Color',colors_p(itrial,:)); hold on;
            %         plot(trial_data(itrial).stride_bins{ipaw}(:,1,1),'*-');
            %errorbar(1:10,TD(imouse).tracks(itrial).stride_bins_mean{ipaw}(:,1),TD(imouse).tracks(itrial).stride_bins_sem{ipaw}(:,1),'.')
            if itrial == 1
                title(paw_lab{ipaw})
            end
            if ipaw == 4 && itrial == length(plot_trials)
                legend(strsplit(num2str(plot_trials)));
            end
            subplot(3,4,ipaw+4);
            plot(TD(imouse).tracks(plot_trials(itrial)).stride_bins_mean{ipaw}(:,2),'o-','Color',colors_p(itrial,:)); hold on;
            
            subplot(3,4,ipaw+8);
            plot(TD(imouse).tracks(plot_trials(itrial)).stride_bins_mean{ipaw}(:,1),TD(imouse).tracks(plot_trials(itrial)).stride_bins_mean{ipaw}(:,2),'o-','Color',colors_p(itrial,:)); hold on;
            
            if ipaw == 1
%                 ylabel(TD(imouse).tracks(itrial).trial_type)
            end
        end
    end
end

%% Paw vs. paw
paw_lab = {'FR','HR','FL','HL'};
figure;
kproc = 0;
plot_trials = [1,5,14,19];

N = length(plot_trials); % Set higher than maximum number of segmented cells
red = [0, 0, 1];
pink = [255, 192, 203]/255;
colors_p = [linspace(red(1),pink(1),N)', linspace(red(2),pink(2),N)', linspace(red(3),pink(3),N)'];

for imouse = 2%:length(TD)
    for itrial = 1:length(plot_trials)
        kproc = kproc+1;
        plot(TD(imouse).tracks(plot_trials(itrial)).stride_bins_mean{1}(:,1),TD(imouse).tracks(plot_trials(itrial)).stride_bins_mean{2}(:,1),'o-','Color',colors_p(itrial,:)); hold on;
%         errorbar(TD(imouse).tracks(plot_trials(itrial)).stride_bins_mean{1}(:,1),TD(imouse).tracks(plot_trials(itrial)).stride_bins_mean{2}(:,1),TD(imouse).tracks(itrial).stride_bins_sem{1}(:,1),'.','Color',colors_p(itrial,:));
%         errorbar(TD(imouse).tracks(plot_trials(itrial)).stride_bins_mean{1}(:,1),TD(imouse).tracks(plot_trials(itrial)).stride_bins_mean{2}(:,1),TD(imouse).tracks(itrial).stride_bins_sem{2}(:,1),'.','Color',colors_p(itrial,:));

        %             if itrial == 1
        %                 title(paw_lab{ipaw})
        %             end
        if itrial == length(plot_trials)
            legend(strsplit(num2str(plot_trials)));
        end
        %             subplot(2,1,2);
        %             plot(TD(imouse).tracks(plot_trials(itrial)).stride_bins_mean{1}(:,2),TD(imouse).tracks(plot_trials(itrial)).stride_bins_mean{2}(:,2),'o-','Color',colors_p(itrial,:)); hold on;
        %
    end
end

%% ds vs. coo for each paw
for imouse = 1:length(TD)
    figure('Name',TD(imouse).mouse);
    for ipaw = 1:4
        subplot(2,2,ipaw)
        scatter(TD(imouse).trial.double_support(:,ipaw),TD(imouse).trial.coo_body(:,ipaw),60,TD(imouse).trial.col_type,'filled');
        hold on;
        text(TD(imouse).trial.double_support(:,ipaw)+0.02,TD(imouse).trial.coo_body(:,ipaw)-0.02,(num2str(TD(imouse).trial.trial_num)),'FontSize',8);
        xlabel('% ds'); ylabel('coo [mm]')
        title(paw_lab{ipaw});
        grid on; axis equal;
    end
end

%% ds vs. coo for sym
for imouse = 1:length(TD)
    figure('Name',TD(imouse).mouse);
    for isym = 1:2
        subplot(1,2,isym)
        scatter(TD(imouse).trial.sym.double_support(:,isym),TD(imouse).trial.sym.coo_body(:,isym),60,TD(imouse).trial.col_type,'filled');
        hold on;
        text(TD(imouse).trial.sym.double_support(:,isym)+0.02,TD(imouse).trial.sym.coo_body(:,isym)-0.02,(num2str(TD(imouse).trial.trial_num)),'FontSize',8);
        xlabel('% ds'); ylabel('coo [mm]')
        title(sym_lab{isym});
        grid on; axis equal;
    end
end
