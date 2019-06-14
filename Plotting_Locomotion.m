
root_path = 'C:\Users\Virginia Casasnovas\Documents\Virginia\Data\IO\';
% root_path = 'C:\Users\Virginia Casasnovas\Documents\Virginia\Data\MINISCOPES\HeadFreeTreadmillMiniscope\';

experiment_folder = '20190527 - dreadds io saline post-cno';

filesep = '\';

addpath(genpath(root_path));

% Behavior
tracking_main = 'TM TRACKING FILES';
tracking_folder = [root_path,tracking_main,filesep,experiment_folder,filesep]; % Folder with tracking files (.mat)

dir_tracking = dir([tracking_folder]);
folder_tracking = dir_tracking([dir_tracking.isdir]);

mouse_folder = {};
for i = 3:length(folder_tracking)
    mouse_folder{i-2} = folder_tracking(i).name;
end

trials = [2,5,12,18];

for imouse = 1:length(mouse_folder)
    for itrial = 1:length(trials)
        tracking_file = dir([tracking_folder,mouse_folder{imouse},filesep,'*_',num2str(trials(itrial)),'.mat']);
        
        load(tracking_file.name);
        
        figure('Name',[mouse_folder{imouse},'; Trial ',num2str(trials(itrial))]);
        plot(squeeze(final_tracks(1,1:4,:))');
        ylabel('x [mm]'); xlabel('Frames');
        title('Tracks');
    end
end

