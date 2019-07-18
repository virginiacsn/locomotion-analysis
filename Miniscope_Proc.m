clear all;

root_path = 'C:\Users\Virginia Casasnovas\Documents\Virginia\Data\MINISCOPES\HeadFreeTreadmillMiniscope\';
addpath(genpath(root_path));
tracking_folder = 'TM TRACKING FILES';
processed_folder = 'TM PROCESSED FILES';
video_folder = 'TM IMAGING FILES';
neuron_folder = 'TM IMAGING PROCESSING FILES';
video_data_folder = 'TM RAW FILES';

experiment_folder = 'tmp';
mouse_folder = 'MC2370';

trial = 2;
trial_folder = ['miniscope\T',num2str(trial)];

trial_video_data_folder = [video_data_folder,'\',experiment_folder,'\',mouse_folder,'\',trial_folder,'\']; % Folder with imaging timestamp file (.dat)
trial_video_folder = [video_data_folder,'\',experiment_folder,'\',mouse_folder];% Folder with registered imaging files (.tiff)
trial_neuron_folder = [neuron_folder,'\',experiment_folder,'\',mouse_folder,'\']; % Folder with edited segmented cell signal files (.mat)
trial_tracking_folder = [tracking_folder,'\',experiment_folder,'\',mouse_folder,'\']; % Folder with tracking files (.mat)
trial_processed_folder = [processed_folder,'\',experiment_folder,'\',mouse_folder,'\']; % Folder with processed tracking files (.mat)

session_path = [root_path,video_data_folder,'\',experiment_folder,'\',mouse_folder];%'C:\Users\Virginia Casasnovas\Documents\Virginia\Data\MINISCOPES\HeadFreeTreadmillMiniscope\TM RAW FILES\erc 0,25 - 0,375\M3';
input_path = [root_path,video_data_folder,'\',experiment_folder];%'C:\Users\Virginia Casasnovas\Documents\Virginia\Data\MINISCOPES\HeadFreeTreadmillMiniscope\TM RAW FILES\erc 0,25 - 0,375';
tm_imaging_folder =  [root_path,video_folder];%'C:\Users\Virginia Casasnovas\Documents\Virginia\Data\MINISCOPES\HeadFreeTreadmillMiniscope\TM IMAGING FILES';

%% Process imaging data
% convert data to the standard 'TM RAW DATA' folder format if not
% already done
if(isempty(dir([session_path, filesep, '*_m.tif'])))
    convert_miniscope_session_data(session_path);
end


dirs = get_directory_tree_from_path(input_path);
tm_imaging_folder = [tm_imaging_folder, filesep, dirs{end}];

options = [];

options.transform_imaging_videos = 1;
options.copy_syncronization_files = 0;
options.aggregate_imaging_files = 0;
options.remove_dark_frames = 0;
options.register_imaging_videos = 1;

options.overwrite_transform_directory = 0;
options.overwrite_registration_template = 0;
options.overwrite_registered_files = 0;
options.registered_file_suffix = '_Reg';
options.registered_final_suffix = '_Ready';
options.registration_filtered_files_suffix = '_Filtered';

options.registration_template_filename = 'registration_template.tif';
options.register_videos_in_directory_independently = 0;
options.delete_registration_filtered_files = 0;
options.delete_registered_files = 0;
options.delete_transformation_files = 0;

options = get_imaging_pipeline_options(options);

% pre-process pipeplie
process_imaging_pipeline(input_path, tm_imaging_folder, options);


% extract rois using mukamel

%% Loading and processing
% Load tiff
video_file = dir([root_path,trial_video_folder,mouse_folder,'_*_',num2str(trial),'_m_Ready.tif']);
video_name = video_file.name;
info = imfinfo(video_name);

image_stack = [];
num_images = length(info);
for k = 1:num_images
    current_image = imread(video_name, k, 'Info', info);
    image_stack(:,:,k) = current_image;
end

if length(sysClock)~=num_images
    fprintf('\nNumber of frames in video file does not match the number of timestamps.\n');
else
    % Segment cells
    if ~exist([root_path,trial_neuron_folder,'neurons_',num2str(trial),'.mat'],'file') % Segmented file does not exist
        options.nPCs = 100;
        options.used_PCs = []; % If this value is empty the algorithm prompts the user for manual selection
        options.mu = 0.5;
        
        cells_sort_file = [root_path,trial_video_folder,video_name];
        muk = get_mukamel_rois(cells_sort_file, [], options.nPCs, options.used_PCs, options.mu, []);
        
        % convert mukamel rois to carey neuron
        cn = convert_mukamel_to_carey_neuron(muk.ica_segments, muk.seg_centroid);
        
        [cn_data, ~, avg] = get_carey_neurons_mean_intensity(cells_sort_file, cn);
        cn.intensity = cn_data;
        
        % Deconvolute to get spikes
        ops.fs = 30;
        ops.recomputeKernel = 0;
        ops.sensorTau = 0.7;
        ops.estimateNeuropil = 1;
        ops.deconvType = 'L0';
        
        [sp, ca] = wrapperDECONV(ops, cn_data);
        cn.sp = sp;
        cn.sp_thres = sp>max(sp(:))*0.1;
        cn.ca = ca;
        
        % Saving mukamel and carey neurons
        if ~exist([root_path,trial_neuron_folder])
            mkdir([root_path,trial_neuron_folder]);
        end
        
        save([root_path,trial_neuron_folder,'neurons_',num2str(trial)],'cn','muk','options');
        
        % Plot segmentation
        figure('Name','Cell contours');
        imshow(image_stack(:,:,200),[]); hold on;
        for i = 1:cn.n_cells
            imcontour(cn.mask{i})
            hold on;
        end
    else
        % Load segmented cells
        load([root_path,trial_neuron_folder,'neurons_',num2str(trial)]);
        
        % Plot segmentation
        figure('Name','Cell contours');
        imshow(image_stack(:,:,200),[]); hold on;
        for i = 1:cn.n_cells_final
            imcontour(cn.mask_final{i})
            hold on;
        end
    end
end

%% Mukamel old script
%TODO: consider concatenating a few videos beforehand
% TODO: save mukamel output file
options.nPCs = 400;
options.used_PCs = 1:350; % if this value is empty the algorithm prompts the user for manual selection
options.mu = 0.5;

file = 'Z:\LocomotionExperiments\RT Self-paced\LocomotorLearning\TM IMAGING FILES\20181018 - water deprived first test\concatvideo\concatenated.tif';


muk = get_mukamel_rois(file, [], options.nPCs, options.used_PCs, options.mu, []);

% convert mukamel rois to carey neuron
% TODO: save carey neuron to output file
cn = convert_mukamel_to_carey_neuron(muk.ica_segments, muk.seg_centroid);

[cn_data, ~, avg] = get_carey_neurons_mean_intensity(cells_sort_file, cn);
cn.intensity = cn_data;

