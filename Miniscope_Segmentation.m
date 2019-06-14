%% Load file
clear all;
root_path = 'C:\Users\Virginia Casasnovas\Documents\Virginia\Data\MINISCOPES\HeadFreeTreadmillMiniscope\';
experiment_folder = '20190508 - miniscopes tied stair';

filesep = '\';
addpath(genpath(root_path));

raw_main = 'TM RAW FILES';
registered_main = 'TM PROCESSED IMAGING FILES';
segmented_main = 'TM SEGMENTED IMAGING FILES';

raw_path = [root_path,raw_main,filesep,experiment_folder,filesep]; 
registered_path = [root_path,registered_main,filesep,experiment_folder,filesep]; 
segmented_path = [root_path,segmented_main,filesep,experiment_folder,filesep]; 

%% Loading and processing
dir_raw = dir(raw_path);
mouse_raw = dir_raw([dir_raw.isdir]);

mouse_folder = {};
for i = 3:length(mouse_raw)
    mouse_folder{i-2} = mouse_raw(i).name;
end

% Params for segmenting cells
opt_mukamel.nPCs = 1000;
opt_mukamel.used_PCs = 1:300; % If this value is empty the algorithm prompts the user for manual selection
opt_mukamel.mu = 0.5;
rng(1);

for imouse = 1:length(mouse_folder)
    cat_file = rdir([registered_path,mouse_folder{imouse},filesep,'*_Cat.tiff']);
    if isempty(cat_file)
        cat_TIFF([registered_path,mouse_folder{imouse},filesep]);
        cat_file = rdir([registered_path,mouse_folder{imouse},filesep,'*_Cat.tiff']);
    end
    cat_file = cat_file.name;
    
    image_stack = imread(cat_file, 1, 'Info', imfinfo(cat_file));
    
    muk = get_mukamel_rois(cat_file, [], opt_mukamel.nPCs, opt_mukamel.used_PCs, opt_mukamel.mu, []);
    cn = convert_mukamel_to_carey_neuron(muk.ica_segments, muk.seg_centroid);
    cn = edit_segmentation_mukamel(cn,image_stack);
    
    hf = figure('Name',[mouse_folder{imouse},'; cell rois']);
    imshow(image_stack(:,:),[]); hold on;
    for i = 1:cn.n_cells
        imcontour(cn.mask{i},'r')
        hold on;
    end
    
    if ~exist([segmented_path,mouse_folder{imouse},filesep])
        mkdir([segmented_path,mouse_folder{imouse},filesep]);
    end
    
    str_cat_path = strsplit(cat_file,filesep);
    mask_file = [segmented_path,mouse_folder{imouse},filesep,str_cat_path{end}(1:end-8),'Mask.mat'];
    save(mask_file,'cn');
end

opts_deconv.fs = 30;
opts_deconv.recomputeKernel = 0;
opts_deconv.sensorTau = 0.7;
opts_deconv.estimateNeuropil = 1;
opts_deconv.deconvType = 'L0';

for imouse = 1:length(mouse_folder)
    mask_file = rdir([segmented_path,mouse_folder{imouse},filesep,'*_Mask.mat']);
    mask_file = mask_file.name;
    
    mask_session = load(mask_file);
    
    ready_files = rdir([registered_path,mouse_folder{imouse},filesep,'*_Ready.tiff']);
    
    for f = 1:length(ready_files)
        
        n_frames_ready = length(imfinfo(ready_files(f).name));
        
        str_ready_path = strsplit(ready_files(f).name,filesep);
        str_ready_file = strsplit(str_ready_path{end},'_');
        trial = str2num(str_ready_file{end-2});
        
        timestamp_data = importdata([raw_path,mouse_folder{imouse},'\miniscope\T',num2str(trial),'\timestamp.dat']); % Loading imaging timestamp file
        timestamp = timestamp_data.data(:,3);
        
        if length(timestamp)~=n_frames_ready
            fprintf('\nIgnoring trial %d: Number of frames in tiff file does not match the number of timestamps.\n',trial);
        else
            cn.time = timestamp/1000;
            
            [cn_data, ~, avg] = get_carey_neurons_mean_intensity(ready_files(f).name, mask_session.cn);
            cn.intensity = detrend(cn_data,'constant');
            [cn.spikes,~] = wrapperDECONV(opts_deconv,cn_data);
            
            segmented_file =  ([segmented_path,mouse_folder{imouse},filesep,str_ready_path{end}(1:end-10),'Seg.mat']);
            save(segmented_file,'cn')
        end
    end
end
