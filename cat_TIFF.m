function[] = cat_TIFF(input_dir)
tiff_files = dir([input_dir,'*_Ready.tiff']);
% tiff_files = tiff_files([tiff_files.isdir]);

str_tiff = strsplit(tiff_files(1).name,'_');
out_tiff_name = [str_tiff{1},'_',str_tiff{2},'_',str_tiff{3},'_',str_tiff{4},'_Cat.tiff'];

cat_tiff = 1:5:length(tiff_files);

for ifile = 1:length(cat_tiff)
    info = imfinfo(tiff_files(cat_tiff(ifile)).name);
    num_frames = length(info);

    for iframe = 1:num_frames
        current_image = imread(tiff_files(cat_tiff(ifile)).name, iframe, 'Info', info);
        
        imwrite(current_image,[input_dir,out_tiff_name],'tiff','WriteMode','append');
    end
end