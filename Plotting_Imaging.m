%% Path tiff
addpath(genpath('C:\Users\Virginia Casasnovas\Documents\Virginia\Data\IMAGING\'));

show_cell_roi_activity( 'concatenated.tif', cn, 70);

%% Load tiff
fname = 'concatenated.tif';
info = imfinfo(fname);
imageStack = [];
numberOfImages = length(info);
for k = 1:numberOfImages
    currentImage = imread(fname, k, 'Info', info);
    imageStack(:,:,k) = currentImage;
end 

%% Plot tiff

% figure;
% for c = 1:cn.n_cells
%     I = zeros(cn.fov_height, cn.fov_width);
%     I(cn.roi{c}) = 1;
%     imcontour(I); hold on;
%     text(cn.centroid{c}(1), cn.centroid{c}(2), num2str(c));
% end

figure;
colormap(gray);

key = '';
i = 0;

while (~strcmp(key, 'space'))&&(i<size(imageStack,3))
    key = get(gcf,'currentkey');
    i = i+1;
    
    imagesc(imageStack(:,:,i)); hold on;
%     imcontour(I)
    title(['Frame: ',num2str(i)]);
    pause(0.0001);
end
% l = medfilt2(g,[3 3]);

%% Plotting FFT
figure;
colormap(gray);

key = '';
i = 0;

while (~strcmp(key, 'space'))&&(i<size(imageStack,3))
    key = get(gcf,'currentkey');
    i = i+1;
    
    imagesc(log(1+abs(fftshift(fft2(imageStack(:,:,i)))))); hold on;
%     imcontour(I)
    title(['Frame: ',num2str(i)]);
    pause(0.0001);
end

%% Load extracted vars
% load('C:\Users\Virginia Casasnovas\Documents\Virginia\Data\IMAGING\carey_neurons.mat');
% load('C:\Users\Virginia Casasnovas\Documents\Virginia\Data\IMAGING\muk.mat');

neurons = 5:10;
muk.ica_sig_dtr = detrend(muk.ica_sig);
cn.intensity_dtr = detrend(cn.intensity);

ysepmuk = mean(std(muk.ica_sig_dtr(neurons,:)'))*10;
ysepcn = mean(std(cn.intensity_dtr(:,neurons)))*10;

fs = 300;
muk.time = (0:size(muk.ica_sig,2)-1)/fs;
cn.time = (0:size(cn.intensity,1)-1)/fs;

figure;
for i = 1:length(neurons)
    plot(muk.time,muk.ica_sig_dtr(neurons(i),:)'+ysepmuk*i);
    hold on;
end
set(gca,'YTickLabel',[]);
xlabel('Time [s]'); ylabel('Intensity');
title('Mukamel Neurons');


%% Carey neurons
wn = (2/300)*30;
[b,a] = butter(2,wn,'low');
datafilt = filtfilt(b,a,cn.intensity);

figure;
for i = 1:length(neurons)
    plot(cn.time,cn.intensity(:,neurons(i))+ysepcn*i);
    hold on;
%     plot(time,datafilt(:,neurons(i))+ysepcn*i,'LineWidth',2);
end
% set(gca,'YTickLabel',[]);
xlabel('Time [s]'); ylabel('Intensity');
title('Carey Neurons');

%% Carey neurons bkg
figure;
imagesc(imageStack(:,:,50));
hold on;
for c = 1:cn.n_cells
    I = zeros(cn.fov_height, cn.fov_width);
    I(cn.roi{c}) = 1;
    imcontour(I); hold on;
    text(cn.centroid{c}(1), cn.centroid{c}(2), num2str(c));
end

figure;
for i = 1:cn.n_cells
    plot(cn.intensity(:,i));
    hold on;
%     plot(time,datafilt(:,neurons(i))+ysepcn*i,'LineWidth',2);
end
% set(gca,'YTickLabel',[]);
xlabel('Time [s]'); ylabel('Intensity');
title('Carey Neurons Bkg'); 

