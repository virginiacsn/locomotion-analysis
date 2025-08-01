function[cn] = edit_segmentation_mukamel(cn,image_stack)
% Function to edit segmentation output from mukamel algorithm. 
%
% Inputs: 
%   - cn: struct from mukamel output. 
%   - image_stack: NxMxF matrix where NxM is image size and F is the number
%   of frames in video.
% Outputs: 
%   Collected in cn struct.
%   - n_cells_final
%   - mask_final
%   - centroid_final
%   - intensity_final
%   - ca_final
%   - sp_final
%   - sp_thres_final

%% Establish bounding box to neurons of interest
centroids = reshape([cn.centroid{:}],[2,cn.n_cells])';
xcnt = centroids(:,1); ycnt = centroids(:,2);

% hcnt = figure('Name','Cell centroids');
% imshow(image_stack(:,:),[]); hold on;
% for i = 1:cn.n_cells
%     plot(cn.centroid{i}(1),cn.centroid{i}(2),'r*');
%     hold on;
% end
% title('Neuron bounding box');
% 
% [xpts,ypts] = getpts(hcnt);

% xcnt_max = max(xpts); ycnt_max = max(ypts);
% xcnt_min = min(xpts); ycnt_min = min(ypts);

xcnt_max = size(image_stack,1)-10; ycnt_max =  size(image_stack,2)-10;
xcnt_min = 10; ycnt_min = 10;

neurons_final = find(xcnt<xcnt_max & ycnt<ycnt_max & xcnt>xcnt_min & ycnt>ycnt_min);

%% Discard small masks
k = 0;
neurons_final_tmp = [];
for i = 1:length(neurons_final)
    area_mask = sum(cn.mask{neurons_final(i)}(:));
    if area_mask > size(image_stack,1)*size(image_stack,2)*0.005
        k = k+1;
        neurons_final_tmp(k) = neurons_final(i);
    end
end
neurons_final = neurons_final_tmp;

centroids_final = centroids(neurons_final,:);
mask_tmp = cn.mask(neurons_final);

%% Identify overlapping masks with centroid distance matrix
dist_mat = zeros(size(neurons_final,1));
for i = 1:size(neurons_final,1)
    for j = i:size(neurons_final,1)
        dist_mat(i,j) = norm(centroids_final(i,:)-centroids_final(j,:),2);
    end
end

[idx,idy] = find(dist_mat<mean(dist_mat(:))*0.3 & dist_mat>0);

if ~isempty(idx)
    figure;
    for i = 1:length(idx)
        if rem(length(idx),2) == 0
            subplot(2,length(idx)/2,i)
        else
            subplot(1,length(idx),i)
        end
        imshow(image_stack(:,:),[]); hold on;
        n1 = plot(centroids_final(idx(i),1),centroids_final(idx(i),2),'b*');
        imcontour(mask_tmp{idx(i)},'b');
        
        n2 = plot(centroids_final(idy(i),1),centroids_final(idy(i),2),'r*');
        imcontour(mask_tmp{idy(i)},'r');
        
        legend([n1,n2],{num2str(idx(i)),num2str(idy(i))})
    end    
    overlap_neurons = input('\nOverlapping neurons to remove? ');
else
    fprintf('\nNo overlapping neurons.\n');
    overlap_neurons = [];
end

% Remove overlapping neurons
neurons_final_tmp = neurons_final;
neurons_final_tmp(overlap_neurons) = [];

cn_fields = fieldnames(cn);
for ifield = 1:length(cn_fields)
    if strcmp(cn_fields{ifield},'n_cells')
        cn.n_cells = length(neurons_final_tmp);
    elseif iscell(cn.(cn_fields{ifield}))
        cn.(cn_fields{ifield}) = cn.(cn_fields{ifield})(neurons_final_tmp);
    end
end
% cn.intensity = detrend(cn.intensity(:,neurons_final_tmp),'constant');
% cn.ca = cn.ca(:,neurons_final_tmp);
% cn.sp = cn.sp(:,neurons_final_tmp);
% cn.sp_thres = cn.sp_thres(:,neurons_final_tmp);

end