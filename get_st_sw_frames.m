function[st_sw_frames] = get_st_sw_frames(tracking_data)

st_sw_frames = {};
for ipaw = 1:4
    if isfield(tracking_data,'ptsraw')
        if size(tracking_data.ptsraw.stance{ipaw},1) >= size(tracking_data.ptsraw.swing{ipaw},1)
            npts = size(tracking_data.ptsraw.swing{ipaw},1);
        else
            npts = size(tracking_data.ptsraw.stance{ipaw},1);
        end
        st_sw_frames{ipaw} = [tracking_data.ptsraw.stance{ipaw}(1:npts,5), tracking_data.ptsraw.swing{ipaw}(1:npts,5)];
    elseif isfield(tracking_data,'stride_fr')
        st_sw_frames{ipaw} = tracking_data.stride_fr{ipaw}(:,1:2);
    elseif isfield(tracking_data,'strides')
        st_sw_frames{ipaw} = tracking_data.strides{ipaw}(:,1:2);
    end
end
end