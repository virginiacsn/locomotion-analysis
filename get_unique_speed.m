function[unique_spd,trial_spd] = get_unique_speed(tracking_data)

spd = cat(1,[tracking_data.speed_L],[tracking_data.speed_R])';
unique_spd = unique(spd,'rows');
trk_trial = [tracking_data.trial_num];

for ispd = 1:size(unique_spd,1)
    trial_spd{ispd} = trk_trial([tracking_data.speed_L]==unique_spd(ispd,1)&[tracking_data.speed_R]==unique_spd(ispd,2));
end
end