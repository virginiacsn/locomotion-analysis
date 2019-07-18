

int_data = TD(1).imaging(1).int_bkg_det(:,1);

npts = size(int_data,1);
ncells = size(int_data,2);


win = 1.5*30;
bins = floor(npts/win);

base = [];
for ibin = 1:bins-1
    base(ibin) = prctile(int_data(1+(ibin-1)*win:ibin*win),10);  
end
base(ibin+1) = prctile(int_data(1+ibin*win:end),10);  
% base(ibin+2) = int_data(end);

base_interp = movingAvg(interp1([1:win:(npts-win)],base,[1:npts]),10);

figure; plot(1:npts,base_interp,'b'), hold on; plot(1:npts,int_data,'r'), plot(1:npts,int_data-base_interp')


% figure; plot(TD(1).tracking(1).time,squeeze(TD(1).tracking(1).final_tracks(1,1:4,:))')
