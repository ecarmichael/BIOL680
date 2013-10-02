%%% Week 3 Functionality script %%%%%%

% This script compares the old LoadCSC function compared to the new
% LoadCSC_ec function which removes invalid data samples.  

fname  = 'R016-2012-10-08-CSC01a.ncs';
%% old section

old_csc = LoadCSC(fname);

new_csc = LoadCSC_ec(fname);

figure
subplot(211)
hld_1 = plot(diff(Range(old_csc)), 'r');
ylabel('Old LoadCSC Sample');
ylim([-0.5 10])
subplot(212)
hld_2 = plot(diff(Range(new_csc)), 'b');
ylabel('New LoadCSC Sample');
ylim([-0.5 10])
set(gcf,'Color',[1 1 1])

