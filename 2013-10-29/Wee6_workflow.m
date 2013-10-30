%%%%%% Week 6 Workflow %%%%%%%
cd('D:\Promoted\R016-2012-10-03');
addpath('D:\Users\mvdmlab\My_Documents\GitHub\BIOL680\2013-09-30');
fname = 'R016-2012-10-03-CSC04a.Ncs';
csc = LoadCSC(fname);
event_times = [2000 1500 2458 3335 1909 3030];  % sample of event times.  
eventLFPplot(csc,event_times,'decimate_signal','yes','LFP_colour', [0 0 1])