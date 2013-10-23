%%%% Week 5 sample script %%%
cd('D:\Promoted\R042-2013-08-18');
csc_name = 'R042-2013-08-18-CSC03a.ncs';
fc = FindFiles('*.t');
S = LoadSpikes(fc);
csc = LoadCSC(csc_name);
[evt] = detectSWR(csc_name);

neuroplot(S, csc,'evt',evt.t)
