%%%%%%%%% Week 6 Sandbox %%%%%%%%%%%%%%%
%% load the data
cd('D:\Promoted\R016-2012-10-03');
addpath('D:\Users\mvdmlab\My_Documents\GitHub\BIOL680\2013-09-30');
% remember to cd to the correct folder here, may need to get this file from the lab database
fname = 'R016-2012-10-03-CSC04a.Ncs';
 
csc = LoadCSC_ec(fname,'plot_out', 0);
 
hdr = getHeader(csc);
Fs = hdr.SamplingFrequency;
% Fs = Fs-100;
csc = LoadCSC(fname);
%% restrict data
cscR = Restrict(csc,3282,3286); % if you don't have this, implement it (it's one line of code!)
figure 
plot(cscR); % note there are various gamma oscillations present, as well as a large negative-going transient
 
%% construct and plot the spectrogram
[S,F,T,P] = spectrogram(Data(cscR),hanning(512),256,1:200,Fs);
imagesc(T,F,10*log10(P)); % converting to dB as usual
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');  
hold on;
tvec = Range(cscR);
data = Data(cscR);
 
lfp_minmax = 25; lfp_cent = 125; % range and mean of LFP plotting
tvec0 = tvec - tvec(1); % align LFP with spectrogram
data = rescale(data,-lfp_minmax,lfp_minmax); data = data+lfp_cent;
 
lfp_h = plot(tvec0,data,'k');
%% construct Altered and plot the spectrogram

[S,F,T,P] = spectrogram(Data(cscR),hanning(512),384,.01:200,Fs);  
imagesc(T,F,10*log10(P)); % converting to dB as usual
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');  
hold on;
tvec = Range(cscR);
data = Data(cscR);
 
lfp_minmax = 25; lfp_cent = 125; % range and mean of LFP plotting
tvec0 = tvec - tvec(1); % align LFP with spectrogram
data = rescale(data,-lfp_minmax,lfp_minmax); data = data+lfp_cent;
 
lfp_h = plot(tvec0,data,'k');

%% construct Altered and plot the spectrogram  (rectwin)  % Rectwin seems to cause the resolution to appear better by havingh a strict cut off at each window. 

[S,F,T,P] = spectrogram(Data(cscR),rectwin(512),384,.01:200,Fs);  
imagesc(T,F,10*log10(P)); % converting to dB as usual
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');  
hold on;
tvec = Range(cscR);
data = Data(cscR);
 
lfp_minmax = 25; lfp_cent = 125; % range and mean of LFP plotting
tvec0 = tvec - tvec(1); % align LFP with spectrogram
data = rescale(data,-lfp_minmax,lfp_minmax); data = data+lfp_cent;

%% Compare a 522 sample and 1024 smaple window spectrogram with even number of time bins
[S,F,T,P] = spectrogram(Data(cscR),hanning(512),384,.01:200,Fs);  
subplot(211)
imagesc(T,F,10*log10(P)); % converting to dB as usual
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');  
title('512')

[S,F,T,P] = spectrogram(Data(cscR),hanning(1024),896,.01:200,Fs); 
subplot(212)
imagesc(T,F,10*log10(P)); % converting to dB as usual
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');  
title('1024')


%% mind the gaps
cscR = Restrict(csc,3300,3340); % section of data with a gap
 
[S,F,T,P] = spectrogram(Data(cscR),rectwin(256),128,1:200,Fs);
 
imagesc(T,F,10*log10(P)); 
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');  
 
hold on;
tvec = Range(cscR); data = Data(cscR);
 
lfp_minmax = 25; lfp_cent = 125; % range and mean of LFP plotting
tvec0 = tvec - tvec(1); % align LFP with spectrogram
data = rescale(data,-lfp_minmax,lfp_minmax); data = data+lfp_cent;
 
lfp_h = plot(tvec0,data,'k');
xlim([tvec0(1) tvec0(end)]);

%% find the events file
fn = FindFile('*Events.nev');
[EVTimeStamps, EventIDs, TTLs, EVExtras, EventStrings, EVHeader] = Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);
evt = getEvents_value;

%% Field Trip stuff
%% remember to cd to your data folder
fc = {'R016-2012-10-03-CSC04a.ncs'};
data = ft_read_neuralynx_interp(fc);
%% creating a trial function in Fieldtrip
cfg = [];
cfg.trialfun = 'ft_trialfun_lineartracktone2';
cfg.trialdef.hdr = data.hdr;
cfg.trialdef.pre = 2.5;
cfg.trialdef.post = 5;
 
cfg.trialdef.eventtype = 'nosepoke'; % could be 'nosepoke', 'reward', 'cue'
cfg.trialdef.location = 'both'; % could be 'left', 'right', 'both'
cfg.trialdef.block = 'both'; % could be 'value', 'risk'
cfg.trialdef.cue = {'c1','c3','c5'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'}
 
[trl, event] = ft_trialfun_lineartracktone2(cfg);
cfg.trl = trl;
 
data_trl = ft_redefinetrial(cfg,data);

%% create an event_triggered spectrogram
cfg              = []; % start with empty cfg
cfg.output       = 'pow';
cfg.channel      = 'R016-2012-10-03-CSC04a';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:100; % frequencies of interest
cfg.t_ftimwin    = 20./cfg.foi;  % window size: fixed at 0.5s
cfg.toi          = -1:0.05:4; % times of interest
 
TFR = ft_freqanalysis(cfg, data_trl);
 
figure
cfg = []; cfg.channel = 'R016-2012-10-03-CSC04a';
ft_singleplotTFR(cfg, TFR);

%% baseline correction
figure
cfg = [];
cfg.baseline     = [-2 0];
cfg.baselinetype = 'relative';
cfg.channel      = 'R016-2012-10-03-CSC04a';
ft_singleplotTFR(cfg, TFR);

%% Statistics
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'R016-2012-10-03-CSC04a';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100;
cfg.keeptrials   = 'yes'; % need this for stats later
cfg.t_ftimwin    = 20./cfg.foi;  % 20 cycles per time window
 
cfg.toi          = -2:0.05:0; % pre-nosepoke baseline
TFR_pre = ft_freqanalysis(cfg, data_trl);
 
cfg.toi          = 0:0.05:2; % post-nosepoke
TFR_post = ft_freqanalysis(cfg, data_trl);
 
TFR_pre.time = TFR_post.time; % time axes should be identical for comparison

%% t-test
cfg = [];
cfg.channel     = 'R016-2012-10-03-CSC04a';
cfg.latency     = 'all';
cfg.trials      = 'all';
cfg.frequency   = 'all';
    cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
cfg.avgoverfreq = 'no';
cfg.parameter   = 'powspctrm';
cfg.method      = 'stats';
cfg.statistic   = 'ttest2';
cfg.alpha       = 0.05;
 
nTrials1 = size(TFR_pre.powspctrm,1); nTrials2 = size(TFR_post.powspctrm,1);
 
cfg.design = cat(2,ones(1,nTrials1),2*ones(1,nTrials2)); % two conditions
cfg.ivar = 1; % dimension of design var which contains the independent variable (group)
 
stat = ft_freqstatistics(cfg,TFR_post,TFR_pre);
 
cfg.parameter = 'stat';
ft_singleplotTFR(cfg,stat); % plot the t-statistic

%% multichannel statistics
fc = FindFiles('*.ncs'); % get filenames of all LFPs recorded
data = ft_read_neuralynx_interp(fc); % load them all -- this will take a while
 
% define layout
cfg = [];
cfg.layout = 'ordered'; cfg.channel = data.label;
layout = ft_prepare_layout(cfg, data);

%% redo the stats
cfg = [];
cfg.trialfun = 'ft_trialfun_lineartracktone2';
cfg.trialdef.hdr = data.hdr;
cfg.trialdef.pre = 2.5;
cfg.trialdef.post = 5;
 
cfg.trialdef.eventtype = 'nosepoke'; % could be 'nosepoke', 'reward', 'cue'
cfg.trialdef.location = 'both'; % could be 'left', 'right', 'both'
cfg.trialdef.block = 'both'; % could be 'value', 'risk'
cfg.trialdef.cue = {'c1','c3','c5'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'}
 
[trl, event] = ft_trialfun_lineartracktone2(cfg);
cfg.trl = trl;
 
data_trl = ft_redefinetrial(cfg,data);
 
%%
cfg              = [];
cfg.output       = 'pow';
%cfg.channel      = 'R016-2012-10-03-CSC04a';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100;
cfg.keeptrials = 'yes'; % should be needed for subsequent statistics...
cfg.t_ftimwin    = 20./cfg.foi;  % 20 cycles per time window
cfg.toi          = -1:0.05:4;
 
TFR = ft_freqanalysis(cfg, data_trl);

%% Plot the outcome
figure
cfg = [];
cfg.baseline     = [-2 0];
cfg.baselinetype = 'relative';
cfg.layout = layout;
 
ft_multiplotTFR(cfg, TFR); % note this is now multiplot rather than singleplot