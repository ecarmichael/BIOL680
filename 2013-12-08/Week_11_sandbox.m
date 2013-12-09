%%%%%%%%%%%%%%%%%%% Week 10 Sandbox %%%%%%%%%%%%%%%%%%%%%%%%%


%% Spike triggered LFP
sd.fc = FindFiles('*.t');
sd.S = LoadSpikes(sd.fc);
 
csc = LoadCSC_ec('R016-2012-10-03-CSC02d.ncs');
 
cscR = Range(csc); cscD = Data(csc);
 
Fs = 2000; dt = 1./Fs;
cscD = locdetrend(cscD,Fs,[1 0.5]); % remove slow drifts in signal (this can mess up the STA)

w = [-1 1]; % time window to compute STA over
tvec = w(1):dt:w(2); % time axis for STA
 
iC = 2; % only do the third neuron for now
% clear sta;
 
spk_t = Data(sd.S{iC});
 
h = waitbar(0,sprintf('Cell %d/%d...',iC,length(sd.S)));
 
for iSpk = length(spk_t):-1:1 % for each spike...
 
   sta_t = spk_t(iSpk)+w(1);
   sta_idx = nearest(cscR,sta_t); % find index of leading window edge
 
   toAdd = cscD(sta_idx:sta_idx+length(tvec)-1); % grab LFP snippet for this window
   % note this way can be dangerous if there are gaps in the data
 
   sta{iC}(iSpk,:) = toAdd'; % build up matrix of [spikes x samples] to average later
 
   waitbar(iSpk/length(spk_t));
end
 
close(h);
%% plot it
plot(tvec,nanmean(sta{3}),'k','LineWidth',2); 
set(gca,'FontSize',14,'XLim',[-0.5 0.5]); xlabel('time (s)'); grid on;

%%  Asymetery in the STA:  The spike perfers to fire on the rise of the 
% theta wave.  This increase in the theta cycle could coincide with the
% increase in the membrane potential for the neuron(s) that happen to spike
% at this time and would suggest that they will fire once more readily when
% a certain membrane potential is achieved after a period of inactivity due
% to the lower membrane potential.  That is a guess and would not be the
% case for all neurons that are theta modulated, as many neurons could
% prefer even the trough of the theta wave.  It also does not fully explain
% why we do not see a burst in the activity during the times when the theta
% oscillation is above the initial possible spiking threshold.  

%% Faster STA
bin_edges = cscR+(dt/2);
len = length(tvec);
 
clear sta;
for iC = 1:3
 
spk_ts = Restrict(sd.S{iC},cscR(1)-w(1),cscR(end)-w(2));
spk_t = Data(spk_ts)+w(1); % times corresponding to start of window
 
[~,spk_bins] = histc(spk_t,bin_edges); % index into data for start of window
 
spk_bins2 = spk_bins(:, ones(len,1));
toadd = repmat(0:len-1,[length(spk_bins) 1]);
 
spk_bins3 = spk_bins2+toadd;
 
sta{iC} = cscD(spk_bins3);

end
%% Plot the avg STA for each spike
hold on
for iC = 1:3;
    subplot(3,1,iC)
    plot(tvec,nanmean(sta{iC}),'k','LineWidth',2); 
    xlabel(['Cell #' num2str(iC)])
end
hold off

% Cells 1 and 3 perfer to fire at different points within the oscillation
% Cell one seems to prefer the peak while cell 3 perfers the trough.  

%% extract the LFPs and Spikes
spike = ft_read_spike('R016-2012-10-03-TT02_2.t'); % needs fixed read_mclust_t.m
fc = {'R016-2012-10-03-CSC02b.ncs','R016-2012-10-03-CSC02d.ncs','R016-2012-10-03-CSC03d.ncs'};
data = ft_read_neuralynx_interp(fc);
data_all = ft_appendspike([],data, spike);

%% Plot
plot(data_all.time{1},data_all.trial{1}(1,:)) % a LFP
hold on;
plot(data_all.time{1},data_all.trial{1}(4,:)*500,'r') % binarized spike train (0 means no spike, 1 means spike)

%% Compute and plot the STA for neuron # 1
cfg              = [];
cfg.timwin       = [-0.5 0.5]; %
cfg.spikechannel = spike.label{1}; % first unit
cfg.channel      = data.label(1:3); % first 3 LFPs
staAll           = ft_spiketriggeredaverage(cfg, data_all);
 
% plot
figure
plot(staAll.time, staAll.avg(:,:)');
legend(data.label); h = title(cfg.spikechannel); set(h,'Interpreter','none');
set(gca,'FontSize',14,'XLim',cfg.timwin,'XTick',cfg.timwin(1):0.1:cfg.timwin(2)); 
xlabel('time (s)'); grid on;

%% split the data into trials
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
 
data_trl = ft_redefinetrial(cfg,data_all);

%% STA for a specific task segment
cfg              = [];
cfg.timwin       = [-0.5 0.5];
cfg.latency      = [-2.5 0];
cfg.spikechannel = spike.label{1}; % first unit
cfg.channel      = data.label(1:3); % first 3 LFPs
staPre           = ft_spiketriggeredaverage(cfg, data_trli);
 
% plot
figure
subplot(2,2,1:2)
plot(staPre.time, staPre.avg(:,:)');
legend(data.label); h = title(cfg.spikechannel); set(h,'Interpreter','none');
set(gca,'FontSize',14,'XLim',cfg.timwin,'XTick',cfg.timwin(1):0.1:cfg.timwin(2)); 
xlabel('time (s)'); grid on; title('Pre-NosePoke');

%post-nosepoke
cfg              = [];
cfg.timwin       = [-0.5 0.5];
cfg.latency      = [0 2.5];
cfg.spikechannel = spike.label{1}; % first unit
cfg.channel      = data.label(1:3); % first 3 LFPs
staPost           = ft_spiketriggeredaverage(cfg, data_trl);
%

% figure
subplot(2,2,3:4)
plot(staPost.time, staPost.avg(:,:)');
legend(data.label); h = title(cfg.spikechannel); set(h,'Interpreter','none');
set(gca,'FontSize',14,'XLim',cfg.timwin,'XTick',cfg.timwin(1):0.1:cfg.timwin(2)); 
xlabel('time (s)'); grid on; title('Post-NosePoke');

%% Pre Vs non-segmented STA
% The reward approach trials in the hippocampal recordings seem to prefer
% to fire just after the theta peak and in the ventral striatum these cells
% spike on the rise of a very weak theta oscillation that is not in phase
% with the hippocampal theta.  In addition there is more noise in the
% segmented sessions since less data is being used for the calculations.  

%% Removing an artifact
cfg              = [];
cfg.timwin       = [-0.002 0.006]; % remove 4 ms around every spike
cfg.spikechannel = spike.label{1};
cfg.channel      = data.label(2); % only remove spike in the second LFP ('02d')
cfg.method       = 'linear'; % remove the replaced segment with interpolation
data_trli        = ft_spiketriggeredinterpolation(cfg, data_trl);

%% Phase locking using FT
cfg              = [];
cfg.method       = 'convol';
cfg.foi          = 1:1:100;
cfg.t_ftimwin    = 5./cfg.foi; % 5 cycles per frequency
cfg.taper        = 'hanning';
cfg.spikechannel = spike.label{1};
cfg.channel      = data.label{1};
stsConvol        = ft_spiketriggeredspectrum(cfg, data_trli);
plot(stsConvol.freq,nanmean(sq(abs(stsConvol.fourierspctrm{1}))))
hist(angle(stsConvol.fourierspctrm{1}(:,:,9)),-pi:pi/18:pi)
%% Pairwise phase consistency
cfg               = [];
cfg.method        = 'ppc0'; % compute the Pairwise Phase Consistency
cfg.spikechannel  = stsConvol.label;
cfg.channel       = stsConvol.lfplabel; % selected LFP channels
cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
cfg.timwin        = 'all'; % compute over all available spikes in the window
cfg.latency       = [-2.5 0];
statSts           = ft_spiketriggeredspectrum_stat(cfg,stsConvol);

% plot the results
figure
plot(statSts.freq,statSts.ppc0')
xlabel('frequency')
ylabel('PPC')

%% pairwise phase consistency with artifact removal
 %%%%%%%%%%%%%% Comeback to this %%%%%%%%%%%%%%%
 
%% PPC using a sliding window
param = 'ppc0'; % set the desired parameter
 
cfg                = [];
cfg.method         = param;
 
cfg.spikechannel  = stsConvol.label;
cfg.channel       = stsConvol.lfplabel; % selected LFP channels
cfg.avgoverchan    = 'unweighted';
cfg.winstepsize    = 0.01; % step size of the window that we slide over time
cfg.timwin         = 0.5; % duration of sliding window
statSts            = ft_spiketriggeredspectrum_stat(cfg,stsConvol);
 
figure
cfg            = [];
cfg.parameter  = param;
cfg.refchannel = statSts.labelcmb{1,1};
cfg.channel    = statSts.labelcmb{1,2};
cfg.xlim       = [-2 3]; cfg.ylim = [2 30];
ft_singleplotTFR(cfg, statSts)

%% Spike-Field Coherence
cfg              = [];
cfg.output       = 'powandcsd';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100; % frequencies to use
cfg.t_ftimwin    = 5./cfg.foi;  % frequency-dependent, 5 cycles per time window
cfg.keeptrials   = 'yes';
cfg.channel      = {'R016-2012-10-03-CSC02b', 'R016-2012-10-03-TT02_2'};
cfg.channelcmb   = {'R016-2012-10-03-CSC02b', 'R016-2012-10-03-TT02_2'}; % channel pairs to compute csd for
 
cfg.toi          = -2:0.05:3;
 
TFR_pre = ft_freqanalysis(cfg, data_trl);
 
cfg            = [];
cfg.method     = 'ppc'; % compute coherence; other measures of connectivity are also available
fd             = ft_connectivityanalysis(cfg,TFR_pre);
 
iC = 1; % which signal pair to plot
lbl = [fd.labelcmb{1,:}]; % get the label of this pair
imagesc(fd.time,fd.freq,sq(fd.ppcspctrm(iC,:,:))); axis xy; colorbar
xlabel('time (s)'); ylabel('Frequency (Hz)'); title(lbl);

%% Phase precession
csc = myLoadCSC('R016-2012-10-03-CSC02d.ncs');
run(FindFile('*keys.m'));
csc = Restrict(csc,ExpKeys.TimeOnTrack(1),ExpKeys.TimeOffTrack(1));
cscR = Range(csc); cscD = Data(csc);
 
% filter in theta range
Fs = 2000;
Wp = [ 6 10] * 2 / Fs;
Ws = [ 4 12] * 2 / Fs;
[N,Wn] = cheb1ord( Wp, Ws, 3, 20); % determine filter parameters
[b_c1,a_c1] = cheby1(N,0.5,Wn); % builds filter
 
csc_filtered = filtfilt(b_c1,a_c1,cscD);
phi = angle(hilbert(csc_filtered));

% get the theta phase for each spike
S = LoadSpikes({'R016-2012-10-03-TT02_1.t'});
spk_t = Data(S{1});
spk_phi = interp1(cscR,phi,spk_t,'nearest');
 
hist(spk_phi,-pi:pi/18:pi)
% plot as a function of position
[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 1 1 1 1 1], 1, 1, []);
Timestamps = Timestamps*10^-6;
toRemove = (X == 0 & Y == 0);
X = X(~toRemove); Y = Y(~toRemove); Timestamps = Timestamps(~toRemove);
 
spk_x = interp1(Timestamps,X,spk_t,'linear');
spk_y = interp1(Timestamps,Y,spk_t,'linear');
 
plot(X,Y,'.','Color',[0.5 0.5 0.5],'MarkerSize',1); axis off; hold on;
h = scatterplotC(spk_x,spk_y,spk_phi,'Scale',[-pi pi],'solid_face',1,'plotchar','.');