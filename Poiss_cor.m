function [xcorr1, xbin1] = Poiss_cor(Spike_data,cell1_id, cell2_id, varargin)
%% This function will take two cells and generate the spike density 
% function of each cell and then use thier firing rates over a time (tvec)
% to generate an inhomogenous poisson process.  These artificial Poisson 
% data and then used to test if any cross-correlation exists between the 
% two artificial sets.  

% Inputs:
% S: spike data
% Cell1_id: a numeric value for the identity of the first neuron to be used
% Cell2_id: a numeric value for the identity of the second neuron to be used
%
% Outputs.  
% xcorr1: the correlation estimates for the two Poisson processes
% xbins1: the bin centers in seconds

%% Initialize Variables

t_start = 3200; % start of the session in time (sec)
t_end = 5650;  % end of the recording session in sec
ccf_binsize = 0.01; % bin size used to calculate the cross correlation estimate.  
binsize = 0.001; % what size you would like the bins to be for the SDF histograms
gau_SD = .050; % in ms
fig = 2;  % this is a boolean that will display the verification figures of the SDF
extract_varargin;
%% for each of the two neurons, restrict the data to [3200 5650] (the time interval when the rat was running on the track)
 spk_t1 = Data(Restrict(Spike_data{cell1_id},3200,5650));
 spk_t2 = Data(Restrict(Spike_data{cell2_id},3200,5650));
%% compute the spike density function for each, making sure that your tvec runs from 3200 to 5650 also, and that you have a 50ms SD for the Gaussian convolution kernel

tbin_edges = t_start:binsize:t_end;
tbin_centers = tbin_edges(1:end-1)+binsize/2;
 
spk_count1 = histc(spk_t1,tbin_edges);
spk_count1 = spk_count1(1:end-1);
spk_count2 = histc(spk_t2,tbin_edges);
spk_count2 = spk_count2(1:end-1);

% apply the guassian filter
gauss_window = 1./binsize; % 1 second window
gauss_SD = gau_SD./binsize; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize
spk_gau1 = conv2(spk_count1,gk,'same');
spk_gau2  = conv2(spk_count2,gk,'same');

if fig ==1;
    main = figure(3);
    subplot(2,1,1)
    bar(tbin_centers,spk_gau1); xlim([3200 5650])
    subplot(2,1,2)
    bar(tbin_centers,spk_gau2); xlim([3200 5650])
end

%% to use these SDFs to generate Poisson spike trains, convert the firing rates given by the SDF to a probability of emitting a spike in a given bin. (As you did above for a 0.47 Hz constant firing rate.)
% % % get the ISIs
% % isi_1 = diff(spk_gau1); % interspike intervals
% % isi_2 =  diff(spk_gau2);
% % isi_edges = 0:binsize:0.25; % bin edges for ISI histogram
% % isi_centers = isi_edges(1:end-1)+dt/2; % for plotting
% %  
% % isih_1 = histc(isi_1,isi_edges);
% % isih_2 = histc(isi_2,isi_edges);
% % 
% % if fig ==1;
% %     ISI_fig = figure;
% %     subplot(2,1,1)
% %     bar(isi_centers,isih_1(1:end-1));
% %     xlim([0 0.25])
% %     subplot(2,1,2)
% %     bar(isi_centers,isih_2(1:end-1));
% %      xlim([0 0.25])
% % end
%% generate Poisson spike trains, making sure to use the same tvec
tvec = tbin_edges(1:end-1);
dt = 0.001;
tvec = t_start:dt:t_end;
tvec = tvec(1:end-1);
pspike_1 = (spk_gau1.*10^-3)'; % probability of generating a spike in bin
rng default; % reset random number generator to reproducible state
spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx_1 = find(spk_poiss < pspike_1); % index of bins with spike
spk_poiss_t_1 = tvec(spk_poiss_idx_1)'; % use idxs to get corresponding spike time
% line([spk_poiss_t_1 spk_poiss_t_1],[-1 -0.5],'Color',[0 1 0]); % note, plots all spikes in one command

pspike_2 = (spk_gau2.*10^-3)'; % probability of generating a spike in bin
rng default; % reset random number generator to reproducible state
spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx_2 = find(spk_poiss < pspike_2); % index of bins with spike
spk_poiss_t_2 = tvec(spk_poiss_idx_2)'; % use idxs to get corresponding spike time
% line([spk_poiss_t_2 spk_poiss_t_2],[-1 -0.5],'Color',[0 0 0]); % note, plots all spikes in one command

%% convert Poisson spike trains to ts objects and compute the ccf

spk_poiss_ts_1 = ts(spk_poiss_t_1); spk_poiss_ts_2 = ts(spk_poiss_t_2); 
spk_ts_1 = ts(spk_t1); spk_ts_2 = ts(spk_t2); 
[xcorr1,xbin1] = ccf(spk_poiss_ts_1,spk_poiss_ts_2,ccf_binsize,1);
[xcorr2,xbin2] = ccf(spk_ts_1,spk_ts_2,ccf_binsize,1);

figure; subplot(211); plot(xbin1,xcorr1);
set(gca,'FontSize',10); xlabel('lag (s)'); ylabel('xcorr');
title(sprintf('Poisson Signal %d-%d',cell1_id,cell2_id));
subplot(212); plot(xbin2,xcorr2);
set(gca,'FontSize',10); xlabel('lag (s)'); ylabel('xcorr');
title(sprintf('Original Signal %d-%d',cell1_id,cell2_id));

end
