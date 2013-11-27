function neuroplot(spikes,csc,varargin)
% function neuroplot(spikes,csc,varargin)
%
% inputs:
%
% spikes: {nCells x 1} cell array of ts objects (spike train)
% csc: {nCSCs x 1} cell array of tsd objects (LFPs)
%
% varargins:
%
% cscColor: [nCSC x 3] RGB values to plot CSCs, default []
% spikeColor: [nCells x 3] RGB values to plot spikes, default []
%
% evt: [nEvents x 1] event times to plot, default []
% evtColor: [nEvents x 3] RGB values to plot events, default []
%
% interactiveMode: boolean to enable/disable arrow key navigation
%
%default_zoom:  [nCells x 2]  This determines the section that will be
%magnified in the plot.  By default it is set to 4000 - 4002s.  

%% Initilize default values
cscColor = [.7 .7 .7];
spikeColor = [1 0 0];
evtColor = [0 0 0];
evt = 0:10:max(Data(spikes{1}));
interactiveMode = 1; % for the interactive figure navigator 1 is on and 0 is off.  
default_zoom = [4000 4002];
extract_varargin;  % this will allow for modifications to the default values

%% extract the spike times from the S cell of ts values
spike_times = cell(1,length(spikes));
for ss = 1:length(spikes)
    spike_times{ss} = Data(spikes{ss,1}); %this should adjust the timescale to being in seconds
end



%% make a raster plot

%this portion will determine if the interactive figure viewer will be used.
if interactiveMode == 1
    f_hdl = figure('KeyPressfcn',@figure_move);
else
    figure
end

hold on
for ss = 1:length(spike_times)
    for ii = 1:length(spike_times{ss})
        line_plot = line([spike_times{ss}(ii) spike_times{ss}(ii)], [(-.5+ss) (0.5+ss)]);
        set(line_plot,'Color', spikeColor)  % this colour is set by the vara
    end
end

%% get the LFP data
data = Data(csc);
time = Range(csc);

%adjust the LFP data so that it is visable in the raster
data_mod = data./500 +(length(spike_times)+5);
csc_plot = plot(time,data_mod);
set(csc_plot,'Color', cscColor)

%% Add the event lines
for ev = 1:length(evt)
    line(evt(ev),'Color',evtColor)
end

%% Make some restrictions and alter the figure.
set(gcf,'Color',[1 1 1])
set(gca,'Visible','off','XLim', default_zoom) %takes the specified event time and plots 1s on either side

% set(gca,'YTick',[1:length(spike_times{ss})])









