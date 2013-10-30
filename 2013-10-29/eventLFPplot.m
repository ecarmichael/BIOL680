function eventLFPplot(csc,event_times,varargin)
%
% INPUTS
%
% csc: [1 x 1] mytsd, LFP signal to be plotted
% event_times: [nEvents x 1] double with event times to align LFP on
%
% varargins (with defaults):
%
% t_window: [2 x 1] double indicating time window to use, e.g. [-1 3] for 1 second before to 3 seconds after event times

%% Variables
Fs = 2000;
t_window = [-1 3];
decimate_signal = 'yes'; % If you want it to deceimate the signal make this =1, or 0 if you want tyo retain the signal
dec_factor = 10; % this is the factor by which you want to decimate the data
LFP_colour = [1 0 0];  %default colour is red
filter_type = 'none'; % filter if needed.
filter_freq = [140 180];  % the desired frequency range of the filter.  Default is within the SWR band
filter_width = 1; % this is who much on either side of the bands you would like to include in the Ws
extract_varargin


%% for each event:
figure;        hold on
if max(event_times)>max(Range(csc))
    error('The event times are outside of the data range')
end
for ii = 1:length(event_times)
    % *extract the corresponding piece of LFP
    event_window_idx = Data(Restrict(csc,(event_times(ii)+t_window(1)),(event_times(ii)+t_window(2))));
    
    % *replace the original time axis with a new one based on the time window asked for
    event_window_tvec = Range(Restrict(csc,(event_times(ii)+t_window(1)),(event_times(ii)+t_window(2))));
    event_window_tvec = event_window_tvec-event_window_tvec(1)+t_window(1);
    
    % (Optional) Decimate the signals
    if strcmp(decimate_signal,'yes')==1 || strcmp(decimate_signal,'Yes')==1 || strcmp(decimate_signal,'YES')==1 ;
        event_window_tvec = downsample(event_window_tvec,dec_factor);  % I chose to downsample this because I am not worreid about frequency issues in the tvec
        event_window_idx = decimate(event_window_idx,dec_factor);
        dec_Fs = Fs/dec_factor; % though not used it is good to keep it in line.
    end
    
    % (Optional) filter the signal  NOTE:  The filter was placed after the
    % decimation since it is assumed based on the anti-aliasing component
    % of the decimation process, that any aliasing issues would be removed
    % prior to filtering below.  It is unclear what aliasing problems would
    % arrise if the signal was filter first and then decimated.  
    
%%%%%%%%%%%%This is the section that would have done the filtering but I could not get it to output anything remotely reasonable.  Perhaps I will try again later.      
% % % % %     if strcmp(filter_type,'butter')==1
% % % % %         Wp = [ filter_freq(1)  filter_freq(2)] * 2 / Fs;
% % % % %         Ws = [ filter_freq(1)-filter_width  filter_freq(2)+filter_width] * 2 / Fs;
% % % % %         [N,Wn] = buttord( Wp, Ws, 3, 20);
% % % % %         [B,A] = butter(N,Wn);
% % % % %         event_window_idx = filtfilt(B, A,event_window_idx);
% % % % %     end
% % % % %     if strcmp(filter_type,'chebyshev')==1
% % % % %         Wp = [ filter_freq(1)  filter_freq(2)] * 2 / Fs;
% % % % %         Ws = [ filter_freq(1)-filter_width  filter_freq(2)+filter_width] * 2 / Fs;
% % % % %        [N,Wp] = cheb1ord( Wp, Ws, 3, 20);
% % % % %        [B,A] = cheby1(N,0.5,Wp);
% % % % %         event_window_idx = filtfilt(B, A,event_window_idx);
% % % % %     end
    
    % *optionally, rescale the LFP
    
    %this is done during the plot
    
    % *add a y-offset to the LFP to plot one above the other
    offset = ii*1000;
    % *plot the current piece of LFP
    h(ii) = plot(event_window_tvec, event_window_idx-mean(event_window_idx)+offset,'Color', LFP_colour);
    % add a vertical line to the plot indicating time zero
    line(0, 0:10:max(event_window_idx-mean(event_window_idx)+offset)+.5*(offset/ii),'Color','k')
    % further niceties: add an option to decimate, to color, to filter before plotting
    
end