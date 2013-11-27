function evt = detectSWR(csc,varargin)
% function evt = detectSWR(csc,varargin)
%
% detect putative sharp wave-ripple events in csc input signal
%
% INPUTS
% csc: [1 x 1 mytsd]
%
% varargins (with defaults):
% ripple_band = [140 180]; % frequency band to use
% threshold = 5; % number of SDs above mean to use for detection
%
% OUTPUTS
% evt: [1 x 1 struct] with fields
% .t: [1 x nEvents double] times (in s) of events
% .pwr: [1 x nEvents double] power of events (in SDs above mean)

%% initialize default variables
ripple_band = [140 180]; % frequency band to use
threshold = 5; % number of SDs above mean to use for detection
Fs = 2000;
time_window = []; % the range of times that you wish to restrict the data to.  
extract_varargin;

%% Load the signal
cscL = LoadCSC(csc);
tvec = Range(cscL);
if isempty(time_window)==1
    round_tvec = round(tvec);
    distanace_from_session_end = 900; % corresponds roughly to 15mins
    start_time = find(round_tvec == (round(tvec(end)))-distanace_from_session_end);
end
cscR = Restrict(cscL,tvec(start_time(1)),tvec(end));
% cscR = Restrict(cscL,tvec(end)-900,tvec(end)); % rough values for 5mins from the end of the session. 
cscD = Data(cscR);
tvec = Range(cscR);
%% filter in the ripple band
% Chebyshev attempt using the ripple range from Logothetis et al., 2012
Wp = [ripple_band] * 2 / Fs;
Ws = [ripple_band(1)-2 ripple_band(2)+2] * 2 / Fs;
[N,Wn] = cheb1ord( Wp, Ws, 3, 20); 
[b_chev,a_chev] = cheby1(N,0.5,Wn);

y = filtfilt(b_chev,a_chev,cscD);
subplot(311)
plot(tvec,cscD,'b',tvec,y,'r');
subplot(312)
hold on
plot_range = 6025:6027;
plot_ind = find(tvec < plot_range(end) & tvec>plot_range(1));
to_plot = plot_ind(1):plot_ind(end);
plot(tvec(to_plot),cscD(to_plot),'b',tvec(to_plot),y(to_plot),'r');

%% convert to power envelope
ripple_power = y.^2;
ripple_power_filtered = medfilt1(ripple_power,51); % this value was selected since the 401 value from the 200ms sample would cause the cumputer to freeze.  
plot(tvec(to_plot),ripple_power_filtered(to_plot),'g');
hold off
%% convert to z-score (SDs from the mean; use nanmean() and nanstd())
ripple_mean = nanmean(ripple_power_filtered);
ripple_sd = nanstd(ripple_power_filtered);
ripple_z = nan*zeros(length(ripple_power_filtered),1);
for ii  = 1:length(ripple_power_filtered)
    ripple_z(ii,1) = (ripple_power_filtered(ii)-ripple_mean)/ripple_sd;
end
subplot(313)
hold on
plot(tvec(to_plot),ripple_z(to_plot),'g');
%% find times when above threshold
ripple_detect = ripple_z>threshold;
plot(tvec(to_plot),ripple_detect(to_plot),'-k');


%% find crossings from below to above threshold and vice versa (can use diff())
ripple_diff = diff(ripple_detect); 
plot(tvec(to_plot),ripple_diff(to_plot),'-m');
ripple_start = ripple_diff >0;
ripple_end = ripple_diff <0;
hold off

%% get center time and power, return
xx = 1; yy = 1;
for ii = 1:length(ripple_start)
    if ripple_start(ii) == 1;
        ripple_start_point(xx) = ii;
        xx = xx+1;
    end
    if ripple_end(ii) == 1;
        ripple_end_point(yy) = ii;
        yy = yy+1;
    end
end
nRipples = length(ripple_start_point); % used for tracking the number of SWR within the desired data range.  
for ii = 1:length(ripple_start_point)
    ripple_times(ii) = round(median(ripple_start_point(ii):ripple_end_point(ii)));
end
evt.t = tvec(ripple_times);
evt.pwr =  ripple_power_filtered(ripple_times);
end