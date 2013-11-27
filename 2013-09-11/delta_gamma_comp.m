%%%%%%%%% Week 7 Assignment %%%%%%%%%%%%

% this script will extract the low-gamma and delta events from a given
% segment of data "csc".


%% Define any variables
delta_threshold = 5;
time_window = [2700 3300]; % restrict to the risk trial only

%% load and restrict the data
cd('D:\Promoted\R016-2012-10-03')
% remember to cd to data first
fname = 'R016-2012-10-03-CSC04a.Ncs';
csc = LoadCSC(fname);

cscR = Restrict(csc,time_window(1),time_window(2)); % risk session only
data = Data(cscR);
tvec = Range(cscR);
Fs = 2000;
% [S,F,T,P] = spectrogram(Data(cscR),hanning(512),384,1:0.25:200,Fs); % spectrogram

%% define frequency ranges of interest
delta = [3 4]; % delta
low_gamma = [45 65];% low-gamma

%% design filters for frequency ranges
%delta filter
Wp = delta * 2 / Fs;
Ws = [delta(1)-2 delta(2)+2] * 2 / Fs;
[N,Wn] = cheb1ord( Wp, Ws, 3, 20);
[delta_b_chev,delta_a_chev] = cheby1(N,0.5,Wn);

% low gamma filter
Wp = low_gamma * 2 / Fs;
Ws = [low_gamma(1)-5 low_gamma(2)+5] * 2 / Fs;
[N,Wn] = cheb1ord( Wp, Ws, 3, 20);
[gamma_b_chev,gamma_a_chev] = cheby1(N,0.5,Wn);
% fvtool(gamma_b_chev,gamma_a_chev)
%% filter the data (remember to use filtfilt!)
delta_filt = filtfilt(delta_b_chev,delta_a_chev,data);
gamma_filt = filtfilt(gamma_b_chev,gamma_a_chev,data);

%% isolate the "good delta" periods
delta_power = delta_filt.^2;
delta_power_filtered = medfilt1(delta_power,101);
delta_mean = nanmean(delta_power_filtered);
delta_sd = nanstd(delta_power_filtered );
delta_z = nan*zeros(length(delta_power_filtered ),1);
for ii  = 1:length(delta_power_filtered )
   delta_z(ii,1) = (delta_power_filtered (ii)-delta_mean)/delta_sd;
end

% Replace values below the thrshold z-score with zeros. 
subthreshold_delta = delta_z<delta_threshold;
good_delta = delta_filt; good_gamma = gamma_filt;
good_delta(subthreshold_delta==1)=0;
good_gamma(subthreshold_delta==1)=0;
%% extract delta phase and low gamma power
% delta
delta_phase = angle(hilbert(good_delta));

% gamma
gamma_pwr = abs(hilbert(good_gamma));

%plot to maintain sanity
[ax_h,h1,h2] = plotyy(good_tvec,gamma_pwr,good_tvec,delta_phase);

%% use averageXbyYbin to plot relationship (ideally with standard deviations)
phi_edges = -pi:pi/8:pi;
[x_avg, x_std] = averageXbyYbin(gamma_pwr,delta_phase,phi_edges);

pow_bin = x_avg;
 
pow_bin(end-1) = pow_bin(end-1)+pow_bin(end); % add counts on last edge to preceding bin
pow_bin = pow_bin(1:end-1); % trim
 
phi_centers = phi_edges(1:end-1)+pi/16; % convert edges to centers

errorbar(phi_centers,pow_bin, x_std);
