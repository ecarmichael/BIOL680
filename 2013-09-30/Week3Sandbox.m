%%%%% Week 3 sandbox %%%%%%%%%%

%%
fs1 = 2000;
tvec = 0:1/fs1:4; % construct time axis, sampled at 2kHz
 
freq1 = 100;
y = sin(2*pi*freq1*tvec); % 100Hz signal
 
ax1 = subplot(211);
stem(tvec,y); title('original');

%% subsample #1
subsample_factor = 4;
 
tvec2 = tvec(subsample_factor:subsample_factor:end); % take every 4th sample
y2 = y(subsample_factor:subsample_factor:end);
 
ax2 = subplot(212);
stem(tvec2,y2,'r'); title('subsampled');
xlabel('time (s)');

xl = [1 1.04];
set(ax1,'XLim',xl); set(ax2,'XLim',xl);

%% Reconstruct the subsampled data set

hold on;
 
y_interp = interp1(tvec2,y2,tvec,'linear');
p1 = plot(tvec,y_interp,'b');
 
y_interp2 = interp1(tvec2,y2,tvec,'spline');
p2 = plot(tvec,y_interp2,'g');
 
legend([p1 p2],{'linear','spline'},'Location','Northeast'); legend boxoff

%% Retry with a subsample of 10

fs1 = 2000;
tvec = 0:1/fs1:4; % construct time axis, sampled at 2kHz

freq1 = 100;
y = sin(2*pi*freq1*tvec); % 100Hz signal

subsample_factor = 10;

tvec2 = tvec(1:subsample_factor:end); % take every 4th sample
y2 = y(1:subsample_factor:end);

% ax2 = subplot(212);
stem(tvec2,y2,'r'); title('subsampled');
xlabel('time (s)');
legend([p1 p2],{'linear','spline'})

%% Retry 
figure(2)
subsample_factor = 4;
 
tvec2 = tvec(2:subsample_factor:end); % best case scenario -- can detect 100Hz signal
y2 = y(2:subsample_factor:end);
 
subplot(212)
stem(tvec2,y2,'r');
set(gca,'XLim',xl);

%% NEw test
%% 2kHz Fs, 100Hz signal with 450Hz signal superimposed
figure(3)
 
fs1 = 2000;
tvec = 0:1/fs1:4;
 
freq1 = 100;
freq2 = 450; % note, above Nyquist frequency for our target subsampled Fs
 
y = sin(2*pi*freq1*tvec) + 0.5.*sin(2*pi*freq2*tvec);
 
subplot(211)
stem(tvec,y)
set(gca,'XLim',xl);
%% ss -- we don't care about the 450Hz signal, but...
subsample_factor = 4;
 
tvec2 = tvec(1:subsample_factor:end);
y2 = y(1:subsample_factor:end);
 
subplot(212)
stem(tvec2,y2,'r');
set(gca,'XLim',xl);

%% using decimate instead of downsample

dec_tvec = decimate(tvec,4);
dec_y = decimate(y,4);
subplot(2,1,1)
stem(tvec2,y2,'r');
set(gca,'XLim',xl);
subplot(2,1,2)
stem(dec_tvec,dec_y,'b')
set(gca,'XLim',xl);


%% Loading function
cd('D:\Promoted\R016-2012-10-08');
fname = 'R016-2012-10-08-CSC03b.Ncs';
[Timestamps, ~, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(fname, [1 1 1 1 1], 1, 1, []);

% Why no spikes int he data file?  They have a frequency component in the
% 1000-5000Hz range and thus would not be sampled properly in this data
% file since it has a freqeuncy of 2000Hz

% The smallest voltage change that can be detected in a 16 bit system
% should be the total range of the voltage +/-1500uV --> 3000 divided by
% the maxium number of values in the 16 bit system (65536).  SO, 3000/65536
% = 0.0458 SO 0.458 bits per volt
BitsperVolt = str2double(Header{15,1}(13:end));

Samples_in_mV = Samples.*1000.*BitsperVolt;

Timestamps_s = Timestamps./1000000;

% Since the number of samples in the time series are less than the maximum
% that can be recorded during this time, there must be gaps.  Proof
plot(diff(Timestamps))

value_ids = (Timestamps_s>ExpKeys.TimeOnTrack(1)&Timestamps_s<ExpKeys.TimeOffTrack(1));
TimestampsValue = Timestamps(value_ids);
for ii = 1:size(Samples,1)
SamplesValue(ii,:) = Samples(ii,value_ids);
end
NumberOfValidSamplesValue = NumberOfValidSamples(value_ids);

%% Find the non






