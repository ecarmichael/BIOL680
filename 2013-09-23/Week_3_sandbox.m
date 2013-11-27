%%%%%%%% BIOL 680 Week 3 sandbox  %%%%%%%%%%%%%%

%% create a sample data set of 100Hz sampled at 2kHz
fs1 = 2000;
tvec = 0:1/fs1:4; % construct time axis, sampled at 2kHz
 
freq1 = 100;
y = sin(2*pi*freq1*tvec); % 100Hz signal
 
ax1 = subplot(211);
stem(tvec,y); title('original');
%%  Naive downsampling
subsample_factor = 4;
 
tvec2 = tvec(subsample_factor:subsample_factor:end); % take every 4th sample
y2 = y(subsample_factor:subsample_factor:end);
 
ax2 = subplot(212);
stem(tvec2,y2,'r'); title('subsampled');
xlabel('time (s)');