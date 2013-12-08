%%%%%%%%%%%%%%%%%%%%% Week 4: Fourier Sandbox %%%%%%%%%%%%%%%%%%%%%

%% plot a simple sinusoid
Fs = 100; % in samples per second (Hz)
t0 = 0; t1 = 1; % start and end times
tvec = t0:1./Fs:t1; % construct time axis

f = 2; % frequency of sine to plot
y = sin(2*pi*f*tvec); % note sin() expects arguments in radians, not degrees (see sind())

stem(tvec,y);
%% Manipulate the phase and amplitude
phi = pi/2;

figure;
y = sin(2*pi*f*tvec + phi); % a phase shift
stem(tvec,y);
hold on;
plot(tvec,cos(2*pi*f*tvec),'LineWidth',2) % notice, cosine is simply phase shifted sine
legend('Sin (Phase Shifted)', 'Cos');
hold off;

a = 2;

figure;
y = a.*sin(2*pi*f*tvec + phi); % amplitude change
stem(tvec,y);

%% frequency modulation (FM) signal
f2 = 10;
m = 2;

subplot(311)
s1 = sin(2*pi*f*tvec);
plot(tvec,s1); title('message');

subplot(312);
s2 = sin(2*pi*f2*tvec);
plot(tvec,s2); title('carrier');

subplot(313);
s3 = sin(2*pi*f2*tvec + m.*sin(2*pi*f*tvec - pi/2));
plot(tvec,s3); title('FM signal');

%% Harmonic series example
mag = [0.1 0 1.3 0.5]; % magnitudes for each term
pha = [-pi/6 0 pi 2*pi/3]; % phases for each term
f = 2; % base frequency

signal_out = zeros(size(tvec));
for ii = 1:numel(mag) % note, the book chapter uses i, not best practice!
    
    this_signal = mag(ii)*cos(2*pi*f*ii*tvec + pha(ii));
    plot(tvec,this_signal,'r:'); hold on;
    signal_out = signal_out + this_signal; % build the sum
    
end
figure;
plot(tvec,signal_out,'LineWidth',2);

%% Why should you not use "i" as a variable.  i is an imaginary unit in MatLab and has a set value.  If you use i as a variable you will overwrite this value.  Same thing goes of pi

%% Decomposing a signal
x = round(rand(1,8)*10); % generate a length 8 vector of integers between 0 and 10
xlen = length(x);

% get magnitudes and phases of Fourier series
X = fft(x);
Xmag = abs(X); % magnitudes, a_n
Xphase = angle(X); % phases, phi_n

n = 0:xlen-1;
t = 0:0.05:xlen-1; % a finer timescale to show the smooth signal later

for iH = xlen-1:-1:0 % reconstruct each harmonic
    s(iH+1,:) = Xmag(iH+1)*cos(2*pi*n*iH/xlen + Xphase(iH+1))/xlen;
    sm(iH+1,:) = Xmag(iH+1)*cos(2*pi*t*iH/xlen + Xphase(iH+1))/xlen;
    % detail: xlen appears here because the fundamental frequency used by fft() depends on this
end

ssum = sum(s);
smsum = sum(sm);0

figure;
plot(n, x, 'go', t, smsum, 'b', n, ssum, 'r*');
legend({'original','sum - all','sum - points only'});

%% why use a backwards for loop.  If a for loop is adding data to a
% variable throughout a loop without preallocating a working space, it will
% have to take all the data from the previous loop and create a new set of
% data with the additional information so that each loop requires more and
% more memory.  By preallocating space, matlab can "fill" a matrix rather
% than rebuilding it with each loop.
% By using a reverse for loop this creates a variable with the maixmal
% length that can be filled with each pass rather than adding on to it.

%% Interpreting the output of an fft()
Fs = 20; % in samples per second (Hz)
t0 = 0; t1 = 1; % start and end times
tvec = t0:1/Fs:t1-(1/Fs); % construct time axis; generate exactly 20 samples
tvec = t0:1/Fs:t1;
f = 2; % signal frequency
y = sin(2*pi*f*tvec); % construct signal, a 2Hz sine wave sampled at 20Hz for 1s

yfft = fft(y,length(y));
yfft_mag = abs(yfft); yfft_ph = angle(yfft);
% yfft_mag = yfft; yfft_ph = angle(yfft);  %contains imaginary components
stem(yfft_mag)

Npoints = length(y);
F = [-Npoints/2:Npoints/2-1]./Npoints; % construct frequency axis

yfft_mag = fftshift(yfft_mag); % align output, see note below
stem(F,yfft_mag);

xlabel('Frequency (Fs^{-1})');
%%
tvec = t0:1/Fs:t1;
nPoints = [length(tvec) 64 256 1024];

for iP = 1:length(nPoints) % repeat fft with different numbers of points
    
    nP = nPoints(iP);
    subplot(2,2,iP);
    
    y = sin(2*pi*f*tvec);
    yfft = fft(y,nP);
    yfft_mag = abs(yfft); yfft_ph = angle(yfft);
    
    F = [-nP/2:nP/2-1]./nP;
    yfft_mag = fftshift(yfft_mag);
    plot(F,yfft_mag,'kx',F,yfft_mag,'k');
    
    title(sprintf('%d point FFT',nP));
    xlabel('Frequency (Fs^{-1})');
    
end

%% If the number of points is less than the length of the signal, the
% resolution is decreased

%% zero padding  test: results are the same when the last 10 values of the
% tvec are zero padded.
figure
tvec = t0:1/Fs:t1;
tvec(end-10) = 0;
nPoints = [length(tvec) 64 256 1024];

for iP = 1:length(nPoints) % repeat fft with different numbers of points
    
    nP = nPoints(iP);
    subplot(2,2,iP);
    
    y = sin(2*pi*f*tvec);
    yfft = fft(y,nP);
    yfft_mag = abs(yfft); yfft_ph = angle(yfft);
    
    F = [-nP/2:nP/2-1]./nP;
    yfft_mag = fftshift(yfft_mag);
    plot(F,yfft_mag,'kx',F,yfft_mag,'k');
    
    title(sprintf('%d point FFT',nP));
    xlabel('Frequency (Fs^{-1})');
    
end

%% Spectral leakage
tvec = t0:1/Fs:t1-(1/Fs);
nRepeats = [1 2 4 8];

nP =  1024;
for iP = 1:length(nRepeats)
    
    subplot(2,2,iP);
    
    y = sin(2*pi*f*tvec);
    y = repmat(y,[1 nRepeats(iP)]); % repeat the signal a number of times
    %     y_ham = y.*hamming(length(y))';  % uses a hamming window instead.
    yfft = fft(y,nP);    yfft_mag = abs(yfft); yfft_ph = angle(yfft);
    
    F = [-nP/2:nP/2-1]./nP;
    yfft_mag = fftshift(yfft_mag);
    plot(F,yfft_mag,'kx',F,yfft_mag,'k');
    
    title(sprintf('%d repeats',nRepeats(iP)));
    xlabel('Frequency (Fs^{-1})');
    
end

%% Windowing
nP = 25;
nPFFT = 1024;

windows = {'rectwin','triang','hamming','hanning','blackman'};
cols = 'rgbcmyk';

for iW = 1:length(windows)
    
    eval(cat(2,'wn = ',windows{iW},'(nP);')); % make sure you understand this
    wn = wn./sum(wn);
    
    subplot(211); % plot the window
    plot(wn,cols(iW),'LineWidth',2); hold on;
    
    subplot(212);
    yfft = fft(wn,nPFFT);
    yfft_mag = abs(yfft); yfft_ph = angle(yfft);
    
    F = [-nPFFT/2:nPFFT/2-1]./nPFFT;
    yfft_mag = fftshift(yfft_mag);
    
    h(iW) = plot(F,yfft_mag,cols(iW),'LineWidth',2); hold on;
    
end

xlabel('Frequency (Fs^{-1})');
legend(h,windows);
set(gca, 'YScale', 'log');  % create a log scale in the Y axis

%% The eval statement creates a set a values (length = nP) using of the
% windowing methods.  Example in the first case it sets "wn" to be equal
% to a 'rectwin(nP)' and then stores this value as the variable "wn'

%%  Verification that the hamming window does change the sidelobes.

figure(3)
tvec = t0:1/Fs:t1-(1/Fs);
nRepeats = [1 2 4 8];

nP =  1024;
for iP = 1:length(nRepeats)
    
    subplot(2,2,iP);
    
    y = sin(2*pi*f*tvec);
    y = repmat(y,[1 nRepeats(iP)]); % repeat the signal a number of times
    y_ham = y.*hamming(length(y))';  % uses a hamming window instead.
    yfft = fft(y_ham,nP);    yfft_mag = abs(yfft); yfft_ph = angle(yfft);
    
    F = [-nP/2:nP/2-1]./nP;
    yfft_mag = fftshift(yfft_mag);
    plot(F,yfft_mag,'kx',F,yfft_mag,'k');
    
    title(sprintf([' Hamming ' '%d repeats'],nRepeats(iP)));
    xlabel('Frequency (Fs^{-1})');
    
end

%% Robust spectral estimation methods
[Pxx,F] = periodogram(y,[],nP,Fs);
plot(F,Pxx); xlabel('Frequency (Hz)');
hold on;
[Pxx,F] = periodogram(y,hanning(length(y)),nP,Fs);
plot(F,Pxx,'r');

%% Pwelching the data for visualization and PSDs
Fs = 20; % in samples per second (Hz)
t0 = 0; t1 = 1;
f = 2;
nRepeats = 4;

tvec = t0:1/Fs:t1-(1/Fs);

nP =  1024;
y = sin(2*pi*f*tvec);
y = repmat(y,[1 nRepeats]);

[Pxx,F] = periodogram(y,rectwin(length(y)),nP,Fs);
plot(F,Pxx);

hold on;
wSize = 40;
[Pxx,F] = pwelch(y,rectwin(wSize),wSize/2,nP,Fs);
plot(F,Pxx,'r'); xlabel('Frequency (Hz)');

%% pitfalls of real world signals.
Fs = 20; % in samples per second (Hz)
t0 = 0; t1 = 1; f = 2;
nP =  1024;
gaps = [5 10 15]; % idx of samples to be removed

tvec = t0:1/Fs:t1;%-(1/Fs);
y = sin(2*pi*f*tvec);

subplot(211)
plot(tvec,y,'k*'); hold on;

yfft = fft(y,nP);
yfft_mag = abs(yfft); yfft_ph = angle(yfft);

F = [-nP/2:nP/2-1]./nP;
yfft_mag = fftshift(yfft_mag);

subplot(212);
plot(F,yfft_mag,'k-x'); hold on;

xlabel('Frequency (Fs^{-1})');

% signal with gaps
y = sin(2*pi*f*tvec);
y2 = y;
y2(gaps) = []; tvec(gaps) = []; % remove

subplot(211);
plot(tvec,y2,'bo'); hold on;

yfft = fft(y2,nP);
yfft_mag = abs(yfft); yfft_ph = angle(yfft);

F = [-nP/2:nP/2-1]./nP;
yfft_mag = fftshift(yfft_mag);

subplot(212);
plot(F,yfft_mag,'b-x');
legend('y','y2 (with gaps)')

%% Real data:
% cd to your R016-2012-10-08 folder
csc = LoadCSC('R016-2012-10-08-CSC04d.ncs');
run(FindFile('*keys.m'));
% restrict to prerecord, leaving some time (10s) before rat actually goes on track
csc_pre = Restrict(csc,0,ExpKeys.TimeOnTrack(1)-10);
csc_preR = Range(csc_pre);
csc_preD = Data(csc_pre);
csc_preDec = Data(csc_pre);
% check if sampling is ok
plot(diff(csc_preR)); % only minimal differences
Fs = 1./mean(diff(csc_preR));

% downsample
dsf = 4;
csc_preD_vs = decimate(csc_preD,dsf);
csc_preR_vs = downsample(csc_preDec,dsf);
Fs = Fs./dsf;

wSize = 1024;
[Pxx,F] = periodogram(csc_preD_vs,hamming(length(csc_preD)),length(csc_preD),Fs);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);

% Pwelch of the above data using a Hamming window
[Pxx,F] = pwelch(csc_preD_vs,hamming(wSize),wSize/2,wSize,Fs);
plot(F,10*log10(Pxx),'r'); xlabel('Frequency (Hz)');
xlim([0 150]);


%% %%%%%%%%%%%%%%%%%%%%%%% Assignment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% part 1: PSD of whitenoise
white = rand(6000,1);
wSize = 1024;
Fs = 2000;
[Pxx,F] = pwelch(white,hamming(wSize),wSize/2,wSize,Fs);
plot(F,Pxx,'r'); xlabel('Frequency (Hz)');
xlim([0 150]); ylim([-5 5]);

% since white noise is entirely random if should not have any dominate
% peaks and thus would yield a flat PSD. It does not follow the
% 1/f profile since the values at different frequencies do not lead to a
% change in the power.

%% part 2: HC PSD vs Vs
csc = LoadCSC('R016-2012-10-08-CSC02b.ncs'); %found via the ExpKeys
csc_pre = Restrict(csc,0,ExpKeys.TimeOnTrack(1)-10);
csc_preR = Range(csc_pre);
csc_preD = Data(csc_pre);
Fs = 1./mean(diff(csc_preR));
% Decimate
dsf = 4;
csc_preD_hc = decimate(csc_preD,dsf);
Fs = Fs./dsf;
% pwelch
wSize = 1024;
[Phc,F] = pwelch(csc_preD_hc,hamming(wSize),wSize/2,wSize,Fs);
hold on
plot(F,10*log10(Phc),'r'); xlabel('Frequency (Hz)');
plot(F,10*log10(Pxx),'b'); xlabel('Frequency (Hz)'); % Vs (assuming "goodGamma" was a Vs recording)
xlim([0 100]); legend('Hippocampus','Ventral Striatum')

%% window sizes
wSize = [256 512 4096 16384];
for ii = 1:length(wSize)
    subplot(2,2,ii)
    hold on
    [Phc,F] = pwelch(csc_preD_hc,hamming(wSize(ii)),wSize(ii)/2,wSize(ii),Fs);
    [Pxx,F] = pwelch(csc_preD_vs,hamming(wSize(ii)),wSize(ii)/2,wSize(ii),Fs);
    plot(F,10*log10(Phc),'r'); xlabel('Frequency (Hz)');
    plot(F,10*log10(Pxx),'b'); xlabel('Frequency (Hz)'); % Vs (assuming "goodGamma" was a Vs recording)
    xlim([0 100]); legend('Hippocampus','Ventral Striatum')
    title(['Sampling Window Size: ' num2str(wSize(ii))]);
end

% smaller window sizes will lead to poor resolution in the frequency range
% and will also attenuate the power at each sampled frequency since the
% power is being taken from a smaller window that will smooth power over a
% larger window.  This is not easy to put into words, but I am trying to
% say that when you have a 256 smaple window this will bin several data
% points together such as 56.000Hz and 56.005Hz and give an answer of say
% 20db, when really 56.000Hz is at 40db and 56.000Hz is much lower.  By
% increasing the sampling window this will allo for the separation of these
% smaller differences in the frequency.  This also does introduce more data
% points and leads to a more noisey looking signal.

%% decimate Vs downsample for the VS data
csc_preR = Range(csc_pre);
Fs = 1./mean(diff(csc_preR));
csc_preD = Data(csc_pre);
csc_preDec = Data(csc_pre);
dsf = 4;
csc_preD_vs = decimate(csc_preD,dsf);
csc_preR_vs = downsample(csc_preDec,dsf);
Fs = Fs./dsf;

wSize = [256 512 4096 16384];
for ii = 1:length(wSize)
    subplot(2,2,ii)
    hold on
    [Pdec,F] = pwelch(csc_preD_vs,hamming(wSize(ii)),wSize(ii)/2,wSize(ii),Fs);
    [Pdown,F] = pwelch(csc_preR_vs,hamming(wSize(ii)),wSize(ii)/2,wSize(ii),Fs);
    plot(F,10*log10(Pdec),'r'); xlabel('Frequency (Hz)');
    plot(F,10*log10(Pdown),'b'); xlabel('Frequency (Hz)'); % Vs (assuming "goodGamma" was a Vs recording)
    xlim([0 150]); legend('Decimated','Downsampled')
    title(['Sampling Window Size: ' num2str(wSize(ii))]);
end

%  The downsampled signal has slight differences in the power. This could
%  be due to the downsampling process creating artifacts
%  throughout the set that manifest as slightly different peaks across all the 
%  frequencies in the PSD.  As would be expected these artifacts are more 
%  apparent at higher frequncies, as can be seen >100Hz. The
%  size of the sampling window does not seem to have a clear effect on the
%  relationship between the downsampled and decimated data.