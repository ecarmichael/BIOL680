%% filter example
NT = 10;
b = [0.9 0.8];
a = [ 1 0.5];
x = zeros(NT,1);
x(1) = 1;
y = filter(b,a,x);
stem(y);

%%
load count.dat;
x = count(:,1);
 a = 1; % a_0 is the (hidden) coefficient on the left side, in front of y(n)
b = [1/4 1/4 1/4 1/4]; % four b's of 1/4 each so we get the mean
 
y = filter(b,a,x); % x is the original signal, y the filtered version
t = 1:length(x);
plot(t,x,'-.',t,y,'-'), grid on
legend('Original','Filtered',2)

%% white noise butterworth filter
% set up time axis
Fs = 500;
tvec = 1:1/Fs:10;
 
% generate white noise
x = rand(length(tvec),1);
nSamples = 512;
% get PSD
[Porig,Forig] = pwelch(x,hanning(nSamples),nSamples/2,length(x),Fs);
 
% design filter
W1 =50/(Fs/2);
W2 = 100/(Fs/2);
[b,a] = butter(4,[W1 W2]);
y = filter(b,a,x);
 
% get PSD
[Pfilt,Ffilt] = pwelch(y,hanning(nSamples),nSamples/2,length(y),Fs);
 
% plot the resulting PSDs
subplot(121)
plot(Forig,10*log10(Porig));
ylim([-150 -20]); xlim([0 250])
subplot(122)
plot(Ffilt,10*log10(Pfilt));
ylim([-150 -20]); xlim([0 250])

%% Better butter
% design filter a better filter
Wp = [ 50 100] * 2 / Fs; % passband - between 50 and 100 Hz
Ws = [ 45 105] * 2 / Fs; % stopband
[N,Wn] = buttord( Wp, Ws, 3, 20); % determine filter parameters
[b2,a2] = butter(N,Wn); % builds filter
y = fvtool(b,a,b2,a2);

%% Chebyshev filter
Wp = [ 50 100] * 2 / Fs; 
Ws = [ 48 102] * 2 / Fs;
[N,Wn] = cheb1ord( Wp, Ws, 3, 20); 
[b_c1,a_c1] = cheby1(N,0.5,Wn);
fvtool(b2,a2,b_c1,a_c1)

%% Chebyshev2 filter
Wp = [ 50 100] * 2 / Fs; 
Ws = [ 48 102] * 2 / Fs;
[N,Wn] = cheb1ord( Wp, Ws, 3, 20); 
[b_c1,a_c1] = cheby2(N,0.5,Wn);
fvtool(b2,a2,b_c1,a_c1)

%% Filtering a realistic signal
Fs = 500; dt = 1./Fs;
t = [0 10];
tvec = t(1):dt:t(2)-dt;
 
s1 = sin(2*pi*80*tvec+pi/6);
s2 = sin(2*pi*40*tvec);
s = s1 + s2;
 
sf = filter(b_c1,a_c1,s);
figure
plot(tvec,s,'k',tvec,sf,'r--'); hold on;
legend({'original','filtered'});
xlim([0 0.2]);

%% Correct for the phase shift by running the filter forwards and backwards
sf = filtfilt(b_c1,a_c1,s);
 
plot(tvec,s,'k',tvec,sf,'r--'); hold on;
legend({'original','filtered'});
xlim([0 0.2]);


%% compare freq responses
Fs = 500; dt = 1./Fs;
t = [0 10];
tvec = t(1):dt:t(2)-dt;
 
x = rand(size(tvec)); % white noise input
[P,F] = pwelch(x,hanning(512),256,2^14,Fs);
 
y1 = filter(b_c1,a_c1,x);
[P1,F1] = pwelch(y1,hanning(512),256,2^14,Fs);
 
    y2 = filtfilt(b_c1,a_c1,x);
[P2,F2] = pwelch(y2,hanning(512),256,2^14,Fs);
 
plot(F,10*log10(P),F,10*log10(P1),F,10*log10(P2));
legend({'original','filter','filtfilt'});

%% Notch Filtering
[b,a] = butter(10, [59 61] * 2 / Fs, 'stop');
fvtool(b,a);

[z,p,k] = butter(10, [59 61] * 2 / Fs, 'stop'); % note, we ask for 3 outputs instead of 2
[sos,g] = zp2sos(z,p,k); % convert to SOS format
h = dfilt.df2sos(sos,g); % create filter object
fvtool(h);

%% notch filter test 
% set up time axis
Fs = 500;
tvec = 1:1/Fs:10;
 
% generate white noise
x = rand(length(tvec),1);
nSamples = 512;

% use the new notch with the filt filt function filtfilt(SOS, G, x)
[z,p,k] = butter(4, [59 61] * 2 / Fs, 'stop'); % note, we ask for 3 outputs instead of 2
[sos,g] = zp2sos(z,p,k); % convert to SOS format
y = filtfilt(sos, g, x) ;
[P,F] = pwelch(y,hanning(512),256,2^14,Fs);
plot(F,10*log10(P));

%% detecting movement artifacts
cd('D:\Promoted\R016-2012-10-08')
csc = LoadCSC('R016-2012-10-08-CSC02b.ncs');
cscR = Restrict(csc,1270,1272);
plot(cscR)

%% guess the noise bands and hope
x = Data(cscR);
tvec = Range(cscR);
 
Fs = 2000;
Wp = [ 180 220] * 2 / Fs;
Ws = [ 178 222] * 2 / Fs;
[N,Wn] = cheb1ord( Wp, Ws, 3, 20); % determine filter parameters
[b_c1,a_c1] = cheby1(N,0.5,Wn); % builds filter
 
fvtool(b_c1,a_c1); % remember to check your filter!
 
y = filtfilt(b_c1,a_c1,x);
subplot(311)
plot(tvec,x,'b',tvec,y,'r');
chew_power = y.^2;
subplot(312)
plot(tvec,x,'b',tvec,chew_power,'r');
subplot(313)
chew_power_filtered = medfilt1(chew_power,101); % filter window is specified in samples, so this is ~50ms
plot(tvec,x,'b',tvec,chew_power_filtered,'r');