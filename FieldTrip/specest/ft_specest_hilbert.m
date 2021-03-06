function [spectrum,freqoi,timeoi] = ft_specest_hilbert_new(dat, time, varargin)

% FT_SPECEST_HILBERT performs a spectral estimation of data by repeatedly
% applying a bandpass filter and then doing a hilbert transform.
%
% Use as
%   [spectrum,freqoi,timeoi] = specest_hilbert(dat,time,...)
% where
%   dat      = matrix of chan*sample
%   time     = vector, containing time in seconds for each sample
%   spectrum = matrix of chan*freqoi*timeoi of fourier coefficients
%   freqoi   = vector of frequencies in spectrum
%   timeoi   = vector of timebins in spectrum
%
% Optional arguments should be specified in key-value pairs and can include:
%   timeoi    = vector, containing time points of interest (in seconds)
%   freqoi    = vector, containing frequencies (in Hz)
%   pad       = number, indicating time-length of data to be padded out to in seconds (used for spectral interpolation, NOT filtering)
%   padtype   = string, indicating type of padding to be used (see ft_preproc_padding, default: zero)
%   width     = number or vector, width of band-pass surrounding each element of freqoi
%   filttype  = string, filter type, 'but' or 'fir' or 'firls'
%   filtorder = number or vector, filter order
%   filtdir   = string, filter direction,  'twopass', 'onepass' or 'onepass-reverse' 
%   verbose   = output progress to console (0 or 1, default 1)
%   polyorder = number, the order of the polynomial to fitted to and removed from the data
%                  prior to the fourier transform (default = 0 -> remove DC-component)
%
% See also FT_FREQANALYSIS, FT_SPECEST_MTMFFT, FT_SPECEST_TFR, FT_SPECEST_MTMCONVOL, FT_SPECEST_WAVELET

% Copyright (C) 2010, Robert Oostenveld
%
% $Id: ft_specest_hilbert.m 7919 2013-04-17 13:44:49Z roevdmei $

% get the optional input arguments
freqoi    = ft_getopt(varargin, 'freqoi');
timeoi    = ft_getopt(varargin, 'timeoi', 'all');
width     = ft_getopt(varargin, 'width', 1);
filttype  = ft_getopt(varargin, 'filttype');    if isempty(filttype),  error('you need to specify filter type'),         end
filtorder = ft_getopt(varargin, 'filtorder');   if isempty(filtorder), error('you need to specify filter order'),        end
filtdir   = ft_getopt(varargin, 'filtdir');     if isempty(filtdir),   error('you need to specify filter direction'),    end
pad       = ft_getopt(varargin, 'pad');
padtype   = ft_getopt(varargin, 'padtype', 'zero');
polyorder = ft_getopt(varargin, 'polyorder', 0);
fbopt     = ft_getopt(varargin, 'feedback');
verbose   = ft_getopt(varargin, 'verbose', true);

if isempty(fbopt),
  fbopt.i = 1;
  fbopt.n = 1;
end

% Set n's
[nchan,ndatsample] = size(dat);

% Remove polynomial fit from the data -> default is demeaning
if polyorder >= 0
  dat = ft_preproc_polyremoval(dat, polyorder, 1, ndatsample);
end

% Determine fsample and set total time-length of data
fsample = 1./mean(diff(time));
dattime = ndatsample / fsample; % total time in seconds of input data

% Zero padding
if round(pad * fsample) < ndatsample
  error('the padding that you specified is shorter than the data');
end
if isempty(pad) % if no padding is specified padding is equal to current data length
  pad = dattime;
end

padlength  = floor((pad-dattime)* fsample)./2;

% set a default sampling for the frequencies-of-interest
if isempty(freqoi),
  freqoi = linspace(2*width, (fsample/3), 50);
end
% check for freqoi = 0 and remove it
if any(freqoi==0)
  freqoi(freqoi==0) = [];
end
nfreqoi = length(freqoi);

% Set timeboi and timeoi
timeoiinput = timeoi;
offset = round(time(1)*fsample);
if isnumeric(timeoi) % if input is a vector
  timeoi   = unique(round(timeoi .* fsample) ./ fsample);
  timeboi  = round(timeoi .* fsample - offset) + 1;
  ntimeboi = length(timeboi);
elseif strcmp(timeoi,'all') % if input was 'all'
  timeboi  = 1:length(time);
  ntimeboi = length(timeboi);
  timeoi   = time;
end

% throw a warning if input timeoi is different from output timeoi
if isnumeric(timeoiinput)
  if numel(timeoiinput) ~= numel(timeoi) % timeoi will not contain double time-bins when requested
    warning('output time-bins are different from input time-bins, multiples of the same bin were requested but not given');
  else
    if any(abs(timeoiinput-timeoi) >= eps*1e6)
      warning('output time-bins are different from input time-bins');
    end
  end
end


% expand width to array if constant width
if numel(width) == 1
  width = ones(1,nfreqoi) * width;
end

% expand filter order to array if constant filterorder
if numel(filtorder) == 1
  filtorder = ones(1,nfreqoi) * filtorder;
end


% create filter frequencies and check validity
filtfreq = [];
invalidind = [];
for ifreqoi = 1:nfreqoi
  tmpfreq = [freqoi(ifreqoi)-width(ifreqoi) freqoi(ifreqoi)+width(ifreqoi)];
  if all((sign(tmpfreq) == 1))
    filtfreq(end+1,:) = tmpfreq;
  else
      invalidind = [invalidind ifreqoi];
      warning(sprintf('frequency %.2f Hz cannot be estimated with resolution %.2f Hz', freqoi(ifreqoi), width(ifreqoi)));
  end
end

freqoi(invalidind) = [];
nfreqoi = size(filtfreq,1);

% preallocate the result and perform the transform
spectrum = complex(nan(nchan, nfreqoi, ntimeboi), nan(nchan, nfreqoi, ntimeboi));
for ifreqoi = 1:nfreqoi
  str = sprintf('frequency %d (%.2f Hz)', ifreqoi,freqoi(ifreqoi));
  [st, cws] = dbstack;
  if length(st)>1 && strcmp(st(2).name, 'ft_freqanalysis') && verbose
    % specest_convol has been called by ft_freqanalysis, meaning that ft_progress has been initialised
    ft_progress(fbopt.i./fbopt.n, ['trial %d, ',str,'\n'], fbopt.i);
  elseif verbose
    fprintf([str, '\n']);
  end
  
  % filter
  flt = ft_preproc_bandpassfilter(dat, fsample, filtfreq(ifreqoi,:), filtorder(ifreqoi), filttype, filtdir); 
  
  % transform and insert
  dum = transpose(hilbert(transpose(ft_preproc_padding(flt, padtype, padlength))));
  spectrum(:,ifreqoi,:) = dum(:,timeboi+padlength);
%   dum = transpose(hilbert(transpose([flt repmat(postpad,[nchan, 1])])));
%   spectrum(:,ifreqoi,:) = dum(:,timeboi);
end
