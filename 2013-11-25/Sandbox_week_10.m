%%%%%%%%%%%%%%%%% Week 10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd 'D:\Promoted\R042-2013-08-18';
 
sd.fc = FindFiles('*.t');
sd.fc = cat(1,sd.fc,FindFiles('*._t')); % also load poorly isolated cells (if you don't have these, get from database)
sd.S = LoadSpikes(sd.fc);

% remember to unzip VT1.zip if not yet done
[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 1 1 1 1 1], 1, 1, []);
Timestamps = Timestamps*10^-6;

% remove zero samples, which happen when no LED is detected
toRemove = (X == 0 & Y == 0);
X = X(~toRemove); Y = Y(~toRemove); Timestamps = Timestamps(~toRemove);

sd.x = tsd(Timestamps,X');
sd.y = tsd(Timestamps,Y');

% the vido sampling frequency seems to be 30fps

%% Extract the main run period
t = [3250 5650];
for iC = 1:length(sd.S)
    sd.S{iC} = Restrict(sd.S{iC},t(1),t(2));
end
sd.x = Restrict(sd.x,t(1),t(2));
sd.y = Restrict(sd.y,t(1),t(2));

%% plot the position data
plot(Data(sd.x),Data(sd.y),'.','Color',[0.5 0.5 0.5],'MarkerSize',1);
axis off; hold on;
%% add in the spiking data from one cell
% get x and y coordinate for times of spike
iC = 5;
spk_x = interp1(Range(sd.x),Data(sd.x),Data(sd.S{iC}),'linear');
spk_y = interp1(Range(sd.y),Data(sd.y),Data(sd.S{iC}),'linear');
 
h = plot(spk_x,spk_y,'.r');

%% bin the spikes in space
SET_xmin = 10; SET_ymin = 10; SET_xmax = 640; SET_ymax = 480;
SET_nxBins = 63; SET_nyBins = 47;
 
spk_binned = ndhist(cat(1,spk_x',spk_y'),[SET_nxBins; SET_nyBins],[SET_xmin; SET_ymin],[SET_xmax; SET_ymax]);
 
imagesc(spk_binned')
axis xy; colorbar

occ_binned = ndhist(cat(1,Data(sd.x)',Data(sd.y)'),[SET_nxBins; SET_nyBins],[SET_xmin; SET_ymin],[SET_xmax; SET_ymax]);
 
% this is a sample count, so need to convert to seconds (1/30s per sample) to get time
VT_Fs = 30;
tc = spk_binned./(occ_binned .* (1./VT_Fs)); % firing rate is spike count divided by time
 
pcolor(tc'); shading flat
axis xy; colorbar; axis off

%% bin at a finer scale

SET_nxBins = 630; SET_nyBins = 470; % more bins
 
kernel = gausskernel([30 30],8); % 2-D gaussian for smoothing: 30 points in each direction, SD of 8 bins
 
% spikes
spk_binned = ndhist(cat(1,spk_x',spk_y'),[SET_nxBins; SET_nyBins],[SET_xmin; SET_ymin],[SET_xmax; SET_ymax]);
spk_binned = conv2(spk_binned,kernel,'same');
 
% occupancy
occ_binned = ndhist(cat(1,Data(sd.x)',Data(sd.y)'),[SET_nxBins; SET_nyBins],[SET_xmin; SET_ymin],[SET_xmax; SET_ymax]);
occ_binned = conv2(occ_binned,kernel,'same');
 
occ_binned(occ_binned < 0.01) = 0;  % this removes the minutely small values that will cause 
tc = spk_binned./(occ_binned .* (1 / VT_Fs));
tc(isinf(tc)) = NaN;  % this ensures that NaN values are not included in the bins which would cause then to appear to have a high firing rate since they are being 
 
figure;
pcolor(tc'); shading flat; axis off
axis xy; colorbar


%% updated tuning curve
kernel = gausskernel([4 4],2); % 2-D gaussian, width 4 bins, SD 2
 
SET_xmin = 10; SET_ymin = 10; SET_xmax = 640; SET_ymax = 480;
SET_nxBins = 63; SET_nyBins = 47;
 
spk_binned = ndhist(cat(1,spk_x',spk_y'),[SET_nxBins; SET_nyBins],[SET_xmin; SET_ymin],[SET_xmax; SET_ymax]);
spk_binned = conv2(spk_binned,kernel,'same'); % smoothing
 
occ_binned = ndhist(cat(1,Data(sd.x)',Data(sd.y)'),[SET_nxBins; SET_nyBins],[SET_xmin; SET_ymin],[SET_xmax; SET_ymax]);
occ_mask = (occ_binned < 5);
occ_binned = conv2(occ_binned,kernel,'same'); % smoothing
 
occ_binned(occ_mask) = 0; % don't include bins with less than 5 samples
 
VT_Fs = 30;
tc = spk_binned./(occ_binned .* (1 / VT_Fs));
tc(isinf(tc)) = NaN;
 
pcolor(tc'); shading flat
axis xy; colorbar; axis off

%% do the same for all cells in the data set
clear tc
nCells = length(sd.S);
for iC = 1:nCells
    spk_x = interp1(Range(sd.x),Data(sd.x),Data(sd.S{iC}),'linear');
    spk_y = interp1(Range(sd.y),Data(sd.y),Data(sd.S{iC}),'linear');
 
    spk_binned = ndhist(cat(1,spk_x',spk_y'),[SET_nxBins; SET_nyBins],[SET_xmin; SET_ymin],[SET_xmax; SET_ymax]);
    spk_binned = conv2(spk_binned,kernel,'same');
 
    tc = spk_binned./(occ_binned .* (1 / VT_Fs));
    tc(isinf(tc)) = NaN;
 
    sd.tc{iC} = tc;
 
end

%% inspect the new tuning curves

ppf = 25; % plots per figure
for iC = 1:length(sd.S)
    nFigure = ceil(iC/ppf);
    figure(nFigure);
 
    subplot(5,5,iC-(nFigure-1)*ppf);
    pcolor(sd.tc{iC}); shading flat; axis off;
 
end

% though I do not have a cell that is completely green, it out likely imply
% that there is an interneuron that is firing at all locations with the
% same probability.  

%% get the (n) spike counts
clear Q;
binsize = 0.25;
tvec = t(1):binsize:t(2);
for iC = length(sd.S):-1:1
 
    spk_t = Data(sd.S{iC});
    Q(iC,:) = histc(spk_t,tvec);
 
end
nActiveNeurons = sum(Q > 0);

%%
% look at it
figure
imagesc(tvec,1:nCells,Q)

%% unwraping the 
clear tc
nBins = numel(occ_binned);
nCells = length(sd.S);
for iC = nCells:-1:1
    tc(:,:,iC) = sd.tc{iC};
end
tc = reshape(tc,[size(tc,1)*size(tc,2) size(tc,3)]);
occUniform = repmat(1/nBins,[nBins 1]);

%% decoding algorithm
len = length(tvec);
p = nan(length(tvec),nBins);
for iB = 1:nBins
    tempProd = nansum(log(repmat(tc(iB,:)',1,len).^Q));
    tempSum = exp(-binsize*nansum(tc(iB,:),2));
    p(:,iB) = exp(tempProd)*tempSum*occUniform(iB);
end

%% renormalize and set the posterior time windows.  
p = p./repmat(sum(p,2),1,nBins);
p(nActiveNeurons < 1,:) = 0;

% rats actual position
xBinEdges = linspace(SET_xmin,SET_xmax,SET_nxBins+1);
yBinEdges = linspace(SET_ymin,SET_ymax,SET_nyBins+1);
 
xTempD = Data(sd.x); xTempR = Range(sd.x);
yTempD = Data(sd.y);
 
gS = find(~isnan(xTempD) & ~isnan(yTempD));
xi = interp1(xTempR(gS),xTempD(gS),tvec,'linear');
yi = interp1(xTempR(gS),yTempD(gS),tvec,'linear');
 
xBinned = (xi-xBinEdges(1))./median(diff(xBinEdges));
yBinned = (yi-yBinEdges(1))./median(diff(yBinEdges));

%% visualize it all 
goodOccInd = find(occ_binned > 0);
for iT = 1:size(p,1)
    cla;
    temp = reshape(p(iT,:),[SET_nxBins SET_nyBins]);
    toPlot = nan(SET_nxBins,SET_nyBins);
    toPlot(goodOccInd) = temp(goodOccInd);
    
    pcolor(toPlot); axis xy; hold on; caxis([0 0.5]);
    shading flat; axis off;
    
    hold on;
    h = figure; set(h,'Position',[100 100 320 240]);
    plot(yBinned(iT),xBinned(iT),'ow','MarkerSize',15);
 
    h = title(sprintf('t %.2f, nCells %d',tvec(iT),nActiveNeurons(iT))); 
    if nActiveNeurons(iT) == 0
        set(h,'Color',[1 0 0]);
    else
        set(h,'Color',[0 0 0]);
    end
   f(iT) = getframe(gcf); % store current frame
   drawnow;
end
% why is there an active neuron yet there is no actual display on the heat
% map? - this could be due to the activity of the cell being very low and
% thus not appearing to be firing on the blue back drop.  If we had set a
% higher threshold for the display of the cell...

%% export it to a movie file.  