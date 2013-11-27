%%%%%%%%%%%%%  BIOL 680 Week 2 Sandbox %%%%%%%%%%%%%%%%%%%%

%% Initialize

addpath('D:\Users\mvdmlab\My_Documents\GitHub\BIOL680\2013-09-23');
cd('D:\Promoted\R042-2013-08-18');
%% load the data (note, may need to unzip position data first)
fc = FindFiles('*.t');
S = LoadSpikes(fc);

[csc,csc_info] = LoadCSC('R042-2013-08-18-CSC03a.ncs');

[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 1 1 1 1 1], 1, 1, [] );

%% Data verification
data = Data(csc);
time = Range(csc);
csc_short = Restrict(tsdX, 5950,6050);

Timestamps_s = Timestamps.*10^-6;
tsdX = tsd(Timestamps_s,X');
tsdY = tsd(Timestamps_s,Y');

X_short = Restrict(tsdX,5950,6050);


%% Ploting with handles
fh_cb = @figure_move; % Create function handle for fig_move function
figure('KeyPressfcn',fh_cb)
plot(time,data)
starting_position = get(gca,'XLim');
set(gcf,'Color',[0 0 0]);  % makes the background black
set(gca,'Color',[0 0 0]);  % makes the current axis black
set(gca,'XColor', [1 1 1]); %makes the X axis while
set(gca,'YColor', [1 1 1]); %makes the Y axis while
set(gcf,'InvertHardCopy','off');  % Ensures that the colour is not automatically inverted for printing on white paper.  
% print(gcf,'-dpng','-r300','R042-2013-08-18-LFPsnippet.png');
hold on; box off;
 
csc_mean = nanmean(Data(csc));
xr = get(gca,'XLim');
 
mean_hdl = plot(xr,[csc_mean csc_mean]);
set(mean_hdl,'LineWidth',[2],'Color','r')
set(gca,'XLim',[5989 5990],'FontSize',[24])
set(gcf,'Name', 'mean_hdl')

%% Figure callback functions
% sqr_fn = @(x) x.^2;
% sqr_fn(2)

fh_cb = @fig_move; % Create function handle for fig_move function



