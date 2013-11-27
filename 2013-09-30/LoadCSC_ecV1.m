function csc = LoadCSC_ecV1(fname,varargin)
% function csc = LoadCSC(fname)
%
% load a Neuralynx .ncs file correctly
%
% INPUTS:
% fname: [1 x n] char, i.e. a filename such as 'R042-2013-08-18-CSC01a.ncs'
% time_units: [1 x n] this is the units that you wish the time vector to be
% converted to.  (defualt: 100000, seconds)
% tUnits: [1 x n] string.  This is the string for the time units ex:
% time_units = 100000, then tUnits = 'sec'
% voltage-units: [1 x n] This is the unit for the voltage (default: 1000,
% milivolts)
%
% OUTPUTS:
%
% csc: [1 x 1] mytsd
%% Define some variables.
run(FindFile('*keys.m'))
% type = 'all';  % this is the type of data you wish to display.
voltage_units = 1000; %mV
time_units = 1000000; % conversion to seconds.
tUnits = 'sec';
extract_varargin; % allows the user to set these variables

if time_units == 1000000
    tUnits = 'sec';
elseif time_units == 10000
    tunits = 'centi_sec';
elseif time_units == 1000
    tunits = 'deci_sec';
else
    tUnits = 'unspecified';  % This should actually cover more units but that would take more time and space, but I am sure that you get the idea.
end
%% Remind the user that the file needs to be in the path.
c_dir = cd;

if strcmp(c_dir((length(c_dir)-14):end),fname(1:15))==0
    error('You are not in the right dirctory.  Please add the data folder to the path or cd to the data folder')
end



%%  Load the data using the old function, and convert the sample voltages and time units.
[Timestamps, ~, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(fname, [1 1 1 1 1], 1, 1, []);

BitsperVolt = str2double(Header{15,1}(13:end));
Samples_in_mV = Samples.*voltage_units.*BitsperVolt; % converts the samples to mV
Samples_in_mV = NaN*ones(size(Samples));
for ii = 1:size(Samples,1)
    Samples_in_mV(ii,:) = Samples(ii,:).*voltage_units.*BitsperVolt;
end
Timestamps_s = Timestamps./time_units; % converts the time to seconds

%% Identify gaps in the signal using the NumberOfValidSamples variable from the data
% find the gaps using the NumberOfValidSamples
invalid_data_ids = (NumberOfValidSamples<512);
invalid_data_ids_plot(invalid_data_ids<1) = NaN;
invalid_data_ids_plot(invalid_data_ids==1) = 1.1;

% %sanity check plot with the ids
% figure; hold on;
% plot(1:length(NumberOfValidSamples), NumberOfValidSamples==512, 'b')
% plot(1:length(NumberOfValidSamples),invalid_data_ids_plot,'r','Marker', '*')
% ylim([-0.5 2])
% xlim([0 length(NumberOfValidSamples)])

%% Reshape the samples and times into a vector
Samples_in_mV_1D =  (reshape(Samples_in_mV,[size(Samples_in_mV,1)*size(Samples_in_mV,2) 1]))';
% New Time Series
Timestamps_s_1D = NaN(size(Samples_in_mV,1),size(Samples_in_mV,2));
for ii =1:length(Timestamps_s)
    Timestamps_s_1D(:,ii) = Timestamps_s(ii);
end
Timestamps_s_1D = (reshape(Timestamps_s_1D,[size(Timestamps_s_1D,1)*size(Timestamps_s_1D,2) 1]))';


if length(Timestamps_s_1D) ~= length(Samples_in_mV_1D)
    error('The Timestamps vector and the Samples vector are not the same length')
end

%% Index, remove and verify the seconds with invalid data points
    index_invalid = (Samples_in_mV_1D == 0);
    Samples_in_mV_1D(index_invalid) = [];
    Timestamps_s_1D(index_invalid) = [];
    
    if length(Timestamps_s_1D) ~= length(Samples_in_mV_1D)
        error('The updated Timestamps vector and the updated Samples vector are not the same length')
    end
%% convert everything in to a tsd using the mytsd function with the built in header component.

csc = mytsd(Timestamps_s_1D, Samples_in_mV_1D, tUnits, Header);

%plot for sanity
inv_inds = sum(NumberOfValidSamples<512);
total_invalid = sum(NumberOfValidSamples(inv_inds));
total_samples = length(Timestamps_s_1D);
percentage = total_invalid/total_samples;
figure
subplot(211)
plot(csc,'b')
subplot(212)
plot(diff(Range(csc)),'r')
xlabel(['Invalid points:  ' num2str(percentage) '%'])
text(length(Range(csc))-length(Range(csc))/4, 0.5*(max(Range(csc))), ['Total Invalid Points:' total_invalid])

end