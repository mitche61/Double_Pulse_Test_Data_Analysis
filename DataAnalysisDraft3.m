%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Analyzing data from double pulse test csv file                         %
%                                                                        %
% To update for a different user:                                        %
% Change finalSpreadSheetFilePath, basePath, resultsPath                 %
% To analyze a new file:                                                 %
% Assign first set of variables to determine which ranges to analyze     %
% to write the resultant data to in the final spreadsheet                %
%                                                                        %
% Michaela Mitchell                                                      % 
% 7/21/2020                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

makeFolderForPlots = 1;
allData = 1;
range1 = 1; %Do range1 and other seperate bc zoom in is different
range2 = 0;
range3 = 0;
range4 = 0;
rowNum = 16;

finalSpreadSheetFilePath = 'C:\Users\mcmit\Desktop\Resonant_Test_Fixture_Data3.xlsx';

%% Read in and Parse Data
filename = "tek0015ALL.csv";
basePath = "G:\.shortcut-targets-by-id\1ud_ZivbIGgY-81vBF_LfsHRn2LZAyKzl\Michaela Mitchell\OIF  Project\Resonant Text Fixture Data\2020-06-16\";
resultsPath = strcat("G:\.shortcut-targets-by-id\1ud_ZivbIGgY-81vBF_LfsHRn2LZAyKzl\Michaela Mitchell\OIF  Project\Resonant Text Fixture Data\Matlab_Plots\",filename);
filePath = strcat(basePath, filename);
data = readtable(filePath, 'ReadVariableNames', false);
temp = table2cell(data);
indexOfNumericData = find(strcmp(temp, 'TIME'), 1);
indexOfSamplePeriod = find(strcmp(temp, 'Sample Interval'), 1);
useableData = data((indexOfNumericData+1):end, 1:4);
% Divide data into time and the three channels
time = str2double(table2array(useableData(:,1)));
Vgs = str2double(table2array(useableData(:,2))); % Gate to source voltage (not important right now)
Vds = str2double(table2array(useableData(:,3))); % Vds
current = str2double(table2array(useableData(:,4))); % Current
% Determine/assign variables
sampleInterval = str2double(table2array(data((indexOfSamplePeriod),2))); % From header file
VDC = 400; %DC Voltage Bus, 400V
Vmin = 10; %Minimum voltage goes to roughly 10V


%% Make new folder for resultant plots and write filename to excel file
if (makeFolderForPlots == 1)
    mkdir(resultsPath)
    rangeWrite = strcat('A',num2str(rowNum),':A',num2str(rowNum));
    xlswrite(finalSpreadSheetFilePath, filename, rangeWrite)
end

%% Plot full data
if (allData == 1)
plotTitle = strcat(filename, ' Current and Vds Vs. Time');
myPlot(plotTitle, time, 'Time (s)', Vds, 'Voltage (V)', current, 'Current (A)', resultsPath);
makeNotesOnWeirdData = input('Please take down any notes regarding how this data is wonky. Press enter when ready to continue.');
end

%% Moving average and new plot for dv/dt to determine switching instances
VdsAveraged = movmean(Vds, 5);
timeAveraged = movmean(time, 5);
dvdt = diff(VdsAveraged)./diff(timeAveraged);
threshold = 0.7*max(dvdt);
j = 1;
for i=1:length(dvdt)
    if (abs(dvdt(i)) >= threshold)
        dvdt(i) = dvdt(i);
        relevantIndicies(j) = i;
        j = j+1;
    else
        dvdt(i) = 0;
    end
end
differenceBetweenRelevantIndicies = diff(relevantIndicies);
j = 1;
for i=1:length(differenceBetweenRelevantIndicies)
    if (abs(differenceBetweenRelevantIndicies(i)) >= 200)
        extraRelevantIndicies(j) = relevantIndicies(i);
        j = j+1;
    end
end

% Sanity check
relevantTimes = zeros(1, length(extraRelevantIndicies));
for i = 1:length(relevantTimes)
    relevantTimes(i) = time(extraRelevantIndicies(i));
end
if (length(relevantTimes) ~= 4)
    relevantTimes(4) = 2.605e-5;
end

figure
yyaxis left
plot(timeAveraged, VdsAveraged)
for i=1:length(relevantTimes)
    xline(relevantTimes(i))
end
sanityCheck1 = input('Press enter if the vertical lines align with the switching instances.');

indiciesPerSeconds = length(time)/(max(time) - min(time));
if (range1 == 1) % Zoom in alot of first turn on
rangeBeyondSpikeIndex = (1.5e-7*indiciesPerSeconds);
else % Zoom in a little for other instances
rangeBeyondSpikeIndex = (4e-7*indiciesPerSeconds);
end

if (range1 == 1)
%% Range 1
% Zoom in on data
section1TurnOnStartTime = relevantTimes(1);
section1TurnOnStartIndex = find(abs(time - (section1TurnOnStartTime))==min(abs(time - (section1TurnOnStartTime))))-rangeBeyondSpikeIndex;
section1TurnOnEndTime = relevantTimes(1);
section1TurnOnEndIndex = find(abs(time - (section1TurnOnEndTime))==min(abs(time - (section1TurnOnEndTime))))+rangeBeyondSpikeIndex;

% Plot zoomed in data 
plotTitle = strcat(filename, ' Section 1 Turn On');
myPlot(plotTitle, time(section1TurnOnStartIndex:section1TurnOnEndIndex), 'Time (s)', Vds(section1TurnOnStartIndex:section1TurnOnEndIndex), 'Voltage (V)', current(section1TurnOnStartIndex:section1TurnOnEndIndex), 'Current (A)');
sanityCheckSection1 = input('Press enter if the plot shows the relevant range for the 1st turn on event.');

% Call function, sanity check bc sometimes it zero pads the output
section1TurnOn = turnOnFunction(section1TurnOnStartIndex, section1TurnOnEndIndex, time, Vds, VDC, Vmin, current, sampleInterval);
if (length(section1TurnOn) ~= 5)
    disp('I got confused, check out this data and fix it manually in the spreadsheet');
end
% Write results to excel file
rangeWrite = strcat('B', num2str(rowNum), ':F', num2str(rowNum));
xlswrite(finalSpreadSheetFilePath, section1TurnOn, rangeWrite)
end

if (range2 == 1)
%% Range 2 Turn Off 1
section2TurnOff1StartTime = relevantTimes(2);
section2TurnOff1StartIndex = find(abs(time - (section2TurnOff1StartTime))==min(abs(time - (section2TurnOff1StartTime))))-rangeBeyondSpikeIndex;
section2TurnOff1EndTime = relevantTimes(2);
section2TurnOff1EndIndex = find(abs(time - (section2TurnOff1EndTime))==min(abs(time - (section2TurnOff1EndTime))))+rangeBeyondSpikeIndex;

plotTitle = strcat(filename, ' Section 2 Turn Off 1');
myPlot(plotTitle, time(section2TurnOff1StartIndex:section2TurnOff1EndIndex), 'Time (s)', Vds(section2TurnOff1StartIndex:section2TurnOff1EndIndex), 'Voltage (V)', current(section2TurnOff1StartIndex:section2TurnOff1EndIndex), 'Current (A)');
sanityCheckSection2TurnOff1 = input('Press enter if the plot shows the relevant range for 1st turn off event.');

section2TurnOff1 = DVDT_time_DIDT_time_EFunc(section2TurnOff1StartIndex, section2TurnOff1EndIndex, time, Vds, VDC, Vmin, current, sampleInterval);
rangeWrite = strcat('G', num2str(rowNum), ':L', num2str(rowNum));
xlswrite(finalSpreadSheetFilePath, section2TurnOff1, rangeWrite)
end

if (range3 == 1)
%% Range 3 Turn On 2
section2TurnOnStartTime = relevantTimes(3);
section2TurnOnStartIndex = find(abs(time - (section2TurnOnStartTime))==min(abs(time - (section2TurnOnStartTime))))-rangeBeyondSpikeIndex;
section2TurnOnEndTime = relevantTimes(3);
section2TurnOnEndIndex = find(abs(time - (section2TurnOnEndTime))==min(abs(time - (section2TurnOnEndTime))))+rangeBeyondSpikeIndex;

plotTitle = strcat(filename, ' Section 2 Turn On');
myPlot(plotTitle, time(section2TurnOnStartIndex:section2TurnOnEndIndex), 'Time (s)', Vds(section2TurnOnStartIndex:section2TurnOnEndIndex), 'Voltage (V)', current(section2TurnOnStartIndex:section2TurnOnEndIndex), 'Current (A)');
sanityCheckSection2TurnOn = input('Press enter if the plot shows the relevant range for 2nd turn on event.');

section2TurnOn = DVDT_time_DIDT_time_EFunc(section2TurnOnStartIndex, section2TurnOnEndIndex, time, Vds, VDC, Vmin, current, sampleInterval);
rangeWrite = strcat('M', num2str(rowNum), ':R', num2str(rowNum));
xlswrite(finalSpreadSheetFilePath, section2TurnOn, rangeWrite)
end

if (range4 == 1)
%% Range 4 Turn Off 2
section2TurnOff2StartTime = relevantTimes(4);
section2TurnOff2StartIndex = find(abs(time - (section2TurnOff2StartTime))==min(abs(time - (section2TurnOff2StartTime))))-rangeBeyondSpikeIndex;
section2TurnOff2EndTime = relevantTimes(4);
section2TurnOff2EndIndex = find(abs(time - (section2TurnOff2EndTime))==min(abs(time - (section2TurnOff2EndTime))))+rangeBeyondSpikeIndex;

plotTitle = strcat(filename, ' Section 2 Turn Off 2');
myPlot(plotTitle, time(section2TurnOff2StartIndex:section2TurnOff2EndIndex), 'Time (s)', Vds(section2TurnOff2StartIndex:section2TurnOff2EndIndex), 'Voltage (V)', current(section2TurnOff2StartIndex:section2TurnOff2EndIndex), 'Current (A)');
sanityCheckSection2TurnOn = input('Press enter if the plot shows the relevant range for 2nd turn off event.');

section2TurnOff2 = DVDT_time_DIDT_time_EFunc(section2TurnOff2StartIndex, section2TurnOff2EndIndex, time, Vds, VDC, Vmin, current, sampleInterval);
rangeWrite = strcat('S', num2str(rowNum), ':X', num2str(rowNum));
xlswrite(finalSpreadSheetFilePath, section2TurnOff2, rangeWrite)
end

%% Function stuff
function DVDT_time_DIDT_time_E = DVDT_time_DIDT_time_EFunc(startIndex, endIndex, time, Vds, VDC, Vmin, current, sampleInterval)
Vds = Vds(startIndex:endIndex);
time = time(startIndex:endIndex);
current = current(startIndex:endIndex);

entireDV = VDC - Vmin; %Evened out max - evened out min

ninetyPercentVds = 0.9*(entireDV);
ninetyPercentVdsIndex = find(abs(Vds - ninetyPercentVds)==min(abs(Vds - ninetyPercentVds)));
tenPercentVds = 0.1*(entireDV);
tenPercentVdsIndex = find(abs(Vds - tenPercentVds)==min(abs(Vds - tenPercentVds)));
DVDTIndicies = [ninetyPercentVdsIndex(1), tenPercentVdsIndex(1)];

timeValueDVDT = (time(max(DVDTIndicies))-time(min(DVDTIndicies)));
DVDT = (Vds(max(DVDTIndicies))-Vds(min(DVDTIndicies)))/timeValueDVDT; 

ninetyPercentCurrent = 0.9*(max(current));
ninetyPercentCurrentIndex = find(abs(current - ninetyPercentCurrent)==min(abs(current - ninetyPercentCurrent)));
tenPercentCurrent = 0.1*(max(current));
tenPercentCurrentIndex = find(abs(current - tenPercentCurrent)==min(abs(current - tenPercentCurrent)));
DIDTIndicies = [ninetyPercentCurrentIndex(1), tenPercentCurrentIndex(1)];

timeValueDIDT = (time(max(DIDTIndicies))-time(min(DIDTIndicies)));
DIDT = (current(max(DIDTIndicies))-current(min(DIDTIndicies)))/timeValueDIDT; 

allRelevantIndicies = [ninetyPercentVdsIndex(1), tenPercentVdsIndex(1), ninetyPercentCurrentIndex(1), tenPercentCurrentIndex(1)];

intgStartIndex = min(allRelevantIndicies);
%intgStartIndex = find(time == 2.601e-5);
intgEndIndex = max(allRelevantIndicies);
hold on
xline(time(intgStartIndex), 'LineWidth', 1)
xline(time(intgEndIndex), 'LineWidth', 1)
sanityCheckIntegration = input('Press enter if youre good to integrate between these bounds for E.');

approximation = max(Vds)*max(current)*(time(intgEndIndex)-time(intgStartIndex))/2;

VdsTimesCurrent = Vds(intgStartIndex:intgEndIndex).*current(intgStartIndex:intgEndIndex);
E = sum(VdsTimesCurrent)*sampleInterval;

DVDT_time_DIDT_time_E = [DVDT, timeValueDVDT, DIDT, timeValueDIDT, approximation, E];
end

%% Everything for section 1
function section1TurnOn = turnOnFunction(startIndex, endIndex, time, Vds, VDC, Vmin, current, sampleInterval)
entireDV = VDC - Vmin; %Evened out max - evened out min

Vds = Vds(startIndex:endIndex); % Parse
time = time(startIndex:endIndex); % Parse
current = current(startIndex:endIndex);

ninetyPercentVds = 0.9*(entireDV);
ninetyPercentVdsIndex = find(abs(Vds - ninetyPercentVds)==min(abs(Vds - ninetyPercentVds)));
tenPercentVds = 0.1*(max(Vds));
tenPercentVdsIndex = find(abs(Vds - tenPercentVds)==min(abs(Vds - tenPercentVds)));

timeValueVds = time(tenPercentVdsIndex)-time(ninetyPercentVdsIndex);
DVDT = (tenPercentVds-ninetyPercentVds)/timeValueVds; 

[Idpk, IdpkIndex] = max(current);
intgStartIndex = find(time == 5.6e-9);

% Find integration end index
currentAfterPeak = current(IdpkIndex:(IdpkIndex+50));
[~, intgEndIndex] = min(currentAfterPeak);
intgEndIndex = intgEndIndex+IdpkIndex;
hold on
xline(time(intgStartIndex(1)), 'LineWidth', 1)
xline(time(intgEndIndex(1)), 'LineWidth', 1)
sanityCheckSection1Integration = input('Press enter if youre good to integrate between these bounds for Q and Eon.');

Q = abs(trapz(current(intgStartIndex:intgEndIndex), time(intgStartIndex:intgEndIndex)));

approximation = max(Vds)*max(current)*(time(intgEndIndex)-time(intgStartIndex))/2;

VdsTimesCurrent = Vds(intgStartIndex:intgEndIndex).*current(intgStartIndex:intgEndIndex);
Eon = sum(VdsTimesCurrent)*sampleInterval;

section1TurnOn = [DVDT, Idpk, Q, approximation, Eon];
end

%% Plotting function
function myPlot(plotTitle, xaxis, xAxisLabel, yaxisleft, ylabelleft, yaxisright, ylabelright, resultsPath)
    fig = figure;
    left_color = [0 0 0];
    right_color = [0 0 1];
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);
    
    % Figure Properties:
        AxisFontSize        = 8;
        ImageSize           = [0 0 3.2 3]; % Width x Height
        PlotLineWidth       = 4;
        BorderGridLineWidth = 1.3;
    % PLOT APPEARANCE
       set(gca,'fontsize', AxisFontSize, ...
                'fontweight', 'bold',...
                'FontName','Times',...
                'LineWidth',BorderGridLineWidth,...
                'XGrid','on', ...
                'YGrid','on');
       set(gcf,'PaperUnits', 'inches',...
                'PaperPosition', ImageSize);
    grid off
    yyaxis left
    plot(xaxis, yaxisleft, 'k')
    ylabel(ylabelleft, 'fontsize',8, 'fontweight', 'bold');
    yyaxis right
    plot(xaxis, yaxisright, 'b')
    ylabel(ylabelright, 'fontsize',8, 'fontweight', 'bold');
    legend (ylabelleft, ylabelright) 
    xlabel(xAxisLabel, 'fontsize',8, 'fontweight', 'bold');
    title(plotTitle, 'fontsize', 9, 'fontweight', 'bold');
    
    %saveas(gcf, plotTitle, 'jpg')
    %movefile(strcat(plotTitle, '.jpg'), resultsPath, 'f');

end


