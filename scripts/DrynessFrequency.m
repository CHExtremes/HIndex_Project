%%
%This code is designed to take files from a specific folder, if they are
%excel files they are converted into csv files. those csv files are then
%used to produce an H-index based on temperature for each station.

%Note: you will need to download the "Weather folder" and add it to your
%directory file to ensure this program works. 
clear 
clc


%% 1
%This script will look at the entire data set for the station and find the
%longest dry period. It will then create and array from 1 to that dry
%period and count the number of times for the whole data set that each dry
%period length occured.
folderName = 'Weather_CSV'; %variable for easy change of folder name
folderInfo = dir(folderName);  %creates a structure array with all the file names in "folderName"
folderLength = length(folderInfo); 
A = struct2cell(folderInfo);
for i = 3:folderLength
   stationNames(1,(i-2)) =  string(A(1,i));
   tableStationNames(1,(i-2))= erase(stationNames(1,(i-2)),".csv");
end
newFolder = strcat(folderName);
clc
%% 2
%this script calculates the number of times each frequency of dry days
%occurs for the entire data set. 
tic
%stationLength = 1;
%calls the path of the current file directory
stationLength = length(stationNames);
folder = strcat(pwd,'/',newFolder); %calls the path of the current file directory
%DIndex = zeros(124,23);
%DIndex(:,1) = (1890:2013)';

%dryFreq1 = (1:150)';
dryFreq = (1:150)';
for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path    
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name 
    dryness = 1; %dryness is the minimum precipitation threshold for a day to be considered 'wet'
    currentDryPeriod = 0;
    tempDryPeriod = 0;
    tempDryFreq = table();
    for j = 1:length(temporaryFile.RAIN)%for each year at this station
        if temporaryFile.RAIN(j) < dryness
            tempDryPeriod = tempDryPeriod + 1;
        else
            if currentDryPeriod < tempDryPeriod
                currentDryPeriod = tempDryPeriod;
            end
            tempDryPeriod = 0;
        end
    end
    maxDP = currentDryPeriod;
    tempDryFreq.MAGNITUDE = transpose(1:maxDP); %creates a column showing the range of magnitudes being calculated
    for h = 1:maxDP %this first loop is the range of dry periods that are being analyzed for their frequency
        exactDP = 0; %variable for the number of times the dry period is exactly h days
        greaterThanDP = 0; %variable for the number of times the dry period is h days or greater
       for j = 1:length(temporaryFile.RAIN)%for the entire precipitation data set
           if temporaryFile.RAIN(j) < dryness%this if statement counts the length of the dry period
               tempDryPeriod = tempDryPeriod + 1;
           else %if the day is had precipitation greater than the dryness threshold, we then determine if the dry period(DP) was greater than or equal to h length
               if tempDryPeriod >= h 
                   greaterThanDP = greaterThanDP +1; 
                   if tempDryPeriod ==h
                        exactDP = exactDP +1;                      
                   end                  
               end
               tempDryPeriod = 0;
           end          
       end
       tempDryFreq.ExactDP(h) = exactDP;
       tempDryFreq.GreaterThanDP(h) = greaterThanDP;
    end
    %this part of the loop will stow the dryness frequency for each station
    for j = 1:length(tempDryFreq.MAGNITUDE)
        if dryFreq(j,1) == tempDryFreq.MAGNITUDE(j)
            %dryFreq1(j,i+1) = tempDryFreq.ExactDP(j);
            dryFreq(j,i+1) = tempDryFreq.GreaterThanDP(j);
        end
    end 
end
%% 3
scatter(tempDryFreq.MAGNITUDE,tempDryFreq.ExactDP, 'filled')
hold on
scatter(tempDryFreq.MAGNITUDE, tempDryFreq.GreaterThanDP, 'filled')
hold off
bar(tempDryFreq.GreaterThanDP)
hold on
bar(tempDryFreq.ExactDP)

%timeDIndex = toc;

%% 4
figure('Name', 'D-Index')
for i = 1:stationLength %for each station in station names
    %this section creates an array of subplots where each station has it
    A = array2table(dryFreq(dryFreq(:,1+i) ~= 0,[ 1 1+i])); %for each station, create a new array with the non-zero values and their corresponding years
    A.Properties.VariableNames{'Var1'} = 'Magnitude'; 
    A.Properties.VariableNames{'Var2'} = 'Frequency';
    subplot(4,6,i)
    l = scatter(A.Magnitude,A.Frequency,5,'filled'); %adds a line to the plot for additional clarity   
    
    hold on %add each station to the same plot
    xlim([-10 150]) 
    xlabel('Frequency')
    ylim([-5 6300])
    ylabel('Magnitude')
    tempNames = split(tableStationNames(i), '_');
    title(tempNames(2,1))
end
%% 5
%for this seciont, I need to update the script from 2 to be able to
%determine the number of times a dry period of a specific lengths occurs
%for each given year at each given station. 
tic
clc
stationLength = 1;%length(stationNames);
folder = strcat(pwd,'/',newFolder); %calls the path of the current file directory
%DIndex = (1890:2013)';

dryFreq = (1:150)';
frequency = 1;
dryness = 1; %dryness is the minimum precipitation threshold for a day to be considered 'wet'
for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path    
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    DIndex = (temporaryFile.YEAR(1):temporaryFile.YEAR(end))';
    %so I want this script to look at the amount of times a dry period of
    %each length listed in dryFreq occurs in a given year at a given
    %station. So the outer loop will, as always, be the station loop. Next
    %will be the frequency loop. Then I'll have a year loop. Then a
    %counting loop. Before the Station loop restarts, I'll need to save the
    %frequency information.
    for j = 1:150 %for the range of the maximum dry period, this script doesn't need to change based on maxDP because the loop only fills DIndex out to the maxDP
        for h = temporaryFile.YEAR(1):temporaryFile.YEAR(end) %for each year
            year = temporaryFile(temporaryFile.YEAR==h,:);
            counter = 0;
            tempFreq = 0;
            %right now this loop doesn't properly count the number of dry
            %days. It currently counts the number of dry periods properly,
            %but doesn't count when there are multiple smaller dry days in
            %a period.
            for k = 1:height(year) %For each year, this loop finds the number of times there is a dry period of length j                
                if year.RAIN(k) <= dryness %counts the number of times in a row that there are days with rain
                    counter = counter + 1;
                    if counter >= j
                        tempFreq = tempFreq+1;
                        counter = 0;
                    end                   
                else
                    if counter >= j
                        tempFreq = tempFreq+1;                        
                    end
                    counter = 0;
                end
            end
            for k = DIndex(1,1):DIndex(end,1)
                if DIndex(k-DIndex(1,1)+1,1) == h 
                    if tempFreq > 0 %DIndex(h-temporaryFile.YEAR(1)+1,j+1)
                        DIndex(k-DIndex(1,1)+1,j+1)= tempFreq; %ensures that the year for DIndex matches the Year for the station
                    end
                end
            end
        end     
    end
    m = find(DIndex);           
end

timeD = toc;

%%
clc
A = DIndex(1:(end-1),1);
[r,c] = size(DIndex);
Results = table;
for i = 11:21
    A = array2table(DIndex(DIndex(:,1+i) ~= 0,[ 1 1+i])); %for each station, create a new array with the non-zero values and their corresponding years
    A.Properties.VariableNames{'Var1'} = 'YEAR'; 
    A.Properties.VariableNames{'Var2'} = 'DIndex';
    subplot(2,6,i-10)
    plot(A.YEAR(1:(end-1)),A.DIndex(1:(end-1)))   
    hold on %add each station to the same plot
    mdl = fitlm(A, 'DIndex ~ YEAR'); %performs a linear regression for the Year and the D Index
    z = plot(mdl); %plots the linear regression and data points
    %below are changes to the colors and markers of the plot for additional
    %clarity I removed 95% error bars (z(3) and z(4) to make the graph less cluttered. We
    %can turn these on later if we want to visually analyze the error
    %margins.
    z(1).Color = '#D95319'; %sets data points to be orange
    z(1).Marker = '.';
    z(1).MarkerSize = 10;
    z(2).Color = 'k';
    z(2).LineWidth = 1;
    %z(3).Color = 'none';
    %z(4).Color = 'none';
    legend('off'); %hides the automatic legend generated by fitlm
    ylim([0 30])
    Results.DIndexSlope(i) = round(table2array(mdl.Coefficients(2,1)),3); %calls the slope given for the linear regression of the data using the fitlm function
    Results.DIndexAverage(i) = round(mean(A.DIndex),1); %calculates the average D index for the station and adds that to a new table
    Results.DIndexRsqr(i) = mdl.Rsquared.Ordinary; %calls the r^2 value from the fitlm function and inputs it into a new table
    Results.DIndexPValue(i) = round(table2array(mdl.Coefficients(2,4)),3); %calls the pValue from the fitlm function and inputs it into a new table
end