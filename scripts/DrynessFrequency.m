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
B = struct2cell(folderInfo);
for i = 3:folderLength
   stationNames(1,(i-2)) =  string(B(1,i));
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
dryFreq2 = (1:150)';
for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path    
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    %creates an array from the starting year to the ending year of the stations available weather data
    YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
    temporaryDIndex = table(YEAR); %creates an column array for the years of the H-Indecies 
    temporaryDIndex.DIndex = zeros(height(temporaryDIndex),1);
    %counter = 0;
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
            dryFreq2(j,i+1) = tempDryFreq.GreaterThanDP(j);
        end
    end 
end
%%
scatter(tempDryFreq.MAGNITUDE,tempDryFreq.ExactDP, 'filled')
hold on
scatter(tempDryFreq.MAGNITUDE, tempDryFreq.GreaterThanDP, 'filled')
hold off
bar(tempDryFreq.GreaterThanDP)
hold on
bar(tempDryFreq.ExactDP)

%timeDIndex = toc;

%% 3
figure('Name', 'D-Index')
for i = 1:stationLength %for each station in station names
    %this section creates an array of subplots where each station has it
    A = array2table(dryFreq2(dryFreq2(:,1+i) ~= 0,[ 1 1+i])); %for each station, create a new array with the non-zero values and their corresponding years
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