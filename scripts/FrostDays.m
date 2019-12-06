%
%This code is designed to take files from a specific folder, if they are
%excel files they are converted into csv files. those csv files are then
%used to produce an H-index based on temperature for each station.

%Note: you will need to download the "Weather folder" and add it to your
%directory file to ensure this program works. 
clear all
clc


%% 1
folderName = 'Weather_CSV'; %variable for easy change of folder name
folderInfo = dir(folderName);  %creates a structure array with all the file names in "folderName"
folderLength = length(folderInfo); 
B = struct2cell(folderInfo);
for i = 3:folderLength
   stationNames(1,(i-2)) =  string(B(1,i));
   tableStationNames(1,(i-2))= erase(stationNames(1,(i-2)),".csv");
end
newFolder = strcat(folderName);

%% 2
clc
tic
stationLength = length(stationNames);
folder = strcat(pwd,'/',newFolder); %calls the path of the current file directory

FrostIndex = zeros(124,23);
FrostIndex(:,1) = (1890:2013)';

for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path    
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    %creates an array from the starting year to the ending year of the stations available weather data
    YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
    temporaryFrostIndex = table(YEAR); %creates an column array for the years Frost days 
    temporaryFrostIndex.FrostIndex = zeros(height(temporaryFrostIndex),1);
    counter = 0;
    for j = temporaryFrostIndex.YEAR(1):temporaryFrostIndex.YEAR(end)%for each year at this station
        year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary matrix for the given year       
        frostTemp = 32;
        %Count the number of times where the daily temp is less than or
        %equal to the frostTemp                
            counter = 0;
            for h = 1:height(year)% for days in this year
                if year.TMIN(h) <= frostTemp %If the value at row h and column TMIN Column are less than frostTemp
                    counter = counter + 1; %increase counter by 1
                end    
            end                    
        temporaryFrostIndex.FrostIndex(j-temporaryFrostIndex.YEAR(1)+1) = counter; %stores the number of Frost Days value for each year in temporaryFrostIndex      
    end
    
    for j = 1:height(temporaryFrostIndex) %for the number of years at the current station
        for h = 1:length(FrostIndex) %for full array of years being analyzed
           if temporaryFrostIndex.YEAR(j) == FrostIndex(h,1) %Checks to make sure that the years are the same for the given station
              FrostIndex(h,i+1)=temporaryFrostIndex.FrostIndex(j); %if the years are the same, then the yearly value for the station is stored in a column specifically for that station
           end
        end
    end
end
Slopes = table;
Slopes.NAME = tableStationNames';
figure('Name', 'C-Index')
for i = 1:stationLength %for each station in station names
    %this section creates an array of subplots where each station has its
    %own plot
    tempNames = split(tableStationNames(i), '_');
    A = FrostIndex(FrostIndex(:,1+i) ~= 0,[ 1 1+i]); %for each station, create a new array with the non-zero values and their corresponding years
    x = A(:,1);%x-axis is the year
    y = A(:,2);%y-axis is the Frost Days
    subplot(4,6,i)
    plot(x,y)
    hold on
    scatter(x, y,15,'filled')  
    s = lsline;
    s.Color = 'k';
    xlabel('Year')
    xlim([1890 2014])
    ylabel('Number of Frost Days')
    ylim([50 200])
    %this section  creates a table of the station names and the slope for
    %each trend line
    p = polyfit(x,y,1);
    Slopes.FrostIndexSlope(i) = round(p(1),3);
    Slopes.FrostIndexAverage(i) = round(mean(A(:,2)),1);
    title(compose(tempNames(2,1)+"\n"+Slopes.FrostIndexAverage(i)+"\n"+Slopes.FrostIndexSlope(i)));
end
timeFrostIndex = toc;