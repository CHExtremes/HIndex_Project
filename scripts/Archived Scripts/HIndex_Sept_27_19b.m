%%
%This code is designed to take files from a specific folder, if they are
%excel files they are converted into csv files. those csv files are then
%used to produce an H-index based on temperature for each station.

%Note: you will need to download the "Weather folder" and add it to your
%directory file to ensure this program works. 
clear all
clc


%% 5
% simplifeid version of blocks 3&4 that also adds an output array that
% contains the names of the stations
clc
folderName = 'Weather'; %variable for easy change of folder name
folderInfo = dir(folderName);  %creates a structure array with all the file names in "folderName"
folderLength = length(folderInfo); 
stationNames = strings([1,(folderLength-2)]); %creates an open array for the station names that are in folderName folder
H = {}; 
TMaxColumn = 4;
TMinColumn = 5;
newFolder = strcat(folderName,'_CSV'); %creates a variable for the folder name that will be used to store the new CSV files.
mkdir(newFolder); %creates the the new folder to store the CSV files. 

for i = 2:folderLength
    H(1,i) = cellstr(folderInfo(i).name); %converts the cells in "folderInfo" into strings
    T1 = endsWith(H(1,i),'.xlsx'); %checks if the file is an xlsx file    
    if (T1 == 1) && (startsWith(H(1,i),'~$') == 0)
    [status,sheets] = xlsfinfo(strcat(folderName,'/',char(H(1,i))));
        for j =1:length(sheets)
            num = xlsread(strcat(folderName,'/',char(H(1,i))),string(sheets(j))); %reads each excel file in folderName and records it as two matricies, one for strings, the other for numbers
            text = {'MONTH','DAY','YEAR','TMAX','TMIN','RAIN','SNOW'};
            file = [text; num2cell(num)];
            writecell(file,strcat(newFolder,'/',string(sheets(j)),'.csv'));%writes the num matrix as a csv file with a name drawn from the title position of the text matrix text matrix
            stationNames(1,j) = strcat(string(sheets(j)),'.csv'); %adds the current station name to the stationNames variable
        end
    end
end
%% 6
folderName = 'Weather_CSV'; %variable for easy change of folder name
folderInfo = dir(folderName);  %creates a structure array with all the file names in "folderName"
folderLength = length(folderInfo); 
B = struct2cell(folderInfo);
for i = 3:folderLength
   stationNames(1,(i-2)) =  string(B(1,i));
   tableStationNames(1,(i-2))= erase(stationNames(1,(i-2)),".csv");
end
newFolder = strcat(folderName);
%% 7
%this will be the first steps of creating the C/H index. it will start by
%anazyling a single station on a yearly basis and produce a bar graph at
%the end
clc

stationLength = length(stationNames);
%stationLength = 1;
folder = strcat(pwd,'/',newFolder); %calls the path of the current file directory
HIndex = zeros(130,23);
HIndex(:,1) = 1890:2019';

for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    %creates an array from the starting year to the ending year of the stations available weather data
    YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
    temporaryHIndex = table(YEAR); %creates an column array for the years of the H-Indecies 
    temporaryHIndex.HIndex = zeros(height(temporaryHIndex),1);
    for j = temporaryHIndex.YEAR(1):temporaryHIndex.YEAR(end)%for each year at this station
        year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary table for the given year        
        currentTemp = round(max(year.TMAX)); %records max temp for the given year         
        %Count the number of times where the daily temp is greater than or
        %equal to that temp
        counter = 0;
        while counter < currentTemp %checks to see if the counter is smaller than the currentTemp. This is to make sure that the value is an H-Index value.
            counter = 0;
            for h = 1:height(year)% for days in this year                
                if year.TMAX(h) >= currentTemp %If the value at row h and column TMaxColumn are greater than currentTemp
                    counter = counter + 1; %increase counter by 1
                end    
            end
            if counter < currentTemp %if the counter is smaller than currentTemp then the H-index is not valid, so we reduce it by one and repeat the loop.
                currentTemp = currentTemp - 1;
            end
        end
        %temporaryHIndex((j-temporaryHIndex(1,1)+1), 2) = currentTemp; %stores the max H-index value for each year in temporaryHIndex      
        temporaryHIndex.HIndex(j-temporaryHIndex.YEAR(1)+1) = currentTemp; 
    end   
    for j = 1:height(temporaryHIndex) %for the number of years at the current station
        for h = 1:length(HIndex) %for full array of years being analyzed
           if temporaryHIndex.YEAR(j) == HIndex(h,1) %Checks to make sure that the years are the same for the given station
              HIndex(h,(i+1))=temporaryHIndex.HIndex(j); %if the years are the same, then the yearly value for the station is stored in a column specifically for that station
           end
        end
    end
end

for i = 1:length(tableStationNames) %for each station in station names
    A = HIndex(HIndex(:,1+i) ~= 0,[ 1 1+i]); %for each station, create a new array with the non-zero values and their corresponding years
    x = A(:,1);%x-axis is the year
    y = A(:,2);%y-axis is the HIndex
    hold on %add each station to the same plot
    plot(x, y)  
    legend(tableStationNames(1,1:i),'Location','South'); %include a legend for the station names
end
hold off
%barh(x,y)

%% 8 
%create a C-Index
clc

%stationLength = 1;
 %calls the path of the current file directory
CIndex = (1890:2019)';
for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path    
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
    %creates an array from the starting year to the ending year of the stations available weather data
    temporaryCIndex =  table(YEAR); %creates an column array for the years of the H-Indecies 
    temporaryCIndex.CIndex = zeros(height(temporaryCIndex),1);
    counter = 0;
    for j = temporaryCIndex.YEAR(1):temporaryCIndex.YEAR(end)%for each year at this station
        year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary matrix for the given year       
        currentTemp = round(min(year.TMIN)); %records max temp for the given year        
        %Count the number of times where the daily temp is greater than or
        %equal to that temp
        while counter < 32-currentTemp %checks to see if the counter is smaller than the currentTemp. This is to make sure that the value is an H-Index value.
            counter = 0;
            for h = 1:length(B)% for days in this year
                if year.TMIN(h) <= currentTemp %If the value at row h and column TMaxColumn are greater than currentTemp
                    counter = counter + 1; %increase counter by 1
                end    
            end
            if counter < 32-currentTemp %if the counter is smaller than currentTemp then the H-index is not valid, so we reduce it by one and repeat the loop.
                currentTemp = currentTemp + 1;
            end
        end
        temporaryCIndex.CIndex(j-temporaryCIndex.YEAR(1)+1) = currentTemp; %stores the max H-index value for each year in temporaryHIndex      
    end
    
    for j = 1:height(temporaryCIndex) %for the number of years at the current station
        for h = 1:length(CIndex) %for full array of years being analyzed
           if temporaryCIndex.YEAR(j) == CIndex(h,1) %Checks to make sure that the years are the same for the given station
              CIndex(h,i+1)=temporaryCIndex.CIndex(j); %if the years are the same, then the yearly value for the station is stored in a column specifically for that station
           end
        end
    end
end
for i = 1:stationLength %for each station in station names
    A = CIndex(CIndex(:,1+i) ~= 0,[ 1 1+i]); %for each station, create a new array with the non-zero values and their corresponding years
    x = A(:,1);%x-axis is the year
    y = A(:,2);%y-axis is the HIndex
    hold on %add each station to the same plot
    plot(x, y)  
    legend(tableStationNames(1,1:i),'Location','South'); %include a legend for the station names
end
