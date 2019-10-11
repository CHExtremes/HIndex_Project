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
    T = endsWith(H(1,i),'.xlsx'); %checks if the file is an xlsx file
    if (T == 1) && (startsWith(H(1,i),'~$') == 0)
        [num, text] = xlsread(strcat(folderName,'/',char(H(1,i)))); %reads each excel file in folderName and records it as two matricies, one for strings, the other for numbers
        for j = 1:length(num) %For each file read, changes the temperature units from C to F
            num(j, TMaxColumn) = (num(j, TMaxColumn)*9/5)+32; %convert tempuratures from C to F for TMax
            num(j, TMinColumn) = (num(j, TMinColumn)*9/5)+32; %convert tempuratures from C to F for TMin
        end
        writematrix(num,strcat(newFolder,'/',char(text(1,2)),'.csv'));%writes the num matrix as a csv file with a name drawn from the title position of the text matrix text matrix
        stationNames(1,(i-2)) = convertCharsToStrings(strcat(char(text(1,2)),'.csv')); %adds the current station name to the stationNames variable
    end
        
end
%% 6
%this will be the first steps of creating the C/H index. it will start by
%anazyling a single station on a yearly basis and produce a bar graph at
%the end
clc

stationLength = length(stationNames);
%stationLength = 1;
folder = strcat(pwd,'/',newFolder); %calls the path of the current file directory
HIndex = (1890:2019)';
for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path
    temporaryFile = readmatrix(fullFileName); %creates a temporary matrix of the the data for the current station name.
    yearColumn = 1;
    dayColumn = 3;
    TMaxColumn = 4;
    TMinColumn = 5;
    %creates an array from the starting year to the ending year of the stations available weather data
    temporaryHIndex = transpose(min(temporaryFile(:,yearColumn)):max(temporaryFile(:,yearColumn))); %creates an column array for the years of the H-Indecies 
    counter = 0;
    for j = temporaryHIndex(1,1):temporaryHIndex(end,1)%for each year at this station
        year = find(temporaryFile==j); %locates the index values for the given year
        B = temporaryFile(year,:); %creates a temporary matrix for the given year
        currentTemp = round(max(B(:,TMaxColumn))); %records max temp for the given year 
        
        %Count the number of times where the daily temp is greater than or
                %equal to that temp
        while counter < currentTemp %checks to see if the counter is smaller than the currentTemp. This is to make sure that the value is an H-Index value.
            counter = 0;
            for h = 1:length(B)% for days in this year
                
                if B(h,TMaxColumn) >= currentTemp %If the value at row h and column TMaxColumn are greater than currentTemp
                    counter = counter + 1; %increase counter by 1
                end    
            end
            if counter < currentTemp %if the counter is smaller than currentTemp then the H-index is not valid, so we reduce it by one and repeat the loop.
                currentTemp = currentTemp - 1;
            end
        end
        temporaryHIndex((j-temporaryHIndex(1,1)+1), 2) = currentTemp; %stores the max H-index value for each year in temporaryHIndex      
    end
    
    for j = 1:length(temporaryHIndex) %for the number of years at the current station
        for h = 1:length(HIndex) %for full array of years being analyzed
           if temporaryHIndex(j,1) == HIndex(h,1) %Checks to make sure that the years are the same for the given station
              HIndex(h,i+1)=temporaryHIndex(j,2); %if the years are the same, then the yearly value for the station is stored in a column specifically for that station
           end
        end
    end
end
x = HIndex(:,1);
y1 = HIndex(:,2);
y2 = HIndex(:,3);
plot(x,y1)
hold on
plot(x,y2)
%barh(x,y)

%% 7 
%create a C-Index
clc


%stationLength = 1;
 %calls the path of the current file directory
CIndex = (1890:2019)';
for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path
    temporaryFile = readmatrix(fullFileName); %creates a temporary matrix of the the data for the current station name.
    yearColumn = 1;
    dayColumn = 3;
    TMaxColumn = 4;
    TMinColumn = 5;
    %creates an array from the starting year to the ending year of the stations available weather data
    temporaryCIndex = transpose(min(temporaryFile(:,yearColumn)):max(temporaryFile(:,yearColumn))); %creates an column array for the years of the H-Indecies 
    counter = 0;
    for j = temporaryCIndex(1,1):temporaryCIndex(end,1)%for each year at this station
        year = find(temporaryFile==j); %locates the index values for the given year
        B = temporaryFile(year,:); %creates a temporary matrix for the given year
        currentTemp = round(min(B(:,TMinColumn))); %records max temp for the given year 
        
        %Count the number of times where the daily temp is greater than or
                %equal to that temp
        while counter < 32-currentTemp %checks to see if the counter is smaller than the currentTemp. This is to make sure that the value is an H-Index value.
            counter = 0;
            for h = 1:length(B)% for days in this year
                if B(h,TMinColumn) <= currentTemp %If the value at row h and column TMaxColumn are greater than currentTemp
                    counter = counter + 1; %increase counter by 1
                end    
            end
            if counter < 32-currentTemp %if the counter is smaller than currentTemp then the H-index is not valid, so we reduce it by one and repeat the loop.
                currentTemp = currentTemp + 1;
            end
        end
        temporaryCIndex((j-temporaryCIndex(1,1)+1), 2) = currentTemp; %stores the max H-index value for each year in temporaryHIndex      
    end
    
    for j = 1:length(temporaryCIndex) %for the number of years at the current station
        for h = 1:length(CIndex) %for full array of years being analyzed
           if temporaryCIndex(j,1) == CIndex(h,1) %Checks to make sure that the years are the same for the given station
              CIndex(h,i+1)=temporaryCIndex(j,2); %if the years are the same, then the yearly value for the station is stored in a column specifically for that station
           end
        end
    end
end
x = CIndex(:,1);
y1 = CIndex(:,2);
y2 = CIndex(:,3);
plot(x,y1)
hold on
plot(x,y2)