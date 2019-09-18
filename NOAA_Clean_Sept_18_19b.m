%%
%The goal of this script is to clean the data from NOAA scripts and ensure
%that it isn't filed with missing data
clear all
clc

%% 1 
%seperates the NOAA data into stations and removes NaN values from the TMIN
%and TMAX data set. Currently the program removes too much and will need to
%be adjusted to substitue Dr. Rahmani's data for any missing data. 
clc
folderName = 'Weather'; %variable for easy change of folder name
folderInfo = dir(folderName);  %creates a structure array with all the file names in "folderName"
folderLength = length(folderInfo); 
stationNames = strings([1,(folderLength-2)]); %creates an open array for the station names that are in folderName folder
H = {};
TMaxColumn = 9;
TMinColumn = 10;
newFolder = strcat(folderName,'_CSV'); %creates a variable for the folder name that will be used to store the new CSV files.
mkdir(newFolder); %creates the the new folder to store the CSV files.

for i = 2:folderLength
    H(1,i) = cellstr(folderInfo(i).name); %converts the cells in "folderInfo" into strings
    T1 = endsWith(H(1,i),'.xlsx'); %checks if the file is an xlsx file
    T2 = endsWith(H(1,i),'.csv');
    %if (T1 == 1) && (startsWith(H(1,i),'~$') == 0)
    %   [num, text] = xlsread(strcat(folderName,'/',char(H(1,i)))); %reads each excel file in folderName and records it as two matricies, one for strings, the other for numbers
    %   writematrix(num,strcat(newFolder,'/',char(text(1,2)),'.csv'));%writes the num matrix as a csv file with a name drawn from the title position of the text matrix text matrix
    %   stationNames(1,(i-2)) = convertCharsToStrings(strcat(char(text(1,2)),'.csv')); %adds the current station name to the stationNames variable
    %end
    if (T2 == 1) && (startsWith(H(1,i),'~$') == 0)
        table = readtable(strcat(folderName,'/',char(H(1,i)))); %reads each excel file in folderName and records it as two matricies, one for strings, the other for numbers
        C=unique(table.NAME); %find the unique station names
        %toDelete = isnan(table.TMIN) == 1;
        %table(toDelete,:) = [];
        for j = 1:length(C)
            temporaryName = convertCharsToStrings(char(C(j,1))); %create a variable for station names
            toDivide = table.NAME == temporaryName; %find where the station name exists in the NOAA table
            newTable = table(toDivide,:); %create a seperate table that only contains a given stations names
            for h = min(newTable.YEAR):max(newTable.YEAR)
                B = newTable(newTable.YEAR == h,:);
                T = isnan(B.TMAX);
                for n = 1:length(T)
                    if T(n) == 1 
                        newTable(newTable.YEAR == h,:) = [];
                    end
                end
            end
            %newTable(isnan(newTable.TMIN) == 1,:) = [];
            %newTable(newTable.YEAR == h,:) = [];
            writetable(newTable,strcat(newFolder,'/',temporaryName,'.csv')); %writes the num matrix as a csv file with a name drawn from the title position of the text matrix text matrix
            stationNames(1,j) = strcat(temporaryName,'.csv'); %adds the current station name to the stationNames variable
            
        end
        
        %The loop below isolates the data from each NOAA station
        
            %writematrix(num,strcat(newFolder,'/',char(C(h,1)),'.csv')); %writes the num matrix as a csv file with a name drawn from the title position of the text matrix text matrix
            %stationNames(1,h) = strcat(char(C(h,1)),'.csv'); %adds the current station name to the stationNames variable
            
        
    end
end
y = unique(newTable.YEAR);
%% 2
clc

stationLength = length(stationNames);
%stationLength = 1;
folder = strcat(pwd,'/',newFolder); %calls the path of the current file directory
for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path
    temporaryFile = readmatrix(fullFileName); %creates a temporary matrix of the the data for the current station name.
    yearColumn = 1;
    dayColumn = 3;
    TMaxColumn = 7;
    TMinColumn = 8;
    range = transpose(min(temporaryFile(:,yearColumn)):max(temporaryFile(:,yearColumn)));
    for j = range(1,1):range(end,1)%for each year at this station
        year = find(temporaryFile==j); %locates the index values for the given year
        B = temporaryFile(year,:); %creates a temporary matrix for the given year
        if B(:,TMaxColumn)
            
        end
    end
end
temporaryFile = temporaryFile(isnan(temporaryFile(:,TMinColumn)) ~= 1,:);