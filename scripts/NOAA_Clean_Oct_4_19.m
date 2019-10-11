%% 1
%The goal of this script is to clean the data from NOAA scripts and ensure
%that it isn't filed with missing data
clear all
clc
%% 2

%seperates the NOAA data into stations and removes NaN values from the TMIN
%and TMAX data set. Currently the program removes too much and will need to
%be adjusted to substitue Dr. Rahmani's data for any missing data. 
clc
folderName = 'Weather'; %variable for easy change of folder name
folderInfo = dir(folderName);  %creates a structure array with all the file names in "folderName"
folderLength = length(folderInfo); 
stationNames = strings([1,(folderLength-2)]); %creates an open array for the station names that are in folderName folder
H = {};
newFolder = strcat(folderName,'_CSV'); %creates a variable for the folder name that will be used to store the new CSV files.
mkdir(newFolder)
for i = 2:folderLength
    H(1,i) = cellstr(folderInfo(i).name);
    T = endsWith(H(1,i),'.csv'); %For scripts that are already .csv files, T2 is 1
    if (T == 1) && (startsWith(H(1,i),'~$') == 0)
        table = readtable(strcat(folderName,'/',char(H(1,i)))); %reads each excel file in folderName and records it as two matricies, one for strings, the other for numbers
        C=unique(table.NAME); %find the unique station names
        for j = 1:length(C)
            temporaryName = string(C(j,1)); %create a variable for station names
            toDivide = table.NAME == temporaryName; %find where the station name exists in the NOAA table
            newTable = table(toDivide,:); %create a seperate table that only contains a given stations names
            %for h = min(newTable.YEAR):2018 %for an array the length of years in the stations recorded history
            %    B = newTable(newTable.YEAR == h,:); %creates a new table, B, that is just the current year, 
            %         newTable(newTable.YEAR == h,:) = [];                               
            %end 
            %newTable(isnan(newTable.TMIN) == 1,:) = [];
            %newTable(newTable.YEAR == h,:) = [];
            writetable(newTable,strcat(newFolder,'/',temporaryName,'.csv')); %writes the num matrix as a csv file with a name drawn from the title position of the text matrix text matrix
            stationNames(1,j+2) = strcat(temporaryName,'.csv'); %adds the current station name to the stationNames variable
        end
    end
end

%% 2 
%this section will be to check if the data from the NOAA files fits within
%our error tolerance (10% per year, 5% per month of missing data)

%probably just create a cell array of all the station names, that will be
%fastest. But might write a script to compare the current station name to
%the other station names in the file.


%a for loop to look at each station name
    % search an array of all the station names in a given file to the
    % current station name provided by the cell array
        %if the name is the same, and the error is within tolerance, then
        %merge the two files, probably by using tables since that seems to
        %be the easiet method. The trouble will be that NOAA column names
        %might be different than our files so that might result in needing
        %to make some changes to adjust for that
    %End loop