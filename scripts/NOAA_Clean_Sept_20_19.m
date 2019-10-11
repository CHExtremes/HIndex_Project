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
newFolder = strcat(folderName,'_CSV'); %creates a variable for the folder name that will be used to store the new CSV files.
mkdir(newFolder); %creates the the new folder to store the CSV files.

for i = 2:folderLength
    H(1,i) = cellstr(folderInfo(i).name); %converts the cells in "folderInfo" into strings
    T1 = endsWith(H(1,i),'.xlsx'); %checks if the file is an xlsx file
    TMaxColumn = 4;
    TMinColumn = 5;

    if (T1 == 1) && (startsWith(H(1,i),'~$') == 0)
       [num, text] = xlsread(strcat(folderName,'/',char(H(1,i)))); %reads each excel file in folderName and records it as two matricies, one for strings, the other for numbers
       for j = 1:length(num) %For each file read, changes the temperature units from C to F
            num(j, TMaxColumn) = (num(j, TMaxColumn)*9/5)+32; %convert tempuratures from C to F for TMax
            num(j, TMinColumn) = (num(j, TMinColumn)*9/5)+32; %convert tempuratures from C to F for TMin
       end
       file = [text; num2cell(num)];
       writecell(file,strcat(newFolder,'/',char(text(1,2)),'.csv'));%writes the num matrix as a csv file with a name drawn from the title position of the text matrix text matrix
       stationNames(1,(i-2)) = convertCharsToStrings(strcat(char(text(1,2)),'.csv')); %adds the current station name to the stationNames variable
    end
end
%% 2
%I split off this section because it will be faster to test it if I don't
%have to run the excel loops more than once. 
TMaxColumn = 9;
TMinColumn = 10;
for i = 2:folderLength
    T2 = endsWith(H(1,i),'.csv'); %For scripts that are already .csv files, T2 is 1
    if (T2 == 1) && (startsWith(H(1,i),'~$') == 0)
        table = readtable(strcat(folderName,'/',char(H(1,i)))); %reads each excel file in folderName and records it as two matricies, one for strings, the other for numbers
        C=unique(table.NAME); %find the unique station names
        %toDelete = isnan(table.TMIN) == 1;
        %table(toDelete,:) = [];
        for j = 1:length(C)
            temporaryName = convertCharsToStrings(char(C(j,1))); %create a variable for station names
            toDivide = table.NAME == temporaryName; %find where the station name exists in the NOAA table
            newTable = table(toDivide,:); %create a seperate table that only contains a given stations names
            %The goal of the below loop is to remove years that have large
            %sections of NaN values without removing years that have single
            %NaN values. 
            for h = min(newTable.YEAR):max(newTable.YEAR) %for an array the length of years in the stations recorded history
                B = newTable(newTable.YEAR == h,:); %creates a new table, B, that is just the current year, h
                if height(B) ~= 365 && height(B) ~= 366 % if the number of elements in the given year is not 365 or 366, remove that year from the data set.
                     newTable(newTable.YEAR == h,:) = [];
                else
                    for n = 1:height(B) %for the given year, Check each TMAX value, if it is NaN, check the values adjacet to it. If those are also NaN, then remove that year from the data set. 
                    T = isnan(B.TMAX(n));
                        if n == 1
                            T2 = isnan(B.TMAX(n+1));
                            if T == 1 && T2 == 1
                                newTable(newTable.YEAR == h,:) = [];
                            end
                        elseif n  == height(B)
                            T3 = isnan(B.TMAX(n-1));
                            if T == 1 && T3 == 1
                                newTable(newTable.YEAR == h,:) = [];
                            end
                        else
                            T2 = isnan(B.TMAX(n+1));
                            T3 = isnan(B.TMAX(n-1));
                            if T == 1 && T2 == 1 && T3 == 1 
                                newTable(newTable.YEAR == h,:) = [];
                            end                   
                        end
                    end
                end
            end 
            %newTable(isnan(newTable.TMIN) == 1,:) = [];
            %newTable(newTable.YEAR == h,:) = [];
            writetable(newTable,strcat(newFolder,'/',temporaryName,'.csv')); %writes the num matrix as a csv file with a name drawn from the title position of the text matrix text matrix
            stationNames(1,j+2) = strcat(temporaryName,'.csv'); %adds the current station name to the stationNames variable
           
        end
        X= min(newTable.YEAR):max(newTable.YEAR);
        %The loop below isolates the data from each NOAA station
        
            %writematrix(num,strcat(newFolder,'/',char(C(h,1)),'.csv')); %writes the num matrix as a csv file with a name drawn from the title position of the text matrix text matrix
            %stationNames(1,h) = strcat(char(C(h,1)),'.csv'); %adds the current station name to the stationNames variable
            
        
    end
end
y = unique(newTable.YEAR);
%% 3
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