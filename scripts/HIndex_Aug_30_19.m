% This first block of code reads the data files given and converts them into
% csv files. 
clc
clear
%Note: Run Speed: xlsread fastest(~15s), readmatrix(1.5minutes), readtable(2.5minute), then
%readcell(4.5min)
[num,text] = xlsread('KS_144972_CO18902009.xlsx');
writematrix(num, 'KS_Manhattan.csv');

%% 2
%okay, that works, now I need to write it so that the program will read
%through a file, locate the spreadsheet files and then convert them into CSV
%files with appropriate names. 
%% 3
% this code puts all of the names of a folder into a cell array
FolderName = 'Weather';
folderInfo = dir(FolderName); %creates a structure array of all the files within a folder 
%Note: this function does not work if your current folder is the one you 
%are performing the dir function on

folderLength = length(folderInfo); %establishes the size of the loop below
H = {}; %creates an empty cell array

for i = 2:folderLength %the array starts at two because the first two values in the cell are filler values
   H(1,i) = cellstr(folderInfo(i).name); %stores each file name into the cell array
end
%% 4
% This code will create csv files of all the excel files in the folder
% described in the step above and will add the name of the station.
for i = 2:folderLength
    T = endsWith(H(1,i),'.xlsx'); %checks if the file is an xlsx file
    if (startsWith(H(1,i),'~$') == 1) 
        H(1,i) = 0;
    end
    if (T == 1)
        [num, text] = xlsread(char(H(1,i))); %reads each excel file in the folder and records it as a matrix
        writematrix(num,strcat(char(text(1,2)),num2str(i),'.csv')); %writes the num matrix as a csv file
    end
        
end
%% 5
% simplifeid version of blocks 3&4 that also adds an output array that
% contains the names of the stations

folderName = 'Weather'; %variable for easy change of folder name
folderInfo = dir(folderName);  %creates a structure array with all the file names in "folderName"
folderLength = length(folderInfo); 
stationNames = strings([1,(folderLength-2)]);
H = {}; 
TMaxColumn = 4;
TMinColumn = 5;

for i = 2:folderLength
    H(1,i) = cellstr(folderInfo(i).name); %converts the cells in "folderInfo" into strings
    T = endsWith(H(1,i),'.xlsx'); %checks if the file is an xlsx file
    if (startsWith(H(1,i),'~$') == 1) 
        H(1,i) = 0;
    end
    if (T == 1)
        [num, text] = xlsread(char(H(1,i)));
        for j = 1:length(num) %For each file read, changes the temperature units from C to F
            num(j, TMaxColumn) = (num(j, TMaxColumn)*9/5)+32; %convert tempuratures from C to F for TMax
            num(j, TMinColumn) = (num(j, TMinColumn)*9/5)+32; %convert tempuratures from C to F for TMin
        end
        writematrix(num,strcat(char(text(1,2)),num2str(i),'.csv'));
        stationNames(1,(i-2)) = convertCharsToStrings(strcat(char(text(1,2)),num2str(i),'.csv'));
    end
        
end
%% 6
%this will be the first steps of creating the C/H index. it will start by
%anazyling a single station on a yearly basis and produce a bar graph at
%the end


%stationLength = length(stationNames);
stationLength = 1;
folder = 'C:/users/mstone2232/Desktop'; %this section is not robust and will need to be adjusted at a later date.
for i = 1:stationLength %for each station
    baseFileName = stationNames(i);
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path
    tempFile = readmatrix(fullFileName);
    yearColumn = 1;
    dayColumn = 3;
    TMaxColumn = 4;
    TMinColumn = 5;
    %creates an array from the starting year to the ending year of the stations available weather data
    tempHIndex = transpose(min(tempFile(:,yearColumn)):max(tempFile(:,yearColumn)));
    counter = 0;
    for j = tempHIndex(1,1):tempHIndex(end,1)
        %Step 1: produce H index for one year       
        %Find the max T for a given year
        
        year = find(tempFile==j); %locates the index values for the given year
        B = tempFile(year,:); %creates a temporary matrix for the given year
        currentTemp = round(max(B(:,TMaxColumn)));  
        
        %Count the number of times where the daily temp is greater than or
                %equal to that temp
        while counter < currentTemp
            counter = 0;
            for h = 1:length(B)
                
                if B(h,TMaxColumn) >= currentTemp
                    counter = counter + 1;
                end
                
            end
            %if the counter is less than the value of the temperature,
                %reduce the temperature by 1 degree F and repeat the loop
            if counter < currentTemp
                currentTemp = currentTemp - 1;
            end
        end
        tempHIndex((j-tempHIndex(1,1)+1), 2) = currentTemp; %stores the max H-index value for each year in tempHIndex      
    end
end
writematrix(tempHIndex,'test.csv');
x = tempHIndex(:,1);
y = tempHIndex(:,2);
%barh(x,y)