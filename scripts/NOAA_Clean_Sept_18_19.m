%%
%The goal of this script is to clean the data from NOAA scripts and ensure
%that it isn't filed with missing data
clear all
clc

%% 1 
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
    T2 = endsWith(H(1,i),'.csv');
    %if (T1 == 1) && (startsWith(H(1,i),'~$') == 0)
    %    [num, text] = xlsread(strcat(folderName,'/',char(H(1,i)))); %reads each excel file in folderName and records it as two matricies, one for strings, the other for numbers
    %    writematrix(num,strcat(newFolder,'/',char(text(1,2)),'.csv'));%writes the num matrix as a csv file with a name drawn from the title position of the text matrix text matrix
    %    stationNames(1,(i-2)) = convertCharsToStrings(strcat(char(text(1,2)),'.csv')); %adds the current station name to the stationNames variable
    %end
    if (T2 == 1) && (startsWith(H(1,i),'~$') == 0)
        [num, text] = xlsread(strcat(folderName,'/',char(H(1,i)))); %reads each excel file in folderName and records it as two matricies, one for strings, the other for numbers
        C = unique(text(:,2));
        %A = CIndex(CIndex(:,1+i) ~= 0,[ 1 1+i]); 
        nameColumn = 2;
        for h = 1:length(C)-2
            temporaryIndex = {};
            for j = 1:length(text) 
                T3 = convertCharsToStrings(char(text(j, nameColumn)));
                T4 = convertCharsToStrings(char(C(h,1)));
                if T3 == T4
                    temporaryIndex(j,:) = text(j,nameColumn);
                else
                    temporaryIndex(j,:) = {'0'};
                end
            end
            A = temporaryIndex(convertCharsToStrings(char(temporaryIndex(1,:)) ~= '1'),:);
            
        end
        %writematrix(num,strcat(newFolder,'/',char(text(1,2)),'.csv')); %writes the num matrix as a csv file with a name drawn from the title position of the text matrix text matrix
        %stationNames(1,(i-2)) = convertCharsToStrings(strcat(char(text(1,2)),'.csv')); %adds the current station name to the stationNames variable
   end
end

