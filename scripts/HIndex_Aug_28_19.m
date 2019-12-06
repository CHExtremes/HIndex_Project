% This firs block of code reads the data files given and converts them into
% csv files. 
clc
clear

num = readmatrix('KS_144972_CO18902009.xlsx');
writematrix(num, 'KS_Manhattan.csv');

%% 2
%okay, that works, now I need to write it so that the program will read
%through a file, locate the spreadsheet files and then convert them into CSV
%files with appropriate names. 
%% 3
% this code puts all of the names of a folder into a cell array
FolderName = 'Weather'
Info = dir(FolderName); %creates a structure array of all the files within a folder 
%Note: this function does not work if your current folder is the one you 
%are performing the dir function on

s = length(Info); %establishes the size of the loop below
H = {}; %creates an empty cell array

for i = 2:s %the array starts at two because the first two values in the cell are filler values
   H(1,i) = cellstr(Info(i).name); %stores each file name into the cell array
end
%% 4
% This code will create csv files of all the excel files in the folder
% described in the step above
for i = 2:s
    T = endsWith(H(1,i),'.xlsx'); %checks if the file is an xlsx file
    if (startsWith(H(1,i),'.~$') == 1) 
        H(1,i) = 0;
    end
    if (T == 1) && (T1 == 0 )
        tempMat = readmatrix(char(H(1,i)));
        writematrix(tempMat,strcat('penguins',num2str(i),'.csv'));
    end
        
end
%% 5
%This is a section to test part of the above code in isolation
if (T == 1) && (T1 == 0 )
        tempMat = readmatrix(char(H(1,7)));
        writematrix(tempMat,strcat('penguins',num2str(7),'.csv'));
else
    print('nope');
end