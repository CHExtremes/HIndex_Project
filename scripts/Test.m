
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
p = 90;
clc
tic
%stationLength = 1;
p = 90;
folder = strcat(pwd,'/',newFolder); %calls the path of the current file directory
percentileHIndex = zeros(121,24);
percentileHIndex(:,1) = (1893:2013)';
baseFileName = stationNames(3); %this is the name of the file excluding file type. 
fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path
temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
for i = min(temporaryFile.YEAR):max(temporaryFile.YEAR)
   counter = 0;
   year = temporaryFile(temporaryFile.YEAR == i,:);
   a = prctile(year.TMAX,p);
   b(i-min(temporaryFile.YEAR)+1,1) = a;
   for j = 1:height(year)
       if year.TMAX(j)>= a
          counter = counter + 1; 
       end
   end  
   percentileHIndex(i-min(temporaryFile.YEAR)+1,2) = counter;
   counter = 0;
end
    
    
    