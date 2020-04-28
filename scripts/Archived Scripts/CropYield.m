%% 1 Universal Constants
folderName2 = 'Crop_Yield'; %variable for easy change of folder name
folderInfo2 = dir(folderName2);  %creates a structure array with all the file names in "folderName"
folderLength2 = length(folderInfo2); 
B = struct2cell(folderInfo2);
for i = 3:folderLength2
   names(1,(i-2)) =  string(B(1,i));
end
newFolder2 = strcat(folderName2);
folder2 = strcat(pwd,'/',newFolder2);
%creates number arrays and tables for the various decades and climate zones
%used in Kansas for this study
decades = cell2table({1911:1920 1921:1930 1931:1940 1941:1950 1951:1960 1961:1970 1971:1980 1981:1990 1991:2000 2001:2010}'); %creates a table for each decade between 1911 and 2010
west = [2,3,4,8,9,14,20,21]; 
EWcentral = [5,10,11,15,16,17];
east = [6,7,12,13,18,19,22,23,24];

north = [2,3,4,5,6,7,11,12];
NScentral = [8,9,10,13,17];
south = [14,15,16,18,19,20,21,22,23,24];
whole = 2:24;
climateZone = {west, EWcentral, east, north, NScentral, south,whole};

baseFileName2 = names; %this is the name of the file excluding file type. 
fullFileName2 = fullfile(folder2, baseFileName2); %creates a variable for the full file path to ensure no errors related to file path
temporaryFile2 = readtable(fullFileName2); %creates a temporary matrix of the the data for the current station name.
%% creates a table for the wheat yields in bushels/acre for 1911 to 2010
wYield = table();
wYield.YEAR = temporaryFile2.Year(temporaryFile2.Year >= 1911 & temporaryFile2.Year <= 2010); 
wYield.YIELD = round(temporaryFile2.Value(temporaryFile2.Year >= 1911 & temporaryFile2.Year <= 2010),1);

base = mean(wYield.YIELD);
yDecades = decades;
%tis loop calculates the average H/L index for each decade from 1911-2010
for j = 1:height(yDecades)          
            yDecades.MEAN(j) = mean(wYield.YIELD(wYield.YEAR <= max(yDecades.Var1(j,:)) & wYield.YEAR >= min(yDecades.Var1(j,:))));            
            yDecades.DIFF(j) = yDecades.MEAN(j)-base;        
        decadeNames (j) = compose(num2str(min(yDecades.Var1(j,:))-1)+"s");         
end
figure('name','Decadal');
bar(categorical(decadeNames),yDecades.MEAN); %creates a bar graph showing the difference between the decadal average and the centenial average of 1911-2010


title("Difference of Grain Yeild decadal and centenial average for 1911-2010")
ylabel("Change in Grain Yield (Bu/acre/year)")
xlabel("Decade")

figure('name', 'Yearly')
bar(categorical(wYield.YEAR),wYield.YIELD);
