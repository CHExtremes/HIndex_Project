%%
%This code is designed to take files from a specific folder, if they are
%excel files they are converted into csv files. those csv files are then
%used to produce an H-index based on temperature for each station.

%Note: you will need to download the "Weather folder" and add it to your
%directory file to ensure this program works. 
clear 
clc

%% 1 Universal Constants
%Weather data
folderName = 'Weather_CSV'; %variable for easy change of folder name
folderInfo = dir(folderName);  %creates a structure array with all the file names in "folderName"
folderLength = length(folderInfo); %creates a variable for the lenght of the folderInfo variable
B = struct2cell(folderInfo); %creates a temporary variable, B, with folderInfor as cell rather than struct
for i = 3:folderLength
   stationNames(1,(i-2)) =  string(B(1,i)); %creates a variable for the station names used for all indices
   tableStationNames(1,(i-2))= erase(stationNames(1,(i-2)),".csv"); %creates an adjusted variable for using the station names for table loops later on
end
newFolder = strcat(folderName); %creates a variable for the folder Name given by folderName

%the Below script repeates what the above does, but for the grain yeilds
%instead of the weather data
%constants for crop yield data
folderName2 = 'Crop_Yield'; %variable for easy change of folder name
folderInfo2 = dir(folderName2);  %creates a structure array with all the file names in "folderName"
folderLength2 = length(folderInfo2); 
B = struct2cell(folderInfo2);
for i = 3:folderLength2
   names(1,(i-2)) =  string(B(1,i));
end
newFolder2 = strcat(folderName2);

folder = strcat(pwd,'/',newFolder); %concatenates the string of the current folder (pwd) and the newFolder name of weather data. This allows for easy identification of the folder later on
folder2 = strcat(pwd,'/',newFolder2); %repeat for grain data

%creates number arrays and tables for the various decades and climate zones
%used in Kansas for this study. These describe the climate zones used in
%grain yield folders
decades = cell2table({1911:1920 1921:1930 1931:1940 1941:1950 1951:1960 1961:1970 1971:1980 1981:1990 1991:2000 2001:2010}'); %creates a table for each decade between 1911 and 2010
west = [10, 20, 30];
EWcentral = [40, 50, 60];
east = [70, 80, 90];

north = [10, 40, 70];
NScentral = [20, 50, 80];
south = [30, 60, 90];
whole = 10:10:90;
agClimateZone = {west, EWcentral, east, north, NScentral, south, whole};
%weather climate zones
west = [2,3,4,8,9,14,20,21]; 
EWcentral = [5,10,11,15,16,17];
east = [6,7,12,13,18,19,22,23,24];

north = [2,3,4,5,6,7,11,12];
NScentral = [8,9,10,13,17];
south = [14,15,16,18,19,20,21,22,23,24];
whole = 2:24;
climateZone = {west, EWcentral, east, north, NScentral, south,whole};
%% 2 Calculates the decadal averages for grain yield
%this scrip calculates the decadal average of grain yield. Currently this
%only does this for the grain yield of wheat by county. I will need to
%adjust the script later to do this for all files in the grain yield folder
%if we decide to continue with the comparision to grain yield
y1 = 1908;
y2 = 2013;
%tempNames = {'totWheat','irrWheat','noIrrWheat','totCorn','irrCorn','noIrrCorn','totSoy','irrSoy','noIrrSoy', 'totSorgh','irrSorgh','noIrrSorgh'};
Yield = struct('totWheat',[],'irrWheat',[],'noIrrWheat',[],'totCorn',[],'irrCorn',[],'noIrrCorn',[],'totSoy',[],'irrSoy',[],'noIrrSoy',[], 'totSorgh',[],'irrSorgh',[],'noIrrSorgh',[]);%{'totWheat','irrWheat','noIrrWheat','totCorn','irrCorn','noIrrCorn','totSoy','irrSoy','noIrrSoy', 'totSorgh','irrSorgh','noIrrSorgh'});
Yield;

T = fieldnames(Yield); %creates a cell array containing the field names in the Yield structure array

%S = struct('type','.','subs','REGION');

for i = 1:length(names)
    baseFileName2 = names(i); %calls the file for grain yield at name i in the grain yeild folder
    fullFileName2 = fullfile(folder2, baseFileName2); %adjust the file name to be redundant and less likely to cause errors
    temporaryFile2 = readtable(fullFileName2); %creates a table for the file with the name given from the above variables
    
    tYield = table(); %creates a new empety table Yield to store the data for grain yield for the given dacade
    tYield.YEAR = temporaryFile2.Year(temporaryFile2.Year >= y1 & temporaryFile2.Year <= y2); %sets the YEAR variable of the Yield table to be equal to the Year variable from temporaryFile2 for values of 1911 to 2010 
    %I need to adjust these scripts to have more flexible year selection
    tYield.REGION = temporaryFile2.AgDistrictCode(temporaryFile2.Year >= y1 & temporaryFile2.Year <= y2); %sets the REGION variable of Yield equal to the AgDistrictCode variable from tempoaryFile2 where the Year variable is between 1911 and 2010 
    tYield.YIELD = round(temporaryFile2.Value(temporaryFile2.Year >= y1 & temporaryFile2.Year <= y2),1);%sets the YIELD variable of Yield equal to the Value variable from tempoaryFile2 where the Year variable is between 1911 and 2010
    %Note, for this section to work, you do need to order the files in the
    %folder to be in the same order as the field names in the Yield
    %structure. This shouldn't be a major issue if you manage the files
    %correctly and means that most of the editing is just adding a number
    %prefix to files to achieve the proper order.
    S = struct('type','.','subs',char(T(i,1))); %creats a indexing structure array with for a dot type and a field designated by the field names of the Yield array at position i of the loop
    Yield = subsasgn(Yield,S,table2cell(tYield)); %converts the tYield table to a cell array and then assigns that cell array to the field described by S above
    %note, for each field of the structure, column 1 is 
    
    base = mean(tYield.YIELD); %base is the baseline average determiend as average of all yields in all regions across the entire period of record
    yDecades = decades; %creates a new table yDecades to show the results of the grain yields by decade and region
    %this loop calculates the average grain yield by decade
    for j = 1:height(yDecades) %for each decade
        for h = 1:7 %for each climate region (W, EW central, E, N, NS central, S, Entire state)
                z = cell2mat(agClimateZone(h)); %creates a double array, z, from the cell at position h in agClimateZone
                tempAvgYield = tYield; %createes a temporay table that is adjusted with each iteration of the loop this allows the Yield variable to be used in later scripts to calculate yearly values
                if h ~= 7 %if h is any region but the entire state
                    tempAvgYield = tempAvgYield(tempAvgYield.REGION == z(1) | tempAvgYield.REGION == z(2) | tempAvgYield.REGION == z(3),:); %limits temp yield to rows where the rows are equal to z at postions 1, 2, and 3 (I think I could do this with a third loop, but its not worth my time right now)
                end
            yDecades.MEAN(j,h) = mean(tempAvgYield.YIELD(tempAvgYield.YEAR <= max(yDecades.Var1(j,:)) & tempAvgYield.YEAR >= min(yDecades.Var1(j,:))));  %calculates the average grain yield of the decade for each region          
            yDecades.DIFF(j,h) = yDecades.MEAN(j,h)-base; %calculates the difference between the decadal average and the base average       
            decadeNames (j) = compose(num2str(min(yDecades.Var1(j,:))-1)+"s");   %creates an array for the period for graphing (i.e. 1920s, 1930s, etc)         
        end
    end
end

%%
figure('name','Decadal');
bar(categorical(decadeNames),wYDecades.MEAN); %creates a bar graph showing the difference between the decadal average and the centenial average of 1911-2010

figure('name', 'Yearly')
bar(categorical(wYield.YEAR),wYield.YIELD); %creates a bar graph showing the total yield vs. time on a yearl basis

%% Calculating All indecies test script
clc
tic
stationLength = length(stationNames); %variable for the number stations
Index = struct('HIndex',[],'LIndex',[],'CIndex',[],'WIndexmm',[],'WIndexin',[],'DIndex',[]);
T = fieldnames(Index); %creates a cell array containing the field names in the Yield structure array
%creats a indexing structure array with for a dot type and a field designated by the field names of the Index array at position i of the loop

%pre-allocates space in the Index structure to reduce runtime
for i = 1:6
    S = struct('type','.','subs',char(T(i,1)));      
    Index = subsasgn(Index,S,zeros(2013-1890+1,24));
    
end
Index.HIndex(:,1) = (1890:2013)';
Index.LIndex(:,1) = (1890:2013)';
dryness = 1; %must have less than 1mm of precip to be considered a "dry" day
%for i=1:6 %number of indices
    %T = fieldnames(Index); %creates a cell array containing the field names in the Yield structure array
   % S = struct('type','.','subs',char(T(i,1)));  %creats a indexing structure array with for a dot type and a field designated by the field names of the Index array at position i of the loop
   
    for j = 1:stationLength
       q = j
       baseFileName = stationNames(j); %this is the name of the file excluding file type. 
       fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path
       temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
       temporaryFile.RAIN = round(temporaryFile.RAIN,0); %rounds the precipitation column to be whole numbers since you can't have fractional days at our time scale
       temporaryFile.RAIN2 = round(temporaryFile.RAIN/25.4,2); %creates a second column for precipitation in inches
       
       %YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
       temporaryIndex = table(transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR)),'VariableNames',"YEAR"); %creates an column array for the years of the H-Indecies        
       temporaryIndex.Index = zeros(height(temporaryIndex),1);
       
       criteria = table();%creates an array of counters and targets that can be used similtaneously
       
       %for the variable table, position 1 is H, 2 is L, 3 is C, 4 is
       %W(mm), 5 is W(in), and 6 is D index.
       criteria.counter = zeros(length(T),1);
       criteria.target = zeros(length(T),1); %target will be what the condition variable is measured against (for example, H index target is the max Tmax value)
                  
       
       for h = temporaryIndex.YEAR(1):temporaryIndex.YEAR(end)%for each year at this station
           year = temporaryFile(temporaryFile.YEAR==h,:); %create a temporary table containing only the rows with the current year
           
           %This next step calculates the target value (i.e the maximum
           %possible value for dry period length, temperature, or
           %precipitation) for each index and stores that variable in the
           %"target" field of the "criteria" table
           criteria.target(1) = max(year.TMAX); %H
           criteria.target(2) = max(year.TMIN); %L
           criteria.target(3) = 32 - min(year.TMIN); %The maximum difference between 0°C and the minimum temperature for a given year
           criteria.target(4) = max(year.RAIN); %W1
           criteria.target(5) = max(year.RAIN2);%W2
           
           tempDryPeriod = 0; %temporary variable used to count the longest dry period
           %this loop identifies the longest dry period for the given year
           
           for D = 1:height(year)%for the current year                
               if year.RAIN(D) < dryness %if the precipitation at day (D) is less than the minimum threshold for dryness
                    tempDryPeriod = tempDryPeriod + 1; %increase the number of dry days by one
               else %else, we encounter a day with precipitation 
                    if tempDryPeriod > criteria.target(6) %if number of days in tempDryPeriod is greater than the number of days in the target dry period
                         criteria.target(6) = tempDryPeriod; %change the value of the current dry period (target at position 6) to be the value of the dry period coutner (tempDryPeriod)
                    end
                    tempDryPeriod = 0; %reset the counter and repeat the loop
               end          
           end  
           
           %criteria.counter([1 2]) < criteria.target([1 2])
           %condition = 0; %to insure the while loop starts, condition is initially set to 0
           %A = sum((criteria.counter([1 2]) < criteria.target([1 2]))) ~= 2;
           %B = logical([1 1]')
           %C = A~=B
           %l = 0
           criteria.target(2:6) = 0;
           while sum((criteria.counter(:) < criteria.target(:))) ~= 0 %while the condition variable is not true
               %l = l +1
               criteria.counter(:) = 0;%resets all of the counter variables to 0 at the start of each iteration of the while loop
               
               for k = 1:height(year) %for the current year
                   
                   %H-Index
                   if year.TMAX(k) >= criteria.target(1) %If the value at row h and column TMaxColumn are greater than currentTemp                    
                        criteria.counter(1) = criteria.counter(1) + 1; %increase counter by 1°F
                   end
%                    %L-Index
%                    if year.TMIN(k) >= criteria.target(2) %If the value at row k in column TMIN are greater than the current target Tmin                    
%                         criteria.counter(2) = criteria.counter(2) + 1; %increase counter by 1°F
%                    end
%                    %C-Index
%                    if year.TMIN(k) <= 32-criteria.target(3) %If the value at row k in column TMIN are less than the current target Tmin                    
%                         criteria.counter(3) = criteria.counter(3) + 1; %increase counter by 1°F
%                    end
%                    %W-index in mm
%                    if year.RAIN(k) >= criteria.target(4) %If the value at row k in column RAIN are less than the current target precipitation                    
%                         criteria.counter(4) = criteria.counter(4) + 1; %increase counter by 1mm
%                    end
%                    %W-index in inches
%                    if year.RAIN2(k) >= criteria.target(5) %If the value at row k in column RAIN2 are less than the current target precipitation                  
%                         criteria.counter(5) = criteria.counter(5) + 0.01; %increase counter by 0.01in
%                    end
%                    %D-Index
%                    if year.RAIN(k) >= dryness %If the value at row k in column TMIN are less than the current target Tmin                    
%                         criteria.counter(6) = criteria.counter(6) + 1; %increase counter by 1
%                    end
                   
               end
               %H
               if criteria.counter(1) < criteria.target(1) %if the counter is smaller than the target Tmax, then the H-index is not valid, so we reduce it by one and repeat the loop.
                  criteria.target(1) = criteria.target(1) - 1; %decrease the target max tmax by 1°F
               end
               %L
%                if criteria.counter(2) < criteria.target(2) %if the counter is smaller than the target Tmin, then the L-index is not valid, so we reduce it by one and repeat the loop.
%                   criteria.target(2) = criteria.target(2) - 1; %decrease the target max tmin by 1°F
%                end
%                %C
%                if criteria.counter(3) < criteria.target(3) %if the counter is smaller than the target Tmin, then the C-index is not valid, so we reduce it by one and repeat the loop.
%                   criteria.target(3) = criteria.target(3) - 1; %decrease the difference between min tmin and 32°F by 1°F
%                end
%                %Wmm
%                if criteria.counter(4) < criteria.target(4) %if the counter is smaller than the target Precip, then the W-index(mm) is not valid, so we reduce it by one and repeat the loop.
%                   criteria.target(4) = criteria.target(4) - 1; %decrease the target precipitation by 1mm
%                end
%                %Win
%                if criteria.counter(5) < criteria.target(5) %if the counter is smaller than the target Precip, then the W-index(mm) is not valid, so we reduce it by one and repeat the loop.
%                   criteria.target(5) = criteria.target(5) - 0.01; %decrease the target precipitation by 0.01in
%                end
%                %D
%                if criteria.counter(6) < criteria.target(6) %if the counter is smaller than the target Precip, then the W-index(mm) is not valid, so we reduce it by one and repeat the loop.
%                   criteria.target(6) = criteria.target(6) - 1; %decrease the target dry period by 1 day
%                end
           end
           
           for k = 1:6
               temporaryIndex.Index(h+1-temporaryIndex.YEAR(1),k) = criteria.target(k); %set the Index value for the current year at the current station equal to the index described by criteria.target(i), basically its just storing the index value for given year in this temporary table              
           end
       end
       %storing the indices for the current station in the Index structure
       %array
       for h = 1:height(temporaryIndex) %for the number of recorded years at the current station
           for k = 1:length(Index.HIndex) %for the maximum number of years across all stations
               if temporaryIndex.YEAR(h) == Index.HIndex(k,1) %if the row in temporaryIndex.YEAR has the same year as the row in the first column of Index.HIndex
                  Index.HIndex(k,(j+1)) = round(temporaryIndex.Index(h,1),0); %stores the first column in the Index field of temporary index (the HIndex values) in the hth column of the HIndex field of the Index structure array 
                  %repeating process done with H index for all the other
                  %indices
                  Index.LIndex(k,(j+1)) = round(temporaryIndex.Index(h,2),0);
                  Index.CIndex(k,(j+1)) = round(temporaryIndex.Index(h,3),0);
                  Index.WIndexmm(k,(j+1)) = round(temporaryIndex.Index(h,4),0);
                  Index.WIndexin(k,(j+1)) = round(temporaryIndex.Index(h,5),0)*100;
                  Index.DIndex(k,(j+1)) = round(temporaryIndex.Index(h,6),0);
                  %there might be a faster/more elegant way to store these
                  %values using subsasgn or subsref, but I couldn't figure
                  %out how that would work.
               end
           end
       end
   end
%end

t = table(1);
S = struct('type','.','subs','newField');
t = subsasgn(t,S,2);
toc
%% 2 Calculating H-index
% This script section calculates the H and L-Index which are a measure of the maximum number of times where Tmax(H-Index) or Tmin(L-Index)  occurs the same number of times as its magnitude in fahrenheit. 

clc
tic %start calculating run time

stationLength = length(stationNames); %variable for the number stations

%for if you want to calculate the L-index(using Tmin) instead of the H-index (using Tmax)
useLIndex = 1; 

%predefining the size of H and L index
HIndex = zeros(2013-1890+1,24);
HIndex(:,1) = (1890:2013)'; %sets the first column as the maximum array of years

if useLIndex == 1
    LIndex = HIndex; %creates a seperate variable to store the results for L-index if it is being used
end

for i = 1:stationLength %for each station
    
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    
    %Rounds the temperature to a whole number, as H and L index cannot be
    %fractions of days
    if useLIndex == 1
        temporaryFile.TMIN = round(temporaryFile.TMIN,0);
    else
        temporaryFile.TMAX = round(temporaryFile.TMAX,0);
    end
    
    YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
    temporaryHIndex = table(YEAR); %creates an column array for the years of the H-Indecies 
    temporaryHIndex.HIndex = zeros(height(temporaryHIndex),1);
    counter = 0;
    for j = temporaryHIndex.YEAR(1):temporaryHIndex.YEAR(end)%for each year at this station
        year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary table for the given year        
        
        %records max temp for the given year
        if useLIndex == 1
            currentTemp = max(year.TMIN);
        else
            currentTemp = max(year.TMAX);
        end
        
        %Count the number of times where the daily temp is greater than or
        %equal to that temp 
        counter = 0;
        while counter < currentTemp %checks to see if the counter is smaller than the currentTemp. This is to make sure that the value is an H-Index value.
            counter = 0;
            for h = 1:height(year)% for days in this year                
                if useLIndex == 1
                    if year.TMIN(h) >= currentTemp %If the value at row h and column TMaxColumn are greater than currentTemp                    
                        counter = counter + 1; %increase counter by 1
                    end
                else
                    if year.TMAX(h) >= currentTemp %If the value at row h and column TMaxColumn are greater than currentTemp                    
                        counter = counter + 1; %increase counter by 1
                    end
                end
            end
            if counter < currentTemp %if the counter is smaller than currentTemp, then the H-index is not valid, so we reduce it by one and repeat the loop.
                currentTemp = currentTemp - 1;
            end
        end              
        temporaryHIndex.HIndex(j-temporaryHIndex.YEAR(1)+1) = currentTemp; %the H index is recoded in the same row as the given year
    end
    
    %this loop records the final H/L index 
    for j = 1:height(temporaryHIndex) %for the number of years at the current station
        if useLIndex == 1
            for h = 1:length(LIndex) %for full array of years being analyzed
               if temporaryHIndex.YEAR(j) == LIndex(h,1) %Checks to make sure that the years are the same for the given station
                  LIndex(h,(i+1))=temporaryHIndex.HIndex(j); %if the years are the same, then the yearly value for the station is stored in a column specifically for that station
               end
            end
        else
            for h = 1:length(HIndex) %for full array of years being analyzed
               if temporaryHIndex.YEAR(j) == HIndex(h,1) %Checks to make sure that the years are the same for the given station
                  HIndex(h,(i+1))=temporaryHIndex.HIndex(j); %if the years are the same, then the yearly value for the station is stored in a column specifically for that station
               end
            end
        end
    end
end
timeHindex = toc; %record the run time as timeHindex
%% 3 graphing H-index time series
%creates a table to store all statistical values for all results for a
%given Index
HResults = table; %pre-define a table for the statistical values calculate later in this section
HResults.NAME = tableStationNames'; %set the rows of this table to be the sation names

%variables to control period analyzed
startYear = 1981;
stopYear = 2010;
useStartYear = 0;
useStopYear = 0;

%creates a new figure named H-Index. Really this should be done for both L
%and  H index using an If statement
figure('Name', 'H-Index')


for i = 1:length(tableStationNames) %for each station in station names 
    %this section creates an array of subplots where each station has its
    %own plot showing the data points in orange, the change in time in
    %blue, and the general trend in black.
    tempNames = split(tableStationNames(i), '_'); %assigns a variable to the stations name and cuts out any unnecissary labels used for organization purposes
    HResults.NAME(i) = tempNames(2,1); %defines station name for the given row of the results table
    
    %for each station, create a new array with the non-zero values and their corresponding years
    if useLIndex == 1
        A = array2table(LIndex(LIndex(:,1+i) ~= 0,[ 1 1+i]));
    else 
        A = array2table(HIndex(HIndex(:,1+i) ~= 0,[ 1 1+i])); 
    end
    
    %creat a new temporary table A to use in graphing
    A.Properties.VariableNames{'Var1'} = 'YEAR'; 
    A.Properties.VariableNames{'Var2'} = 'HIndex';  
    
    %create a new temporary table B that is table A from startYear to
    %stopYear
    B = A(A.YEAR>=startYear & A.YEAR <=stopYear,:);  
    
    %creats a set of subplots with a width of 6 and a height of 4
    subplot(4,6,i)
    %reduces the white space between graphs for the current graph
    pos=get(gca,'Position');
    set(gca,'Position',[pos(1,1) pos(1,2) pos(1,3)+.03 pos(1,4)]) %subplot position
    
   
    lA = plot(A.YEAR,A.HIndex); %adds a line, lA, to the plot to show the change in H/L over time
    lA.LineWidth = .5; %sets the line width to the minimum of 0.5pt
    
    hold on %add each addition graph to the same subplot    
    mdlA = fitlm(A, 'HIndex ~ YEAR'); %performs a linear regression for the Year and the H Index    
    zA = plot(mdlA); %plots the linear regression and data point    
    
    xA = A.YEAR;
    yA = A.HIndex;
    
    meanB = B;
    meanB.HIndex(:) = mean(B.HIndex);
    zB = plot(meanB.YEAR, meanB.HIndex);
    zB.Color = 'r';  
    zB.LineWidth = 1.15;
    %below are changes to the colors and markers of the plot for additional
    %clarity I removed 95% error bars (z(3) and z(4) to make the graph less cluttered. We
    %can turn these on later if we want to visually analyze the error
    %margins.
    zA(1).Color = 'none'; %point color changed to make them invisible for less busy graphs%'#D95319'; %sets data points to be orange
    zA(1).Marker = '.';
    zA(1).MarkerSize = 10;
    zA(2).Color = 'g';
    zA(2).LineWidth = 1.15;
    zA(3).Color = 'none';
    zA(4).Color = 'none'; 
    zA(2).LineStyle = '-';
   
    legend('off'); %hides the automatic legend generated by fitlm
    
    %if statement to set axis on left and bottom edge of subplot matrix
    %instead of on each plot
    if i == 19 || i == 20 || i == 21 || i == 22 || i == 23 
        xlabel('Year', 'FontSize', 11)
        set(gca,'Xtick',[1921 1951 1981 2011]);
        %xticks(gca, 1890:2010, 4)
    else
        xlabel('')
        xticks('')
    end
    
    %if statement that allows greater control on the years examined
    if useStartYear == 1 && useStopYear == 0
        xlim([startYear 2014])
    elseif useStartYear == 1 && useStopYear == 1
        xlim([startYear (stopYear+1)])
    elseif useStartYear == 0 && useStopYear == 1
        xlim([1890 (stopYear+1)])
    else     
        xlim([1890 2014])
    end
    
    % %sets boundries for the y axis for all graphs to be equal
    if useLIndex == 1 
        ylim([55 72])
    else
        ylim([79 95])
    end
    
    %If statement to make sure that only far left subplots have axis labels
    %and ticks marks
    if i == 1 || i == 7 || i == 13 || i == 19
        if useLIndex == 1
            ylabel('L - Index (days/year)', 'FontSize', 11)
        else
            ylabel('H - Index (days/year)', 'FontSize', 11)
        end        
        yticks('auto')
    else
        ylabel('')
        yticks('')
    end
    %this section  creates a table statistically important values for both
    %the POR and the specified period
    
    
    %this code runs trends analysis on x and y as independant and depedant variables
    %it tests the hypothesis of no correlation against the alternative
    %hypothesis of a nonzero correlation. so if p value is smaller than 0.05,
    %we reject the hypothesis.
    %analysis for the whole period of record
    [tau,pD1]=corr(xA,yA,'type','kendall'); %kendall method
    tau_pA(i,1:3)=[i,tau,pD1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
    [rho,p2]=corr(xA,yA,'type','spearman');%spearman method
    rho_pA(i,1:3)=[i,rho,p2];
    [r,p3]=corr(xA,yA);%pearson (linear) method
    r_pA(i,1:3)=[i,r,p3]; %pearson(Least square method) method corrcoef(x,y); 
    
    %Allocates a number to each station based on there Climate Division in
    %the state of Kansas
    
    if HResults.NAME(i,1) == "Saint Francis" 
        HResults.climateDivision(i) = 1;
    elseif HResults.NAME(i,1) ==  "Oberlin"  
        HResults.climateDivision(i) = 1;
    elseif HResults.NAME(i,1) ==   "Colby"
        HResults.climateDivision(i) = 1;      
    elseif HResults.NAME(i) == "Phillipsburg"
        HResults.climateDivision(i) = 2;
    elseif HResults.NAME(i) == "Minneapolis" 
        HResults.climateDivision(i) = 2;    
    elseif HResults.NAME(i) == "Horton" 
        HResults.climateDivision(i) = 3;
    elseif HResults.NAME(i) == "Atchison" 
        HResults.climateDivision(i) = 3;
    elseif HResults.NAME(i) == "Manhattan"
        HResults.climateDivision(i) = 3;    
    elseif HResults.NAME(i) == "Tribune"
        HResults.climateDivision(i) = 4;
    elseif HResults.NAME(i) == "Wakeeney"
        HResults.climateDivision(i) = 4;    
    elseif HResults.NAME(i) == "Hays"
        HResults.climateDivision(i) = 5;
    elseif HResults.NAME(i) == "McPherson"
        HResults.climateDivision(i) = 5;    
    elseif HResults.NAME(i) == "Ottawa"
        HResults.climateDivision(i) = 6;
    elseif HResults.NAME(i) == "Lakin" 
        HResults.climateDivision(i) = 7;
    elseif HResults.NAME(i) == "Elkhart" 
        HResults.climateDivision(i) = 7;
    elseif HResults.NAME(i) == "Ashland"
        HResults.climateDivision(i) = 7;    
    elseif HResults.NAME(i) == "Larned"
        HResults.climateDivision(i) = 8;
    elseif HResults.NAME(i) == "MedicineLodge"
        HResults.climateDivision(i) = 8;    
    elseif HResults.NAME(i) == "Winfield"
        HResults.climateDivision(i) = 9;             
    elseif HResults.NAME(i) == "Sedan" 
        HResults.climateDivision(i) = 9;         
    elseif HResults.NAME(i) == "Independence" 
        HResults.climateDivision(i) = 9;            
    elseif HResults.NAME(i) == "FortScott" 
        HResults.climateDivision(i) = 9;         
    elseif HResults.NAME(i) == "Columbus"
        HResults.climateDivision(i) = 9;         
    end
   
    
    
    HResults.slopePOR(i) = round(table2array(mdlA.Coefficients(2,1)),3);%calls the slope given for the linear regression of the data using the fitlm function
    HResults.adjustedSlopePOR(i) = HResults.slopePOR(i)*100; %this converts the slope from days/year to days/century
    HResults.rPOR(i) = r_pA(i,2);
    HResults.rSqrPOR(i) = mdlA.Rsquared.Ordinary; %calls the r^2 value from the fitlm function and inputs it into a new table
    HResults.rPValuePOR(i) = round(table2array(mdlA.Coefficients(2,4)),3); %calls the pValue from the fitlm function and inputs it into a new table
    HResults.rhoPOR(i) = rho_pA(i,2);
    HResults.rhoPValuePOR(i) = rho_pA(i,3);
    HResults.tauPOR(i) = tau_pA(i,2);
    HResults.tauPValuePOR(i)= tau_pA(i,3);
    HResults.minPOR(i) = min(A.HIndex);
    HResults.maxPOR(i) = max(A.HIndex);
    HResults.meanPOR(i) = round(mean(A.HIndex),1); %calculates the average H index for the station and adds that to a new table
    HResults.medianPOR(i) = median(A.HIndex);
    
    l = yline(HResults.meanPOR(i));
    l.Color = ('m');
    l.LineWidth = 1.25;
    l.LineStyle = '--';
    
    HResults.minSP(i) = min(B.HIndex);
    HResults.maxSP(i) = max(B.HIndex);
    HResults.meanSP(i) = round(mean(B.HIndex),1);
    HResults.medianSP(i) = median(B.HIndex);    
  
    title(compose(tempNames(2,1)+"\n"+num2str(HResults.slopePOR(i)*100,"%#.1f")),'FontSize', 11);
end


%use the below script when you want to automatically make tiff files for
%the given graphs. There is an issue with this function in that it doesn't
%expand the window before saving the file, so it compresses subplots to an
%unreadible level. I'll need to find a way to fix this later if we still
%want to use it over manually saving graphs.
%set(gcf,'PaperPositionMode','auto') %set the print area same as paper
%print('-dtiff','-r600', '12.3.19_H_Index_Subplots')
%% 4.0.1 Calculating time series vs. yearly average grain yield
%This section compares average grain yield to average H-index for every
%year from y1 to y2 
y1 = 1908; %smallest year for climate data
y2 = 2013; %largest year for climate data
r = 1+y2-y1; %the size of the range between y1 and y2 for creating row indices later

if useLIndex == 1
    avgLIndex = table(); %creates a new table for calculating the climate region averages of L Index   
else
    avgHIndex = table(); %creates a new table for calculating the climate region averages of H Index
end


%This loop calculates the regional averages for the climate index based on
%the climate zones for Kansas.
for i = 1:length(climateZone)
    if useLIndex == 1       
        avgLIndex.YEAR = LIndex(:,1);       
        LIndex(LIndex==0) = NaN;
        avgLIndex.LIndex(:,i) = mean(LIndex(:,cell2mat(climateZone(i))),2,'omitnan'); 
    else
        
        avgHIndex.YEAR = HIndex(:,1);
        HIndex(HIndex==0) = NaN;
        avgHIndex.HIndex(:,i) = mean(HIndex(:,cell2mat(climateZone(i))),2,'omitnan');       
    end
end 

if useLIndex == 1
    avgLIndex = avgLIndex(avgLIndex.YEAR >= y1 & avgLIndex.YEAR <= y2,:);
else
    avgHIndex = avgHIndex(avgHIndex.YEAR >= y1 & avgHIndex.YEAR <= y2,:);
end
%This loop creates a structure array that contains the regional averages
%for each crop 

avgYield = struct('totWheat',[],'irrWheat',[],'noIrrWheat',[],'totCorn',[],'irrCorn',[],'noIrrCorn',[],'totSoy',[],'irrSoy',[],'noIrrSoy',[], 'totSorgh',[],'irrSorgh',[],'noIrrSorgh',[]); 
% avgWYield.YEAR = zeros(r,1); %Pre-allocating to reduce errors
% avgWYield.YIELD = zeros(r,7);
for h = 1:length(names)
    
    T = fieldnames(avgYield); %creates a cell array containing the names of all the crop s in the avgYield structure array
    S = struct('type','.','subs',T(h,1)); %creates the indexing structure array for each crop in T for the h itteration of the loop
    B = table(zeros(r,1),zeros(r,1),'VariableNames',{'YEAR','YIELD'}); %creates a temporary table used to alter the table's format before being stowed in avgYield
    
    for i = 1:7 %1:length(climateZone) %for each climate zone
       
        z = cell2mat(agClimateZone(i)); %creates a vector, z, containing all of the regions numbers that are in the climate zone i
       
        %step below creates a temporary table, tempYield, based on data from the crop
        %identified in S. The table has three columns, labeled YEAR,
        %REGION, and YIELD. It was easier for me to work with table than
        %cells for the next steps.
               
        tempAvgYield = cell2table(subsref(Yield,S),'VariableNames',{'YEAR','REGION','YIELD'});
        %the below if statement alster the table to only contain data from 
        %the region i. In this loop i == 7 is the statewhide average, so we
        %don't alter tempYield.
        
        if i ~= 7
            %climate zone in kansas are, in ever case, defined by 3 climate
            %region markers. therfore, as long as i is not 7, z will be a
            %1x3 matrix of the regional identifiers fro different climate
            %zones.
            tempAvgYield = tempAvgYield(tempAvgYield.REGION == z(1) | tempAvgYield.REGION == z(2) | tempAvgYield.REGION == z(3),:); %
        end
        
      
        for j = y1:y2 %for the maximum range of years at any station defined by y1 as minimum year and y2 as maximum year

             B.YEAR(1-y1+j) = j; %creates a year column for temp table B where the year is the range of years for climate data                         
             B.YIELD(1-y1+j,i) = mean(tempAvgYield.YIELD(tempAvgYield.YEAR == j)); %where the years of j an tempAvgYield are equal, set the cell at that position equal to the average yeild of all crops during that year.                 

        end    
               
    end
    B(isnan(B.YIELD(:,1)),:) = []; %removes all rows where there are
    %NaN values in the YIELD field in any column
    avgYield = subsasgn(avgYield,S,table2cell(B));
end

%avgYield(isnan(avgYield.YIELD(:,1)),:) = [];
%% 4.0.2.1 Sorting data and Pettitt test
%this script will need at least three loops. First we will have an outer
clc
%sorting loop
CP = table();
%fileName='annual_pettitt_sort_sample.xlsx';
%avgLIndex = readtable(fileName);
cropCP = struct('totWheat',[],'irrWheat',[],'noIrrWheat',[],'totCorn',[],'irrCorn',[],'noIrrCorn',[],'totSoy',[],'irrSoy',[],'noIrrSoy',[], 'totSorgh',[],'irrSorgh',[],'noIrrSorgh',[]); 
for h = 1:length(names) %for the number of crops described in names   
    %this will define what part of avgYield(and thus what crop) the loop is
    %running
    T = fieldnames(avgYield); %creates a cell array containing the names of all the crop s in the avgYield structure array
    S = struct('type','.','subs',T(h,1)); %creates the indexing structure array for each crop in T for the h itteration of the loop
    
    for j= 1:length(climateZone) %for the number of climate zones in Kansas
        %creates a temporary variable to more easily sort by index value.
        CP.Region = ["west", "EWcentral", "east", "north", "NScentral", "south", "whole"]'; %sets the 
        
        B = cell2table(subsref(avgYield,S),'VariableNames',{'YEAR','Index'});           
       

%         if sum(ismember("IndexRank" ,B.Properties.VariableNames)) ==1 || sum(ismember("U" ,B.Properties.VariableNames)) ==1 
%             B = removevars(B, ["U", "IndexRank"]);
%         end
        B.Index = B.Index(:,j);
        %create a rank column based on year in ascending order
        for i = 1:height(B)          
            B.YearRank(i) = i;              
        end
        %sort temp table B by ascending index magnitude
        B = sortrows(B,2);
        %store the altered order of the Year Rank in avg""Index

      
        %avgYield.IndexRank(:,j) = B.YearRank;
        

        %Pettitt loop
        for k = 1:height(B) %for the number of rows in avgLIndex
            B.U(k)=2*sum(B.YearRank(1:k))-k*(height(B)+1); %calculates the rank statistic
            %storing rank statistic                      
            if abs(B.U(k))==max(abs(B.U)) %determines the postion of the max rank statistic
                c=k;
            end
        end

        %determines the year the change point occured at
        CP.YEAR(j) = c+min(B.YEAR)-1;
        K=max(abs(B.U)); %calculates the statistical change point test statistic
        K_c=(-log(0.05)*((height(B)^3)+(height(B)^2))/6)^.5; %calculating the critical value at the give station
        CP.ChangePoint(j) = K; %stores the U value as the change point
        if K>K_c %if K is greater than K_c, than the change point is statistically significant P<0.05 and we calculate 1-P
            CP.CL(j)=1;
            CP.P(j) = 1 - exp((-6*K^2)/(height(B)^3+height(B)^2));
        else

            CP.CL(j) = 0;
            CP.P(j) = 1 - exp((-6*K^2)/(height(B)^3+height(B)^2)); %still calculating 1-P just in case results are significant for P<0.10
        end     
      
        %avgYield.U(:,j) =B.U;
        cropCP = subsasgn(cropCP,S,table2cell(CP)); %store results of CP for this iteration of the loop in the structure array cropCP
        %reset the values of B and K to reduce errors
        B = table();
        K = 0;
    end
end
% 4.0.2.2 Climite calculations for Pettitt
%useLIndex = 0;
CP = table();
for j= 1:length(climateZone)
    %creates a temporary variable to more easily sort by index value.
    CP.Region = ["west", "EWcentral", "east", "north", "NScentral", "south", "whole"]'; %sets the 
    
    if useLIndex == 1
        B = renamevars(avgLIndex,"LIndex","Index");

    else
        B = renamevars(avgHIndex,"HIndex","Index");
    end
   


    if sum(ismember("IndexRank" ,B.Properties.VariableNames)) ==1 || sum(ismember("U" ,B.Properties.VariableNames)) ==1 
        B = removevars(B, ["U", "IndexRank"]);
    end
    B.Index = B.Index(:,j);
    %create a rank column based on year in ascending order
    for i = 1:height(B)
       B.YearRank(i) = i; 
    end
    %sort temp table B by ascending index magnitude
    B = sortrows(B,2);
    %store the altered order of the Year Rank in avg""Index


    if useLIndex == 1
        avgLIndex.IndexRank(:,j) = B.YearRank;
    else
        avgHIndex.IndexRank(:,j) = B.YearRank;
    end
   
    %Pettitt loop
    for k = 1:height(B) %for the number of rows in avgLIndex
        B.U(k)=2*sum(B.YearRank(1:k))-k*(height(B)+1); %calculates the rank statistic
        %storing rank statistic                      
        if abs(B.U(k))==max(abs(B.U)) %determines the postion of the max rank statistic
            c=k;
        end
    end

    %determines the year the change point occured at
    CP.YEAR(j) = c+min(B.YEAR)-1;
    K=max(abs(B.U)); %calculates the statistical change point test statistic
    K_c=(-log(0.05)*((height(B)^3)+(height(B)^2))/6)^.5; %calculating the critical value at the give station
    CP.ChangePoint(j) = K;
    if K>K_c
        CP.CL(j)=1;
        CP.P(j) = 1 - exp((-6*K^2)/(height(B)^3+height(B)^2));
    else

        CP.CL(j) = 0;
        CP.P(j) = 1 - exp((-6*K^2)/(height(B)^3+height(B)^2));
    end     

    if useLIndex == 1
        avgLIndex.U(:,j)= B.U;
        climateLCP = CP;
    else
        avgHIndex.U(:,j) = B.U;
        climateHCP = CP;
    end
    
    B = table();
    K = 0;
end

        
%% 4.0.2 Graphing time series vs. Yearly Average Grain Yield
tic
locations = ["western", "east-west central", "eastern", "northern", "north-south central", "southern", "statewide"];

for h = 1:length(names)
    
    T = fieldnames(avgYield); %creates a cell array containing the names of all the crop s in the avgYield structure array
    S = struct('type','.','subs',T(h,1)); %creates the indexing structure array for each crop in T for the h itteration of the loop
    
    for j = 1:2 %for the two different climate gradients (North to South and East to West)
        if j == 1 %j==1 represents the EW gradient
            fig = figure('Name', compose("TIndex vs " + string(cell2mat(T(h,1)))+"Yield, EW"),'Position',[10 50 1200 1000]); %name the figure "Tindex vs (name of crop) Yield, EW"
            climZone = [1 2 3 7]; %defines the climate regions that are EW, as well as the averge (7)
        else %repeat of above, but for NS           
            climZone = [4 5 6 7];
            fig = figure('Name', compose("TIndex vs " + string(cell2mat(T(h,1)))+"Yield, NS"),'Position',[10 50 1200 1000]);
        end
        
        %useLIndex = 1;
        % cropClimateStats = table(); %variables for staring statistical data from
        % the following graphs. This isn't complete yet.
        % cropClimateStats.Region = locations';
        for i = climZone %for the climate zone, either east to west or north to south
            B = cell2table(subsref(cropCP,S));
            changeP1 = B.Var2(i); %this is the change point for the crop yield
            changeP2 = changeP1; %this is the location of the change point for the climate, use the below script for the actual climate change point
        %     if useLIndex == 1
        %         changeP2 = climateLCP.YEAR(i);
        %     else
        %         changeP2 = climateHCP.YEAR(i);
        %     end
        
            if i ~= 7 %for all regions but the avager
                if j == 1
                    subplot(2,2,i) %the subplot position for the climate regions 
                else
                    subplot(2,2,i-3)
                end
            else
                subplot(2,2,4) %the subplot position for the statewide average
            end
            
            A = cell2table(subsref(avgYield,S),'VariableNames',{'YEAR','YIELD'}); %creates a temporary table for the crop yield stored in the avgYield structure array
            
            yyaxis left
            b = bar(A.YEAR,A.YIELD(:,i),1);
            %Adjusts ylimit to be the same for all graphs of a given crop
            %type
            if h == 1 || h==2 || h==3
                ylim([0 70]);
            end

            if h == 4||h==5||h==6
                ylim([0 220])
            end

            if h == 7||h==8||h==9
                ylim([0 70]);
            end

            if h == 10||h==11||h==12
                ylim([0 125]);
            end
            
            hold on

            mdlC = fitlm(A.YEAR,A.YIELD(:,i));
            zC = plot(mdlC);

            mdlC1 = fitlm(A.YEAR(A.YEAR<=changeP1),A.YIELD(A.YEAR<=changeP1,i));
            zC1 = plot(mdlC1);

            mdlC2 = fitlm(A.YEAR(A.YEAR>changeP1),A.YIELD(A.YEAR>changeP1,i));
            zC2 = plot(mdlC2);

            zC(1).Color = 'none'; 
            zC(1).Marker = '.';
            zC(1).MarkerSize = 10;
            zC(2).Color = 'm';
            zC(2).LineStyle = '--';
            zC(2).LineWidth = 1.15;
            zC(3).Color = 'none';
            zC(4).Color = 'none'; 

            zC1(1).Color = 'none'; 
            zC1(1).Marker = '.';
            zC1(1).MarkerSize = 10;
            zC1(2).Color = 'm';
            zC1(2).LineWidth = 1.15;
            zC1(3).Color = 'none';
            zC1(4).Color = 'none'; 

            zC2(1).Color = 'none'; 
            zC2(1).Marker = '.';
            zC2(1).MarkerSize = 10;
            zC2(2).Color = 'm';
            zC2(2).LineWidth = 1.15;
            zC2(3).Color = 'none';
            zC2(4).Color = 'none'; 

            ylabel(compose('Average' + string(cell2mat(T(h,1))) + ' Yield (Bu/acre)'))
           
            yyaxis right
            
            if useLIndex ==  1
                p = plot(avgLIndex.YEAR, avgLIndex.LIndex(:,i));
                ylim([57 68])

            else
                p = plot(avgHIndex.YEAR, avgHIndex.HIndex(:,i));
                ylim([80 94])
            end

            xlim([min(A.YEAR)-1 max(A.YEAR)+1])%xlim([1925 2015])
            hold on

            p.LineWidth = 1.15;
            p.Color = 'r';
            if j  == 1 %only for the first loop
                if useLIndex == 1 %modeling the linear trend for L index across the whole period, before the change point, and after the change point
                    mdlD = fitlm(avgLIndex.YEAR,avgLIndex.LIndex(:,i));                
                    mdlD1 = fitlm(avgLIndex.YEAR(avgLIndex.YEAR<=changeP2),avgLIndex.LIndex(avgLIndex.YEAR<=changeP2,i));                 
                    mdlD2 = fitlm(avgLIndex.YEAR(avgLIndex.YEAR>changeP2),avgLIndex.LIndex(avgLIndex.YEAR>changeP2,i));   

                else %modeling the linear trend for H index across the whole period, before the change point, and after the change point
                    mdlD = fitlm(avgHIndex.YEAR,avgHIndex.HIndex(:,i));                    
                    mdlD1 = fitlm(avgHIndex.YEAR(avgHIndex.YEAR<=changeP2),avgHIndex.HIndex(avgHIndex.YEAR<=changeP2,i));                  
                    mdlD2 = fitlm(avgHIndex.YEAR(avgHIndex.YEAR>changeP2),avgHIndex.HIndex(avgHIndex.YEAR>changeP2,i));   

                end
            end
            
            zD = plot(mdlD);
            zD1 = plot(mdlD1);
            zD2 = plot(mdlD2);
            

            zD(1).Color = 'none'; 
            zD(1).Marker = '.';
            zD(1).MarkerSize = 10;
            zD(2).Color = 'k';
            zD(2).LineStyle = '--';
            zD(2).LineWidth = 1.15;
            zD(3).Color = 'none';
            zD(4).Color = 'none'; 

            zD1(1).Color = 'none'; 
            zD1(1).Marker = '.';
            zD1(1).MarkerSize = 10;
            zD1(2).Color = 'k';
            zD1(2).LineWidth = 1.15;
            zD1(3).Color = 'none';
            zD1(4).Color = 'none';

            zD2(1).Color = 'none'; 
            zD2(1).Marker = '.';
            zD2(1).MarkerSize = 10;
            zD2(2).Color = 'k';
            zD2(2).LineWidth = 1.15;
            zD2(3).Color = 'none';
            zD2(4).Color = 'none';

            if useLIndex == 1
                ylabel("Average L - Index (days/year)")
                title(compose(locations(i) + " average L-Index " + string(cell2mat(T(h,1)))+ " Yield vs. year"))
%                 if i == 1 || 4
%                     lgd = legend([p b zD(2) zC(2)],{'L-Index','Grain Yield','Trend of H-Index','Trend of Grain Yield'});
%                 end
            else
                ylabel("Average H - Index (days/year)")
                title(compose(locations(i) +" average H-Index " + string(cell2mat(T(h,1)))+ " Yield vs. year"))
%                 if i == 1 || 4
%                     lgd = legend([p b zD(2) zC(2)],{'L-Index','Grain Yield','Trend of H-Index','Trend of Grain Yield'});
%                 end
            end 
            legend off
            xlabel('Year')
        end
        
        
        %defines the current date
        t = datetime('now');
        t.Format = 'dd.MM.yyyy';
        t = string(t);       
        
        %creates a legend and defines the name of the file that will be
        %saved to location defined later by fname
        if useLIndex == 1
            lgd = legend([p b zD(2) zC(2)],{'L-Index','Grain Yield','Trend of L-Index','Trend of Grain Yield'});
            if j == 1
                fileName = compose(t + "_LIndex_" +  string(cell2mat(T(h,1)))+"EW_countyYield_subplot");
            else
                fileName = compose(t + "_LIndex_" +  string(cell2mat(T(h,1)))+"NS_countyYield_subplot");
            end
        else
            lgd = legend([p b zD(2) zC(2)],{'H-Index','Grain Yield','Trend of H-Index','Trend of Grain Yield'});
            if j == 1
                fileName = compose(t + "_HIndex_" +  string(cell2mat(T(h,1)))+"EW_countyYield_subplot.tiff");
            else
                fileName = compose(t + "_HIndex_" +  string(cell2mat(T(h,1)))+"NS_countyYield_subplot.tiff");
            end
            
        end
        lgd.Position(1) = 0.01;
        lgd.Position(2) = 0.47;
        
        
        
        
        if h == 1 || h==2 || h==3
            fname = 'F:\github\HIndex\Outputs\Current Figures\Crop Yield Results\Wheat';
        end
        
        if h == 4||h==5||h==6
            fname = 'F:\github\HIndex\Outputs\Current Figures\Crop Yield Results\Corn';
        end
        
        if h == 7||h==8||h==9
            fname = 'F:\github\HIndex\Outputs\Current Figures\Crop Yield Results\Soybean';
        end
        
        if h == 10||h==11||h==12
            fname = 'F:\github\HIndex\Outputs\Current Figures\Crop Yield Results\Sorghum';
        end
       
        saveas(fig,fullfile(fname,fileName))
        
        
        
    end
end

toc
%%
    figure('Name', "TIndex vs Yield, NS")
    for i = [4 5 6 7]  
        changeP1 = CP.YEAR(i);
        changeP2 = changeP1;
    %     if useLIndex == 1
    %         changeP2 = climateLCP.YEAR(i);
    %     else
    %         changeP2 = climateHCP.YEAR(i);
    %     end
        subplot(2,2,i-3)    

        yyaxis left
        b = bar(avgYield.YEAR,avgYield.YIELD(:,i),1);
        hold on
        ylim([0 60])

        mdlC = fitlm(avgYield.YEAR,avgYield.YIELD(:,i));
        zC = plot(mdlC);

        mdlC1 = fitlm(avgYield.YEAR(avgYield.YEAR<=changeP1),avgYield.YIELD(avgYield.YEAR<=changeP1,i));
        zC1 = plot(mdlC1);

        mdlC2 = fitlm(avgYield.YEAR(avgYield.YEAR>changeP1),avgYield.YIELD(avgYield.YEAR>changeP1,i));
        zC2 = plot(mdlC2);

        zC(1).Color = 'none'; 
        zC(1).Marker = '.';
        zC(1).MarkerSize = 10;
        zC(2).Color = 'm';
        zC(2).LineStyle = '--';
        zC(2).LineWidth = 1.15;
        zC(3).Color = 'none';
        zC(4).Color = 'none'; 

        zC1(1).Color = 'none'; 
        zC1(1).Marker = '.';
        zC1(1).MarkerSize = 10;
        zC1(2).Color = 'm';
        zC1(2).LineWidth = 1.15;
        zC1(3).Color = 'none';
        zC1(4).Color = 'none'; 

        zC2(1).Color = 'none'; 
        zC2(1).Marker = '.';
        zC2(1).MarkerSize = 10;
        zC2(2).Color = 'm';
        zC2(2).LineWidth = 1.15;
        zC2(3).Color = 'none';
        zC2(4).Color = 'none'; 
        ylabel('Average Wheat Yield (Bu/acre)')

        yyaxis right
         if useLIndex ==  1
             p = plot(avgLIndex.YEAR, avgLIndex.LIndex(:,i));
             ylim([57 68])
         else
             p = plot(avgHIndex.YEAR, avgHIndex.HIndex(:,i));
             ylim([80 94])
         end

        xlim([1925 2015])
        hold on

        p.LineWidth = 1.15;
        p.Color = 'r';
        if useLIndex == 1
            mdlD = fitlm(avgLIndex.YEAR,avgLIndex.LIndex(:,i));   
            zD = plot(mdlD);
            mdlD1 = fitlm(avgLIndex.YEAR(avgLIndex.YEAR<=changeP2),avgLIndex.LIndex(avgLIndex.YEAR<=changeP2,i));   
            zD1 = plot(mdlD1);
            mdlD2 = fitlm(avgLIndex.YEAR(avgLIndex.YEAR>changeP2),avgLIndex.LIndex(avgLIndex.YEAR>changeP2,i));   
            zD2 = plot(mdlD2);
        else
            mdlD = fitlm(avgHIndex.YEAR,avgHIndex.HIndex(:,i));   
            zD = plot(mdlD); 
             mdlD1 = fitlm(avgHIndex.YEAR(avgHIndex.YEAR<=changeP2),avgHIndex.HIndex(avgHIndex.YEAR<=changeP2,i));   
            zD1 = plot(mdlD1);
            mdlD2 = fitlm(avgHIndex.YEAR(avgHIndex.YEAR>changeP2),avgHIndex.HIndex(avgHIndex.YEAR>changeP2,i));   
            zD2 = plot(mdlD2);
        end

        zD(1).Color = 'none'; 
        zD(1).Marker = '.';
        zD(1).MarkerSize = 10;
        zD(2).Color = 'k';
        zD(2).LineStyle = "--";
        zD(2).LineWidth = 1.15;
        zD(3).Color = 'none';
        zD(4).Color = 'none'; 

        zD1(1).Color = 'none'; 
        zD1(1).Marker = '.';
        zD1(1).MarkerSize = 10;
        zD1(2).Color = 'k';
        zD1(2).LineWidth = 1.15;
        zD1(3).Color = 'none';
        zD1(4).Color = 'none';

        zD2(1).Color = 'none'; 
        zD2(1).Marker = '.';
        zD2(1).MarkerSize = 10;
        zD2(2).Color = 'k';
        zD2(2).LineWidth = 1.15;
        zD2(3).Color = 'none';
        zD2(4).Color = 'none'; 

        if useLIndex == 1
            ylabel("Average L - Index (days/year)")
            title(compose(locations(i) + " average L-Index and average wheat yield vs. year"))
        else
            ylabel("Average H - Index (days/year)")
            title(compose(locations(i) +" average H-Index and average wheat yield vs. year"))
        end 
        legend off
       
    end

%% 4.0.2 Graphing H-index decadal analysis
%creates a bar graph with the average H index for each decade between 1911
%and 2010
% useLIndex = 1;

hDecades = decades;
hDecades.MEAN(2,1) = mean(HIndex(HIndex(:,1) <= max(hDecades.Var1(2,:)) & HIndex(:,1) >= min(hDecades.Var1(2,:)),cell2mat(climateZone(1))), 'all');
%calculates the average H/L index for the entire period between 1911-2010
if useLIndex ==1
    base = mean(LIndex(LIndex(:,1) <= hDecades.Var1(10,10) & LIndex(:,1) >= hDecades.Var1(1,1),2:end),'all'); 
else
    base = mean(HIndex(HIndex(:,1) <= hDecades.Var1(10,10) & HIndex(:,1) >= hDecades.Var1(1,1),2:end),'all');
end

%tis loop calculates the average H/L index for each decade from 1911-2010
for j = 1:height(hDecades)
    
    for h = 1:7      
        if useLIndex == 1
            hDecades.MEAN(j,h) = mean(LIndex(LIndex(:,1) <= max(hDecades.Var1(j,:)) & LIndex(:,1) >= min(hDecades.Var1(j,:)),cell2mat(climateZone(h))), 'all');            
            hDecades.DIFF(j,h) = hDecades.MEAN(j,h)-base;
        else
            hDecades.MEAN(j,h) = mean(HIndex(HIndex(:,1) <= max(hDecades.Var1(j,:)) & HIndex(:,1) >= min(hDecades.Var1(j,:)),cell2mat(climateZone(h))), 'all');            
            hDecades.DIFF(j,h) = hDecades.MEAN(j,h)-base;
        end
        
        decadeNames (j) = compose(num2str(min(hDecades.Var1(j,:))-1)+"s");        
    end
end
figure('name','Hbar');
bar(categorical(decadeNames),hDecades.DIFF); %creates a bar graph showing the difference between the decadal average and the centenial average of 1911-2010

if useLIndex == 1
    title("Difference of L-index decadal and centenial average for 1911-2010")
    ylabel("Change in L-index (days/year)")
    legend("West","NE Central", "East", "North", "NS Central", "South", "Entire State", 'location', 'northwest')
else
    title("Difference of H-index decadal and centenial average for 1911-2010")
    ylabel("Change in H-index (days/year)")
    legend("West","NE Central", "East", "North", "NS Central", "South", "Entire State")
end
xlabel("Decade")

%% 4.1 graphing decadal analysis subplots
figure('name','Hbar');
for i = 1:2
    subplot(2,1,i)
    if i == 1
        %creates a bar graph showing the difference between the decadal average and the centenial average of 1911-2010
        bar(categorical(decadeNames),hDecades.DIFF(:,[1,2,3,7])); %graphs climate zones of west, EW central, and East        
        legend("West","NE Central", "East", "Entire State")
    else
        bar(categorical(decadeNames),hDecades.DIFF(:,[4,5,6,7])); %graphs climate zones for north, NS central, and South
        legend("North", "NS Central", "South", "Entire State")
    end
    
    if useLIndex == 1
        title("Difference of L-index decadal and centenial average for 1911-2010")
        ylabel("Change in L-index (days/year)")
    else
        title("Difference of H-index decadal and centenial average for 1911-2010")
        ylabel("Change in H-index (days/year)")
    end
    xlabel("Decade")
    %ylim([-20 20])
end
%%
figure('name','Hbar vs WBar');
for i = 1:2
    subplot(2,1,i)
    if i == 1
        %creates a bar graph showing the difference between the decadal average and the centenial average of 1911-2010
        yyaxis left
        bar(categorical(decadeNames),hDecades.MEAN(:,[1,2,3,7])); %graphs climate zones of west, EW central, and East        
        
        yyaxis right
        bar(categorical(decadeNames),wYDecades.MEAN(:,[1,2,3,7]));
        legend("West","NE Central", "East", "Entire State")
    else
        yyaxis left
        bar(categorical(decadeNames),hDecades.MEAN(:,[4,5,6,7])); %graphs climate zones for north, NS central, and South
        
        yyaxis right
         bar(categorical(decadeNames),wYDecades.MEAN(:,[4,5,6,7]));
        legend("North", "NS Central", "South", "Entire State")
    end
    
    if useLIndex == 1
        title("Difference of L-index decadal and centenial average for 1911-2010")
        ylabel("Change in L-index (days/year)")
    else
        title("Difference of H-index decadal and centenial average for 1911-2010")
        ylabel("Change in H-index (days/year)")
    end
    xlabel("Decade")
    %ylim([-20 20])
end
%% 5 Calculating C-index
%create a C-Index
clc
tic
stationLength = length(stationNames);
%stationLength = 1;
%calls the path of the current file directory
CIndex = zeros(2013-1890+1,23);
CIndex(:,1) = (1890:2013)';

for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path    
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    temporaryFile.TMIN = round(temporaryFile.TMIN,0);
    %creates an array from the starting year to the ending year of the stations available weather data
    YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
    temporaryCIndex = table(YEAR); %creates an column array for the years of the H-Indecies 
    temporaryCIndex.CIndex = zeros(height(temporaryCIndex),1);
    counter = 0;
    for j = temporaryCIndex.YEAR(1):temporaryCIndex.YEAR(end)%for each year at this station
        year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary matrix for the given year       
        currentTemp = min(year.TMIN); %records min temp for the given year        
        %Count the number of times where the daily temp is greater than or
        %equal to that temp
        counter = 0;
        while counter < 32-currentTemp %checks to see if the counter is smaller than the currentTemp. This is to make sure that the value is an H-Index value.
            counter = 0;
            for h = 1:height(year)% for days in this year
                if year.TMIN(h) <= currentTemp %If the value at row h and column TMaxColumn are greater than currentTemp
                    counter = counter + 1; %increase counter by 1
                end    
            end
            if counter < 32-currentTemp %if the counter is smaller than currentTemp then the H-index is not valid, so we reduce it by one and repeat the loop.
                currentTemp = currentTemp + 1;
            end
        end
        temporaryCIndex.CIndex(j-temporaryCIndex.YEAR(1)+1) = 32-currentTemp; %stores the max H-index value for each year in temporaryHIndex      
    end
    
    for j = 1:height(temporaryCIndex) %for the number of years at the current station
        for h = 1:length(CIndex) %for full array of years being analyzed
           if temporaryCIndex.YEAR(j) == CIndex(h,1) %Checks to make sure that the years are the same for the given station
              CIndex(h,i+1)=temporaryCIndex.CIndex(j); %if the years are the same, then the yearly value for the station is stored in a column specifically for that station
           end
        end
    end
end
timeCIndex = toc
%% 7 Graphing C-index time series
CResults = table;
CResults.NAME = tableStationNames';
%for now, I'm only controlling the year from the H-Index section
 startYear = 1981;
 stopYear = 2013;
 useStartYear = 0;
 useStopYear = 0;
figure('Name', 'C-Index')
%periods = 
for i = 1:length(tableStationNames) %for each station in station names 
    %this section creates an array of subplots where each station has its
    %own plot showing the data points in orange, the change in time in
    %blue, and the general trend in black.
    tempNames = split(tableStationNames(i), '_'); %assigns a variable to the stations name and cuts out any unnecissary labels used for organization purposes
    CResults.NAME(i) = tempNames(2,1);
    
    A = array2table(CIndex(CIndex(:,1+i) ~= 0,[ 1 1+i])); %for each station, create a new array with the non-zero values and their corresponding years
    A.Properties.VariableNames{'Var1'} = 'YEAR'; 
    A.Properties.VariableNames{'Var2'} = 'CIndex';  
    B = A(A.YEAR>=startYear & A.YEAR <=stopYear,:);
    
    subplot(4,6,i) %Creates a system of subplots in a 4x6 grid
    pos=get(gca,'Position');
    set(gca,'Position',[pos(1,1) pos(1,2) pos(1,3)+.03 pos(1,4)]) %subplot position
    
    lA = plot(A.YEAR,A.CIndex); %adds a line to the plot for additional clarity
    hold on %add each station to the same plot
    mdlA = fitlm(A, 'CIndex ~ YEAR'); %performs a linear regression for the Year and the H Index    
    zA = plot(mdlA); %plots the linear regression and data point
    
    xA = A.YEAR;
    yA = A.CIndex;
    
    meanB = B;
    meanB.CIndex(:) = mean(B.CIndex);
    zB = plot(meanB.YEAR, meanB.CIndex);
    zB.Color = 'r';  
      
    %below are changes to the colors and markers of the plot for additional
    %clarity I removed 95% error bars (z(3) and z(4) to make the graph less cluttered. We
    %can turn these on later if we want to visually analyze the error
    %margins.
    zA(1).Color = 'none';%'#D95319'; %sets data points to be orange
    zA(1).Marker = '.';
    zA(1).MarkerSize = 10;
    zA(2).Color = 'k';
    zA(2).LineWidth = 1;
    zA(3).Color = 'none';
    zA(4).Color = 'none';
    
    legend('off'); %hides the automatic legend generated by fitlm
    
 
    %adjusting xtick and label to only show on edges
    if i == 19 || i == 20 || i == 21 || i == 22 || i == 23 
        xlabel('Year', 'FontSize', 11)
        set(gca,'Xtick',[1921 1951 1981 2011]);
    else
        xlabel('')
        xticks('')
    end
    
    %adjusting xlim based on year selected
    if useStartYear == 1 && useStopYear == 0
        xlim([startYear 2014])
    elseif useStartYear == 1 && useStopYear == 1
        xlim([startYear (stopYear+1)])
    elseif useStartYear == 0 && useStopYear == 1
        xlim([1890 (stopYear+1)])
    else     
        xlim([1890 2014])
    end
    %adjusting y tick and label to only show on edges
    if i == 1 || i == 7 || i == 13 || i == 19
        ylabel('C - Index (days/year)', 'FontSize', 11)
        yticks('auto')
    else
        ylabel('')
        yticks('')
    end  
    ylim([10 30])
    
    %this section  creates a table statistically important values for both
    %the POR and the specified period
    
    %p = polyfit(x,y,1);
       %this code runs trends analysis on x and y as independant and depedant variables
    %it tests the hypothesis of no correlation against the alternative
    %hypothesis of a nonzero correlation. so if p value is smaller than 0.05,
    %we reject the hypothesis.
    %analysis for the whole period of record
    
    if CResults.NAME(i,1) == "Saint Francis"
        CResults.climateDivision(i) = 1;
    elseif CResults.NAME(i,1) ==  "Oberlin"  
        CResults.climateDivision(i) = 1;
    elseif CResults.NAME(i,1) ==   "Colby"
        CResults.climateDivision(i) = 1;      
    elseif CResults.NAME(i) == "Phillipsburg"
        CResults.climateDivision(i) = 2;
    elseif CResults.NAME(i) == "Minneapolis" 
        CResults.climateDivision(i) = 2;    
    elseif CResults.NAME(i) == "Horton" 
        CResults.climateDivision(i) = 3;
    elseif CResults.NAME(i) == "Atchison" 
        CResults.climateDivision(i) = 3;
    elseif CResults.NAME(i) == "Manhattan"
        CResults.climateDivision(i) = 3;    
    elseif CResults.NAME(i) == "Tribune"
        CResults.climateDivision(i) = 4;
    elseif CResults.NAME(i) == "Wakeeney"
        CResults.climateDivision(i) = 4;    
    elseif CResults.NAME(i) == "Hays"
        CResults.climateDivision(i) = 5;
    elseif CResults.NAME(i) == "McPherson"
        CResults.climateDivision(i) = 5;    
    elseif CResults.NAME(i) == "Ottawa"
        CResults.climateDivision(i) = 6;
    elseif CResults.NAME(i) == "Lakin" 
        CResults.climateDivision(i) = 7;
    elseif CResults.NAME(i) == "Elkhart" 
        CResults.climateDivision(i) = 7;
    elseif CResults.NAME(i) == "Ashland"
        CResults.climateDivision(i) = 7;    
    elseif CResults.NAME(i) == "Larned"
        CResults.climateDivision(i) = 8;
    elseif CResults.NAME(i) == "MedicineLodge"
        CResults.climateDivision(i) = 8;    
    elseif CResults.NAME(i) == "Winfield"
        CResults.climateDivision(i) = 9;             
    elseif CResults.NAME(i) == "Sedan" 
        CResults.climateDivision(i) = 9;         
    elseif CResults.NAME(i) == "Independence" 
        CResults.climateDivision(i) = 9;            
    elseif CResults.NAME(i) == "FortScott" 
        CResults.climateDivision(i) = 9;         
    elseif CResults.NAME(i) == "Columbus"
        CResults.climateDivision(i) = 9;         
    end
    
    [tau,pD1]=corr(xA,yA,'type','kendall'); %kendall method
    tau_pA(i,1:3)=[i,tau,pD1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
    [rho,p2]=corr(xA,yA,'type','spearman');%spearman method
    rho_pA(i,1:3)=[i,rho,p2];
    [r,p3]=corr(xA,yA);%pearson (linear) method
    r_pA(i,1:3)=[i,r,p3]; %pearson(Least square method) method corrcoef(x,y);    
    
    CResults.slopePOR(i) = round(table2array(mdlA.Coefficients(2,1)),3);%calls the slope given for the linear regression of the data using the fitlm function
    CResults.adjustedSlopePOR(i) = CResults.slopePOR(i)*100;
    CResults.rPOR(i) = r_pA(i,2);
    CResults.rSqrPOR(i) = mdlA.Rsquared.Ordinary; %calls the r^2 value from the fitlm function and inputs it into a new table
    CResults.rPValuePOR(i) = round(table2array(mdlA.Coefficients(2,4)),3); %calls the pValue from the fitlm function and inputs it into a new table
    CResults.rhoPOR(i) = rho_pA(i,2);
    CResults.rhoPValuePOR(i) = rho_pA(i,3);
    CResults.tauPOR(i) = tau_pA(i,2);
    CResults.tauPValuePOR(i)= tau_pA(i,3);
    CResults.minPOR(i) = min(A.CIndex);
    CResults.maxPOR(i) = max(A.CIndex);
    CResults.meanPOR(i) = round(mean(A.CIndex),1); %calculates the average H index for the station and adds that to a new table
    CResults.medianPOR(i) = median(A.CIndex);
    
    CResults.minSP(i) = min(B.CIndex);
    CResults.maxSP(i) = max(B.CIndex);
    CResults.meanSP(i) = round(mean(B.CIndex),1);
    CResults.medianSP(i) = median(B.CIndex);
    title(compose(tempNames(2,1)+"\n"+num2str(100*CResults.slopePOR(i),"%#.1f")),'FontSize', 11);
end
%use the below script when you want to automatically make tiff files for
%the given graphs. See H-Index for issues with this function
%set(gcf,'PaperPositionMode','auto') %set the print area same as paper
%print('-dtiff','-r600', 'C_Index_Subplots')

%% 7.0.1 Calculating time series vs. yearly average grain yield
%This section compares average grain yield to average H-index for every
%year from y1 to y2 
y1 = 1926;
y2 = 2007;
  

avgCIndex = table();
avgYield = table();
avgYield.YEAR = zeros(82,1);
avgYield.YIELD = zeros(82,7);

for i = 1:length(climateZone)

    z = cell2mat(agClimateZone(i));
    tempAvgYield = wYield;
    if i ~= 7
        tempAvgYield = tempAvgYield(tempAvgYield.REGION == z(1) | tempAvgYield.REGION == z(2) | tempAvgYield.REGION == z(3),:);
    end    
     
    for j = min(tempAvgYield.YEAR):max(tempAvgYield.YEAR)
        avgYield.YEAR(1-min(tempAvgYield.YEAR)+j) = j;
        avgYield.YIELD(1-min(tempAvgYield.YEAR)+j,i) = mean(tempAvgYield.YIELD(tempAvgYield.YEAR == j));
    end    


    avgCIndex.YEAR = CIndex(:,1);
    avgCIndex.CIndex(:,i) = mean(CIndex(:,cell2mat(climateZone(i))),2);       
  
end

avgCIndex = avgCIndex(avgCIndex.YEAR >= y1 & avgCIndex.YEAR <= y2,:);

% 7.0.2.1 Sorting data and Pettitt test C-Index
%this script will need at least three loops. First we will have an outer
clc
%sorting loop
CP = table();
%fileName='annual_pettitt_sort_sample.xlsx';
%avgLIndex = readtable(fileName);

cropClimateStats = table();
cropClimateStats.Region = locations';
cropClimateStats.cRPValue = zeros(7,3);
cropClimateStats.cSlope = zeros(7,3);
cropClimateStats.cCCPTauPValue = zeros(7,3);
cropClimateStats.cropRPValue = zeros(7,3);
cropClimateStats.cropSlope = zeros(7,3);
cropClimateStats.cropTauPValue = zeros(7,3);
cropClimateStats.cCropCPTauPValue = zeros(7,3);


for e = 1:2
    for j= 1:length(climateZone)
        %creates a temporary variable to more easily sort by index value.
        CP.Region = ["west", "EWcentral", "east", "north", "NScentral", "south", "whole"]';
        if e == 1          
            B = renamevars(avgCIndex,"CIndex","Index");            
        else
            B = renamevars(avgYield,"YIELD","Index");
        end
        
        if j > 1
            B = removevars(B, ["U", "IndexRank"]);
        end
        B.Index = B.Index(:,j);
        %create a rank column based on year in ascending order
        for i = 1:height(B)
           B.YearRank(i) = i; 
        end
        %sort temp table B by ascending index magnitude
        B = sortrows(B,2);
        %store the altered order of the Year Rank in avg""Index
        
        if e == 1         
            avgCIndex.IndexRank(:,j) = B.YearRank;
        else
            avgYield.IndexRank(:,j) = B.YearRank;
        end


        %Pettitt loop
        for k = 1:height(B) %for the number of rows in avgLIndex
            B.U(k)=2*sum(B.YearRank(1:k))-k*(height(B)+1); %calculates the rank statistic
            %storing rank statistic                      
            if abs(B.U(k))==max(abs(B.U)) %determines the postion of the max rank statistic
                c=k;
            end
        end
        
        %determines the year the change point occured at
        CP.YEAR(j) = c+min(B.YEAR)-1;
        K=max(abs(B.U)); %calculates the statistical change point test statistic
        K_c=(-log(0.05)*((height(B)^3)+(height(B)^2))/6)^.5; %calculating the critical value at the give station
        CP.ChangePoint(j) = K;
        if K>K_c
            CP.CL(j)=1;
            CP.P(j) = 1 - exp((-6*K^2)/(height(B)^3+height(B)^2));
        else

            CP.CL(j) = 0;
            CP.P(j) = 1 - exp((-6*K^2)/(height(B)^3+height(B)^2));
        end     
        
        if e ==1
            
            avgCIndex.U(:,j) = B.U;
            climateCCP = CP;
            %calculating statistical significance of slopes based on 
            [tau,pD]=corr(B.YEAR,B.Index,'type','kendall'); %kendall method
            tau_pD(j,1:3)=[j,tau,pD];
            [tau,pD1]=corr(B.YEAR(B.YEAR<=climateCCP.YEAR(j)),B.Index(B.YEAR<=climateCCP.YEAR(j)),'type','kendall'); %kendall method
            tau_pD1(j,1:3)=[j,tau,pD1];
            [tau,pD2]=corr(B.YEAR(B.YEAR>climateCCP.YEAR(j)),B.Index(B.YEAR>climateCCP.YEAR(j)),'type','kendall'); %kendall method
            tau_pD2(j,1:3)=[j,tau,pD2];
            
            cropClimateStats.cRPValue(j) = 1;
            
            cropClimateStats.cSlope(j) = 1;
            
            cropClimateStats.cCCPTauPValue(j,1) = pD;
            cropClimateStats.cCCPTauPValue(j,2) = pD1;
            cropClimateStats.cCCPTauPValue(j,3) = pD2;
            
            
          
        else
            avgYield.U(:,j) =B.U;
            cropCP = CP;
            
            [tau,pC]=corr(B.YEAR,B.Index,'type','kendall'); %kendall method
            tau_pC(j,1:3)=[j,tau,pC];
            [tau,pC1]=corr(B.YEAR(B.YEAR<=cropCP.YEAR(j)),B.Index(B.YEAR<=cropCP.YEAR(j)),'type','kendall'); %kendall method
            tau_pC1(j,1:3)=[j,tau,pC1];
            [tau,pC2]=corr(B.YEAR(B.YEAR>cropCP.YEAR(j)),B.Index(B.YEAR>cropCP.YEAR(j)),'type','kendall'); %kendall method
            tau_pC2(j,1:3)=[j,tau,pC2];
            
            [tau,pD]=corr(B.YEAR,B.Index,'type','kendall'); %kendall method
            tau_pD(j,1:3)=[j,tau,pD];
            [tau,pD1]=corr(B.YEAR(B.YEAR<=cropCP.YEAR(j)),avgCIndex.CIndex(avgCIndex.YEAR<=cropCP.YEAR(j),j),'type','kendall'); %kendall method
            tau_pD1(j,1:3)=[j,tau,pD1];
            [tau,pD2]=corr(B.YEAR(B.YEAR>cropCP.YEAR(j)),avgCIndex.CIndex(avgCIndex.YEAR>cropCP.YEAR(j),j),'type','kendall'); %kendall method
            tau_pD2(j,1:3)=[j,tau,pD2];
            
            cropClimateStats.cropTauPValue(j,1) = pC;
            cropClimateStats.cropTauPValue(j,2) = pC1;
            cropClimateStats.cropTauPValue(j,3) = pC2;
            
            cropClimateStats.cCropCPTauPValue(j,1) = pD;
            cropClimateStats.cCropCPTauPValue(j,2) = pD1;
            cropClimateStats.cCropCPTauPValue(j,3) = pD2;
        end
        
        B = table();
        
    end
    B = table();
end
%% 7.0.2 Graphing time series vs. Yearly Average Grain Yield
locations = ["western", "east-west central", "eastern", "northern", "north-south central", "southern", "statewide"];
figure('Name', "CIndex vs Yield, EW")


for i = [1 2 3 7]
    changeP1 = cropCP.YEAR(i);
    changeP2 = climateCCP.YEAR(i);
    %changeP2 = changeP1;
    
    if i ~= 7
        subplot(2,2,i)
    else
        subplot(2,2,4)
    end

    yyaxis left
    b = bar(avgYield.YEAR,avgYield.YIELD(:,i),1);
    hold on

    mdlC = fitlm(avgYield.YEAR,avgYield.YIELD(:,i));
    zC = plot(mdlC);
    
    mdlC1 = fitlm(avgYield.YEAR(avgYield.YEAR<=changeP1),avgYield.YIELD(avgYield.YEAR<=changeP1,i));
    zC1 = plot(mdlC1);
    
    mdlC2 = fitlm(avgYield.YEAR(avgYield.YEAR>changeP1),avgYield.YIELD(avgYield.YEAR>changeP1,i));
    zC2 = plot(mdlC2);

    zC(1).Color = 'none'; 
    zC(1).Marker = '.';
    zC(1).MarkerSize = 10;
    zC(2).Color = 'm';
    zC(2).LineStyle = '--';
    zC(2).LineWidth = 1.15;
    zC(3).Color = 'none';
    zC(4).Color = 'none'; 
    
    zC1(1).Color = 'none'; 
    zC1(1).Marker = '.';
    zC1(1).MarkerSize = 10;
    zC1(2).Color = 'm';
    zC1(2).LineWidth = 1.15;
    zC1(3).Color = 'none';
    zC1(4).Color = 'none'; 
    
    zC2(1).Color = 'none'; 
    zC2(1).Marker = '.';
    zC2(1).MarkerSize = 10;
    zC2(2).Color = 'm';
    zC2(2).LineWidth = 1.15;
    zC2(3).Color = 'none';
    zC2(4).Color = 'none'; 
    
    ylabel('Average Wheat Yield (Bu/acre)')
    ylabel('Average Wheat Yield (Bu/acre)')
        
    yyaxis right
    p = plot(avgCIndex.YEAR, avgCIndex.CIndex(:,i));
    ylim([10 30])
  
      
    xlim([1925 2010])
    hold on

    p.LineWidth = 1.15;
    p.Color = 'r';
 
    mdlD = fitlm(avgCIndex.YEAR,avgCIndex.CIndex(:,i));   
    zD = plot(mdlD);
    mdlD1 = fitlm(avgCIndex.YEAR(avgCIndex.YEAR<=changeP2),avgCIndex.CIndex(avgCIndex.YEAR<=changeP2,i));   
    zD1 = plot(mdlD1);
    mdlD2 = fitlm(avgCIndex.YEAR(avgCIndex.YEAR>changeP2),avgCIndex.CIndex(avgCIndex.YEAR>changeP2,i));   
    zD2 = plot(mdlD2); 

    zD(1).Color = 'none'; 
    zD(1).Marker = '.';
    zD(1).MarkerSize = 10;
    zD(2).Color = 'k';
    zD(2).LineStyle = '--';
    zD(2).LineWidth = 1.15;
    zD(3).Color = 'none';
    zD(4).Color = 'none'; 
    
    zD1(1).Color = 'none'; 
    zD1(1).Marker = '.';
    zD1(1).MarkerSize = 10;
    zD1(2).Color = 'k';
    zD1(2).LineWidth = 1.15;
    zD1(3).Color = 'none';
    zD1(4).Color = 'none';
    
    zD2(1).Color = 'none'; 
    zD2(1).Marker = '.';
    zD2(1).MarkerSize = 10;
    zD2(2).Color = 'k';
    zD2(2).LineWidth = 1.15;
    zD2(3).Color = 'none';
    zD2(4).Color = 'none';
 


    ylabel("Average C - Index (days/year)")
    title(compose(locations(i) +" average C-Index and average wheat yield vs. year"))
     
    legend off
    xlabel('Year')
    
end


figure('Name', "CIndex vs Yield, NS")
for i = [4 5 6 7]  
    changeP1 = cropCP.YEAR(i);
    changeP2 = climateCCP.YEAR(i);
    %changeP2 = changeP1;
    
    subplot(2,2,i-3)    
    
    yyaxis left
    b = bar(avgYield.YEAR,avgYield.YIELD(:,i),1);
    hold on
    ylim([0 60])
  
    mdlC = fitlm(avgYield.YEAR,avgYield.YIELD(:,i));
    zC = plot(mdlC);
    
    mdlC1 = fitlm(avgYield.YEAR(avgYield.YEAR<=changeP1),avgYield.YIELD(avgYield.YEAR<=changeP1,i));
    zC1 = plot(mdlC1);
    
    mdlC2 = fitlm(avgYield.YEAR(avgYield.YEAR>changeP1),avgYield.YIELD(avgYield.YEAR>changeP1,i));
    zC2 = plot(mdlC2);

    zC(1).Color = 'none'; 
    zC(1).Marker = '.';
    zC(1).MarkerSize = 10;
    zC(2).Color = 'm';
    zC(2).LineStyle = '--';
    zC(2).LineWidth = 1.15;
    zC(3).Color = 'none';
    zC(4).Color = 'none'; 
    
    zC1(1).Color = 'none'; 
    zC1(1).Marker = '.';
    zC1(1).MarkerSize = 10;
    zC1(2).Color = 'm';
    zC1(2).LineWidth = 1.15;
    zC1(3).Color = 'none';
    zC1(4).Color = 'none'; 
    
    zC2(1).Color = 'none'; 
    zC2(1).Marker = '.';
    zC2(1).MarkerSize = 10;
    zC2(2).Color = 'm';
    zC2(2).LineWidth = 1.15;
    zC2(3).Color = 'none';
    zC2(4).Color = 'none'; 
    ylabel('Average Wheat Yield (Bu/acre)')
    
    yyaxis right
     
    p = plot(avgCIndex.YEAR, avgCIndex.CIndex(:,i));
    ylim([10 30])
       
    xlim([1925 2010])
    hold on
    
    p.LineWidth = 1.15;
    p.Color = 'r';
   
    mdlD = fitlm(avgCIndex.YEAR,avgCIndex.CIndex(:,i));   
    zD = plot(mdlD);
    mdlD1 = fitlm(avgCIndex.YEAR(avgCIndex.YEAR<=changeP2),avgCIndex.CIndex(avgCIndex.YEAR<=changeP2,i));   
    zD1 = plot(mdlD1);
    mdlD2 = fitlm(avgCIndex.YEAR(avgCIndex.YEAR>changeP2),avgCIndex.CIndex(avgCIndex.YEAR>changeP2,i));   
    zD2 = plot(mdlD2); 

    zD(1).Color = 'none'; 
    zD(1).Marker = '.';
    zD(1).MarkerSize = 10;
    zD(2).Color = 'k';
    zD(2).LineStyle = '--';
    zD(2).LineWidth = 1.15;
    zD(3).Color = 'none';
    zD(4).Color = 'none'; 
    
    zD1(1).Color = 'none'; 
    zD1(1).Marker = '.';
    zD1(1).MarkerSize = 10;
    zD1(2).Color = 'k';
    zD1(2).LineWidth = 1.15;
    zD1(3).Color = 'none';
    zD1(4).Color = 'none';
    
    zD2(1).Color = 'none'; 
    zD2(1).Marker = '.';
    zD2(1).MarkerSize = 10;
    zD2(2).Color = 'k';
    zD2(2).LineWidth = 1.15;
    zD2(3).Color = 'none';
    zD2(4).Color = 'none'; 
       
    ylabel("Average C - Index (days/year)")
    title(compose(locations(i) +" average C-Index and average wheat yield vs. year"))
 
    legend off
    xlabel('Year')
end
%% 8 Graphing C-index decadal bar graph analysis
cDecades = decades;
cDecades.MEAN(2,1) = mean(CIndex(CIndex(:,1) <= max(cDecades.Var1(2,:)) & CIndex(:,1) >= min(cDecades.Var1(2,:)),cell2mat(climateZone(1))), 'all');
%calculates the average C index for the entire period between 1911-2010
base = mean(CIndex(CIndex(:,1) <= cDecades.Var1(10,10) & CIndex(:,1) >= cDecades.Var1(1,1),2:end),'all');

%this loop calculates the average H/L index for each decade from 1911-2010
for j = 1:height(cDecades)
    for h = 1:7      
        cDecades.MEAN(j,h) = mean(CIndex(CIndex(:,1) <= max(cDecades.Var1(j,:)) & CIndex(:,1) >= min(cDecades.Var1(j,:)),cell2mat(climateZone(h))), 'all');            
        cDecades.DIFF(j,h) = cDecades.MEAN(j,h)-base;
        decadeNames (j) = compose(num2str(min(cDecades.Var1(j,:))-1)+"s");        
    end
end
figure('name','Cbar');
bar(categorical(decadeNames),cDecades.DIFF); %creates a bar graph showing the difference between the decadal average and the centenial average of 1911-2010

title("Difference of C-index decadal and centenial average for 1911-2010")
ylabel("Change in C-index (days/year)")

xlabel("Decade")
legend("West","NE Central", "East", "North", "NS Central", "South", "Entire State")
%% 8.1 C-index decadal bar graph subplots
figure('name','Cbar');
for i = 1:2
    subplot(2,1,i)
    if i == 1
        %creates a bar graph showing the difference between the decadal average and the centenial average of 1911-2010
        bar(categorical(decadeNames),cDecades.DIFF(:,[1,2,3,7])); %graphs climate zones of west, EW central, and East        
        legend("West","NE Central", "East", "Entire State")
    else
        bar(categorical(decadeNames),cDecades.DIFF(:,[4,5,6,7])); %graphs climate zones for north, NS central, and South
        legend("North", "NS Central", "South", "Entire State")
    end
    
    title("Difference of C-index decadal and centenial average for 1911-2010")
    ylabel("Change in C-index (days/year)")
    xlabel("Decade")
end
%% 9 Calculationg W-index
%WIndex will be the Index of wetness. This code will look at each month and
%see if it got any rain at all. 
%start with precip max
clc
tic
%stationLength = 1;
%calls the path of the current file directory
stationLength = length(stationNames);
WIndex = zeros(2013-1890+1,23);
WIndex(:,1) = (1890:2013)';
WIndex2 = WIndex;
%useInches = 0;
for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path    
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    temporaryFile.RAIN = round(temporaryFile.RAIN,0); %rounds the result to be a whole number, since we don't have measure values for mm to the tenthousandths place.
        
    temporaryFile.RAIN2 = round(temporaryFile.RAIN/25.4,2); %converts precipitations values from mm to inches
    %creates an array from the starting year to the ending year of the stations available weather data
    YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
    temporaryWIndex = table(YEAR); %creates an column array for the years of the H-Indecies 
    temporaryWIndex.WIndex = zeros(height(temporaryWIndex),1);
    temporaryWIndex.WIndex2 = zeros(height(temporaryWIndex),1);
    counter = 0;
    counter2 = 0;
    for j = temporaryWIndex.YEAR(1):temporaryWIndex.YEAR(end)%for each year at this station
        year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary matrix for the given year       
        currentPrecip = max(year.RAIN); %records max precip for the given year        
        %Count the number of times where the daily precip is greater than or
        %equal to the maximum precip
        while counter < currentPrecip %checks to see if the counter is smaller than the currentTemp. This is to make sure that the value is an H-Index value.
            counter = 0;
            for h = 1:height(year)% for days in this year
                if year.RAIN(h) >= currentPrecip %If the value at row h and column TMaxColumn are greater than currentTemp
                   counter = counter + 1;
                end    
            end
            
            if counter < currentPrecip %if the counter is smaller than currentTemp then the H-index is not valid, so we reduce it by one and repeat the loop.              
               currentPrecip = currentPrecip - 1 ;             
            end
        end
        
        temporaryWIndex.WIndex(j-temporaryWIndex.YEAR(1)+1) = currentPrecip;        
        
        
        currentPrecip = max(year.RAIN2);
        while counter2 < currentPrecip %checks to see if the counter is smaller than the currentTemp. This is to make sure that the value is an H-Index value.
            counter2 = 0;
            for h = 1:height(year)% for days in this year
                if year.RAIN2(h) >= currentPrecip %If the value at row h and column TMaxColumn are greater than currentTemp
                   
                   counter2 = counter2 + 0.01;
              
                end    
            end
            
            if counter2 < currentPrecip %if the counter is smaller than currentTemp then the H-index is not valid, so we reduce it by one and repeat the loop.             
               currentPrecip = currentPrecip - 0.01;              
            end            
        end
        
        temporaryWIndex.WIndex2(j-temporaryWIndex.YEAR(1)+1) = currentPrecip*100; %stores the max H-index value for each year in temporaryHIndex

    end
    
    for j = 1:height(temporaryWIndex) %for the number of years at the current station
        for h = 1:length(WIndex) %for full array of years being analyzed
           if temporaryWIndex.YEAR(j) == WIndex(h,1) %Checks to make sure that the years are the same for the given station
              WIndex(h,i+1)=temporaryWIndex.WIndex(j); %if the years are the same, then the yearly value for the station is stored in a column specifically for that station
              WIndex2(h,i+1)=temporaryWIndex.WIndex2(j);
           end
        end
    end  
end

timeWIndex = toc;
%% 10 Graphing W-index time series
WResults = table;
WResults.NAME = tableStationNames';

% startYear = 1981;
% stopYear = 2010;
% useStartYear = 0;
% useStopYear = 0;
figure('Name', 'W-Index')


for i = 1:length(tableStationNames) %for each station in station names 
    %this section creates an array of subplots where each station has its
    %own plot showing the data points in orange, the change in time in
    %blue, and the general trend in black.
    tempNames = split(tableStationNames(i), '_'); %assigns a variable to the stations name and cuts out any unnecissary labels used for organization purposes
    WResults.NAME(i) = tempNames(2,1); 
  
    A = array2table(WIndex(WIndex(:,1+i) ~= 0,[ 1 1+i])); %for each station, create a new array with the non-zero values and their corresponding years    
    A.Properties.VariableNames{'Var1'} = 'YEAR'; 
    A.Properties.VariableNames{'Var2'} = 'WIndex'; %W Index in mm
    
    A2 = array2table(WIndex2(WIndex2(:,1+i) ~= 0,[ 1 1+i])); %for each station, create a new array with the non-zero values and their corresponding years   
    A2.Properties.VariableNames{'Var1'} = 'YEAR'; 
    A2.Properties.VariableNames{'Var2'} = 'WIndex'; %W Index in inches
    
    B = A(A.YEAR>=startYear & A.YEAR <=stopYear,:);
    B2 = A2(A2.YEAR>=startYear & A2.YEAR <=stopYear,:); 
    
    subplot(4,6,i)
    pos=get(gca,'Position');
    set(gca,'Position',[pos(1,1) pos(1,2) pos(1,3)+.03 pos(1,4)]) %subplot position
    
   
    lA = plot(A.YEAR,A.WIndex); %adds a line to the plot for the yearly W index in mm
    hold on %add each station to the same plot
    lA2 = plot(A2.YEAR,A2.WIndex); %adds a line to the plot for the yearly W index in inches
    hold on   
    
    mdlA = fitlm(A, 'WIndex ~ YEAR'); %performs a linear regression for the Year and the H Index    
    mdlA2 = fitlm(A2, 'WIndex ~ YEAR');
    
    zA = plot(mdlA); %plots the linear regression and data point    
    zA2 = plot(mdlA2);
    
    xA = A.YEAR;
    yA = A.WIndex;
    xA2 = A2.YEAR;
    yA2 = A2.WIndex;
    
    meanB = B;
    meanB.WIndex(:) = mean(B.WIndex);
    zB = plot(meanB.YEAR, meanB.WIndex);
    zB.Color = 'r'; 
    zB.LineWidth = 1.25;
    
    meanB2 = B2;
    meanB2.WIndex(:) = mean(B2.WIndex);
    zB2 = plot(meanB2.YEAR, meanB2.WIndex);
    zB2.Color = 'c';
    zB2.LineWidth = 1.25;
    %below are cWanges to the colors and markers of the plot for additional
    %clarity I removed 95% error bars (z(3) and z(4) to make the graph less cluttered. We
    %can turn these on later if we want to visually analyze the error
    %margins.
    zA(1).Color = 'none'; %point color changed to make them invisible for less busy graphs%'#D95319'; %sets data points to be orange
    zA(1).Marker = '.';
    zA(1).MarkerSize = 10;
    zA(2).Color = 'k';
    zA(2).LineWidth = 1;
    zA(3).Color = 'none';
    zA(4).Color = 'none'; 
    
    zA2(1).Color = 'none'; %point color changed to make them invisible for less busy graphs%'#D95319'; %sets data points to be orange
    zA2(1).Marker = '.';
    zA2(1).MarkerSize = 10;
    zA2(2).Color = 'k';
    zA2(2).LineWidth = 1;
    zA2(3).Color = 'none';
    zA2(4).Color = 'none';  
    
    legend('off'); %hides the automatic legend generated by fitlm
    
    %if statement to set axis on left and bottom edge of subplot matrix
    %instead of on each plot
    
    if i == 19 || i == 20 || i == 21 || i == 22 || i == 23 
        xlabel('Year', 'FontSize', 11)
        set(gca,'Xtick',[1921 1951 1981 2011]);
        %xticks(gca, 1890:2010, 4)
    else
        xlabel('')
        xticks('')
    end
    
    %if statement that allows greater control on the years examined
    if useStartYear == 1 && useStopYear == 0
        xlim([startYear 2014])
    elseif useStartYear == 1 && useStopYear == 1
        xlim([startYear (stopYear+1)])
    elseif useStartYear == 0 && useStopYear == 1
        xlim([1890 (stopYear+1)])
    else     
        xlim([1890 2014])
    end
    
    % %sets boundries for the y axis for all graphs to be equal
    
    ylim([5 46])
    
    %If statement to make sure that only far left subplots have axis labels
    %and ticks marks
    if i == 1 || i == 7 || i == 13 || i == 19
         
        ylabel('W - Index (days/year)', 'FontSize', 11)
    
        yticks('auto')
    else
        ylabel('')
        yticks('')
    end
    %this section  creates a table statistically important values for both
    %the POR and the specified period
    
    
    %this code runs trends analysis on x and y as independant and depedant variables
    %it tests the hypothesis of no correlation against the alternative
    %hypothesis of a nonzero correlation. so if p value is smaller than 0.05,
    %we reject the hypothesis.
    %analysis for the whole period of record
    [tau,pD1]=corr(xA,yA,'type','kendall'); %kendall method
    tau_pA(i,1:3)=[i,tau,pD1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
    [rho,p2]=corr(xA,yA,'type','spearman');%spearman method
    rho_pA(i,1:3)=[i,rho,p2];
    [r,p3]=corr(xA,yA);%pearson (linear) method
    r_pA(i,1:3)=[i,r,p3]; %pearson(Least square method) method corrcoef(x,y);
    
    [tau,pD1]=corr(xA2,yA2,'type','kendall'); %kendall method
    tau_pA2(i,1:3)=[i,tau,pD1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
    [rho,p2]=corr(xA2,yA2,'type','spearman');%spearman method
    rho_pA2(i,1:3)=[i,rho,p2];
    [r,p3]=corr(xA2,yA2);%pearson (linear) method
    r_pA2(i,1:3)=[i,r,p3]; %pearson(Least square method) method corrcoef(x,y);
    
    %Allocates a number to each station based on there Climate Division in
    %the state of Kansas
    
    if WResults.NAME(i,1) == "Saint Francis"
        WResults.climateDivision(i) = 1;
    elseif WResults.NAME(i,1) ==  "Oberlin"  
        WResults.climateDivision(i) = 1;
    elseif WResults.NAME(i,1) ==   "Colby"
        WResults.climateDivision(i) = 1;      
    elseif WResults.NAME(i) == "Phillipsburg"
        WResults.climateDivision(i) = 2;
    elseif WResults.NAME(i) == "Minneapolis" 
        WResults.climateDivision(i) = 2;    
    elseif WResults.NAME(i) == "Horton" 
        WResults.climateDivision(i) = 3;
    elseif WResults.NAME(i) == "Atchison" 
        WResults.climateDivision(i) = 3;
    elseif WResults.NAME(i) == "Manhattan"
        WResults.climateDivision(i) = 3;    
    elseif WResults.NAME(i) == "Tribune"
        WResults.climateDivision(i) = 4;
    elseif WResults.NAME(i) == "Wakeeney"
        WResults.climateDivision(i) = 4;    
    elseif WResults.NAME(i) == "Hays"
        WResults.climateDivision(i) = 5;
    elseif WResults.NAME(i) == "McPherson"
        WResults.climateDivision(i) = 5;    
    elseif WResults.NAME(i) == "Ottawa"
        WResults.climateDivision(i) = 6;
    elseif WResults.NAME(i) == "Lakin" 
        WResults.climateDivision(i) = 7;
    elseif WResults.NAME(i) == "Elkhart" 
        WResults.climateDivision(i) = 7;
    elseif WResults.NAME(i) == "Ashland"
        WResults.climateDivision(i) = 7;    
    elseif WResults.NAME(i) == "Larned"
        WResults.climateDivision(i) = 8;
    elseif WResults.NAME(i) == "MedicineLodge"
        WResults.climateDivision(i) = 8;    
    elseif WResults.NAME(i) == "Winfield"
        WResults.climateDivision(i) = 9;             
    elseif WResults.NAME(i) == "Sedan" 
        WResults.climateDivision(i) = 9;         
    elseif WResults.NAME(i) == "Independence" 
        WResults.climateDivision(i) = 9;            
    elseif WResults.NAME(i) == "FortScott" 
        WResults.climateDivision(i) = 9;         
    elseif WResults.NAME(i) == "Columbus"
        WResults.climateDivision(i) = 9;         
    end
    
   
    
    
    WResults.slopePOR(i) = round(table2array(mdlA.Coefficients(2,1)),3);%calls the slope given for the linear regression of the data using the fitlm function
    WResults.adjustedSlopePOR(i) = WResults.slopePOR(i)*100; %this converts the slope from days/year to days/century
    WResults.rPOR(i) = r_pA(i,2);
    WResults.rSqrPOR(i) = mdlA.Rsquared.Ordinary; %calls the r^2 value from the fitlm function and inputs it into a new table
    WResults.rPValuePOR(i) = round(table2array(mdlA.Coefficients(2,4)),3); %calls the pValue from the fitlm function and inputs it into a new table
    WResults.rhoPOR(i) = rho_pA(i,2);
    WResults.rhoPValuePOR(i) = rho_pA(i,3);
    WResults.tauPOR(i) = tau_pA(i,2);
    WResults.tauPValuePOR(i)= tau_pA(i,3);
    WResults.minPOR(i) = min(A.WIndex);
    WResults.maxPOR(i) = max(A.WIndex);
    WResults.meanPOR(i) = round(mean(A.WIndex),1); %calculates the average H index for the station and adds that to a new table
    WResults.medianPOR(i) = median(A.WIndex);
    
    WResults.minSP(i) = min(B.WIndex);
    WResults.maxSP(i) = max(B.WIndex);
    WResults.meanSP(i) = round(mean(B.WIndex),1);
    WResults.medianSP(i) = median(B.WIndex);
    
    WResults.slopePOR2(i) = round(table2array(mdlA2.Coefficients(2,1)),3);%calls the slope given for the linear regression of the data using the fitlm function
    WResults.adjustedSlopePOR2(i) = WResults.slopePOR(i)*100; %this converts the slope from days/year to days/century
    WResults.rPOR2(i) = r_pA2(i,2);
    WResults.rSqrPOR2(i) = mdlA2.Rsquared.Ordinary; %calls the r^2 value from the fitlm function and inputs it into a new table
    WResults.rPValuePOR2(i) = round(table2array(mdlA2.Coefficients(2,4)),3); %calls the pValue from the fitlm function and inputs it into a new table
    WResults.rhoPOR2(i) = rho_pA2(i,2);
    WResults.rhoPValuePOR2(i) = rho_pA2(i,3);
    WResults.tauPOR2(i) = tau_pA2(i,2);
    WResults.tauPValuePOR2(i)= tau_pA2(i,3);
    WResults.minPOR2(i) = min(A2.WIndex);
    WResults.maxPOR2(i) = max(A2.WIndex);
    WResults.meanPOR2(i) = round(mean(A2.WIndex),1); %calculates the average H index for the station and adds that to a new table
    WResults.medianPOR2(i) = median(A2.WIndex);
    
    WResults.minSP2(i) = min(B2.WIndex);
    WResults.maxSP2(i) = max(B2.WIndex);
    WResults.meanSP2(i) = round(mean(B2.WIndex),1);
    WResults.medianSP2(i) = median(B2.WIndex);
  
    title(compose(tempNames(2,1)+"\n"+num2str(WResults.slopePOR(i)*100+"\n"+num2str(WResults.slopePOR2(i)*100),"%#.1f")),'FontSize', 11);
end


%use the below script when you want to automatically make tiff files for
%the given graphs. There is an issue with this function in that it doesn't
%expand the window before saving the file, so it compresses subplots to an
%unreadible level. I'll need to find a way to fix this later if we still
%want to use it over manually saving graphs.
%set(gcf,'PaperPositionMode','auto') %set the print area same as paper
%print('-dtiff','-r600', '12.3.19_H_Index_Subplots')

%% 11.0.1 Calculating time series vs. yearly average grain yield
%This section compares average grain yield to average H-index for every
%year from y1 to y2 
y1 = 1926;
y2 = 2007;
  

avgWIndex = table();
avgYield = table();
avgYield.YEAR = zeros(82,1);
avgYield.YIELD = zeros(82,7);

for i = 1:length(climateZone)

    z = cell2mat(agClimateZone(i));
    tempAvgYield = wYield;
    if i ~= 7
        tempAvgYield = tempAvgYield(tempAvgYield.REGION == z(1) | tempAvgYield.REGION == z(2) | tempAvgYield.REGION == z(3),:);
    end    
     
    for j = min(tempAvgYield.YEAR):max(tempAvgYield.YEAR)
        avgYield.YEAR(1-min(tempAvgYield.YEAR)+j) = j;
        avgYield.YIELD(1-min(tempAvgYield.YEAR)+j,i) = mean(tempAvgYield.YIELD(tempAvgYield.YEAR == j));
    end    


    avgWIndex.YEAR = WIndex(:,1);
    avgWIndex.WIndex(:,i) = mean(WIndex(:,cell2mat(climateZone(i))),2);       
  
end

avgWIndex = avgWIndex(avgWIndex.YEAR >= y1 & avgWIndex.YEAR <= y2,:);

%% 11.0.2 Graphing time series vs. Yearly Average Grain Yield
locations = ["western", "east-west central", "eastern", "northern", "north-south central", "southern", "statewide"];
figure('Name', "WIndex vs Yield, EW")

changeP1 = 1968;
for i = [1 2 3 7]
    if i ~= 7
        subplot(2,2,i)
    else
        subplot(2,2,4)
    end

    yyaxis left
    b = bar(avgYield.YEAR,avgYield.YIELD(:,i),1);
    hold on

    mdlC = fitlm(avgYield.YEAR,avgYield.YIELD(:,i));
    zC = plot(mdlC);
    
    mdlC1 = fitlm(avgYield.YEAR(avgYield.YEAR<=changeP1),avgYield.YIELD(avgYield.YEAR<=changeP1,i));
    zC1 = plot(mdlC1);
    
    mdlC2 = fitlm(avgYield.YEAR(avgYield.YEAR>changeP1),avgYield.YIELD(avgYield.YEAR>changeP1,i));
    zC2 = plot(mdlC2);

    zC(1).Color = 'none'; 
    zC(1).Marker = '.';
    zC(1).MarkerSize = 10;
    zC(2).Color = 'c';
    zC(2).LineWidth = 1.15;
    zC(3).Color = 'none';
    zC(4).Color = 'none'; 
    
    zC1(1).Color = 'none'; 
    zC1(1).Marker = '.';
    zC1(1).MarkerSize = 10;
    zC1(2).Color = 'm';
    zC1(2).LineWidth = 1.15;
    zC1(3).Color = 'none';
    zC1(4).Color = 'none'; 
    
    zC2(1).Color = 'none'; 
    zC2(1).Marker = '.';
    zC2(1).MarkerSize = 10;
    zC2(2).Color = 'Y';
    zC2(2).LineWidth = 1.15;
    zC2(3).Color = 'none';
    zC2(4).Color = 'none'; 
    
    ylabel('Average Wheat Yield (Bu/acre)')
    ylabel('Average Wheat Yield (Bu/acre)')
        
    yyaxis right
    p = plot(avgWIndex.YEAR, avgWIndex.WIndex(:,i));
    ylim([8 21])
        
    xlim([1925 2010])
    hold on

    p.LineWidth = 1.15;
    p.Color = 'r';
 
    mdlD = fitlm(avgWIndex.YEAR,avgWIndex.WIndex(:,i));   
    zD = plot(mdlD); 

    zD(1).Color = 'none'; 
    zD(1).Marker = '.';
    zD(1).MarkerSize = 10;
    zD(2).Color = 'g';
    zD(2).LineWidth = 1.15;
    zD(3).Color = 'none';
    zD(4).Color = 'none'; 


    ylabel("Average W - Index (days/year)")
    title(compose(locations(i) +" average W-Index and average wheat yield vs. year"))
     
    legend off
    xlabel('Year')
end


figure('Name', "WIndex vs Yield, NS")
for i = [4 5 6 7]  
    
    subplot(2,2,i-3)    
    
    yyaxis left
    b = bar(avgYield.YEAR,avgYield.YIELD(:,i),1);
    hold on
    ylim([0 60])
  
    mdlC = fitlm(avgYield.YEAR,avgYield.YIELD(:,i));
    zC = plot(mdlC);
    
    mdlC1 = fitlm(avgYield.YEAR(avgYield.YEAR<=changeP1),avgYield.YIELD(avgYield.YEAR<=changeP1,i));
    zC1 = plot(mdlC1);
    
    mdlC2 = fitlm(avgYield.YEAR(avgYield.YEAR>changeP1),avgYield.YIELD(avgYield.YEAR>changeP1,i));
    zC2 = plot(mdlC2);

    zC(1).Color = 'none'; 
    zC(1).Marker = '.';
    zC(1).MarkerSize = 10;
    zC(2).Color = 'c';
    zC(2).LineWidth = 1.15;
    zC(3).Color = 'none';
    zC(4).Color = 'none'; 
    
    zC1(1).Color = 'none'; 
    zC1(1).Marker = '.';
    zC1(1).MarkerSize = 10;
    zC1(2).Color = 'm';
    zC1(2).LineWidth = 1.15;
    zC1(3).Color = 'none';
    zC1(4).Color = 'none'; 
    
    zC2(1).Color = 'none'; 
    zC2(1).Marker = '.';
    zC2(1).MarkerSize = 10;
    zC2(2).Color = 'Y';
    zC2(2).LineWidth = 1.15;
    zC2(3).Color = 'none';
    zC2(4).Color = 'none'; 
    ylabel('Average Wheat Yield (Bu/acre)')
    
    yyaxis right
     
    p = plot(avgWIndex.YEAR, avgWIndex.WIndex(:,i));
    ylim([8 21])
       
    xlim([1925 2010])
    hold on
    
    p.LineWidth = 1.15;
    p.Color = 'r';
   
    mdlD = fitlm(avgWIndex.YEAR,avgWIndex.WIndex(:,i));   
    zD = plot(mdlD); 
   

    zD(1).Color = 'none'; 
    zD(1).Marker = '.';
    zD(1).MarkerSize = 10;
    zD(2).Color = 'g';
    zD(2).LineWidth = 1.15;
    zD(3).Color = 'none';
    zD(4).Color = 'none'; 
   
    ylabel("Average W - Index (days/year)")
    title(compose(locations(i) +" average W-Index and average wheat yield vs. year"))
 
    legend off
    xlabel('Year')
end
%% 11 Graphing W-index decadal analysis
%creates a bar graph with the average H index for each decade between 1911
%and 2010
% useLIndex = 1;

wDecades = decades;
wDecades.MEAN(2,1) = mean(WIndex(WIndex(:,1) <= max(wDecades.Var1(2,:)) & WIndex(:,1) >= min(wDecades.Var1(2,:)),cell2mat(climateZone(1))), 'all');
%calculates the average H/L index for the entire period between 1911-2010
base = mean(WIndex(WIndex(:,1) <= wDecades.Var1(10,10) & WIndex(:,1) >= wDecades.Var1(1,1),2:end),'all');

%tis loop calculates the average H/L index for each decade from 1911-2010
for j = 1:height(wDecades)
    for h = 1:7      

        wDecades.MEAN(j,h) = mean(WIndex(WIndex(:,1) <= max(wDecades.Var1(j,:)) & WIndex(:,1) >= min(wDecades.Var1(j,:)),cell2mat(climateZone(h))), 'all');            
        wDecades.DIFF(j,h) = wDecades.MEAN(j,h)-base;       
        
        decadeNames (j) = compose(num2str(min(wDecades.Var1(j,:))-1)+"s");        
    end
end
figure('name','Wbar');
bar(categorical(decadeNames),wDecades.DIFF); %creates a bar graph showing the difference between the decadal average and the centenial average of 1911-2010

title("Difference of W-index decadal and centenial average for 1911-2010")
ylabel("Change in W-index (days/year)")
legend("West","NE Central", "East", "North", "NS Central", "South", "Entire State")

xlabel("Decade")

%% 11.1 graphing decadal analysis subplots
figure('name','Wbar');
for i = 1:2
    subplot(2,1,i)
    if i == 1
        %creates a bar graph showing the difference between the decadal average and the centenial average of 1911-2010
        bar(categorical(decadeNames),wDecades.DIFF(:,[1,2,3,7])); %graphs climate zones of west, EW central, and East        
        legend("West","NE Central", "East", "Entire State")
    else
        bar(categorical(decadeNames),wDecades.DIFF(:,[4,5,6,7])); %graphs climate zones for north, NS central, and South
        legend("North", "NS Central", "South", "Entire State")
    end
    
    if useLIndex == 1
        title("Difference of L-index decadal and centenial average for 1911-2010")
        ylabel("Change in L-index (days/year)")
    else
        title("Difference of W-index decadal and centenial average for 1911-2010")
        ylabel("Change in W-index (days/year)")
    end
    xlabel("Decade")
end

%%  12 Calculating D-index
%Writes a dry index, D-Index. This could be done by counting x number of
%periods that had x number of days with no rain. I don't think that will be
%very informative since prolonged dry periods(drought) are more likely to
%be influential than knowning that there are more, shorter dry periods.
clc
tic
%stationLength = 1;
%calls the path of the current file directory


stationLength = length(stationNames);
DIndex = zeros(2013-1890+1,23);
DIndex(:,1) = (1890:2013)';
DIndex2 = DIndex;

for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path    
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    
    temporaryFile.RAIN = round(temporaryFile.RAIN,0); %rounds the result to be a whole number, since we don't have measure values for mm to the tenthousandths place.
        
    temporaryFile.RAIN2 = round(temporaryFile.RAIN/25.4,2);
    
    %creates an array from the starting year to the ending year of the stations available weather data   
    YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
    
    temporaryDIndex = table(YEAR); %creates an column array for the years of the H-Indecies 
    temporaryDIndex.DIndex = zeros(height(temporaryDIndex),1);
    counter = 0;
    %dryness is the minimum precipitation threshold for a day to be considered 'wet'
    dryness = 1; 
    dryness2 = 0.04; %dryness limit in inches (converting there are 0.0394 inches per mm)
 
    for j = temporaryDIndex.YEAR(1):temporaryDIndex.YEAR(end)%for each year at this station
            year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary matrix for the given year       
            currentDryPeriod = 0;
            tempDryPeriod = 0;
            counter = 0;
            %for dryness in mm
            for D = 1:height(year) %this loop identifies the longest dry period for the given year               
                if year.RAIN(D) < dryness
                    tempDryPeriod = tempDryPeriod + 1;
                else
                     if tempDryPeriod > currentDryPeriod
                         currentDryPeriod = tempDryPeriod;
                     end
                     tempDryPeriod = 0;
                end          
            end
            
            while counter < currentDryPeriod
                counter = 0;
                tempDryPeriod = 0;
                for h = 1:height(year)% for days in this year
                    
                    if year.RAIN(h) < dryness
                        tempDryPeriod = tempDryPeriod + 1;
                    else
                        if tempDryPeriod >= currentDryPeriod
                            counter = counter + 1;
                        end
                        tempDryPeriod = 0;
                    end                    
                end
                if counter < currentDryPeriod
                    currentDryPeriod = currentDryPeriod-1;
                end
            end
            temporaryDIndex.DIndex(j-temporaryDIndex.YEAR(1)+1) = currentDryPeriod;       
    
            %for Dryness in inches
            currentDryPeriod = 0;
            tempDryPeriod = 0;
            counter2 = 0;
            for D = 1:height(year) %this loop identifies the longest dry period for the given year               
                if year.RAIN2(D) < dryness2
                    tempDryPeriod = tempDryPeriod + 1;
                else
                     if tempDryPeriod > currentDryPeriod
                         currentDryPeriod = tempDryPeriod;
                     end
                    tempDryPeriod = 0;
                end          
            end
            
            while counter2 < currentDryPeriod
                counter2 = 0;
                tempDryPeriod = 0;
                for h = 1:height(year)% for days in this year
                    
                    if year.RAIN2(h) < dryness2
                        tempDryPeriod = tempDryPeriod + 1;
                    else
                        if tempDryPeriod >= currentDryPeriod
                            counter2 = counter2 + 1;
                        end
                        tempDryPeriod = 0;
                    end                    
                end
                if counter2 < currentDryPeriod
                    currentDryPeriod = currentDryPeriod-1;
                end
            end
            temporaryDIndex.DIndex2(j-temporaryDIndex.YEAR(1)+1) = currentDryPeriod;                 
    end
    
    
    for j = 1:height(temporaryDIndex) %for the number of years at the current station
        for h = 1:length(DIndex) %for full array of years being analyzed
           if temporaryDIndex.YEAR(j) == DIndex(h,1) %Checks to make sure that the years are the same for the given station
              DIndex(h,i+1)=temporaryDIndex.DIndex(j); %if the years are the same, then the yearly value for the station is stored in a column specifically for that station
              DIndex2(h,i+1)=temporaryDIndex.DIndex2(j);
           end
        end
    end  
end
timeDIndex = toc;
%% 13 graphing D-index time series

DResults = table;
DResults.NAME = tableStationNames';

% startYear = 1981;
% stopYear = 2010;
% useStartYear = 0;
% useStopYear = 0;
figure('Name', 'D-Index')


for i = 1:length(tableStationNames) %for each station in station names 
    %this section creates an array of subplots where each station has its
    %own plot showing the data points in orange, the change in time in
    %blue, and the general trend in black.
    tempNames = split(tableStationNames(i), '_'); %assigns a variable to the stations name and cuts out any unnecissary labels used for organization purposes
    DResults.NAME(i) = tempNames(2,1); 
  
    A = array2table(DIndex(DIndex(:,1+i) ~= 0,[ 1 1+i])); %for each station, create a new array with the non-zero values and their corresponding years
    
    A.Properties.VariableNames{'Var1'} = 'YEAR'; 
    A.Properties.VariableNames{'Var2'} = 'DIndex';  
    B = A(A.YEAR>=startYear & A.YEAR <=stopYear,:);  
    
    subplot(4,6,i)
    pos=get(gca,'Position');
    set(gca,'Position',[pos(1,1) pos(1,2) pos(1,3)+.03 pos(1,4)]) %subplot position
    
   
    lA = plot(A.YEAR,A.DIndex); %adds a line to the plot for additional clarity
    hold on %add each station to the same plot    
    mdlA = fitlm(A, 'DIndex ~ YEAR'); %performs a linear regression for the Year and the H Index    
    zA = plot(mdlA); %plots the linear regression and data point    
    
    xA = A.YEAR;
    yA = A.DIndex;
    
    meanB = B;
    meanB.DIndex(:) = mean(B.DIndex);
    zB = plot(meanB.YEAR, meanB.DIndex);
    zB.Color = 'r';  
    %below are cWanges to the colors and markers of the plot for additional
    %clarity I removed 95% error bars (z(3) and z(4) to make the graph less cluttered. We
    %can turn these on later if we want to visually analyze the error
    %margins.
    zA(1).Color = 'none'; %point color changed to make them invisible for less busy graphs%'#D95319'; %sets data points to be orange
    zA(1).Marker = '.';
    zA(1).MarkerSize = 10;
    zA(2).Color = 'k';
    zA(2).LineWidth = 1;
    zA(3).Color = 'none';
    zA(4).Color = 'none';   
    
    legend('off'); %hides the automatic legend generated by fitlm
    
    %if statement to set axis on left and bottom edge of subplot matrix
    %instead of on each plot
    if i == 19 || i == 20 || i == 21 || i == 22 || i == 23 
        xlabel('Year', 'FontSize', 11)
        set(gca,'Xtick',[1921 1951 1981 2011]);
        %xticks(gca, 1890:2010, 4)
    else
        xlabel('')
        xticks('')
    end
    
    %if statement that allows greater control on the years examined
    if useStartYear == 1 && useStopYear == 0
        xlim([startYear 2014])
    elseif useStartYear == 1 && useStopYear == 1
        xlim([startYear (stopYear+1)])
    elseif useStartYear == 0 && useStopYear == 1
        xlim([1890 (stopYear+1)])
    else     
        xlim([1890 2014])
    end
    
    % %sets boundries for the y axis for all graphs to be equal
    if useInches == 1 
        ylim([15 46])
    else
        ylim([6 13])
    end
    
    %If statement to make sure that only far left subplots have axis labels
    %and ticks marks
    if i == 1 || i == 7 || i == 13 || i == 19
         
        ylabel('D - Index (days/year)', 'FontSize', 11)
    
        yticks('auto')
    else
        ylabel('')
        yticks('')
    end
    %this section  creates a table statistically important values for both
    %the POR and the specified period
    
    
    %this code runs trends analysis on x and y as independant and depedant variables
    %it tests the hypothesis of no correlation against the alternative
    %hypothesis of a nonzero correlation. so if p value is smaller than 0.05,
    %we reject the hypothesis.
    %analysis for the whole period of record
    [tau,pD1]=corr(xA,yA,'type','kendall'); %kendall method
    tau_pA(i,1:3)=[i,tau,pD1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
    [rho,p2]=corr(xA,yA,'type','spearman');%spearman method
    rho_pA(i,1:3)=[i,rho,p2];
    [r,p3]=corr(xA,yA);%pearson (linear) method
    r_pA(i,1:3)=[i,r,p3]; %pearson(Least square method) method corrcoef(x,y); 
    
    %Allocates a number to each station based on there Climate Division in
    %the state of Kansas
    
    if DResults.NAME(i,1) == "Saint Francis"
        DResults.climateDivision(i) = 1;
    elseif DResults.NAME(i,1) ==  "Oberlin"  
        DResults.climateDivision(i) = 1;
    elseif DResults.NAME(i,1) ==   "Colby"
        DResults.climateDivision(i) = 1;      
    elseif DResults.NAME(i) == "Phillipsburg"
        DResults.climateDivision(i) = 2;
    elseif DResults.NAME(i) == "Minneapolis" 
        DResults.climateDivision(i) = 2;    
    elseif DResults.NAME(i) == "Horton" 
        DResults.climateDivision(i) = 3;
    elseif DResults.NAME(i) == "Atchison" 
        DResults.climateDivision(i) = 3;
    elseif DResults.NAME(i) == "Manhattan"
        DResults.climateDivision(i) = 3;    
    elseif DResults.NAME(i) == "Tribune"
        DResults.climateDivision(i) = 4;
    elseif DResults.NAME(i) == "Wakeeney"
        DResults.climateDivision(i) = 4;    
    elseif DResults.NAME(i) == "Hays"
        DResults.climateDivision(i) = 5;
    elseif DResults.NAME(i) == "McPherson"
        DResults.climateDivision(i) = 5;    
    elseif DResults.NAME(i) == "Ottawa"
        DResults.climateDivision(i) = 6;
    elseif DResults.NAME(i) == "Lakin" 
        DResults.climateDivision(i) = 7;
    elseif DResults.NAME(i) == "Elkhart" 
        DResults.climateDivision(i) = 7;
    elseif DResults.NAME(i) == "Ashland"
        DResults.climateDivision(i) = 7;    
    elseif DResults.NAME(i) == "Larned"
        DResults.climateDivision(i) = 8;
    elseif DResults.NAME(i) == "MedicineLodge"
        DResults.climateDivision(i) = 8;    
    elseif DResults.NAME(i) == "Winfield"
        DResults.climateDivision(i) = 9;             
    elseif DResults.NAME(i) == "Sedan" 
        DResults.climateDivision(i) = 9;         
    elseif DResults.NAME(i) == "Independence" 
        DResults.climateDivision(i) = 9;            
    elseif DResults.NAME(i) == "FortScott" 
        DResults.climateDivision(i) = 9;         
    elseif DResults.NAME(i) == "Columbus"
        DResults.climateDivision(i) = 9;         
    end
   
    
    
    DResults.slopePOR(i) = round(table2array(mdlA.Coefficients(2,1)),3);%calls the slope given for the linear regression of the data using the fitlm function
    DResults.adjustedSlopePOR(i) = DResults.slopePOR(i)*100; %this converts the slope from days/year to days/century
    DResults.rPOR(i) = r_pA(i,2);
    DResults.rSqrPOR(i) = mdlA.Rsquared.Ordinary; %calls the r^2 value from the fitlm function and inputs it into a new table
    DResults.rPValuePOR(i) = round(table2array(mdlA.Coefficients(2,4)),3); %calls the pValue from the fitlm function and inputs it into a new table
    DResults.rhoPOR(i) = rho_pA(i,2);
    DResults.rhoPValuePOR(i) = rho_pA(i,3);
    DResults.tauPOR(i) = tau_pA(i,2);
    DResults.tauPValuePOR(i)= tau_pA(i,3);
    DResults.minPOR(i) = min(A.DIndex);
    DResults.maxPOR(i) = max(A.DIndex);
    DResults.meanPOR(i) = round(mean(A.DIndex),1); %calculates the average H index for the station and adds that to a new table
    DResults.medianPOR(i) = median(A.DIndex);
    
    DResults.minSP(i) = min(B.DIndex);
    DResults.maxSP(i) = max(B.DIndex);
    DResults.meanSP(i) = round(mean(B.DIndex),1);
    DResults.medianSP(i) = median(B.DIndex);    
  
    title(compose(tempNames(2,1)+"\n"+num2str(DResults.slopePOR(i)*100,"%#.1f")),'FontSize', 11);
end


%use the below script when you want to automatically make tiff files for
%the given graphs. There is an issue with this function in that it doesn't
%expand the window before saving the file, so it compresses subplots to an
%unreadible level. I'll need to find a way to fix this later if we still
%want to use it over manually saving graphs.
%set(gcf,'PaperPositionMode','auto') %set the print area same as paper
%print('-dtiff','-r600', '12.3.19_H_Index_Subplots')

%% 14.0.1 Calculating time series vs. yearly average grain yield
%This section compares average grain yield to average H-index for every
%year from y1 to y2 
y1 = 1926;
y2 = 2007;
  

avgDIndex = table();
avgYield = table();
avgYield.YEAR = zeros(82,1);
avgYield.YIELD = zeros(82,7);

for i = 1:length(climateZone)

    z = cell2mat(agClimateZone(i));
    tempAvgYield = wYield;
    if i ~= 7
        tempAvgYield = tempAvgYield(tempAvgYield.REGION == z(1) | tempAvgYield.REGION == z(2) | tempAvgYield.REGION == z(3),:);
    end    
     
    for j = min(tempAvgYield.YEAR):max(tempAvgYield.YEAR)
        avgYield.YEAR(1-min(tempAvgYield.YEAR)+j) = j;
        avgYield.YIELD(1-min(tempAvgYield.YEAR)+j,i) = mean(tempAvgYield.YIELD(tempAvgYield.YEAR == j));
    end    


    avgDIndex.YEAR = DIndex(:,1);
    avgDIndex.DIndex(:,i) = mean(DIndex(:,cell2mat(climateZone(i))),2);       
  
end

avgDIndex = avgDIndex(avgDIndex.YEAR >= y1 & avgDIndex.YEAR <= y2,:);

%% 14.0.2 Graphing time series vs. Yearly Average Grain Yield
locations = ["western", "east-west central", "eastern", "northern", "north-south central", "southern", "statewide"];
figure('Name', "CIndex vs Yield, EW")

changeP1 = 1968;
for i = [1 2 3 7]
    if i ~= 7
        subplot(2,2,i)
    else
        subplot(2,2,4)
    end

    yyaxis left
    b = bar(avgYield.YEAR,avgYield.YIELD(:,i),1);
    hold on

    mdlC = fitlm(avgYield.YEAR,avgYield.YIELD(:,i));
    zC = plot(mdlC);
    
    mdlC1 = fitlm(avgYield.YEAR(avgYield.YEAR<=changeP1),avgYield.YIELD(avgYield.YEAR<=changeP1,i));
    zC1 = plot(mdlC1);
    
    mdlC2 = fitlm(avgYield.YEAR(avgYield.YEAR>changeP1),avgYield.YIELD(avgYield.YEAR>changeP1,i));
    zC2 = plot(mdlC2);

    zC(1).Color = 'none'; 
    zC(1).Marker = '.';
    zC(1).MarkerSize = 10;
    zC(2).Color = 'c';
    zC(2).LineWidth = 1.15;
    zC(3).Color = 'none';
    zC(4).Color = 'none'; 
    
    zC1(1).Color = 'none'; 
    zC1(1).Marker = '.';
    zC1(1).MarkerSize = 10;
    zC1(2).Color = 'm';
    zC1(2).LineWidth = 1.15;
    zC1(3).Color = 'none';
    zC1(4).Color = 'none'; 
    
    zC2(1).Color = 'none'; 
    zC2(1).Marker = '.';
    zC2(1).MarkerSize = 10;
    zC2(2).Color = 'Y';
    zC2(2).LineWidth = 1.15;
    zC2(3).Color = 'none';
    zC2(4).Color = 'none'; 
    
    ylabel('Average Wheat Yield (Bu/acre)')
    ylabel('Average Wheat Yield (Bu/acre)')
        
    yyaxis right
    p = plot(avgDIndex.YEAR, avgDIndex.DIndex(:,i));
    ylim([7 12])
  
      
    xlim([1925 2010])
    hold on

    p.LineWidth = 1.15;
    p.Color = 'r';
 
    mdlD = fitlm(avgDIndex.YEAR,avgDIndex.DIndex(:,i));   
    zD = plot(mdlD); 

    zD(1).Color = 'none'; 
    zD(1).Marker = '.';
    zD(1).MarkerSize = 10;
    zD(2).Color = 'g';
    zD(2).LineWidth = 1.15;
    zD(3).Color = 'none';
    zD(4).Color = 'none'; 


    ylabel("Average D - Index (days/year)")
    title(compose(locations(i) +" average D-Index and average wheat yield vs. year"))
     
    legend off
    xlabel('Year')
end


figure('Name', "DIndex vs Yield, NS")
for i = [4 5 6 7]  
    
    subplot(2,2,i-3)    
    
    yyaxis left
    b = bar(avgYield.YEAR,avgYield.YIELD(:,i),1);
    hold on
    ylim([0 60])
  
    mdlC = fitlm(avgYield.YEAR,avgYield.YIELD(:,i));
    zC = plot(mdlC);
    
    mdlC1 = fitlm(avgYield.YEAR(avgYield.YEAR<=changeP1),avgYield.YIELD(avgYield.YEAR<=changeP1,i));
    zC1 = plot(mdlC1);
    
    mdlC2 = fitlm(avgYield.YEAR(avgYield.YEAR>changeP1),avgYield.YIELD(avgYield.YEAR>changeP1,i));
    zC2 = plot(mdlC2);

    zC(1).Color = 'none'; 
    zC(1).Marker = '.';
    zC(1).MarkerSize = 10;
    zC(2).Color = 'c';
    zC(2).LineWidth = 1.15;
    zC(3).Color = 'none';
    zC(4).Color = 'none'; 
    
    zC1(1).Color = 'none'; 
    zC1(1).Marker = '.';
    zC1(1).MarkerSize = 10;
    zC1(2).Color = 'm';
    zC1(2).LineWidth = 1.15;
    zC1(3).Color = 'none';
    zC1(4).Color = 'none'; 
    
    zC2(1).Color = 'none'; 
    zC2(1).Marker = '.';
    zC2(1).MarkerSize = 10;
    zC2(2).Color = 'Y';
    zC2(2).LineWidth = 1.15;
    zC2(3).Color = 'none';
    zC2(4).Color = 'none'; 
    ylabel('Average Wheat Yield (Bu/acre)')
    
    yyaxis right
     
    p = plot(avgDIndex.YEAR, avgDIndex.DIndex(:,i));
    ylim([7 12])
       
    xlim([1925 2010])
    hold on
    
    p.LineWidth = 1.15;
    p.Color = 'r';
   
    mdlD = fitlm(avgDIndex.YEAR,avgDIndex.DIndex(:,i));   
    zD = plot(mdlD); 
   

    zD(1).Color = 'none'; 
    zD(1).Marker = '.';
    zD(1).MarkerSize = 10;
    zD(2).Color = 'g';
    zD(2).LineWidth = 1.15;
    zD(3).Color = 'none';
    zD(4).Color = 'none'; 
   
    ylabel("Average D - Index (days/year)")
    title(compose(locations(i) +" average D-Index and average wheat yield vs. year"))
 
    legend off
    xlabel('Year')
end
%% 14 Graphing D-index decadal analysis
%creates a bar graph with the average H index for each decade between 1911
%and 2010
% useLIndex = 1;

dDecades = decades;
dDecades.MEAN(2,1) = mean(DIndex(DIndex(:,1) <= max(dDecades.Var1(2,:)) & DIndex(:,1) >= min(dDecades.Var1(2,:)),cell2mat(climateZone(1))), 'all');
%calculates the average H/L index for the entire period between 1911-2010
base = mean(DIndex(DIndex(:,1) <= dDecades.Var1(10,10) & DIndex(:,1) >= dDecades.Var1(1,1),2:end),'all');

%tis loop calculates the average H/L index for each decade from 1911-2010
for j = 1:height(dDecades)
    for h = 1:7      

        dDecades.MEAN(j,h) = mean(DIndex(DIndex(:,1) <= max(dDecades.Var1(j,:)) & DIndex(:,1) >= min(dDecades.Var1(j,:)),cell2mat(climateZone(h))), 'all');            
        dDecades.DIFF(j,h) = dDecades.MEAN(j,h)-base;       
        
        decadeNames (j) = compose(num2str(min(dDecades.Var1(j,:))-1)+"s");        
    end
end
figure('name','Dbar');
bar(categorical(decadeNames),dDecades.DIFF); %creates a bar graph showing the difference between the decadal average and the centenial average of 1911-2010

title("Difference of D-index decadal and centenial average for 1911-2010")
ylabel("Change in D-index (days/year)")
legend("West","NE Central", "East", "North", "NS Central", "South", "Entire State")

xlabel("Decade")

%% 14.1 graphing decadal analysis subplots
figure('name','Dbar');
for i = 1:2
    subplot(2,1,i)
    if i == 1
        %creates a bar graph showing the difference between the decadal average and the centenial average of 1911-2010
        bar(categorical(decadeNames),dDecades.DIFF(:,[1,2,3,7])); %graphs climate zones of west, EW central, and East        
        legend("West","NE Central", "East", "Entire State")
    else
        bar(categorical(decadeNames),dDecades.DIFF(:,[4,5,6,7])); %graphs climate zones for north, NS central, and South
        legend("North", "NS Central", "South", "Entire State")
    end
    
    if useLIndex == 1
        title("Difference of L-index decadal and centenial average for 1911-2010")
        ylabel("Change in L-index (days/year)")
    else
        title("Difference of D-index decadal and centenial average for 1911-2010")
        ylabel("Change in D-index (days/year)")
    end
    xlabel("Decade")
end