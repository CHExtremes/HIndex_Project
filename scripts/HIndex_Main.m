%%
%This code is designed to take files from a specific folder, if they are
%excel files they are converted into csv files. those csv files are then
%used to produce an H-index based on temperature for each station.

%INPUTS: 
%The inputs for this script are weather data (in .csv format). Our
%weather csv's files have columns labled MONTH, DAY, YEAR, TMAX, TMIN, RAIN,and
%SNOW (from left to right). This script also uses grain yields inputs (sourced from USDA-NASS) on a
%monthly time-scale at the county level. Crop yield columns have headers of Program, Year, Period, Week Ending,	
%Geo Level,	State,	State ANSI,	Ag District, Ag District Code, County, County ANSI, Zip Code, Region,
%watershed_code, Watershed,	Commodity,	Data Item,	Domain,	Domain
%Category, Value,	CV (%). However, we only use the Year, Ag District
%Code, and Value column in this script. 

%To ensure the script runs, make sure you have these data downloaded. 
%Update "folderName" to match the name of the folder you used 
%to store the weather data and "folderName2" to match the name of the
%folder you used to store the crop data. 

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

%% 3 Calculating All indecies test script
%The following section calculates h index(an index where the number of 
%occurance is equal to the magnitude of the occurance based climate 
%indices for temperature (H, L, C) and precipitaiton (W(mm),W(in), D). All
%of the sub sections follow a simlar format, so most of the comments
%describing the scripts function will be in H-index with any relavant updates put in
%subsequent index scripts.

%This section takes acround 1000 seconds to run to completion.

clc
tic %start calculating run time

stationLength = length(stationNames); %variable for the number stations

%NOTE: I know this section is clunky. That's because it's is an amalgomation of the index scripts rather than a single
%loop that calculates all the scripts. This is because it takes about 10
%times longer to run each index when they are grouped in a single loop, for some reason.

%H and L index are measures of high maximum and high minimum temperature,
%respectively. They correspond to a number of days in a year with a
%temperature that is equal to or greater than that number of days. H-index
%in Kansas ranges from 70-90 generally and L index ranges from 50-65
%generally. 
for L = 1:2
    %for if you want to calculate the L-index(using Tmin) instead of the H-index (using Tmax)
    useLIndex = 0;
    if L == 1
        HIndex = zeros(2013-1890+1,24); %pre-allocating space for the maximum number of stations and years
        HIndex(:,1) = (1890:2013)'; %sets the first column as the maximum array of years
    end
    
    if L == 2
        LIndex = zeros(2013-1890+1,24); %pre-allocating space for the maximum number of stations and years
        LIndex(:,1) = (1890:2013)'; %sets the first column as the maximum array of years  
        useLIndex = 1;
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

        
        temporaryHIndex = table(transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR)),'VariableNames',"YEAR"); %creates an column array for the years of the H-Indecies 
        temporaryHIndex.HIndex = zeros(height(temporaryHIndex),1); %pre-allocating a column for the H/L index of the current station
        counter = 0;
        for j = temporaryHIndex.YEAR(1):temporaryHIndex.YEAR(end)%for each year at this station
            year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary table for the given year        

            %records max temp for the given year, the script will then work
            %backwards from here to find the maximum H/L index
            if useLIndex == 1
                currentTemp = max(year.TMIN); 
            else
                currentTemp = max(year.TMAX);
            end

            %This while loop counts the number of times where the daily temp is greater than or
            %equal to currentTemp. H/L index occurs when the number of days
            %for which the daily temp is greater than the current temp is
            %equal to the value of that current temp. 
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
end

timeHindex = toc; %record the run time as timeHindex

%C-Index is a measure of couldness. It measures of the difference between a temperature and
%32F(0C) and the number of times that difference occurs during
%that given year. CIndex ranges from 15-30, meaning there are 15-30 days in
%a given year with a temperature that is 15-30 degress fahrenheit below
%32F.

clc
tic

CIndex = zeros(2013-1890+1,23);
CIndex(:,1) = (1890:2013)';

for i = 1:stationLength %for each station
    baseFileName = stationNames(i); 
    fullFileName = fullfile(folder, baseFileName);    
    temporaryFile = readtable(fullFileName); 
    temporaryFile.TMIN = round(temporaryFile.TMIN,0); %C-Index is based on Tmin
    
    
    
    temporaryCIndex = table(transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR)),'VariableNames',"YEAR"); 
    temporaryCIndex.CIndex = zeros(height(temporaryCIndex),1);
    counter = 0;
    
    for j = temporaryCIndex.YEAR(1):temporaryCIndex.YEAR(end)
        year = temporaryFile(temporaryFile.YEAR==j,:);      
        currentTemp = min(year.TMIN); %records min temp for the given year 
        
        
        %this while loop counts the number of times where the daily temp is less than or
        %equal to the current temperature and continues until that counter
        %is greater than 32 minus the current temperature.
        counter = 0;
        while counter < 32-currentTemp 
            counter = 0;
            for h = 1:height(year)
                if year.TMIN(h) <= currentTemp 
                    counter = counter + 1; 
                end    
            end
            if counter < 32-currentTemp %if the counter is smaller than 32-currentTemp then the C-index is not valid, so we increase the current temperature it by one and repeat the loop.
                currentTemp = currentTemp + 1;
            end
        end
        %stores C Index for given year
        temporaryCIndex.CIndex(j-temporaryCIndex.YEAR(1)+1) = 32-currentTemp;     
    end
    %stores the results of CIndex for the given station
    for j = 1:height(temporaryCIndex) 
        for h = 1:length(CIndex) 
           if temporaryCIndex.YEAR(j) == CIndex(h,1) 
              CIndex(h,i+1)=temporaryCIndex.CIndex(j); 
           end
        end
    end
end
timeCIndex = toc
%WIndex will be the Index of wetness. This code examines the number of days
%were the volume of rainfall is the same as that number of days for a given
%year. W-index in Kansas ranges from 9-17mm typically. The second varible,
%W2, measures W-index in inches rather than milimeters. W2 in Kansas ranges from
%15-32, meaning that there are 15-32 days in a year with 0.15-0.32inches of
%rainfall.

clc
tic

stationLength = length(stationNames);
WIndex = zeros(2013-1890+1,23);
WIndex(:,1) = (1890:2013)';
WIndex2 = WIndex;

for i = 1:stationLength
    baseFileName = stationNames(i); 
    fullFileName = fullfile(folder, baseFileName); 
    temporaryFile = readtable(fullFileName); 
    temporaryFile.RAIN = round(temporaryFile.RAIN,0); %rounds the result to be a whole number, since we don't have measure values for mm to the tenthousandths place.
        
    temporaryFile.RAIN2 = round(temporaryFile.RAIN/25.4,2); %converts precipitations values from mm to inches
    %creates an array from the starting year to the ending year of the stations available weather data
  
    temporaryWIndex = table(transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR)),'VariableNames',"YEAR"); %creates an column array for the years of the H-Indecies 
    temporaryWIndex.WIndex = zeros(height(temporaryWIndex),1);
    temporaryWIndex.WIndex2 = zeros(height(temporaryWIndex),1);
    counter = 0;
    counter2 = 0;
    
    for j = temporaryWIndex.YEAR(1):temporaryWIndex.YEAR(end)
        year = temporaryFile(temporaryFile.YEAR==j,:);        
        currentPrecip = max(year.RAIN); 
        
        %Count the number of times where the daily precip is greater than or
        %equal to the maximum precip. If the counter is greater than the
        %currentPrecip, than the W-index is valid for that year, otherwise
        %the loop repeats
        %calculating W-index in mm
        while counter < currentPrecip 
            counter = 0;
            for h = 1:height(year)
                
                if year.RAIN(h) >= currentPrecip 
                   counter = counter + 1;
                end    
                
            end
            
            if counter < currentPrecip %if the counter is smaller than currentPrecip then the W-index is not valid, so we reduce currentPrecip by one and repeat the loop.              
               currentPrecip = currentPrecip - 1 ;             
            end
            
        end
        
        temporaryWIndex.WIndex(j-temporaryWIndex.YEAR(1)+1) = currentPrecip;        
        
        %Calculating W-Index in inches. Increments are in hundreths of an
        %inch rather than a single milimeter. This section takes the most
        %time in the whole section. 
        currentPrecip = max(year.RAIN2);
        while counter2 < currentPrecip 
            counter2 = 0;
            
            for h = 1:height(year)
                if year.RAIN2(h) >= currentPrecip                   
                   counter2 = counter2 + 0.01;
                end    
            end   
            
            if counter2 < currentPrecip           
               currentPrecip = currentPrecip - 0.01;              
            end            
        end        
        temporaryWIndex.WIndex2(j-temporaryWIndex.YEAR(1)+1) = currentPrecip*100;
        
    end
    
    for j = 1:height(temporaryWIndex) 
        for h = 1:length(WIndex)
           if temporaryWIndex.YEAR(j) == WIndex(h,1) 
              WIndex(h,i+1)=temporaryWIndex.WIndex(j);
              WIndex2(h,i+1)=temporaryWIndex.WIndex2(j);
           end
        end
    end 
    
end

timeWIndex = toc;
%Writes a dry index, D-Index. D-index is the number of periods, D, in a given
%year with D days that have precipitation less than some dryness threshold
%(currently we use dryness<1mm). D-index should be between 6-15 periods
%with 6-15, or more, dry days in a row. There is also a calculation for
%D-Index in inches, but it is the exact same as D-index in mm. 
clc
tic

DIndex = zeros(2013-1890+1,23);
DIndex(:,1) = (1890:2013)';
DIndex2 = DIndex;

for i = 1:stationLength 
    baseFileName = stationNames(i); 
    fullFileName = fullfile(folder, baseFileName);  
    temporaryFile = readtable(fullFileName);    
    temporaryFile.RAIN = round(temporaryFile.RAIN,0); %rounds the result to be a whole number, since we don't have measure values for mm to the tenthousandths place.    
    temporaryFile.RAIN2 = round(temporaryFile.RAIN/25.4,2);
       
    temporaryDIndex = table(transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR)),'VariableNames',"YEAR"); %creates an column array for the years of the H-Indecies 
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
            %this loop identifies the longest dry period for the given year        
            for D = 1:height(year)%for a given year        
                if year.RAIN(D) < dryness %if the precipitation at day D is less than the dryness threshold
                    tempDryPeriod = tempDryPeriod + 1; %increpent the temporary dry period
                else
                     if tempDryPeriod > currentDryPeriod %if the temporary dry period is greater than the current dry period
                         currentDryPeriod = tempDryPeriod; %change the current dry period to be equal to that temporary dry period
                     end
                     tempDryPeriod = 0; %reset the temporary dry period coutner and repeat the loop
                end          
            end
            
            %calculating D-Index using a threshold in mm
            while counter < currentDryPeriod
                counter = 0;
                tempDryPeriod = 0;
                for h = 1:height(year)
                    
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
            %maximum dry period based on a threshold in inches
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
            %caclulatind D-Index in inches
            while counter2 < currentDryPeriod
                counter2 = 0;
                tempDryPeriod = 0;
                for h = 1:height(year)
                    
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
    
    %Storing D-index for all statiosn
    for j = 1:height(temporaryDIndex) 
        for h = 1:length(DIndex) 
           if temporaryDIndex.YEAR(j) == DIndex(h,1) 
              DIndex(h,i+1)=temporaryDIndex.DIndex(j); 
              DIndex2(h,i+1)=temporaryDIndex.DIndex2(j);
           end
        end
    end  
end
timeDIndex = toc;
%storing all of the indices in a structure array to be able to more easily
%call them later
Index = struct('HIndex',HIndex,'LIndex',LIndex,'CIndex',CIndex,'WIndexmm',WIndex,'WIndexin',WIndex2,'DIndex',DIndex);

%% 4 graphing indices vs time
%creates a table to store all statistical values for all results for a
%given Index
T = fieldnames(Index); %creates a cell array with the names of the different indices
indexStats = struct('HIndex',[],'LIndex',[],'CIndex',[],'WIndexmm',[],'WIndexin',[],'DIndex',[],'varNames',[]); %creates an empty cell array to store the statistical results for the graphs of each index
limits = {[79 95],[55 72],[10 30],[5 26],[15 46],[6 13]}; %cell array containing the various y limits for the different indices
for k = 1:6
    indexName = cell2mat(T(k)); %gives the name of the current index as a character array
    S = struct('type','.','subs',char(T(k,1)));
    Results = table; %pre-define a table for the statistical values calculate later in this section
    Results.NAME = tableStationNames'; %set the rows of this table to be the sation names
    
    %variables to control period analyzed
    startYear = 1981;
    stopYear = 2010;
    useStartYear = 0;
    useStopYear = 0;

    %creates a new figure named H-Index. Really this should be done for both L
    %and  H index using an If statement
    figure('Name', indexName)


    for i = 1:length(tableStationNames) %for each station in station names 
        %this section creates an array of subplots where each station has its
        %own plot showing the data points in orange, the change in time in
        %blue, and the general trend in black.
        tempNames = split(tableStationNames(i), '_'); %assigns a variable to the stations name and cuts out any unnecissary labels used for organization purposes
        Results.NAME(i) = tempNames(2,1); %defines station name for the given row of the results table

        
        C = subsref(Index,S);                  
        %for each station, create a new array with the non-zero values and their corresponding years
        A = array2table(C(C(:,1+i) ~= 0,[ 1 1+i])); 
        
        %creat a new temporary table A to use in graphing
        A.Properties.VariableNames{'Var1'} = 'YEAR'; 
        A.Properties.VariableNames{'Var2'} = 'Index';  
        
        %create a new temporary table B that is table A from startYear to
        %stopYear
        B = A(A.YEAR>=startYear & A.YEAR <=stopYear,:); %creates a temporary table thats ranges from startYear to stopYear in length

        %creats a set of subplots with a width of 6 and a height of 4
        subplot(4,6,i)
        %reduces the white space between graphs for the current graph
        pos=get(gca,'Position');
        set(gca,'Position',[pos(1,1) pos(1,2) pos(1,3)+.03 pos(1,4)]) %subplot position


        lA = plot(A.YEAR,A.Index); %adds a line, lA, to the plot to show the change in H/L over time
        lA.LineWidth = .5; %sets the line width to the minimum of 0.5pt

        hold on %add each addition graph to the same subplot    
        mdlA = fitlm(A, 'Index ~ YEAR'); %performs a linear regression for the Year and the H Index    
        zA = plot(mdlA); %plots the linear regression and data point    

        xA = A.YEAR; %creates a tempory variable for graphing purposes
        yA = A.Index;

        meanB = B;
        meanB.Index(:) = mean(B.Index);
        %zB = plot(meanB.YEAR, meanB.Index);
        %zB.Color = 'r';  
        %zB.LineWidth = 1.15;
        %below are changes to the colors and markers of the plot for additional
        %clarity I removed 95% error bars (z(3) and z(4) to make the graph less cluttered. We
        %can turn these on later if we want to visually analyze the error
        %margins.
        zA(1).Color = 'none'; %point color changed to make them invisible for less busy graphs%'#D95319'; %sets data points to be orange
        zA(1).Marker = '.';
        zA(1).MarkerSize = 10;
        zA(2).Color = 'r';
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

        %if statement that allows greater control on the years examined,
        %this usually isn't used, but if you want to take a closer look at
        %certain years it helps to have this built in already
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
        
        ylim(cell2mat(limits(k)))
        
        %If statement to make sure that only far left subplots have axis labels
        %and ticks marks
        if i == 1 || i == 7 || i == 13 || i == 19                          
            ylabel(compose(string(indexName) +" (days/year)"), 'FontSize', 11)
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
        %the state of Kansas. I'm sure there is a faster and more elegant way
        %to do this, but the best I could come up with was a bunch of if
        %statements.

        if Results.NAME(i,1) == "Saint Francis" 
            Results.climateDivision(i) = 1;
        elseif Results.NAME(i,1) ==  "Oberlin"  
            Results.climateDivision(i) = 1;
        elseif Results.NAME(i,1) ==   "Colby"
            Results.climateDivision(i) = 1;      
        elseif Results.NAME(i) == "Phillipsburg"
            Results.climateDivision(i) = 2;
        elseif Results.NAME(i) == "Minneapolis" 
            Results.climateDivision(i) = 2;    
        elseif Results.NAME(i) == "Horton" 
            Results.climateDivision(i) = 3;
        elseif Results.NAME(i) == "Atchison" 
            Results.climateDivision(i) = 3;
        elseif Results.NAME(i) == "Manhattan"
            Results.climateDivision(i) = 3;    
        elseif Results.NAME(i) == "Tribune"
            Results.climateDivision(i) = 4;
        elseif Results.NAME(i) == "Wakeeney"
            Results.climateDivision(i) = 4;    
        elseif Results.NAME(i) == "Hays"
            Results.climateDivision(i) = 5;
        elseif Results.NAME(i) == "McPherson"
            Results.climateDivision(i) = 5;    
        elseif Results.NAME(i) == "Ottawa"
            Results.climateDivision(i) = 6;
        elseif Results.NAME(i) == "Lakin" 
            Results.climateDivision(i) = 7;
        elseif Results.NAME(i) == "Elkhart" 
            Results.climateDivision(i) = 7;
        elseif Results.NAME(i) == "Ashland"
            Results.climateDivision(i) = 7;    
        elseif Results.NAME(i) == "Larned"
            Results.climateDivision(i) = 8;
        elseif Results.NAME(i) == "MedicineLodge"
            Results.climateDivision(i) = 8;    
        elseif Results.NAME(i) == "Winfield"
            Results.climateDivision(i) = 9;             
        elseif Results.NAME(i) == "Sedan" 
            Results.climateDivision(i) = 9;         
        elseif Results.NAME(i) == "Independence" 
            Results.climateDivision(i) = 9;            
        elseif Results.NAME(i) == "FortScott" 
            Results.climateDivision(i) = 9;         
        elseif Results.NAME(i) == "Columbus"
            Results.climateDivision(i) = 9;         
        end

        %this section stores various statistical numbers in a dedicated
        %column for the Results table r is pearson correlation coefficient, rho is
        %spearman correlation, tau is for Mann-Kindall correlation. POR is
        %for values that are calculated over the entire period of record.
        Results.slopePOR(i) = round(table2array(mdlA.Coefficients(2,1)),3);%calls the slope given for the linear regression of the data using the fitlm function
        Results.adjustedSlopePOR(i) = Results.slopePOR(i)*100; %this converts the slope from days/year to days/century
        Results.rPOR(i) = r_pA(i,2);
        Results.rSqrPOR(i) = mdlA.Rsquared.Ordinary; %calls the r^2 value from the fitlm function and inputs it into a new table
        Results.rPValuePOR(i) = round(table2array(mdlA.Coefficients(2,4)),3); %calls the pValue from the fitlm function and inputs it into a new table
        Results.rhoPOR(i) = rho_pA(i,2);
        Results.rhoPValuePOR(i) = rho_pA(i,3);
        Results.tauPOR(i) = tau_pA(i,2);
        Results.tauPValuePOR(i)= tau_pA(i,3);
        Results.minPOR(i) = min(A.Index);
        Results.maxPOR(i) = max(A.Index);
        Results.meanPOR(i) = round(mean(A.Index),1); %calculates the average H index for the station and adds that to a new table
        Results.medianPOR(i) = median(A.Index);
        %plots a horizontal line for the average value over the entire
        %period of record. Color is magenta, style is dashed line, line
        %width slightly thicker than standard.
        l = yline(Results.meanPOR(i));
        l.Color = ('k');
        l.LineWidth = 1.25;
        l.LineStyle = '--';
        %stores results for a specific period (SP) defined in variable B by
        %the startYear and stopYear
        Results.minSP(i) = min(B.Index);
        Results.maxSP(i) = max(B.Index);
        Results.meanSP(i) = round(mean(B.Index),1);
        Results.medianSP(i) = median(B.Index);    
        %gives each subplot a title with the station name and the slope of
        %trend for the entire period of record.
        title(compose(tempNames(2,1)+"\n"+num2str(Results.slopePOR(i)*100,"%#.1f")),'FontSize', 11);
    end
    
    S2 = struct('type','.','subs',char(T(k,1))); %describes the current index
    indexStats = subsasgn(indexStats,S2,table2cell(Results)); %stores the table from Results, a cell of the indexStats table
    indexStats.varNames = Results.Properties.VariableNames; %stores the names of all the variables used in the Results table.
end
%Might create a save script for these graphs too

%use the below script when you want to automatically make tiff files for
%the given graphs. There is an issue with this function in that it doesn't
%expand the window before saving the file, so it compresses subplots to an
%unreadible level. I'll need to find a way to fix this later if we still
%want to use it over manually saving graphs.
%set(gcf,'PaperPositionMode','auto') %set the print area same as paper
%print('-dtiff','-r600', '12.3.19_H_Index_Subplots')
%% 5.1 Calculating time series vs. yearly average grain yield
%This section compares average grain yield to average H-index for every
%year from y1 to y2 
y1 = 1908; %smallest year for climate data
y2 = 2013; %largest year for climate data
r = 1+y2-y1; %the size of the range between y1 and y2 for creating row indices later

varNames = {'HIndex','LIndex','CIndex','WIndexmm','WIndexin','DIndex'}; %creating a cell array containing the names of all the indices
avgClimateIndex = table([1893:2013]','VariableNames',{'YEAR'}); %creates a year column ranging from year y1 to year y2

for j = 1:length(varNames) %for each index used
    S = struct('type','.','subs',varNames(j)); %structure variable for subreferencing. uses dot index and field name based on varNames at position j 
    B = subsref(Index,S); %creates temporary array B equal to the Index at with name J
    C = zeros(124,1); %pre-allocates the size of C to reduce errors
    for i = 1:length(climateZone) %for each aggregate climate region (E, W, N, S etc.)
        B(B==0) = NaN; %sets 0 values in the array equal to NaN for easier removal later and more accurate calculation of mean in the next line.
        C(:,i) = mean(B(:,cell2mat(climateZone(i))),2,'omitnan'); %all rows for column i of temporary array C are set equal to the mean of index j for all the stations in climate zone i      
    end 
    C(isnan(C(:,1)),:) = [];
    avgClimateIndex = subsasgn(avgClimateIndex,S,C); %set the variable(defined by S) in table avgClimateIndex equal to the temporary table C
    
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
%% 5.2.1 Sorting data and Pettitt test

clc
%creat CP as a temporary table to store information and cropCP as a
%permanent variable to store the results of each iterration of CP in the
%for loop.
CP = table();
cropCP = struct('totWheat',[],'irrWheat',[],'noIrrWheat',[],'totCorn',[],'irrCorn',[],'noIrrCorn',[],'totSoy',[],'irrSoy',[],'noIrrSoy',[], 'totSorgh',[],'irrSorgh',[],'noIrrSorgh',[]); 

for h = 1:length(names) %for the number of crops described in names   
    %this will define what part of avgYield(and thus what crop) the loop is
    %running
    T = fieldnames(avgYield); %creates a cell array containing the names of all the crop s in the avgYield structure array
    S = struct('type','.','subs',T(h,1)); %creates the indexing structure array for each crop in T for the h itteration of the loop
    
    for j= 1:length(climateZone) %for the number of climate zones in Kansas
        %creates a temporary variable to more easily sort by index value.
        CP.Region = ["west", "EWcentral", "east", "north", "NScentral", "south", "whole"]'; %sets the 
        
        B = cell2table(subsref(avgYield,S),'VariableNames',{'YEAR','Index'});   %creates temporary table B that is used to calculate the         
        B.Index = B.Index(:,j); %sets the Index variable of Table B equal to the climate region at position j
        
        %create a rank column based on year in ascending order
        for i = 1:height(B)          
            B.YearRank(i) = i;              
        end
        
        %sort temp table B by ascending index magnitude
        B = sortrows(B,2);
        %store the altered order of the Year Rank in avg""Index           

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
%% 5.2.2 Climite calculations for Pettitt

CP = table();
T = avgClimateIndex.Properties.VariableNames; %creates a cell array of the Index Names used in avgClimateIndex
climateCP = table();
for l = 2:width(avgClimateIndex) %for each index used
    S = struct('type','.','subs',char(T(l))); %define a structure array with type 'dot' and a field that is the current index at position l
    for j= 1:length(climateZone) %For each climate zone in Kansas
        %creates a temporary variable to more easily sort by index value.
        CP.Region = ["west", "EWcentral", "east", "north", "NScentral", "south", "whole"]'; %sets the Region variable of the CP table to be equal to the names of the different climate regions for Kansas
        B = table(avgClimateIndex.YEAR,subsref(avgClimateIndex,S),'VariableNames',{'YEAR','Index'}); %creates a temporary table, B, that is the same as the Index described by S in the avgClimateIndex table        

        B.Index = B.Index(:,j); %adjusts B index to only be equal to column j
       
        %create a rank column based on year in ascending order
        for i = 1:height(B)
           B.YearRank(i) = i; 
        end
        
        %sort temp table B by ascending index magnitude
        B = sortrows(B,2);

        %Pettitt loop
        for k = 1:height(B) %for the number of rows in avgLIndex
            B.U(k)=2*sum(B.YearRank(1:k))-k*(height(B)+1); %calculates the rank statistic
            %storing rank statistic                      
            if abs(B.U(k))==max(abs(B.U)) %determines the postion of the max rank statistic
                c=k; %c is equal to the row where U is a maxium, this is used later to calculate the year where the max u occurs
            end
        end

        %determines the year the change point occured at
        CP.YEAR(j) = c+min(B.YEAR)-1;
        K=max(abs(B.U)); %calculates the statistical change point test statistic
        K_c=(-log(0.05)*((height(B)^3)+(height(B)^2))/6)^.5; %calculating the critical value at the give station
        CP.ChangePoint(j) = K;
        if K>K_c %if K is larger than K_c, the change point is statistically significant for P<0.05
            CP.CL(j)=1;
            CP.P(j) = 1 - exp((-6*K^2)/(height(B)^3+height(B)^2)); % calculates the exact significane of 1-P
        else
            CP.CL(j) = 0;
            CP.P(j) = 1 - exp((-6*K^2)/(height(B)^3+height(B)^2)); %calculating 1-P incase a large number of points are significant for P<0.10
        end     

        climateCP = subsasgn(climateCP,S,table2array(CP)); %stores the change points for each region in variable named after the climate index described by S for the climateCP
        B = table(); %reset variables for next iteration of the loop
        K = 0;
    end
    
end

        
%% 5.2.3 Graphing time series vs. Yearly Average Grain Yield
%expected run time ~200 seconds
tic
locations = ["western", "east-west central", "eastern", "northern", "north-south central", "southern", "statewide"]; %names of each climate region
climLimits = {[79 95],[55 72],[10 30],[5 26],[15 46],[6 13]}; %y limits of climate indices order with position is H, L, C, W(mm), W(in), D
cropLimits = {[0 70], [0 220], [0 70], [0 125]};  %y limites of crop yields order with position is winter wheat, corn, soybeans, sorghum

cropClimateStats = table(); %variable for storing statistical data from
% the following graphs. This isn't complete yet.
cropClimateStats.Region = locations';

for k = 2:width(avgClimateIndex) %for each index
    T2 = avgClimateIndex.Properties.VariableNames; %T2 is a cell array of the variable names (climate indices) in the avgClimateIndex table
    S2 = struct('type','.','subs',T2(k));
        
    B = table(avgClimateIndex.YEAR, subsref(avgClimateIndex,S2),'VariableNames',{'YEAR','Index'});
             
           
    for h = 1:length(names) %for each crop 

        T = fieldnames(avgYield); %creates a cell array containing the names of all the crop s in the avgYield structure array
        S = struct('type','.','subs',T(h,1)); %creates the indexing structure array for each crop in T for the h itteration of the loop
        %Remove this, since changing to only a single subplot
       
        fig = figure('Name', compose(string(cell2mat(T2(k))) + " vs " + string(cell2mat(T(h,1)))+"Yield"),'Position',[10 50 1500 1000]); %name the figure "Tindex vs (name of crop) Yield, EW"
        
        
        for i = 1:length(locations) %for each climate zone
            E = cell2table(subsref(cropCP,S)); %temporary table 
            changeP1 = E.Var2(i); %this is the change point for the crop yield
            changeP2 = changeP1; %this is the year of the change point for the climate lines, this is currently set to be the same change point as the crop yield change point
            
            %Use thise scripts if you want the change point for the climate
            %line to be based on the climate change point rather than the
            %crop change point
            %indexCP = subsref(climateCP,S2);
%            changeP2 = str2num(index ylim(cell2mat(cropLimits(1)));CP(i,2));

            %places graphs for the E to W regions on the left and the N to
            %S regions on the right with the statewide average across the
            %bottom
            if i == 2
                subplot(4,2,3)
            elseif i == 3
                subplot(4,2,5)
            elseif i == 4
                subplot(4,2,2)
            elseif i == 5
                subplot(4,2,4)
            else
                subplot(4,2,i)
            end                  
                     
            A = cell2table(subsref(avgYield,S),'VariableNames',{'YEAR','YIELD'}); %creates a temporary table for the crop yield stored in the avgYield structure array
                                                             
            yyaxis left %graphs the yield results using the left hand y axis
            b = bar(A.YEAR,A.YIELD(:,i),1); %graphing the crop yields as bars
            %b.FaceColor = 	[0.165 0.042 0.042];%'#EDB120';
            %Adjusts ylimit to be the same for all graphs of a given crop
            %type
            if h == 1 || h==2 || h==3
                ylim(cell2mat(cropLimits(1)));
            end

            if h == 4||h==5||h==6
                 ylim(cell2mat(cropLimits(2)));
            end

            if h == 7||h==8||h==9
                 ylim(cell2mat(cropLimits(3)));
            end

            if h == 10||h==11||h==12
                 ylim(cell2mat(cropLimits(4)));
            end
            
            hold on
           
            %plotting linear trends for yield for entire period
            mdlC = fitlm(A.YEAR,A.YIELD(:,i));
            zC = plot(mdlC);
            %before change point
            mdlC1 = fitlm(A.YEAR(A.YEAR<=changeP1),A.YIELD(A.YEAR<=changeP1,i));
            zC1 = plot(mdlC1);
            %after change point
            mdlC2 = fitlm(A.YEAR(A.YEAR>changeP1),A.YIELD(A.YEAR>changeP1,i));
            zC2 = plot(mdlC2);
            
            %if k == 2 %only calculate these values once (during the first climate index iteration)
                %calculating and storing the P-Value for Mann-Kindal
                %non-parametric regression of crop data
                %calculating tau for the current crop over the entire period  
                [tau,pD]=corr(A.YEAR,A.YIELD(:,i),'type','kendall','rows','complete'); %kendall method
                D(i,1) = pD;
                %before the change point
                [tau,pD]=corr(A.YEAR(A.YEAR<=changeP1),A.YIELD(A.YEAR<=changeP1,i),'type','kendall','rows','complete'); %kendall method
                D(i,2) = pD;
                %after the change point
                [tau,pD]=corr(A.YEAR(A.YEAR>changeP1),A.YIELD(A.YEAR>changeP1,i),'type','kendall','rows','complete'); %kendall method
                D(i,3) = pD;

                %storing the slope and corresponding P-value for the linear models
                D(i,4) = table2array(mdlC.Coefficients(2,1));
                D(i,5) = table2array(mdlC.Coefficients(2,4));
                D(i,6) = table2array(mdlC1.Coefficients(2,1));
                D(i,7) = table2array(mdlC1.Coefficients(2,4));
                D(i,8) = table2array(mdlC2.Coefficients(2,1));
                D(i,9) = table2array(mdlC2.Coefficients(2,4));
            %end
            %adjusting color and line type for grain yield trend lines
            yieldColor = 	'c';
            
            zC(1).Color = 'none'; 
            zC(1).Marker = '.';
            zC(1).MarkerSize = 10;
            zC(2).Color = yieldColor;
            zC(2).LineStyle = '--';
            zC(2).LineWidth = 1.15;
            zC(3).Color = 'none';
            zC(4).Color = 'none'; 

            zC1(1).Color = 'none'; 
            zC1(1).Marker = '.';
            zC1(1).MarkerSize = 10;
            zC1(2).Color = yieldColor;
            zC1(2).LineWidth = 1.15;
            zC1(3).Color = 'none';
            zC1(4).Color = 'none'; 

            zC2(1).Color = 'none'; 
            zC2(1).Marker = '.';
            zC2(1).MarkerSize = 10;
            zC2(2).Color = yieldColor;
            zC2(2).LineWidth = 1.15;
            zC2(3).Color = 'none';
            zC2(4).Color = 'none'; 

            ylabel(compose('Average' + string(cell2mat(T(h,1))) + ' Yield (Bu/acre)')) %y label for the grain yield (left hand y axis)

            yyaxis right %set right hand y axis for climate 

            %creating a temporary table, B, for climate 
            B = table(avgClimateIndex.YEAR, subsref(avgClimateIndex,S2),'VariableNames',{'YEAR','Index'});
                        
            %temporary table for the climate index described by S2 from the avgClimateIndex table
            p = plot(B.YEAR, B.Index(:,i)); %plot year versus climate in region i, only 
         
            %I'll need to adjust the ylim based on the climate index being
            %used. I can probably copy this from the previous section of
            %script for ylimits
            
            ylim(cell2mat(climLimits(k-1))); %sets the y limits for the climate index (right hand y axis)
            xlim([min(A.YEAR)-1 max(A.YEAR)+1]) %set the x-axis limites to be the minimum and maximum years of available crop yield data
            hold on

            p.LineWidth = 1.15; %increases the thinkness of the climate trend line for greater visibility          

       
     
            %only calculate the climate trend lines during the first crop
            %iteration since it doesn't change until the climate index is
            %changed
            %if h == 1
                %the below statements create 3 linear models, 1 for the entire
                %period (mdlD), 1 for the period after the change point (mdlD1) and
                %1 after the change point (mdlD2)
                mdlD = fitlm(B.YEAR,B.Index(:,i));                
                mdlD1 = fitlm(B.YEAR(B.YEAR<=changeP2),B.Index(B.YEAR<=changeP2,i));                 
                mdlD2 = fitlm(B.YEAR(B.YEAR>changeP2),B.Index(B.YEAR>changeP2,i));            
                
                %calculating and storing the P-Value for Mann-Kindal
                %non-parametric regression of climate data
                [tau,pD]=corr(B.YEAR,B.Index(:,i),'type','kendall','rows','complete'); %kendall method
                F(i,1) = pD;
                %before the change point
                [tau,pD2]=corr(B.YEAR(B.YEAR<=changeP1),B.Index(B.YEAR<=changeP1,i),'type','kendall','rows','complete'); %kendall method
                F(i,2) = pD2;
                %after the change point
                [tau,pD3]=corr(B.YEAR(B.YEAR>changeP1),B.Index(B.YEAR>changeP1,i),'type','kendall','rows','complete'); %kendall method
                F(i,3) = pD3;
                
                %storing P-Value for linear climate models
                F(i,4) = table2array(mdlD.Coefficients(2,1));
                F(i,5) = table2array(mdlD.Coefficients(2,4));
                F(i,6) = table2array(mdlD1.Coefficients(2,1));
                F(i,7) = table2array(mdlD1.Coefficients(2,4));
                F(i,8) = table2array(mdlD2.Coefficients(2,1));
                F(i,9) = table2array(mdlD2.Coefficients(2,4));
            %end
            
            %plotting all of the trendlines for linear models created for
            %the climate indices
            zD = plot(mdlD);
            zD1 = plot(mdlD1);
            zD2 = plot(mdlD2);
            
            
            %adjusting the colors, size and line type of the index
            %trend lines
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
               
            ylabel(compose("Average " + string(cell2mat(T2(k))) + " (days/year)")) %labels the right hand axis with the name of the cliamte index defined by T2 at position k
            title(compose(locations(i) + " average " + string(cell2mat(T2(k)))+" and "+ string(cell2mat(T(h,1)))+ " Yield vs. year")) %gives each subplot a title with the climate region, climate index used, and crop used
            
            legend off %removes the legends automatically made by fitlm graphs
            xlabel('Year')
                    
        end
        %stores the statistical results for each crop in D 
        cropClimateStats = subsasgn(cropClimateStats,S,D);
        %stores the statistical results for each Index in F. The POR values
        %don't change with the change in crop, but the Change point does.
        %so this value must be stores with each shift of the crop loop.
        S3 = struct('type','.','subs',strcat(S2.subs,'_vs_',S.subs));
        cropClimateStats = subsasgn(cropClimateStats,S3,F);
        
        %defines the current date
        t = datetime('now');
        t.Format = 'dd.MM.yyyy';
        t = string(t);       

        %creates a legend and defines the name of the file that will be
        %saved to location defined later by fname
        %creats a legend for the climate line, the crop bars, the overall
        %trend line for climate and the overall trend line for crop yield
        lgd = legend([p b zD(2) zC(2)],{string(cell2mat(T2(k))),compose(string(cell2mat(T(h,1)))+" Yield"),compose("Trend of "+ string(cell2mat(T2(k)))),compose("Trend of "+string(cell2mat(T(h,1)))+ "Yield")});           
        %creates a single variable to use as the name for the saved figure
        fileName = compose(t + "_"+ string(cell2mat(T2(k))) +"_" +  string(cell2mat(T(h,1)))+"_countyYield_subplot.tiff");          
      
        %places the legend in the bottom right corner of the plot
        lgd.Position(1) = 0.6;
        lgd.Position(2) = 0.15;
        %increase the text size of the legend for better readibility
        lgd.FontSize = 12;
        %adjusts the save location based on the crop type
        %NOTE: this step could be adjusted so that the file path is based
        %on the current file location so that it doesn't need to be edited
        %by the user if they want to save the files in a different
        %location
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
        
        %NOTE: add the "close" statement to the loop and it seems to help
       
        
        %saveas(fig,fullfile(fname,fileName)) %saves the current figure in the fname folder with the name defined by "fileName"       
        %close
        %NOTE: it might be a good idea to have a script that moves figures
        %currently in the save destination to a new folder for greater
        %clarity. That would probably take a bit of work to write though,
        %so maybe not
    end   
    %storing the statistical results in F for each climate 
   
end
O = stack(cropClimateStats,2:length(cropClimateStats.Properties.VariableNames));
fileName = compose(t + "_"+ "Yield_and_Climate_vs_Year"+"_statistics.csv"); 
%defines the location of the saved file
pathName = 'F:\github\HIndex\Outputs\Current Figures\Crop Yield Results';
%saves the table cropClimateCovariance with the name defined by fileName
%and the location defined by pathName

%writetable(O, fullfile(pathName,fileName))
close all
toc
%% 5.3 Calculating Index VS. Grain Yield
tic
%cropClimateCovariance = table();
locations = ["western", "east-west central", "eastern", "northern", "north-south central", "southern", "statewide"]; %names of each climate region
X = table(); 
for k = 2:width(avgClimateIndex) %for each index, currently set to run for H and W-index(mm)
    T2 = avgClimateIndex.Properties.VariableNames; %T2 is a cell array of the variable names (climate indices) in the avgClimateIndex table
    S2 = struct('type','.','subs',T2(k));
       
           
    for h = 1:length(names) %for each crop  

        T = fieldnames(avgYield); %creates a cell array containing the names of all the crop s in the avgYield structure array
        S = struct('type','.','subs',T(h,1));
        for i = 1:length(locations) %for each climate zone
            A = cell2table(subsref(avgYield,S),'VariableNames',{'YEAR','YIELD'}); %creates a temporary table for the crop yield stored in the avgYield structure array
            B = table(avgClimateIndex.YEAR, subsref(avgClimateIndex,S2),'VariableNames',{'YEAR','Index'});           
            
            E = cell2table(subsref(cropCP,S)); %temporary table 
            changeP1 = E.Var2(i); %this is the change point for the crop yield
            changeP2 = changeP1; %this is the year of the change point for the climate lines, this is currently set to be the same change point as the crop yield change point
            
                       
            %this loop moves the values in yield values in table A to the same
            %years where those yields occur in table B. This makes it easir
            %to use fitlm, since it allows the two data sets to be easily
            %limited to the same sized vector
            for y = 1:length(A.YEAR) %for each year in A
                for j = 1:length(B.YEAR) %for each year in B
                    if B.YEAR(j) == A.YEAR(y) %if the year position j of table B is the same as the year at position y of table A
                       %create a new column in table B name YIELD and set 
                       %all of the columns at row j in B to be equal to all 
                       %the columns in the YIELD variable of row y in table
                       %A
                        B.YIELD(j,:) = A.YIELD(y,:); 
                    end
                end             
            end
           B.YIELD(B.YIELD==0) = NaN; %sets all of the values in the YIELD variable 
           %that are equal to 0 to be equal to NaN instead. This makes
           %ignoring these values easier in the next step
           C = B(~isnan(B.YIELD(:,i)),:);
           
           %creates a linear model with the climate index as the x axis and
           %the crop yield as the y axis. The model only compares years for
           %which both there is a crop yield and climate information
           mdlL = fitlm(C.Index(:,i),C.YIELD(:,i));
           mdlL1 = fitlm(C.Index(C.YEAR<=changeP1,i),C.YIELD(C.YEAR<=changeP1,i));
           %this if statement deals with there being no data for irrigated winter wheat in eastern KS
           if i == 3 && h == 2
               mdlL2 = mdlL; 
           else             
               mdlL2 = fitlm(C.Index(C.YEAR>=changeP1,i),C.YIELD(C.YEAR>=changeP1,i));
           end

           plot(mdlL)
           hold on
           plot(mdlL1)
           hold on
           plot(mdlL2)
           hold off
           
           %total period
           x(i,1) = table2array(mdlL.Coefficients(2,1)); %slope
           x(i,2) = table2array(mdlL.Coefficients(2,4)); %P-Value
           x(i,3) = mdlL.Rsquared.Ordinary; %Rsquared
           %before CP
           x(i,4) = table2array(mdlL1.Coefficients(2,1)); %slope
           x(i,5) = table2array(mdlL1.Coefficients(2,4)); %P-Value
           x(i,6) = mdlL1.Rsquared.Ordinary; %Rsquared
           %after CP
           x(i,7) = table2array(mdlL2.Coefficients(2,1)); %slope
           x(i,8) = table2array(mdlL2.Coefficients(2,4)); %P-Value
           x(i,9) = mdlL2.Rsquared.Ordinary; %Rsquared
           
           x(i,10) = min(C.YEAR); %minimum year
           x(i,11) = changeP1; %change poi nt
           x(i,12) = max(C.YEAR); %maximum year
        end
        %create a new structure variable, S3, with type 'dot' and a subs
        %that of the form "climate_index_name vs crop_name"
        S3 = struct('type','.','subs',strcat(S2.subs,'_vs_',S.subs));
        
        %stores statistics from x in a column samed by S3 in the table
        %cropClimateCovariance
        X = subsasgn(X,S3,x);
        
        
    end
    X.Properties.RowNames = locations;
   
end
X = stack(X,1:length(X.Properties.VariableNames));
X = addvars(X,X.Properties.RowNames,'before',X.Properties.VariableNames(2));
%X = sortrows(X);
cropClimateCovariance = renamevars(X,X.Properties.VariableNames,["Index_vs_Crop","Region","Statistics"]);

%defines the current date
t = datetime('now');
t.Format = 'dd.MM.yyyy';
t = string(t);
%defines the name of the saved file
fileName = compose(t + "_"+ "climate index" +"_vs_" +  "crop yield"+"_statistics.csv"); 
%defines the location of the saved file
pathName = 'F:\github\HIndex\Outputs\Current Figures\Crop Yield Results';

%saves the table cropClimateCovariance with the name defined by fileName
%and the location defined by pathName
%writetable(cropClimateCovariance, fullfile(pathName,fileName))

toc
%% 4.1.0 Graphing H-index decadal analysis
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