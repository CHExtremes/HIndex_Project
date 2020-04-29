%% 
%The purpose of this code is to analyze the centenial trend for the yearling average in the top p
%percentile of TMAX and bottom P percentile of TMIN.

clear
clc

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
%This script counts the number of days per year that are above the p
%percentile
clc
tic
stationLength = length(stationNames);
%stationLength = 1;
p = 75;
folder = strcat(pwd,'/',newFolder); %calls the path of the current file directory
percentileFrequency = zeros(124,23);  %creates an empty matrix for the frequency of days above percentile p
percentileFrequency(:,1) = (1890:2013)'; %adds a year column to the empty matrix for the whole data set
percentileMagnitude = zeros(124,23); %creates an empty matrix for the magnitude of p percentile per year
percentileMagnitude(:,1) = (1890:2013)'; %adds a year column to the empty matrix for the whole data set

for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    
    YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR)); %creates an array from the starting year to the ending year of the stations available weather data
    tPercentileFrequency = table(YEAR); %creates an empty variable for tracking the frequency of days above percentile p for a given station 
    tPercentileFrequency.Frequency = zeros(height(tPercentileFrequency),1); %pre-allocates a collumn named "frequency" to the temporary table
    tPercentileMagnitude = table(YEAR); %creats an empty variable for tracking the magnitude of the temperature for the p percentile
    tPercentileMagnitude.Magnitude = zeros(height(tPercentileMagnitude),1); %pre-allocates a collumn named "magnitude" to the temporary table
    counter = 0;
    
    referenceSet = temporaryFile; %creates a variable "referenceSet" which referes to the current 30 year climtelogical period used as a reference for percentile analysis
    toDelete =  temporaryFile.YEAR <1981 | temporaryFile.YEAR > 2010 ; %
    referenceSet(toDelete,:) = [];
    
    currentTemp = prctile(referenceSet.TMAX,p); %currentTemp is the reference percentile used for the frequency analysis
    for j = tPercentileFrequency.YEAR(1):tPercentileFrequency.YEAR(end)%for each year at this station
        year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary table for the given year        
        tPercentileMagnitude.Magnitude(j-tPercentileFrequency.YEAR(1)+1) = prctile(year.TMAX,p);  %calculates what the p percentile is for this year and stores it in tPercentileMagnitude          
        counter = 0;  
        %This loop counts the number of days per year where TMAX is greater
        %than or equal to the reference temperature, currentTemp
            for h = 1:height(year)% for days in this year                
                if year.TMAX(h) >= currentTemp %If the value at row h and column TMaxColumn are greater than currentTemp
                    counter = counter + 1; %increase counter by 1
                end    
            end                    
        tPercentileFrequency.Frequency(j-tPercentileFrequency.YEAR(1)+1) = counter; %stores the frequency of days above currentTemp for the current year 
    end
    %this loop transfers the information from the temporary percentile
    %matricies to the final percentile matrix
    for j = 1:height(tPercentileFrequency) %for the number of years at the current station
        for h = 1:length(percentileFrequency) %for full array of years being analyzed
           if tPercentileFrequency.YEAR(j) == percentileFrequency(h,1) %Checks to make sure that the years are the same for the given station
              percentileFrequency(h,(i+1))= tPercentileFrequency.Frequency(j); %if the years are the same, then the yearly value for the station is stored in a column specifically for that station
              percentileMagnitude(h,(i+1))= tPercentileMagnitude.Magnitude(j);
           end
        end
    end
end
%% 3
%This section graphs the results for section 2
Slopes = table;
Slopes.NAME = tableStationNames';
figure('Name', 'Frequency')
%graph the frequency of the percentile
for i = 1:length(tableStationNames) %for each station in station names 
    %this section creates an array of subplots where each station has it
    %sown plot
    tempNames = split(tableStationNames(i), '_');
    A = percentileFrequency(percentileFrequency(:,1+i) ~= 0,[ 1 1+i]); %for each station, create a new array with the non-zero values and their corresponding years
    x = A(:,1);%x-axis is the year
    y = A(:,2);%y-axis is the HIndex
    subplot(4,6,i) %Creates a system of subplots in a 4x6 grid
    plot(x, y)  %creates a line connecting the scatter plot for clarity
    hold on %add each station to the same plot
    scatter(x,y,15,'filled')%creates a scatter plot so lsline can be used to graph the linear trendline
    s = lsline; 
    s.Color = 'k'; %changes the color of the trendline to black   
    xlabel('Year')
    xlim([1890 2014])
    ylabel('Number of Days')  
    ylim([30 160])
    %this section  creates a table of the station names and the slope and average for
    %each trend line
    p = polyfit(x,y,1);
    Slopes.percentileFrequencySlope(i) = round(p(1),3);
    Slopes.percentileFrequencyAverage(i) = round(mean(A(:,2)),1);
    title(compose(tempNames(2,1)+"\n"+Slopes.percentileFrequencyAverage(i)+"\n"+Slopes.percentileFrequencySlope(i)));
end
%graphs the magnitudes for the percentile
figure('Name', 'Percentile Magnitude')
for i = 1:length(tableStationNames) %for each station in station names 
    %this section creates an array of subplots where each station has it
    %sown plot
    tempNames = split(tableStationNames(i), '_');
    A = percentileMagnitude(percentileMagnitude(:,1+i) ~= 0,[ 1 1+i]); %for each station, create a new array with the non-zero values and their corresponding years
    x = A(:,1);%x-axis is the year
    y = A(:,2);%y-axis is the HIndex
    subplot(4,6,i) %Creates a system of subplots in a 4x6 grid
    plot(x, y)  %creates a line connecting the scatter plot for clarity
    hold on %add each station to the same plot
    scatter(x,y,15,'filled')%creates a scatter plot so lsline can be used to graph the linear trendline
    s = lsline; 
    s.Color = 'k'; %changes the color of the trendline to black   
    xlabel('Year')
    xlim([1890 2014])
    ylabel('Temperature (°F)')  
    ylim([80 95])
    %this section  creates a table of the station names and the slope for
    %each trend line
    p = polyfit(x,y,1);
    Slopes.percentileMagnitudeSlope(i) = round(p(1),3);
    Slopes.percentileMagnitudeAverage(i) = round(mean(A(:,2)),1);
    title(compose(tempNames(2,1)+"\n"+Slopes.percentileMagnitudeAverage(i)+"\n"+Slopes.percentileMagnitudeSlope(i)));
end   
timePercentileTMAX = toc;
%% 4 
%This script calculates the magnitude of temperature and the frequency of
%days below the p percentile
clc
tic
stationLength = length(stationNames);
%stationLength = 1;
p = 5;
folder = strcat(pwd,'/',newFolder); %calls the path of the current file directory
percentileFrequency = zeros(124,23);  %creates an empty matrix for the frequency of days above percentile p
percentileFrequency(:,1) = (1890:2013)'; %adds a year column to the empty matrix for the whole data set
percentileMagnitude = zeros(124,23); %creates an empty matrix for the magnitude of p percentile per year
percentileMagnitude(:,1) = (1890:2013)'; %adds a year column to the empty matrix for the whole data set

for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    
    YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR)); %creates an array from the starting year to the ending year of the stations available weather data
    tPercentileFrequency = table(YEAR); %creates an empty variable for tracking the frequency of days above percentile p for a given station 
    tPercentileFrequency.Frequency = zeros(height(tPercentileFrequency),1); %pre-allocates a collumn named "frequency" to the temporary table
    tPercentileMagnitude = table(YEAR); %creats an empty variable for tracking the magnitude of the temperature for the p percentile
    tPercentileMagnitude.Magnitude = zeros(height(tPercentileMagnitude),1); %pre-allocates a collumn named "magnitude" to the temporary table
    counter = 0;
    
    referenceSet = temporaryFile; %creates a variable "referenceSet" which referes to the current 30 year climtelogical period used as a reference for percentile analysis
    toDelete =  temporaryFile.YEAR <1981 | temporaryFile.YEAR > 2010 ; %
    referenceSet(toDelete,:) = [];
    
    currentTemp = prctile(referenceSet.TMIN,p); %currentTemp is the reference percentile used for the frequency analysis
    for j = tPercentileFrequency.YEAR(1):tPercentileFrequency.YEAR(end)%for each year at this station
        year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary table for the given year        
        tPercentileMagnitude.Magnitude(j-tPercentileFrequency.YEAR(1)+1) = prctile(year.TMIN,p);  %calculates what the p percentile is for this year and stores it in tPercentileMagnitude          
        counter = 0;  
        %This loop counts the number of days per year where TMAX is greater
        %than or equal to the reference temperature, currentTemp
            for h = 1:height(year)% for days in this year                
                if year.TMIN(h) <= currentTemp %If the value at row h and column TMaxColumn are greater than currentTemp
                    counter = counter + 1; %increase counter by 1
                end    
            end                    
        tPercentileFrequency.Frequency(j-tPercentileFrequency.YEAR(1)+1) = counter; %stores the frequency of days above currentTemp for the current year 
    end
    %this loop transfers the information from the temporary percentile
    %matricies to the final percentile matrix
    for j = 1:height(tPercentileFrequency) %for the number of years at the current station
        for h = 1:length(percentileFrequency) %for full array of years being analyzed
           if tPercentileFrequency.YEAR(j) == percentileFrequency(h,1) %Checks to make sure that the years are the same for the given station
              percentileFrequency(h,(i+1))= tPercentileFrequency.Frequency(j); %if the years are the same, then the yearly value for the station is stored in a column specifically for that station
              percentileMagnitude(h,(i+1))= tPercentileMagnitude.Magnitude(j);
           end
        end
    end
end
%% 5
%This section graphs the results for section 2
Slopes = table;
Slopes.NAME = tableStationNames';
figure('Name', 'Frequency')
%graph the frequency of the percentile
for i = 1:length(tableStationNames) %for each station in station names 
    %this section creates an array of subplots where each station has it
    %sown plot
    tempNames = split(tableStationNames(i), '_');
    A = percentileFrequency(percentileFrequency(:,1+i) ~= 0,[ 1 1+i]); %for each station, create a new array with the non-zero values and their corresponding years
    x = A(:,1);%x-axis is the year
    y = A(:,2);%y-axis is the HIndex
    subplot(4,6,i) %Creates a system of subplots in a 4x6 grid
    plot(x, y)  %creates a line connecting the scatter plot for clarity
    hold on %add each station to the same plot
    scatter(x,y,15,'filled')%creates a scatter plot so lsline can be used to graph the linear trendline
    s = lsline; 
    s.Color = 'k'; %changes the color of the trendline to black   
    xlabel('Year')
    xlim([1890 2014])
    ylabel('Number of Days')  
    %ylim([0 30])
    %this section  creates a table of the station names and the slope and average for
    %each trend line
    p = polyfit(x,y,1);
    Slopes.percentileFrequencySlope(i) = round(p(1),3);
    Slopes.percentileFrequencyAverage(i) = round(mean(A(:,2)),1);
    title(compose(tempNames(2,1)+"\n"+Slopes.percentileFrequencyAverage(i)+"\n"+Slopes.percentileFrequencySlope(i)));
end
%graphs the magnitudes for the percentile
figure('Name', 'Percentile Magnitude')
for i = 1:length(tableStationNames) %for each station in station names 
    %this section creates an array of subplots where each station has it
    %sown plot
    tempNames = split(tableStationNames(i), '_');
    A = percentileMagnitude(percentileMagnitude(:,1+i) ~= 0,[ 1 1+i]); %for each station, create a new array with the non-zero values and their corresponding years
    x = A(:,1);%x-axis is the year
    y = A(:,2);%y-axis is the HIndex
    subplot(4,6,i) %Creates a system of subplots in a 4x6 grid
    plot(x, y)  %creates a line connecting the scatter plot for clarity
    hold on %add each station to the same plot
    scatter(x,y,15,'filled')%creates a scatter plot so lsline can be used to graph the linear trendline
    s = lsline; 
    s.Color = 'k'; %changes the color of the trendline to black   
    xlabel('Year')
    xlim([1890 2014])
    ylabel('Temperature (°F)')  
    %ylim([-10 25])
    %this section  creates a table of the station names and the slope for
    %each trend line
    p = polyfit(x,y,1);
    Slopes.percentileMagnitudeSlope(i) = round(p(1),3);
    Slopes.percentileMagnitudeAverage(i) = round(mean(A(:,2)),1);
    title(compose(tempNames(2,1)+"\n"+Slopes.percentileMagnitudeAverage(i)+"\n"+Slopes.percentileMagnitudeSlope(i)));
end   
timePercentileTMIN = toc;
