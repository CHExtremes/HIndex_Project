%%
%This code is designed to take files from a specific folder, if they are
%excel files they are converted into csv files. those csv files are then
%used to produce an H-index based on temperature for each station.

%Note: you will need to download the "Weather folder" and add it to your
%directory file to ensure this program works. 
clear all
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
%this will be the first steps of creating the C/H index. it will start by
%anazyling a single station on a yearly basis and produce a bar graph at
%the end
clc
tic
stationLength = length(stationNames);
%stationLength = 1;
folder = strcat(pwd,'/',newFolder); %calls the path of the current file directory
HIndex = zeros(124,23);
HIndex(:,1) = (1890:2013)';

for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    %creates an array from the starting year to the ending year of the stations available weather data
    YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
    temporaryHIndex = table(YEAR); %creates an column array for the years of the H-Indecies 
    temporaryHIndex.HIndex = zeros(height(temporaryHIndex),1);
    counter = 0;
    for j = temporaryHIndex.YEAR(1):temporaryHIndex.YEAR(end)%for each year at this station
        year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary table for the given year        
        currentTemp = round(max(year.TMAX)); %records max temp for the given year         
        %Count the number of times where the daily temp is greater than or
        %equal to that temp
        counter = 0;
        while counter < currentTemp %checks to see if the counter is smaller than the currentTemp. This is to make sure that the value is an H-Index value.
            counter = 0;
            for h = 1:height(year)% for days in this year                
                if year.TMAX(h) >= currentTemp %If the value at row h and column TMaxColumn are greater than currentTemp
                    counter = counter + 1; %increase counter by 1
                end    
            end
            if counter < currentTemp %if the counter is smaller than currentTemp then the H-index is not valid, so we reduce it by one and repeat the loop.
                currentTemp = currentTemp - 1;
            end
        end              
        temporaryHIndex.HIndex(j-temporaryHIndex.YEAR(1)+1) = currentTemp; 
    end
    
    for j = 1:height(temporaryHIndex) %for the number of years at the current station
        for h = 1:length(HIndex) %for full array of years being analyzed
           if temporaryHIndex.YEAR(j) == HIndex(h,1) %Checks to make sure that the years are the same for the given station
              HIndex(h,(i+1))=temporaryHIndex.HIndex(j); %if the years are the same, then the yearly value for the station is stored in a column specifically for that station
           end
        end
    end
end

Slopes = table;
Slopes.NAME = tableStationNames';
figure('Name', 'H-Index')
for i = 1:length(tableStationNames) %for each station in station names 
    %this section creates an array of subplots where each station has it
    %sown plot
    tempNames = split(tableStationNames(i), '_');
    A = HIndex(HIndex(:,1+i) ~= 0,[ 1 1+i]); %for each station, create a new array with the non-zero values and their corresponding years
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
    ylabel('H-Index')  
    ylim([79 95])
    %this section  creates a table of the station names and the slope for
    %each trend line
    p = polyfit(x,y,1);
    Slopes.HIndexSlope(i) = round(p(1),3);
    Slopes.HIndexAverage(i) = round(mean(A(:,2)),1);
    title(compose(tempNames(2,1)+"\n"+Slopes.HIndexAverage(i)+"\n"+Slopes.HIndexSlope(i)));
end


    
timeHindex = toc

%% 3 
%create a C-Index
clc
tic
stationLength = length(stationNames);
%stationLength = 1;
 %calls the path of the current file directory
CIndex = zeros(124,23);
CIndex(:,1) = (1890:2013)';

for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path    
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    %creates an array from the starting year to the ending year of the stations available weather data
    YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
    temporaryCIndex = table(YEAR); %creates an column array for the years of the H-Indecies 
    temporaryCIndex.CIndex = zeros(height(temporaryCIndex),1);
    counter = 0;
    for j = temporaryCIndex.YEAR(1):temporaryCIndex.YEAR(end)%for each year at this station
        year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary matrix for the given year       
        currentTemp = round(min(year.TMIN)); %records max temp for the given year        
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
        temporaryCIndex.CIndex(j-temporaryCIndex.YEAR(1)+1) = currentTemp; %stores the max H-index value for each year in temporaryHIndex      
    end
    
    for j = 1:height(temporaryCIndex) %for the number of years at the current station
        for h = 1:length(CIndex) %for full array of years being analyzed
           if temporaryCIndex.YEAR(j) == CIndex(h,1) %Checks to make sure that the years are the same for the given station
              CIndex(h,i+1)=temporaryCIndex.CIndex(j); %if the years are the same, then the yearly value for the station is stored in a column specifically for that station
           end
        end
    end
end
figure('Name', 'C-Index')
for i = 1:stationLength %for each station in station names
    %this section creates an array of subplots where each station has it
    %sown plot
    tempNames = split(tableStationNames(i), '_');
    A = CIndex(CIndex(:,1+i) ~= 0,[ 1 1+i]); %for each station, create a new array with the non-zero values and their corresponding years
    x = A(:,1);%x-axis is the year
    y = A(:,2);%y-axis is the HIndex
    subplot(4,6,i)
    plot(x,y)
    hold on
    scatter(x, y,15,'filled')  
    s = lsline;
    s.Color = 'k';
    xlabel('Year')
    xlim([1890 2014])
    ylabel('C-Index')
    ylim([0 25])
    %this section  creates a table of the station names and the slope for
    %each trend line
    p = polyfit(x,y,1);
    Slopes.CIndexSlope(i) = round(p(1),3);
    Slopes.CIndexAverage(i) = round(mean(A(:,2)),1);
    title(compose(tempNames(2,1)+"\n"+Slopes.CIndexAverage(i)+"\n"+Slopes.CIndexSlope(i)));
end
timeCIndex = toc

%% 4
%WIndex will be the Index of wetness. This code will look at each month and
%see if it got any rain at all. 
%start with precip max
clc
tic
%stationLength = 1;
 %calls the path of the current file directory
stationLength = length(stationNames);
WIndex = zeros(124,23);
WIndex(:,1) = (1890:2013)';
for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path    
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    %creates an array from the starting year to the ending year of the stations available weather data
    YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
    temporaryWIndex = table(YEAR); %creates an column array for the years of the H-Indecies 
    temporaryWIndex.WIndex = zeros(height(temporaryWIndex),1);
    counter = 0;
    for j = temporaryWIndex.YEAR(1):temporaryWIndex.YEAR(end)%for each year at this station
        year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary matrix for the given year       
        currentPrecip = round(max(year.RAIN)); %records max temp for the given year        
        %Count the number of times where the daily temp is greater than or
        %equal to that temp
        while counter < currentPrecip %checks to see if the counter is smaller than the currentTemp. This is to make sure that the value is an H-Index value.
            counter = 0;
            for h = 1:height(year)% for days in this year
                if year.RAIN(h) >= currentPrecip %If the value at row h and column TMaxColumn are greater than currentTemp
                    counter = counter + 1; %increase counter by 1
                end    
            end
            if counter < currentPrecip %if the counter is smaller than currentTemp then the H-index is not valid, so we reduce it by one and repeat the loop.
                currentPrecip = currentPrecip - 1 ;
            end
        end
        temporaryWIndex.WIndex(j-temporaryWIndex.YEAR(1)+1) = currentPrecip; %stores the max H-index value for each year in temporaryHIndex      
    end
    
    for j = 1:height(temporaryWIndex) %for the number of years at the current station
        for h = 1:length(WIndex) %for full array of years being analyzed
           if temporaryWIndex.YEAR(j) == WIndex(h,1) %Checks to make sure that the years are the same for the given station
              WIndex(h,i+1)=temporaryWIndex.WIndex(j); %if the years are the same, then the yearly value for the station is stored in a column specifically for that station
           end
        end
    end  
end

figure('Name', 'W-Index')
for i = 1:stationLength %for each station in station names
    %this section creates an array of subplots where each station has it
    %sown plot
    tempNames = split(tableStationNames(i), '_');
    A = WIndex(WIndex(:,1+i) ~= 0,[ 1 1+i]); %for each station, create a new array with the non-zero values and their corresponding years
    x = A(:,1);%x-axis is the year
    y = A(:,2);%y-axis is the HIndex
    subplot(4,6,i)
    plot(x,y)
    hold on
    scatter(x, y,15,'filled')  
    s = lsline;
    s.Color = 'k';
    title(tempNames(2,1));
    xlabel('Year')
    xlim([1890 2014])
    ylabel('W-Index')
    ylim([5 26])
    %this section  creates a table of the station names and the slope for
    %each trend line
    p = polyfit(x,y,1);
    Slopes.WIndexSlope(i) = round(p(1),3);
    Slopes.WIndexAverage(i) = round(mean(A(:,2)),1);
    title(compose(tempNames(2,1)+"\n"+Slopes.WIndexAverage(i)+"\n"+Slopes.WIndexSlope(i)));
end
timeWIndex = toc
%% 5
%Writes a dry index, D-Index. This could be done by counting x number of
%periods that had x number of days with no rain. I don't think that will be
%very informative since prolonged dry periods(drought) are more likely to
%be influential than knowning that there are more, shorter dry periods.
clc
tic
%stationLength = 1;
 %calls the path of the current file directory
stationLength = length(stationNames);
DIndex = zeros(124,23);
DIndex(:,1) = (1890:2013)';
for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path    
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    %creates an array from the starting year to the ending year of the stations available weather data
    YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
    temporaryDIndex = table(YEAR); %creates an column array for the years of the H-Indecies 
    temporaryDIndex.DIndex = zeros(height(temporaryDIndex),1);
    counter = 0;
    dryness = 1; %dryness is the minimum precipitation threshold for a day to be considered 'wet'
    for j = temporaryDIndex.YEAR(1):temporaryDIndex.YEAR(end)%for each year at this station
            year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary matrix for the given year       
            currentDryPeriod = 0;
            tempDryPeriod = 0;
            counter = 0;
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
    end
    
    for j = 1:height(temporaryDIndex) %for the number of years at the current station
        for h = 1:length(DIndex) %for full array of years being analyzed
           if temporaryDIndex.YEAR(j) == DIndex(h,1) %Checks to make sure that the years are the same for the given station
              DIndex(h,i+1)=temporaryDIndex.DIndex(j); %if the years are the same, then the yearly value for the station is stored in a column specifically for that station
           end
        end
    end  
end
timeDIndex = toc;

figure('Name', 'D-Index')
for i = 1:stationLength %for each station in station names
    %this section creates an array of subplots where each station has it
   
    tempNames = split(tableStationNames(i), '_');
    A = DIndex(DIndex(:,1+i) ~= 0,[ 1 1+i]); %for each station, create a new array with the non-zero values and their corresponding years
    x = A(:,1);%x-axis is the year
    y = A(:,2);%y-axis is the DIndex
    subplot(4,6,i)
    plot(x,y)
    hold on
    scatter(x, y,15,'filled')  
    s = lsline;
    s.Color = 'k';
    title(tempNames(2,1));
    xlabel('Year')
    xlim([1890 2014])
    ylabel('D-Index')
    ylim([5 15])
    %this section  creates a table of the station names and the slope for
    %each trend line
    p = polyfit(x,y,1);
    Slopes.DIndexSlope(i) = round(p(1),3);
    Slopes.DIndexAverage(i) = round(mean(A(:,2)),1);
    title(compose(tempNames(2,1)+"\n"+Slopes.DIndexAverage(i)+"\n"+Slopes.DIndexSlope(i)));
end

