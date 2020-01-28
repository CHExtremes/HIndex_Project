%%
%This code is designed to take files from a specific folder, if they are
%excel files they are converted into csv files. those csv files are then
%used to produce an H-index based on temperature for each station.

%Note: you will need to download the "Weather folder" and add it to your
%directory file to ensure this program works. 
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
%this will be the first steps of creating the C/H index. it will start by
%anazyling a single station on a yearly basis and produce a bar graph at
%the end
clc
tic
stationLength = length(stationNames);
%stationLength = 1;
folder = strcat(pwd,'/',newFolder); %calls the path of the current file directory
startYear = 1981;
stopYear = 2013;
useStartYear = 1;
useStopYear = 0;
HIndex = zeros(stopYear-startYear+1,23);
HIndex(:,1) = (startYear:stopYear)';

for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    %creates an array from the starting year to the ending year of the stations available weather data
    if useStartYear == 1 && useStopYear == 0
        YEAR = transpose(startYear:max(temporaryFile.YEAR));
    elseif useStartYear == 1 && useStopYear == 1
        YEAR = transpose(startYear:stopYear);
    elseif useStartYear == 0 && useStopYear == 1
        YEAR = transpose(min(temporaryFile.YEAR):stopYear);
    else     
        YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
    end
    temporaryHIndex = table(YEAR); %creates an column array for the years of the H-Indecies 
    temporaryHIndex.HIndex = zeros(height(temporaryHIndex),1);
    counter = 0;
    for j = temporaryHIndex.YEAR(1):temporaryHIndex.YEAR(end)%for each year at this station
        year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary table for the given year        
        currentTemp = max(year.TMAX); %records max temp for the given year         
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
timeHindex = toc

Results = table;
Results.NAME = tableStationNames';
figure('Name', 'H-Index')
%periods = 
for i = 1:length(tableStationNames) %for each station in station names 
    %this section creates an array of subplots where each station has its
    %own plot showing the data points in orange, the change in time in
    %blue, and the general trend in black.
    tempNames = split(tableStationNames(i), '_'); %assigns a variable to the stations name and cuts out any unnecissary labels used for organization purposes
    A = array2table(HIndex(HIndex(:,1+i) ~= 0,[ 1 1+i])); %for each station, create a new array with the non-zero values and their corresponding years
    A.Properties.VariableNames{'Var1'} = 'YEAR'; 
    A.Properties.VariableNames{'Var2'} = 'HIndex';    
    subplot(4,6,i) %Creates a system of subplots in a 4x6 grid
    l = plot(A.YEAR,A.HIndex); %adds a line to the plot for additional clarity
    hold on %add each station to the same plot
    mdl = fitlm(A, 'HIndex ~ YEAR'); %performs a linear regression for the Year and the H Index
    z = plot(mdl); %plots the linear regression and data points
    x = A.YEAR;
    y = A.HIndex;
    %this code runs trends analysis on x and y as independant and depedant variables
    %it tests the hypothesis of no correlation against the alternative
    %hypothesis of a nonzero correlation. so if p value is smaller than 0.05,
    %we reject the hypothesis.
    [tau,p1]=corr(x,y,'type','kendall'); %kendall method
    tau_p(i,1:3)=[i,tau,p1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
    [rho,p2]=corr(x,y,'type','spearman');%spearman method
    rho_p(i,1:3)=[i,rho,p2];
    [r,p3]=corr(x,y);%pearson (linear) method
    r_p(i,1:3)=[i,r,p3]; %pearson(Least square method) method corrcoef(x,y);    
    %below are changes to the colors and markers of the plot for additional
    %clarity I removed 95% error bars (z(3) and z(4) to make the graph less cluttered. We
    %can turn these on later if we want to visually analyze the error
    %margins.
    z(1).Color = '#D95319'; %sets data points to be orange
    z(1).Marker = '.';
    z(1).MarkerSize = 10;
    z(2).Color = 'k';
    z(2).LineWidth = 1;
    z(3).Color = 'none';
    z(4).Color = 'none';
    legend('off'); %hides the automatic legend generated by fitlm
    
    xlabel('Year', 'FontSize', 11)
    xlim([startYear (stopYear+1)])
    ylabel('H - Index', 'FontSize', 11)  
    ylim([79 95])
    
    %this section  creates a table of the station names and the slope for
    %each trend line
    %p = polyfit(x,y,1);
    Results.HIndexSlope(i) = round(table2array(mdl.Coefficients(2,1)),3);%calls the slope given for the linear regression of the data using the fitlm function
    Results.HIndexAverage(i) = round(mean(A.HIndex),1); %calculates the average H index for the station and adds that to a new table
    Results.HIndexR(i) = r_p(i,2);
    Results.HIndexRsqr(i) = mdl.Rsquared.Ordinary; %calls the r^2 value from the fitlm function and inputs it into a new table
    Results.HIndexRPValue(i) = round(table2array(mdl.Coefficients(2,4)),3); %calls the pValue from the fitlm function and inputs it into a new table
    Results.HIndexRho(i) = rho_p(i,2);
    Results.HIndexRhoPValue(i) = rho_p(i,3);
    Results.HIndexTau(i) = tau_p(i,2);
    Results.HIndexTauPValue(i)= tau_p(i,3);
    title(compose(tempNames(2,1)+"\n"+Results.HIndexAverage(i)+"\n"+Results.HIndexSlope(i)),'FontSize', 11);
end
%use the below script when you want to automatically make tiff files for
%the given graphs. There is an issue with this function in that it doesn't
%expand the window before saving the file, so it compresses subplots to an
%unreadible level. I'll need to find a way to fix this later if we still
%want to use it over manually saving graphs.
%set(gcf,'PaperPositionMode','auto') %set the print area same as paper
%print('-dtiff','-r600', '12.3.19_H_Index_Subplots')



%% 3 
%create a C-Index
clc
tic
stationLength = length(stationNames);
%stationLength = 1;
%calls the path of the current file directory
CIndex = zeros(stopYear-startYear+1,23);
CIndex(:,1) = (startYear:stopYear)';

for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path    
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    %creates an array from the starting year to the ending year of the stations available weather data
    if useStartYear == 1 && useStopYear == 0
        YEAR = transpose(startYear:max(temporaryFile.YEAR));
    elseif useStartYear == 1 && useStopYear == 1
        YEAR = transpose(startYear:stopYear);
    elseif useStartYear == 0 && useStopYear == 1
        YEAR = transpose(min(temporaryFile.YEAR):stopYear);
    else     
        YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
    end
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

figure('Name', 'C-Index')
for i = 1:stationLength %for each station in station names
    %this section creates an array of subplots where each station has it
    %sown plot
    tempNames = split(tableStationNames(i), '_');
    A = array2table(CIndex(CIndex(:,1+i) ~= 0,[ 1 1+i])); %for each station, create a new array with the non-zero values and their corresponding years
    A.Properties.VariableNames{'Var1'} = 'YEAR'; 
    A.Properties.VariableNames{'Var2'} = 'CIndex';
    subplot(4,6,i)
    l = plot(A.YEAR,A.CIndex); %adds a line to the plot for additional clarity
    hold on %add each station to the same plot
    mdl = fitlm(A, 'CIndex ~ YEAR'); %performs a linear regression for the Year and the C Index
    z = plot(mdl); %plots the linear regression and data points
    x = A.YEAR;
    y = A.CIndex;
    %this code runs trends analysis on x and y as independant and depedant variables
    %it tests the hypothesis of no correlation against the alternative
    %hypothesis of a nonzero correlation. so if p value is smaller than 0.05,
    %we reject the hypothesis.
    [tau,p1]=corr(x,y,'type','kendall'); %kendall method
    tau_p(i,1:3)=[i,tau,p1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
    [rho,p2]=corr(x,y,'type','spearman');%spearman method
    rho_p(i,1:3)=[i,rho,p2];
    [r,p3]=corr(x,y);%pearson (linear) method
    r_p(i,1:3)=[i,r,p3]; %pearson(Least square method) method corrcoef(x,y);
    %below are changes to the colors and markers of the plot for additional
    %clarity I removed 95% error bars (z(3) and z(4) to make the graph less cluttered. We
    %can turn these on later if we want to visually analyze the error
    %margins.
    z(1).Color = '#D95319'; %sets data points to be orange
    z(1).Marker = '.';
    z(1).MarkerSize = 10;
    z(2).Color = 'k';
    z(2).LineWidth = 1;
    z(3).Color = 'none';
    z(4).Color = 'none';
    legend('off'); %hides the automatic legend generated by fitlm  
    
    xlabel('Year','FontSize', 11)
    xlim([startYear (stopYear+1)])
    ylabel('C - Index','FontSize', 11)
    ylim([10 30])
    %this section  creates a table of the station names and the slope for
    %each trend line
   
    Results.CIndexSlope(i) = round(table2array(mdl.Coefficients(2,1)),3); %calls the slope given for the linear regression of the data using the fitlm function
    Results.CIndexAverage(i) = round(mean(A.CIndex),1); %calculates the average C index for the station and adds that to a new table
    Results.CIndexR(i) = r_p(i,2);
    Results.CIndexRsqr(i) = mdl.Rsquared.Ordinary; %calls the r^2 value from the fitlm function and inputs it into a new table
    Results.CIndexRPValue(i) = round(table2array(mdl.Coefficients(2,4)),3); %calls the pValue from the fitlm function and inputs it into a new table
    Results.CIndexRho(i) = rho_p(i,2);
    Results.CIndexRhoPValue(i) = rho_p(i,3);
    Results.CIndexTau(i) = tau_p(i,2);
    Results.CIndexTauPValue(i)= tau_p(i,3);
    title(compose(tempNames(2,1)+"\n"+Results.CIndexAverage(i)+"\n"+Results.CIndexSlope(i)),'FontSize', 11);
end
%use the below script when you want to automatically make tiff files for
%the given graphs. See H-Index for issues with this function
%set(gcf,'PaperPositionMode','auto') %set the print area same as paper
%print('-dtiff','-r600', 'C_Index_Subplots')


%% 4
%WIndex will be the Index of wetness. This code will look at each month and
%see if it got any rain at all. 
%start with precip max
clc
tic
%stationLength = 1;
 %calls the path of the current file directory
stationLength = length(stationNames);
WIndex = zeros(stopYear-startYear+1,23);
WIndex(:,1) = (startYear:stopYear)';
for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path    
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    %creates an array from the starting year to the ending year of the stations available weather data
    if useStartYear == 1 && useStopYear == 0
        YEAR = transpose(startYear:max(temporaryFile.YEAR));
    elseif useStartYear == 1 && useStopYear == 1
        YEAR = transpose(startYear:stopYear);
    elseif useStartYear == 0 && useStopYear == 1
        YEAR = transpose(min(temporaryFile.YEAR):stopYear);
    else     
        YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
    end
    temporaryWIndex = table(YEAR); %creates an column array for the years of the H-Indecies 
    temporaryWIndex.WIndex = zeros(height(temporaryWIndex),1);
    counter = 0;
    for j = temporaryWIndex.YEAR(1):temporaryWIndex.YEAR(end)%for each year at this station
        year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary matrix for the given year       
        currentPrecip = max(year.RAIN); %records max temp for the given year        
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
timeWIndex = toc;

figure('Name', 'W-Index')
for i = 1:stationLength %for each station in station names
    %this section creates an array of subplots where each station has it
    %sown plot
    tempNames = split(tableStationNames(i), '_');
    A = array2table(WIndex(WIndex(:,1+i) ~= 0,[ 1 1+i])); %for each station, create a new array with the non-zero values and their corresponding years
    A.Properties.VariableNames{'Var1'} = 'YEAR'; 
    A.Properties.VariableNames{'Var2'} = 'WIndex';
    subplot(4,6,i)
    l = plot(A.YEAR,A.WIndex); %adds a line to the plot for additional clarity
    hold on %add each station to the same plot
    mdl = fitlm(A, 'WIndex ~ YEAR'); %performs a linear regression for the Year and the C Index
    z = plot(mdl); %plots the linear regression and data points
    x = A.YEAR;
    y = A.WIndex;
    %this code runs trends analysis on x and y as independant and depedant variables
    %it tests the hypothesis of no correlation against the alternative
    %hypothesis of a nonzero correlation. so if p value is smaller than 0.05,
    %we reject the hypothesis.
    [tau,p1]=corr(x,y,'type','kendall'); %kendall method
    tau_p(i,1:3)=[i,tau,p1];%j is stations number in  the loop; tou is kendall tou value; and p1 is the p-value for the test.
    [rho,p2]=corr(x,y,'type','spearman');%spearman method
    rho_p(i,1:3)=[i,rho,p2];
    [r,p3]=corr(x,y);%pearson (linear) method
    r_p(i,1:3)=[i,r,p3]; %pearson(Least square method) method corrcoef(x,y);
    %below are changes to the colors and markers of the plot for additional
    %clarity I removed 95% error bars (z(3) and z(4) to make the graph less cluttered. We
    %can turn these on later if we want to visually analyze the error
    %margins.
    z(1).Color = '#D95319'; %sets data points to be orange
    z(1).Marker = '.';
    z(1).MarkerSize = 10;
    z(2).Color = 'k';
    z(2).LineWidth = 1;
    z(3).Color = 'none';
    z(4).Color = 'none';
    legend('off'); %hides the automatic legend generated by fitlm
    
    xlabel('Year', 'FontSize', 11)
    xlim([startYear (stopYear+1)])
    ylabel('W - Index', 'FontSize', 11)
    ylim([5 26])
    %this section  creates a table of the station names and the slope for
    %each trend line
    
    Results.WIndexSlope(i) = round(table2array(mdl.Coefficients(2,1)),3); %calls the slope given for the linear regression of the data using the fitlm function
    Results.WIndexAverage(i) = round(mean(A.WIndex),1); %calculates the average W index for the station and adds that to a new table
    Results.WIndexR(i) = r_p(i,2);
    Results.WIndexRsqr(i) = mdl.Rsquared.Ordinary; %calls the r^2 value from the fitlm function and inputs it into a new table
    Results.WIndexPValue(i) = round(table2array(mdl.Coefficients(2,4)),3); %calls the pValue from the fitlm function and inputs it into a new table
    Results.WIndexRho(i) = rho_p(i,2);
    Results.WIndexRhoPValue(i) = rho_p(i,3);
    Results.WIndexTau(i) = tau_p(i,2);
    Results.WIndexTauPValue(i)= tau_p(i,3);
    title(compose(tempNames(2,1)+"\n"+Results.WIndexAverage(i)+"\n"+Results.WIndexSlope(i)),'FontSize', 11);
end
%use the below script when you want to automatically make tiff files for
%the given graphs. See H-Index for issues with this function
%set(gcf,'PaperPositionMode','auto') %set the print area same as paper
%print('-dtiff','-r600', 'W_Index_Subplots')

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
DIndex = zeros(stopYear-startYear+1,23);
DIndex(:,1) = (startYear:stopYear)';
for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path    
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    %creates an array from the starting year to the ending year of the stations available weather data
    if useStartYear == 1 && useStopYear == 0
        YEAR = transpose(startYear:max(temporaryFile.YEAR));
    elseif useStartYear == 1 && useStopYear == 1
        YEAR = transpose(startYear:stopYear);
    elseif useStartYear == 0 && useStopYear == 1
        YEAR = transpose(min(temporaryFile.YEAR):stopYear);
    else     
        YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
    end
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
    A = array2table(DIndex(DIndex(:,1+i) ~= 0,[ 1 1+i])); %for each station, create a new array with the non-zero values and their corresponding years
    A.Properties.VariableNames{'Var1'} = 'YEAR'; 
    A.Properties.VariableNames{'Var2'} = 'DIndex';
    subplot(4,6,i)
    l = plot(A.YEAR,A.DIndex); %adds a line to the plot for additional clarity
    hold on %add each station to the same plot
    mdl = fitlm(A, 'DIndex ~ YEAR'); %performs a linear regression for the Year and the D Index
    z = plot(mdl); %plots the linear regression and data points
    %below are changes to the colors and markers of the plot for additional
    %clarity I removed 95% error bars (z(3) and z(4) to make the graph less cluttered. We
    %can turn these on later if we want to visually analyze the error
    %margins.
    z(1).Color = '#D95319'; %sets data points to be orange
    z(1).Marker = '.';
    z(1).MarkerSize = 10;
    z(2).Color = 'k';
    z(2).LineWidth = 1;
    z(3).Color = 'none';
    z(4).Color = 'none';
    legend('off'); %hides the automatic legend generated by fitlm  
    
    xlabel('Year', 'FontSize', 11)
    xlim([startYear (stopYear+1)])
    ylabel('D - Index', 'FontSize', 11)
    ylim([5 15])
    %this section  creates a table of the station names and the slope for
    %each trend line
    
    
    Results.DIndexSlope(i) = round(table2array(mdl.Coefficients(2,1)),3); %calls the slope given for the linear regression of the data using the fitlm function
    Results.DIndexAverage(i) = round(mean(A.DIndex),1); %calculates the average D index for the station and adds that to a new table
    Results.DIndexRsqr(i) = mdl.Rsquared.Ordinary; %calls the r^2 value from the fitlm function and inputs it into a new table
    Results.DIndexPValue(i) = round(table2array(mdl.Coefficients(2,4)),3); %calls the pValue from the fitlm function and inputs it into a new table
    title(compose(tempNames(2,1)+"\n"+Results.DIndexAverage(i)+"\n"+Results.DIndexSlope(i)), 'FontSize', 11);
end
%use the below script when you want to automatically make tiff files for
%the given graphs. See H-Index for issues with this function
%set(gcf,'PaperPositionMode','auto') %set the print area same as paper
%print('-dtiff','-r600', 'D_Index_Subplots')
