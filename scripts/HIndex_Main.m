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
HIndex = zeros(2013-1890+1,23);
HIndex(:,1) = (1890:2013)';

% % This section of script can be used to make the entire H-Index data set
% adjustable by date. However, I've removed it for now because I can adjust
% date during the graphin section more effectively
% if useStartYear == 1 && useStopYear == 1
%     HIndex = zeros(stopYear-startYear+1,23);
%     HIndex(:,1) = (startYear:stopYear)';
% elseif useStartYear == 1 && useStopeYear == 0
%     HIndex = zeros(2013-startYear+1,23);
%     HIndex(:,1) = (startYear:2013)';
% elseif useStartYear == 0 && useStopYear == 1
%     HIndex = zeros(1890-startYear+1,23);
%     HIndex(:,1) = (1890:stopYear)';
% else
%     HIndex = zeros(2013-1890+1,23);
%     HIndex(:,1) = (1890:2013)';
% end

for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    temporaryFile.TMAX = round(temporaryFile.TMAX,0);
    %creates an array from the starting year to the ending year of the stations available weather data
%     commented out this section of script for the same reason as above.
%     if useStartYear == 1 && useStopYear == 0
%         YEAR = transpose(startYear:max(temporaryFile.YEAR));
%     elseif useStartYear == 1 && useStopYear == 1
%         YEAR = transpose(startYear:stopYear);
%     elseif useStartYear == 0 && useStopYear == 1
%         YEAR = transpose(min(temporaryFile.YEAR):stopYear);
%     else     
%         YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
%     end
    YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
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
%%
%creates a table to store all statistical values for all results for a
%given Index
HResults = table;
HResults.NAME = tableStationNames';
%
startYear = 1981;
stopYear = 2010;
useStartYear = 0;
useStopYear = 0;
figure('Name', 'H-Index')
%size of the plot
m = 4;
n = 6;
%position of plots

for i = 1:length(tableStationNames) %for each station in station names 
    %this section creates an array of subplots where each station has its
    %own plot showing the data points in orange, the change in time in
    %blue, and the general trend in black.
    tempNames = split(tableStationNames(i), '_'); %assigns a variable to the stations name and cuts out any unnecissary labels used for organization purposes
    A = array2table(HIndex(HIndex(:,1+i) ~= 0,[ 1 1+i])); %for each station, create a new array with the non-zero values and their corresponding years
    A.Properties.VariableNames{'Var1'} = 'YEAR'; 
    A.Properties.VariableNames{'Var2'} = 'HIndex';  
    
    B = A(A.YEAR>=startYear & A.YEAR <=stopYear,:);
    subplot(4,6,i)
    pos=get(gca,'Position');
    set(gca,'Position',[pos(1,1) pos(1,2) pos(1,3)+.03 pos(1,4)]) %subplot position
    
    %subplot('Position',[(j-1)*1/width (height-i)*1/numRecsDown 1/width 1/height])
    lA = plot(A.YEAR,A.HIndex); %adds a line to the plot for additional clarity
    hold on %add each station to the same plot
    
    mdlA = fitlm(A, 'HIndex ~ YEAR'); %performs a linear regression for the Year and the H Index    
    zA = plot(mdlA); %plots the linear regression and data point
    
    mdlB = fitlm(B, 'HIndex ~ YEAR');
    zB = plot(mdlB);
    
    xA = A.YEAR;
    yA = A.HIndex;
    xB = B.YEAR;
    yB = B.HIndex;  
      
    %below are changes to the colors and markers of the plot for additional
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
    zB(1).Color = 'none';
    zB(1).Marker = '.';
    zB(1).MarkerSize = 10;
    zB(2).Color = 'r';
    zB(2).LineWidth = 1;
    zB(3).Color = 'none';
    zB(4).Color = 'none';
    
    legend('off'); %hides the automatic legend generated by fitlm
    
    if i == 19 || i == 20 || i == 21 || i == 22 || i == 23 
        xlabel('Year', 'FontSize', 11)
        set(gca,'Xtick',[1921 1951 1981 2011]);
        %xticks(gca, 1890:2010, 4)
    else
        xlabel('')
        xticks('')
    end
    
    if useStartYear == 1 && useStopYear == 0
        xlim([startYear 2014])
    elseif useStartYear == 1 && useStopYear == 1
        xlim([startYear (stopYear+1)])
    elseif useStartYear == 0 && useStopYear == 1
        xlim([1890 (stopYear+1)])
    else     
        xlim([1890 2014])
    end
    
    ylim([79 95]) %sets boundries for the y axis for all graphs
    %ylim([55 72])
    if i == 1 || i == 7 || i == 13 || i == 19
        ylabel('H - Index (days/year)', 'FontSize', 11)
        yticks('auto')
    else
        ylabel('')
        yticks('')
    end
    %this section  creates a table statistically important values for both
    %the POR and the specified period
    
    %p = polyfit(x,y,1);
       %this code runs trends analysis on x and y as independant and depedant variables
    %it tests the hypothesis of no correlation against the alternative
    %hypothesis of a nonzero correlation. so if p value is smaller than 0.05,
    %we reject the hypothesis.
    %analysis for the whole period of record
    [tau,p1]=corr(xA,yA,'type','kendall'); %kendall method
    tau_pA(i,1:3)=[i,tau,p1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
    [rho,p2]=corr(xA,yA,'type','spearman');%spearman method
    rho_pA(i,1:3)=[i,rho,p2];
    [r,p3]=corr(xA,yA);%pearson (linear) method
    r_pA(i,1:3)=[i,r,p3]; %pearson(Least square method) method corrcoef(x,y);    
    
    HResults.slopePOR(i) = round(table2array(mdlA.Coefficients(2,1)),3);%calls the slope given for the linear regression of the data using the fitlm function
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
    
     %Statistical analysis for a specified period denoted by either a B or SP
     %(Specified Period)
    [tau,p1]=corr(xB,yB,'type','kendall'); %kendall method
    tau_pB(i,1:3)=[i,tau,p1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
    [rho,p2]=corr(xB,yB,'type','spearman');%spearman method
    rho_pB(i,1:3)=[i,rho,p2];
    [r,p3]=corr(xB,yB);%pearson (linear) method
    r_pB(i,1:3)=[i,r,p3]; %pearson(Least square method) method corrcoef(x,y); 
    
    HResults.slopeSP(i) = round(table2array(mdlB.Coefficients(2,1)),3);%calls the slope given for the linear regression of the data using the fitlm functioN
    HResults.rSP(i) = r_pB(i,2);
    HResults.rSqrSP(i) = mdlB.Rsquared.Ordinary; %calls the r^2 value from the fitlm function and inputs it into a new table
    HResults.rPValueSP(i) = round(table2array(mdlB.Coefficients(2,4)),3); %calls the pValue from the fitlm function and inputs it into a new table
    HResults.rhoSP(i) = rho_pB(i,2);
    HResults.rhoPValueSP(i) = rho_pB(i,3);
    HResults.tauSP(i) = tau_pB(i,2);
    HResults.tauPValueSP(i)= tau_pB(i,3);
    HResults.minSP(i) = min(B.HIndex);
    HResults.maxSP(i) = max(B.HIndex);
    HResults.meanSP(i) = round(mean(B.HIndex),1);
    HResults.medianSP(i) = median(B.HIndex);
    title(compose(tempNames(2,1)+"\n"+num2str(HResults.slopePOR(i)*100,"%#.1f")+"\n"+num2str(HResults.slopeSP(i)*100,"%#.1f")),'FontSize', 11);
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
%%
CResults = table;
CResults.NAME = tableStationNames';
%for now, I'm only controlling the year from the H-Index section
% startYear = 1981;
% stopYear = 2013;
% useStartYear = 0;
% useStopYear = 0;
figure('Name', 'C-Index')
%periods = 
for i = 1:length(tableStationNames) %for each station in station names 
    %this section creates an array of subplots where each station has its
    %own plot showing the data points in orange, the change in time in
    %blue, and the general trend in black.
    tempNames = split(tableStationNames(i), '_'); %assigns a variable to the stations name and cuts out any unnecissary labels used for organization purposes
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
    mdlB = fitlm(B, 'CIndex ~ YEAR');
    zB = plot(mdlB);
    xA = A.YEAR;
    yA = A.CIndex;
    xB = B.YEAR;
    yB = B.CIndex;  
      
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
    zB(1).Color = 'none';%'#D95319'; %sets data points to be orange
    zB(1).Marker = '.';
    zB(1).MarkerSize = 10;
    zB(2).Color = 'r';
    zB(2).LineWidth = 1;
    zB(3).Color = 'none';
    zB(4).Color = 'none';
    
    legend('off'); %hides the automatic legend generated by fitlm
    %determines if the scale of the graph needs to be changed based on the
    %size of the year
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
    [tau,p1]=corr(xA,yA,'type','kendall'); %kendall method
    tau_pA(i,1:3)=[i,tau,p1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
    [rho,p2]=corr(xA,yA,'type','spearman');%spearman method
    rho_pA(i,1:3)=[i,rho,p2];
    [r,p3]=corr(xA,yA);%pearson (linear) method
    r_pA(i,1:3)=[i,r,p3]; %pearson(Least square method) method corrcoef(x,y);    
    
    CResults.slopePOR(i) = round(table2array(mdlA.Coefficients(2,1)),3);%calls the slope given for the linear regression of the data using the fitlm function
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
    
     %Statistical analysis for a specified period denoted by either a B or SP
     %(Specified Period)
    [tau,p1]=corr(xB,yB,'type','kendall'); %kendall method
    tau_pB(i,1:3)=[i,tau,p1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
    [rho,p2]=corr(xB,yB,'type','spearman');%spearman method
    rho_pB(i,1:3)=[i,rho,p2];
    [r,p3]=corr(xB,yB);%pearson (linear) method
    r_pB(i,1:3)=[i,r,p3]; %pearson(Least square method) method corrcoef(x,y); 
    
    CResults.slopeSP(i) = round(table2array(mdlB.Coefficients(2,1)),3);%calls the slope given for the linear regression of the data using the fitlm functioN
    CResults.rSP(i) = r_pB(i,2);
    CResults.rSqrSP(i) = mdlB.Rsquared.Ordinary; %calls the r^2 value from the fitlm function and inputs it into a new table
    CResults.rPValueSP(i) = round(table2array(mdlB.Coefficients(2,4)),3); %calls the pValue from the fitlm function and inputs it into a new table
    CResults.rhoSP(i) = rho_pB(i,2);
    CResults.rhoPValueSP(i) = rho_pB(i,3);
    CResults.tauSP(i) = tau_pB(i,2);
    CResults.tauPValueSP(i)= tau_pB(i,3);
    CResults.minSP(i) = min(B.CIndex);
    CResults.maxSP(i) = max(B.CIndex);
    CResults.meanSP(i) = round(mean(B.CIndex),1);
    CResults.medianSP(i) = median(B.CIndex);
    title(compose(tempNames(2,1)+"\n"+num2str(100*CResults.slopePOR(i),"%#.1f")+"\n"+num2str(100*CResults.slopeSP(i),"%#.1f")),'FontSize', 11);
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
WIndex = zeros(2013-1890+1,23);
WIndex(:,1) = (1890:2013)';
useInches = 0;
for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path    
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    temporaryFile.RAIN = round(temporaryFile.RAIN,0); %rounds the result to be a whole number, since we don't have measure values for mm to the tenthousandths place.
    if useInches == 1        
        temporaryFile.RAIN = round(temporaryFile.RAIN/25.4,2); %converts precipitations values from mm to inches
    end
    %creates an array from the starting year to the ending year of the stations available weather data
    YEAR = transpose(min(temporaryFile.YEAR):max(temporaryFile.YEAR));
    temporaryWIndex = table(YEAR); %creates an column array for the years of the H-Indecies 
    temporaryWIndex.WIndex = zeros(height(temporaryWIndex),1);
    counter = 1;
    for j = temporaryWIndex.YEAR(1):temporaryWIndex.YEAR(end)%for each year at this station
        year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary matrix for the given year       
        currentPrecip = max(year.RAIN); %records max precip for the given year        
        %Count the number of times where the daily precip is greater than or
        %equal to the maximum precip
        while counter < currentPrecip %checks to see if the counter is smaller than the currentTemp. This is to make sure that the value is an H-Index value.
            counter = 0;
            for h = 1:height(year)% for days in this year
                if year.RAIN(h) >= currentPrecip %If the value at row h and column TMaxColumn are greater than currentTemp
                   if useInches == 1
                       counter = counter + 0.01;
                   else
                       counter = counter + 1;
                   end %increase counter by 1
                end    
            end
            if counter < currentPrecip %if the counter is smaller than currentTemp then the H-index is not valid, so we reduce it by one and repeat the loop.
                if useInches == 1
                    currentPrecip = currentPrecip - 0.01;
                else
                    currentPrecip = currentPrecip - 1 ;
                end
            end
        end
        if useInches == 1
            temporaryWIndex.WIndex(j-temporaryWIndex.YEAR(1)+1) = currentPrecip*100; %stores the max H-index value for each year in temporaryHIndex
        else
            temporaryWIndex.WIndex(j-temporaryWIndex.YEAR(1)+1) = currentPrecip;
        end
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
%%
WResults = table;
WResults.NAME = tableStationNames';
%for now, I'm only controlling the year from the H-Index section
% startYear = 1981;
% stopYear = 2013;
% useStartYear = 0;
% useStopYear = 0;
figure('Name', 'W-Index')
%formatSpec = "%f";
for i = 1:length(tableStationNames) %for each station in station names 
    %this section creates an array of subplots where each station has its
    %own plot showing the data points in orange, the change in time in
    %blue, and the general trend in black.
    tempNames = split(tableStationNames(i), '_'); %assigns a variable to the stations name and cuts out any unnecissary labels used for organization purposes
    A = array2table(WIndex(WIndex(:,1+i) ~= 0,[ 1 1+i])); %for each station, create a new array with the non-zero values and their corresponding years
    A.Properties.VariableNames{'Var1'} = 'YEAR'; 
    A.Properties.VariableNames{'Var2'} = 'WIndex';
    B = A(A.YEAR>=startYear & A.YEAR <=stopYear,:);
    
    subplot(4,6,i) %Creates a system of subplots in a 4x6 grid
    pos=get(gca,'Position');
    set(gca,'Position',[pos(1,1) pos(1,2) pos(1,3)+.03 pos(1,4)]) %subplot position
    
    lA = plot(A.YEAR,A.WIndex); %adds a line to the plot for additional clarity
    hold on %add each station to the same plot
    mdlA = fitlm(A, 'WIndex ~ YEAR'); %performs a linear regression for the Year and the H Index    
    zA = plot(mdlA); %plots the linear regression and data point
    mdlB = fitlm(B, 'WIndex ~ YEAR');
    zB = plot(mdlB);
    xA = A.YEAR;
    yA = A.WIndex;
    xB = B.YEAR;
    yB = B.WIndex;  
      
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
    zB(1).Color = 'none';%'#D95319'; %sets data points to be orange
    zB(1).Marker = '.';
    zB(1).MarkerSize = 10;
    zB(2).Color = 'r';
    zB(2).LineWidth = 1;
    zB(3).Color = 'none';
    zB(4).Color = 'none';
    
    legend('off'); %hides the automatic legend generated by fitlm
    
    %determines if the scale of the graph needs to be changed based on the
    %size of the year
    
    %x tick and label adjustments
   if i == 19 || i == 20 || i == 21 || i == 22 || i == 23 
        xlabel('Year', 'FontSize', 11)
        set(gca,'Xtick',[1921 1951 1981 2011]);
    else
        xlabel('')
        xticks('')
   end
    
    if useStartYear == 1 && useStopYear == 0
        xlim([startYear 2014])
    elseif useStartYear == 1 && useStopYear == 1
        xlim([startYear (stopYear+1)])
    elseif useStartYear == 0 && useStopYear == 1
        xlim([1890 (stopYear+1)])
    else     
        xlim([1890 2014])
    end
    %y tick and label adjustments
    if i == 1 || i == 7 || i == 13 || i == 19
        ylabel('W - Index (days/year)', 'FontSize', 11)
        yticks('auto')
    else
        ylabel('')
        yticks('')
    end
    if useInches == 1
       ylim([15 46]) 
    else
       ylim([5 26])
    end
    
    %this section  creates a table statistically important values for both
    %the POR and the specified period
    
    %p = polyfit(x,y,1);
       %this code runs trends analysis on x and y as independant and depedant variables
    %it tests the hypothesis of no correlation against the alternative
    %hypothesis of a nonzero correlation. so if p value is smaller than 0.05,
    %we reject the hypothesis.
    %analysis for the whole period of record
    [tau,p1]=corr(xA,yA,'type','kendall'); %kendall method
    tau_pA(i,1:3)=[i,tau,p1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
    [rho,p2]=corr(xA,yA,'type','spearman');%spearman method
    rho_pA(i,1:3)=[i,rho,p2];
    [r,p3]=corr(xA,yA);%pearson (linear) method
    r_pA(i,1:3)=[i,r,p3]; %pearson(Least square method) method corrcoef(x,y);    
    
    if useInches == 1
        WResults.slopePOR(i) = round(table2array(mdlA.Coefficients(2,1)),3);
    else
        WResults.slopePOR(i) = round(table2array(mdlA.Coefficients(2,1)),3);%calls the slope given for the linear regression of the data using the fitlm functioN
    end%calls the slope given for the linear regression of the data using the fitlm function
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
    
     %Statistical analysis for a specified period denoted by either a B or SP
     %(Specified Period)
    [tau,p1]=corr(xB,yB,'type','kendall'); %kendall method
    tau_pB(i,1:3)=[i,tau,p1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
    [rho,p2]=corr(xB,yB,'type','spearman');%spearman method
    rho_pB(i,1:3)=[i,rho,p2];
    [r,p3]=corr(xB,yB);%pearson (linear) method
    r_pB(i,1:3)=[i,r,p3]; %pearson(Least square method) method corrcoef(x,y); 
    
    if useInches == 1
        WResults.slopeSP(i) = round(table2array(mdlB.Coefficients(2,1)),3, 'significant');
    else
        WResults.slopeSP(i) = round(table2array(mdlB.Coefficients(2,1)),3, 'significant');%calls the slope given for the linear regression of the data using the fitlm functioN
    end
    WResults.rSP(i) = r_pB(i,2);
    WResults.rSqrSP(i) = mdlB.Rsquared.Ordinary; %calls the r^2 value from the fitlm function and inputs it into a new table
    WResults.rPValueSP(i) = round(table2array(mdlB.Coefficients(2,4)),3); %calls the pValue from the fitlm function and inputs it into a new table
    WResults.rhoSP(i) = rho_pB(i,2);
    WResults.rhoPValueSP(i) = rho_pB(i,3);
    WResults.tauSP(i) = tau_pB(i,2);
    WResults.tauPValueSP(i)= tau_pB(i,3);
    WResults.minSP(i) = min(B.WIndex);
    WResults.maxSP(i) = max(B.WIndex);
    WResults.meanSP(i) = round(mean(B.WIndex),1);
    WResults.medianSP(i) = median(B.WIndex);
    %num2str(100*WResults.slopePOR(23), "%#g")
    title(compose(tempNames(2,1)+"\n"+num2str(100*WResults.slopePOR(i),"%#.1f")+"\n"+num2str(100*WResults.slopeSP(i),"%#.1f")),'FontSize', 11);
end
%use the below script when you want to automatically make tiff files for
%the given graphs. See H-Index for issues with this function
%set(gcf,'PaperPositionMode','auto') %set the print area same as paper
%print('-dtiff','-r600', 'W_Index_Subplots')

