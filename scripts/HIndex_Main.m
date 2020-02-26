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
%% 2
%this will be the first steps of creating the C/H index. it will start by
%anazyling a single station on a yearly basis and produce a bar graph at
%the end
clc
tic
stationLength = length(stationNames);
%stationLength = 1;

%for if you want to calculate the L-index(using Tmin) instead of the H-index (using Tmax)
useLIndex = 0; 

folder = strcat(pwd,'/',newFolder); %calls the path of the current file directory
HIndex = zeros(2013-1890+1,23);
HIndex(:,1) = (1890:2013)';

if useLIndex == 1
    LIndex = HIndex; %creates a seperate variable to store the results for L-index if it is being used
end

for i = 1:stationLength %for each station
    baseFileName = stationNames(i); %this is the name of the file excluding file type. 
    fullFileName = fullfile(folder, baseFileName); %creates a variable for the full file path to ensure no errors related to file path
    temporaryFile = readtable(fullFileName); %creates a temporary matrix of the the data for the current station name.
    
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
            if counter < currentTemp %if the counter is smaller than currentTemp then the H-index is not valid, so we reduce it by one and repeat the loop.
                currentTemp = currentTemp - 1;
            end
        end              
        temporaryHIndex.HIndex(j-temporaryHIndex.YEAR(1)+1) = currentTemp; 
    end
    
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
timeHindex = toc
%%
%creates a table to store all statistical values for all results for a
%given Index
HResults = table;
HResults.NAME = tableStationNames';

startYear = 1981;
stopYear = 2010;
useStartYear = 0;
useStopYear = 0;
figure('Name', 'H-Index')


for i = 1:length(tableStationNames) %for each station in station names 
    %this section creates an array of subplots where each station has its
    %own plot showing the data points in orange, the change in time in
    %blue, and the general trend in black.
    tempNames = split(tableStationNames(i), '_'); %assigns a variable to the stations name and cuts out any unnecissary labels used for organization purposes
    HResults.NAME(i) = tempNames(2,1);
    
    if useLIndex == 1
        A = array2table(LIndex(LIndex(:,1+i) ~= 0,[ 1 1+i]));
    else 
        A = array2table(HIndex(HIndex(:,1+i) ~= 0,[ 1 1+i])); %for each station, create a new array with the non-zero values and their corresponding years
    end
    
    A.Properties.VariableNames{'Var1'} = 'YEAR'; 
    A.Properties.VariableNames{'Var2'} = 'HIndex';  
    B = A(A.YEAR>=startYear & A.YEAR <=stopYear,:);  
    
    subplot(4,6,i)
    pos=get(gca,'Position');
    set(gca,'Position',[pos(1,1) pos(1,2) pos(1,3)+.03 pos(1,4)]) %subplot position
    
   
    lA = plot(A.YEAR,A.HIndex); %adds a line to the plot for additional clarity
    hold on %add each station to the same plot    
    mdlA = fitlm(A, 'HIndex ~ YEAR'); %performs a linear regression for the Year and the H Index    
    zA = plot(mdlA); %plots the linear regression and data point    
    
    xA = A.YEAR;
    yA = A.HIndex;
    
    meanB = B;
    meanB.HIndex(:) = mean(B.HIndex);
    zB = plot(meanB.YEAR, meanB.HIndex);
    zB.Color = 'r';  
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
    [tau,p1]=corr(xA,yA,'type','kendall'); %kendall method
    tau_pA(i,1:3)=[i,tau,p1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
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
%%
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

%%
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
end

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
    
    [tau,p1]=corr(xA,yA,'type','kendall'); %kendall method
    tau_pA(i,1:3)=[i,tau,p1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
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

%%
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
%%
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
