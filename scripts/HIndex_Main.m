%%
%This code is designed to take files from a specific folder, if they are
%excel files they are converted into csv files. those csv files are then
%used to produce an H-index based on temperature for each station.

%Note: you will need to download the "Weather folder" and add it to your
%directory file to ensure this program works. 
clear 
clc


%% 1 Universal Constants
folderName = 'Weather_CSV'; %variable for easy change of folder name
folderInfo = dir(folderName);  %creates a structure array with all the file names in "folderName"
folderLength = length(folderInfo); 
B = struct2cell(folderInfo);
for i = 3:folderLength
   stationNames(1,(i-2)) =  string(B(1,i));
   tableStationNames(1,(i-2))= erase(stationNames(1,(i-2)),".csv");
end
newFolder = strcat(folderName);

%constants for crop yield data
folderName2 = 'Crop_Yield'; %variable for easy change of folder name
folderInfo2 = dir(folderName2);  %creates a structure array with all the file names in "folderName"
folderLength2 = length(folderInfo2); 
B = struct2cell(folderInfo2);
for i = 3:folderLength2
   names(1,(i-2)) =  string(B(1,i));
end
newFolder2 = strcat(folderName2);

%creates number arrays and tables for the various decades and climate zones
%used in Kansas for this study
decades = cell2table({1911:1920 1921:1930 1931:1940 1941:1950 1951:1960 1961:1970 1971:1980 1981:1990 1991:2000 2001:2010}'); %creates a table for each decade between 1911 and 2010
west = [10, 20, 30];
EWcentral = [40, 50, 60];
east = [70, 80, 90];

north = [10, 40, 70];
NScentral = [20, 50, 80];
south = [30, 60, 90];
whole = 10:10:90;
agClimateZone = {west, EWcentral, east, north, NScentral, south, whole};
%% 2 Calculates the decadal averages for grain yield
folder2 = strcat(pwd,'/',newFolder2);

for i = 5:5%i = 1:length(names)
    baseFileName2 = names(i); 
    fullFileName2 = fullfile(folder2, baseFileName2); 
    temporaryFile2 = readtable(fullFileName2);
    %temporaryFile2.Year = cell2mat(temporaryFile2.Year);
    %temporaryFile2.Value = cell2mat(temporaryFile2.Value);

    Yield = table();
    Yield.YEAR = temporaryFile2.Year(temporaryFile2.Year >= 1911 & temporaryFile2.Year <= 2010); 
    Yield.REGION = temporaryFile2.AgDistrictCode(temporaryFile2.Year >= 1911 & temporaryFile2.Year <= 2010); 
    Yield.YIELD = round(temporaryFile2.Value(temporaryFile2.Year >= 1911 & temporaryFile2.Year <= 2010),1);

    base = mean(Yield.YIELD);
    yDecades = decades;
    %tis loop calculates the average H/L index for each decade from 1911-2010
    for j = 1:height(yDecades) 
        for h = 1:7
                z = cell2mat(agClimateZone(h));
                tempYield = Yield;
                if h ~= 7
                    tempYield = tempYield(tempYield.REGION == z(1) | tempYield.REGION == z(2) | tempYield.REGION == z(3),:);
                end
            yDecades.MEAN(j,h) = mean(tempYield.YIELD(tempYield.YEAR <= max(yDecades.Var1(j,:)) & tempYield.YEAR >= min(yDecades.Var1(j,:))));            
            yDecades.DIFF(j,h) = yDecades.MEAN(j,h)-base;        
            decadeNames (j) = compose(num2str(min(yDecades.Var1(j,:))-1)+"s");            
        end
    end
    if i == 5
        wYield = Yield;
        wYDecades = yDecades;
    end
end


%figure('name','Decadal');
%bar(categorical(decadeNames),yDecades.MEAN); %creates a bar graph showing the difference between the decadal average and the centenial average of 1911-2010

%figure('name', 'Yearly')
%bar(categorical(wYield.YEAR),wYield.YIELD);

%% 2 Calculating H-index
%this will be the first steps of creating the C/H index. it will start by
%anazyling a single station on a yearly basis and produce a bar graph at
%the end
clc
tic
stationLength = length(stationNames);
%stationLength = 1;

west = [2,3,4,8,9,14,20,21]; 
EWcentral = [5,10,11,15,16,17];
east = [6,7,12,13,18,19,22,23,24];

north = [2,3,4,5,6,7,11,12];
NScentral = [8,9,10,13,17];
south = [14,15,16,18,19,20,21,22,23,24];
whole = 2:24;
climateZone = {west, EWcentral, east, north, NScentral, south,whole};
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
%% 3 graphing H-index time series
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
    lA.LineWidth = .5;
    
    hold on %add each station to the same plot    
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
%year from 1911-2010 
y1 = 1911;
y2 = 2010;

if useLIndex == 1
    avgLIndex = table();
    avgWYield = table();
    avgWYield.YEAR = zeros(82,1);
    avgWYield.YIELD = zeros(82,7);
else
    avgHIndex = table();
    avgWYield = table();
    avgWYield.YEAR = zeros(82,1);
    avgWYield.YIELD = zeros(82,7);
end

for i = 1:length(climateZone)

    z = cell2mat(agClimateZone(i));
    tempWYield = wYield;
    if i ~= 7
        tempWYield = tempWYield(tempWYield.REGION == z(1) | tempWYield.REGION == z(2) | tempWYield.REGION == z(3),:);
    end
    
     
    for j = min(tempWYield.YEAR):max(tempWYield.YEAR)
        avgWYield.YEAR(1-min(tempWYield.YEAR)+j) = j;
        avgWYield.YIELD(1-min(tempWYield.YEAR)+j,i) = mean(tempWYield.YIELD(tempWYield.YEAR == j));
    end
    


    if useLIndex == 1       
        avgLIndex.YEAR = LIndex(:,1);
        avgLIndex.LIndex(:,i) = mean(LIndex(:,cell2mat(climateZone(i))),2); 
    else
        
        avgHIndex.YEAR = HIndex(:,1);
        avgHIndex.HIndex(:,i) = mean(HIndex(:,cell2mat(climateZone(i))),2);       
    end
end 

if useLIndex == 1
    avgLIndex = avgLIndex(avgLIndex.YEAR >= y1 & avgLIndex.YEAR <= y2,:);
else
    avgHIndex = avgHIndex(avgHIndex.YEAR >= y1 & avgHIndex.YEAR <= y2,:);
end
%% 4.0.3 Graphing time series vs. Yearly Average Grain Yield
locations = ["western", "east-west central", "eastern", "northern", "north-south central", "southern", "statewide"];
figure('Name', "TIndex vs Yield, EW")
for i = [1 2 3 7]
    if i ~= 7
        subplot(2,2,i)
    else
        subplot(2,2,4)
    end

    yyaxis left
    b = bar(avgWYield.YEAR,avgWYield.YIELD(:,i),1);
    hold on

    mdlC = fitlm(avgWYield.YEAR,avgWYield.YIELD(:,i));
    zC = plot(mdlC);

    zC(1).Color = 'none'; 
    zC(1).Marker = '.';
    zC(1).MarkerSize = 10;
    zC(2).Color = 'c';
    zC(2).LineWidth = 1.15;
    zC(3).Color = 'none';
    zC(4).Color = 'none'; 
    
    ylabel('Average Wheat Yield (Bu/acre)')
        
    yyaxis right
    if useLIndex ==  1
        p = plot(avgLIndex.YEAR, avgLIndex.LIndex(:,i));
        ylim([57 68])
    else
        p = plot(avgHIndex.YEAR, avgHIndex.HIndex(:,i));
        ylim([80 94])
    end
      
    xlim([1910 2011])
    hold on

    p.LineWidth = 1.15;
    p.Color = 'r';
    if useLIndex == 1
        mdlD = fitlm(avgLIndex.YEAR,avgLIndex.LIndex(:,i));   
        zD = plot(mdlD);
    else
        mdlD = fitlm(avgHIndex.YEAR,avgHIndex.HIndex(:,i));   
        zD = plot(mdlD); 
    end

    zD(1).Color = 'none'; 
    zD(1).Marker = '.';
    zD(1).MarkerSize = 10;
    zD(2).Color = 'g';
    zD(2).LineWidth = 1.15;
    zD(3).Color = 'none';
    zD(4).Color = 'none'; 

    if useLIndex == 1
        ylabel("Average L - Index (days/year)")
        title(compose(locations(i) + " average L-Index and average wheat yield vs. year"))
    else
        ylabel("Average H - Index (days/year)")
        title(compose(locations(i) +" average H-Index and average wheat yield vs. year"))
    end 
    legend off
    xlabel('Year')
end


figure('Name', "TIndex vs Yield, NS")
for i = [4 5 6 7]  
    
    subplot(2,2,i-3)    
    
    yyaxis left
    b = bar(avgWYield.YEAR,avgWYield.YIELD(:,i),1);
    hold on
    ylim([0 60])
  
    mdlC = fitlm(avgWYield.YEAR,avgWYield.YIELD(:,i));
    zC = plot(mdlC);

    zC(1).Color = 'none'; 
    zC(1).Marker = '.';
    zC(1).MarkerSize = 10;
    zC(2).Color = 'c';
    zC(2).LineWidth = 1.15;
    zC(3).Color = 'none';
    zC(4).Color = 'none'; 
    ylabel('Average Wheat Yield (Bu/acre)')
    
    yyaxis right
     if useLIndex ==  1
         p = plot(avgLIndex.YEAR, avgLIndex.LIndex(:,i));
         ylim([57 68])
     else
         p = plot(avgHIndex.YEAR, avgHIndex.HIndex(:,i));
         ylim([80 94])
     end
    
    xlim([1910 2011])
    hold on
    
    p.LineWidth = 1.15;
    p.Color = 'r';
    if useLIndex == 1
        mdlD = fitlm(avgLIndex.YEAR,avgLIndex.LIndex(:,i));      
        zD = plot(mdlD);
    else
        mdlD = fitlm(avgHIndex.YEAR,avgHIndex.HIndex(:,i));   
        zD = plot(mdlD); 
    end

    zD(1).Color = 'none'; 
    zD(1).Marker = '.';
    zD(1).MarkerSize = 10;
    zD(2).Color = 'g';
    zD(2).LineWidth = 1.15;
    zD(3).Color = 'none';
    zD(4).Color = 'none'; 

    if useLIndex == 1
        ylabel("Average L - Index (days/year)")
        title(compose(locations(i) + " average L-Index and average wheat yield vs. year"))
    else
        ylabel("Average H - Index (days/year)")
        title(compose(locations(i) +" average H-Index and average wheat yield vs. year"))
    end 
    legend off
    xlabel('Year')
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
    [tau,p1]=corr(xA,yA,'type','kendall'); %kendall method
    tau_pA(i,1:3)=[i,tau,p1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
    [rho,p2]=corr(xA,yA,'type','spearman');%spearman method
    rho_pA(i,1:3)=[i,rho,p2];
    [r,p3]=corr(xA,yA);%pearson (linear) method
    r_pA(i,1:3)=[i,r,p3]; %pearson(Least square method) method corrcoef(x,y);
    
    [tau,p1]=corr(xA2,yA2,'type','kendall'); %kendall method
    tau_pA2(i,1:3)=[i,tau,p1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
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
    dryness2 = 0.04; %dryness limit in inches (converting there are 0.0394 inches per mm)
    for j = temporaryDIndex.YEAR(1):temporaryDIndex.YEAR(end)%for each year at this station
            year = temporaryFile(temporaryFile.YEAR==j,:); %locates the index values for the given year and creates a temporary matrix for the given year       
            currentDryPeriod = 0;
            tempDryPeriod = 0;
            counter = 0;
            counter2 = 0;
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
            
            while counter2 < currentDryPeriod
                counter2 = 0;
                tempDryPeriod = 0;
                for h = 1:height(year)% for days in this year
                    
                    if year.RAIN(h) < dryness
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
    [tau,p1]=corr(xA,yA,'type','kendall'); %kendall method
    tau_pA(i,1:3)=[i,tau,p1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
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