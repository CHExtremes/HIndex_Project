%% Trend Analysis
%this code runs trends analysis on x and y as independant and depedant variables
%it tests the hypothesis of no correlation against the alternative
%hypothesis of a nonzero correlation. so if p value is smaller than 0.05,
%we reject the hypothesis.
[tou,p1]=corr(x,y,'type','kendall'); %kendall method
tou_p(j,1:3)=[j,tou,p1];%j is stations number in the loop; tou is kendall tou value; and p1 is the p-value for the test.
[rho,p2]=corr(x,y,'type','spearman');%spearman method
rho_p(j,1:3)=[j,rho,p2];
[r,p3]=corr(x,y);%pearson (linear) method
r_p(j,1:3)=[j,r,p3]; %pearson(Least square method) method corrcoef(x,y);