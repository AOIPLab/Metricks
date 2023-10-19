clear all;
close all;
clc;

% Found Here: https://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
% Adapted to accept our data by Mina Gaffney 
% Adapted on 11/03/2020
% Commented out Rob plotting 04/05/2021
% v2 - 02/17/2022
% Joe edited to make aspect ratio 1:1 & have axes limits be +/- 150 
% Joe edited for clean axes
% Joe removed Rob method to clean up a bit

% Load in a spreadsheet
[inname1, inpath1] = uigetfile({'*.xlsx';'*.csv';'*.xls'}, 'Please Select Input PCD Spreadsheet','MultiSelect','on');
MaxDensity = readmatrix(fullfile(inpath1,inname1));

[inname2, inpath2] = uigetfile({'*.xlsx';'*.csv';'*.xls'}, 'Please Select Input CDC Spreadsheet','MultiSelect','on');
IsoDensity95 = readmatrix(fullfile(inpath2,inname2));

globalscale = questdlg('Would you like output axes to be scaled globally?', ...
      'confirmation', ...
      'YES', 'NO', 'NO');
  if strcmpi(globalscale, 'NO')
      scaling = 0;
  elseif strcmpi(globalscale, 'YES')
      scaling = 1; 
  end

Outputfolder = fullfile(inpath1,'ErrorEllipseOutput');
mkdir(Outputfolder);

NumGraders = (size(MaxDensity,2)-1)/2; %Num graders is 1/2 the number of columns in the spreadsheet minus the first row (subject ID)

XdataStart = 2; 
XdataEnd = NumGraders +1; 
YdataStart = XdataEnd +1; 
YdataEnd = XdataEnd + NumGraders;

MaxX = max([MaxDensity(:,XdataStart:XdataEnd),IsoDensity95(:,XdataStart:XdataEnd)],[],'all');
MaxY = max([MaxDensity(:,YdataStart:YdataEnd),IsoDensity95(:,YdataStart:YdataEnd)],[],'all');
MinX = min([MaxDensity(:,XdataStart:XdataEnd),IsoDensity95(:,XdataStart:XdataEnd)],[],'all');
MinY = min([MaxDensity(:,YdataStart:YdataEnd),IsoDensity95(:,YdataStart:YdataEnd)],[],'all');

for i = 1:size(MaxDensity,1)
      
    SubjectID = MaxDensity(i,1);
    HoldingMaxDensity(:,1) = MaxDensity(i,XdataStart:XdataEnd)';
    HoldingMaxDensity(:,2) = MaxDensity(i,YdataStart:YdataEnd)';
    HoldingIsoDensity95(:,1) = IsoDensity95(i,XdataStart:XdataEnd)';
    HoldingIsoDensity95(:,2) = IsoDensity95(i,YdataStart:YdataEnd)';
    
% Calculate the eigenvectors and eigenvalues
covariance = cov(HoldingMaxDensity);
[eigenvec, eigenval ] = eig(covariance);
covariance2 = cov(HoldingIsoDensity95);
[eigenvec2, eigenval2 ] = eig(covariance2);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);
[largest_eigenvec_ind_c2, r2] = find(eigenval2 == max(max(eigenval2)));
largest_eigenvec2 = eigenvec2(:, largest_eigenvec_ind_c2);


% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));
largest_eigenval2 = max(max(eigenval2));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
end
if(largest_eigenvec_ind_c2 == 1)
    smallest_eigenval2 = max(eigenval2(:,2));
    smallest_eigenvec2 = eigenvec2(:,2);
else
    smallest_eigenval2 = max(eigenval2(:,1));
    smallest_eigenvec2 = eigenvec2(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));
angle2 = atan2(largest_eigenvec2(2), largest_eigenvec2(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end
if(angle2 < 0)
    angle2 = angle2 + 2*pi;
end

% Get the coordinates of the data mean
avg = mean(HoldingMaxDensity);
% avg = mean(HoldingMaxDensity)+ HoldingMeanDif;
avg2 = mean(HoldingIsoDensity95);

% Get the confidence interval error ellipse
% chi square table https://people.richland.edu/james/lecture/m170/tbl-chi.html
% For 2 degrees of freedom 1-0.05 = 0.95 --> 5.991 --> sqrt(5.991) =
% For 99% CI use sqrt(9.210) = 3.0348
% For 95% CI use sqrt(5.991) = 2.4477
% For 90% CI use sqrt(4.605) = 2.1459

chisquare_val = 2.4477; 
theta_grid = linspace(0,2*pi);
phi = angle;
X0=avg(1);
Y0=avg(2);
a=chisquare_val*sqrt(largest_eigenval); %major axis
b=chisquare_val*sqrt(smallest_eigenval); %minor axis
chisquare_val2 = 2.4477; %square root of 5.991
theta_grid2 = linspace(0,2*pi);
phi2 = angle2;
X02=avg2(1);
Y02=avg2(2);
a2=chisquare_val2*sqrt(largest_eigenval2); %major axis
b2=chisquare_val2*sqrt(smallest_eigenval2); %minor axis

Area(i,1) = SubjectID; 
Area(i,2) = pi*a*b;
Area(i,3) = pi*a2*b2; 

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );
ellipse_x_r2  = a2*cos( theta_grid2 );
ellipse_y_r2  = b2*sin( theta_grid2 );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
R2 = [ cos(phi2) sin(phi2); -sin(phi2) cos(phi2) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
r_ellipse2 = [ellipse_x_r2;ellipse_y_r2]' * R2;

% Draw the error ellipse
p1 = plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-b');
hold on;
p2 = plot(r_ellipse2(:,1) + X02,r_ellipse2(:,2) + Y02,'--m');

% Plot the original data
p3 = plot(HoldingMaxDensity(:,1), HoldingMaxDensity(:,2), 'ob','MarkerSize',4);
set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
p4 = plot(HoldingIsoDensity95(:,1), HoldingIsoDensity95(:,2), '.m','MarkerSize',6);
set(get(get(p4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
if scaling == 1;
    xlim([-150, 150]);
    ylim([-150, 150]);
end

% Set the axis labels
%hXLabel = xlabel('x'); %Comment out for clean figure
%hYLabel = ylabel('y'); %Comment out for clean figure
%set(gca,'xticklabel',[])
%set(gca,'yticklabel',[])
daspect([1 1 1]); %1:1 aspect ratio
hold off;
saveas(gcf,fullfile(Outputfolder, num2str(SubjectID)));
saveas(gcf,fullfile(Outputfolder, [num2str(SubjectID), '.tif']));
close; 


% Plot max density separately 
p5 = plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-b');
hold on;
p6 = plot(HoldingMaxDensity(:,1), HoldingMaxDensity(:,2), 'ob','MarkerSize',4);
set(get(get(p6,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

if scaling == 1;
    xlim([-150, 150]);
    ylim([-150, 150]);
% xlim([MinX-100, MaxX+100]);
% ylim([MinY-100, MaxY+100]);
end

% Set the axis labels
%hXLabel = xlabel('x'); %Comment out for clean figure
%hYLabel = ylabel('y'); %Comment out for clean figure
%set(gca,'xticklabel',[])
%set(gca,'yticklabel',[])
daspect([1 1 1]); %1:1 aspect ratio
hold off;
saveas(gcf,fullfile(Outputfolder, [num2str(SubjectID), '_95CI_PCD.fig']));
saveas(gcf,fullfile(Outputfolder, [num2str(SubjectID), '_95CI_PCD.tif']));
close; 


% Plot 95% Iso density
p7 = plot(r_ellipse2(:,1) + X02,r_ellipse2(:,2) + Y02,'--m');
hold on;
p8 = plot(HoldingIsoDensity95(:,1), HoldingIsoDensity95(:,2), '.m','MarkerSize',6);
set(get(get(p8,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

if scaling == 1;
    xlim([-150, 150]);
    ylim([-150, 150]);
end
% Set the axis labels
%hXLabel = xlabel('x'); %Comment out for clean figure
%hYLabel = ylabel('y'); %Comment out for clean figure
%set(gca,'xticklabel',[])
%set(gca,'yticklabel',[])

daspect([1 1 1]); %1:1 aspect ratio
hold off;
saveas(gcf,fullfile(Outputfolder, [num2str(SubjectID), '_95CI_CDC.fig']));
saveas(gcf,fullfile(Outputfolder, [num2str(SubjectID), '_95CI_CDC.tif']));
close;

clear HoldingMaxDensity; clear HoldingIsoDensity95;
end

AreaHeader = {'Subject ID', 'Max Density Area', 'Iso 95 Area'};
FinalAreaOut = [AreaHeader; num2cell(Area)]; 
writecell(FinalAreaOut, fullfile(inpath1,'ErrorEllipseAreas.csv'));