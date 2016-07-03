%FluxMapCalc-10-22

%INSTRUCTIONS:  FLIR output data is read in through .csv files for each
%frame of data.  Analysis parameters can either be input on each execution
%or entered directly into the code.  Per execution input is recommended for
%limited use studies.  Entering parameters directly in the code is
%recommended for multiple studies.  Code data will be utilized by entering
%'n' into command prompt when asked for manual input

%WARNINGS:
%Using code with ROI smaller than half the FLIR image is not recommended.
%    Random temperature error will not be filtered out for small ROIs
%Code will only work for ROIs with corner angles <180°

%Required Sub-routines
%   ReadTemperatureData_3_24_14
%   ReadTimeStampData_3_24_14
%   CompressTemperatureData
%   FilterTemperatureData
%   FilterTemperatureDifferential
%   RectTemperatureData

%Record of Revision********************************************************
%Date      Programmer          Description
%May 2012  C Rippe             Orginal Code
%Oct 2012  C Rippe             Code changed to ignore outter pixels and
%                              solving algorithm improved
%Mar 2013  C Rippe             1: Reorganized parameter input code to call
%                              FLIR file in individual code cell.  This
%                              allows for user to load matlab workspaces
%                              instead of re-importing data
%                              2: Added code to calculate temperature of
%                              different plate based on calculated flux
%Apr 2014  C Rippe             1: Changed FLIR file input code to skip
%                              files that are not found to allow for frame
%                              skipping when exporting frame from examinIR
%                              2: Reincorporatd outter pixels into
%                              calculation to allow for more accurate
%                              mapping when using DFLUX routine in Abaqus
%**************************************************************************
clear; clc;
%BASIC OUTPUT VARIABLES
% X -           "Horizontal" location of elements in flux matrices
% Y -           "Vertical" location of elements is flux matrices
% Flux -        Incident cold surface heat flux at each element
% CondFlux -    Conduction flux into each element
% ConvFlux -    Convection flux into each element
% RadFlux -     Radiation flux into each element
% StoreFulx -   Thermal Storage energy differential for each element

% CODE INPUT PARAMETERS****************************************************
File = 'SS304_FluxMap_40kW';  %Root name of data file
%Files should be named:
%(Root Name)_(Frame Number)
Frames = 291;                           %Total number of included frames
%Pixel locations of RoI [Horizontal Coordinate, Vertical Coordinate]
TopLeft = [1 1];
TopRight = [590 1];
BottomLeft = [1 480];
BottomRight = [640 480];
Compression = 0;

Width = .6096;                             %x-direction length (m)
Height = .6096;                            %y-direction length (m)
Thickness = .00079375;                   %plate thickness (m)
Properties = 'SS304';                    %plate material

Tinf = 298;                             %Air Temperature - K
h = 9.1;                                 %convection coefficient - W/(m2-K)
eps = 0.95;                             %emmissivity of the plate surface
filter = 'y';                           %Filter temperature matrix

% NO INPUT PARAMETERS BEYOND THIS POINT************************************

% ASK FOR INPUT PARAMETERS FROM USER FOR LOCATING FILES AND CHOOSING PROPERTIES
userin = input('Manually enter FLIR file paramters? (y/n) = ','s');
if userin == 'y'
    
    disp('File name structure should be [Root Name]_[Frame Number]')
    File = input('Root Name = ','s');
    Beginning = input('First Frame Number = ');
    Frames = input('Total Number of Frames = ');
    disp('Pixel locations should be entered in the form [ColumnID RowID]')
    disp('To use entire image, enter a single 0 for the TopLeft pixel')
    TopLeft = input('Top-Left Pixel of ROI = ');
    if length(TopLeft) ~= 0
        BottomLeft = input('Bottom-Left Pixel of ROI = ');
        TopRight = input('Top-Right Pixel of ROI = ');
        BottomRight = input('Bottom-Right Pixel of ROI = ');
    else
        TopLeft = [1 1];
        TopRight = [460 1];
        BottomLeft = [1 480];
        BottomRight = [640 480];
    end
end

userin = input('Manually enter plate parameters? (y/n) = ','s');
if userin == 'y'
    Width = input('Plate Width (m) = ');
    Height = input('Plate Height (m) = ');
    Thickness = input('Plate Thickness (m) = ');
    disp('Material Choices: AA6061 AA5083 SS304')
    Properties = input('Material = ','s');
end


userin = input('Manually enter analysis paramters? (y/n) = ','s');
if userin == 'y'
    Tinf = input('Ambient Temperature (°C) = ')+273;
    h = input('Convection Coefficient (W/m^2-K) = ');
    eps = input('Surface Emmissivity = ');
    filter = input('Filter Temperature Data? (y/n) = ','s');
    disp('Compression value of 0 will use original data')
    disp('Compression value of 1 will average 3x3 pixel blocks')
    disp('Compression value of 2 will average 5x5 pixel blocks, etc...')
    Compression = input('Compression Value = ');
end

% Analysis Constants
DataStart = [7 0];                      %Location in data start in CSV file
sigma = 5.67e-8;                        %steffan-boltzman constant - W/(m2-K4)

%Input material properties for chosen material
if strcmp(Properties,'AA6061') == 1
    k = 186;                                %plate conductivity - W/(m-K)
    cp = 1010;                               %plate specific heat - J/(kg-K)
    ro = 2690;                               %plate density = kg/m3
elseif strcmp(Properties,'AA5083') == 1
    k = 186;                                %plate conductivity - W/(m-K)
    cp = 1010;                               %plate specific heat - J/(kg-K)
    ro = 2660;                               %plate density = kg/m3
elseif strcmp(Properties,'SS304') == 1
    k = 16.2;                               %plate conductivity - W/(m-K)
    cp = 515;                               %plate specific heat - J/(kg-K)
    ro = 7900;                              %plate density - kg/m3
else
    error('Invalid material ID. Material Choices are:  AA6061 AA5083 SS304')
end

% READ ExaminIR GENERATED CSV FILES
%Read full temperature matrix
TemperatureFull = ReadTemperatureData_3_24_14(File,Frames,DataStart);
%Read Timestamp data from files
TimeData = ReadTimeStampData_3_24_14(File,Frames);

% PRECALUCLATE NECESSARY VARIABLES

TopLeft = fliplr(TopLeft);
TopRight = fliplr(TopRight);
BottomLeft = fliplr(BottomLeft);
BottomRight = fliplr(BottomRight);

FullSize = size(TemperatureFull);
FullWidthPixels = FullSize(2);
FullHeightPixels = FullSize(1);

clc;
if (TopLeft(1) ~= 0)
    TemperatureRect = RectTemperatureData(TemperatureFull,TopLeft,TopRight,BottomLeft,BottomRight);
else
    TemperatureRect = TemperatureFull;
end
    
if (strcmp(filter,'y')==1)
    TemperatureFiltered = FilterTemperatureData(TemperatureRect);
else
    TemperatureFiltered = TemperatureRect;
end

if Compression ~= 0
    TemperatureCompressed = CompressTemperatureData(TemperatureFiltered,Compression,Frames);
else
    TemperatureCompressed = TemperatureFiltered;
end


TemperatureSize = size(TemperatureCompressed);    %Variable for size of temperature matrix
Temperature = TemperatureCompressed+273;          %Convert Celcius temperature to Kelvin

%Input Temperature field parameters
WidthPixels = TemperatureSize(2);       %Number of pixels in x-direction
HeightPixels = TemperatureSize(1);      %Number of pixels in y-direction
Timesteps = Frames;                     %Total number of temperature datasets
Pixels = WidthPixels*HeightPixels;      %total number of pixels
DeltaX = Width/(WidthPixels-1);         %meters
DeltaY = Height/(HeightPixels-1);       %meters

%Allocate space for different pixel properties
Flux = zeros(HeightPixels,WidthPixels,Timesteps-1);     %Incident Flux
CondFlux = Flux;                                        %Conduction Flux
ConvFlux = Flux;                                        %Convection Flux
RadFlux = Flux;                                         %Radiation Flux
StoreFlux = Flux;                                       %Storage Equivalent Flux
TimeTDiff = Flux;
SpaceTDiff = Flux;

%Define display variables
XFull = linspace(0,Width,FullWidthPixels);
YFull = linspace(0,Height,FullHeightPixels);
X = linspace(0,Width,WidthPixels);
Y = linspace(0,Height,HeightPixels);

figure(1);
tempdisplay = flipud(TemperatureFull(:,:,end));
subplot(2,1,1); pcolor(XFull,YFull,tempdisplay); shading flat; colorbar;
title('Final Temperature Image'); xlabel('Approximate Location (m)'); ylabel('Approximate Location (m)');
colormin = min(min(Temperature(:,:,end),[],2),[],1)-273;
colormax = max(max(Temperature(:,:,end),[],2),[],1)-273;
axis([0 Width 0 Height -1 1 colormin colormax]);
tempdisplay = flipud(Temperature(:,:,end));
subplot(2,1,2); pcolor(X,Y,tempdisplay-273); shading flat; colorbar;
title('Rectified Temperature Image'); xlabel('Location (m)'); ylabel('Location (m)');
axis([0 Width 0 Height -1 1 colormin colormax]);

disp('Calculation paused for review. Press any key to continue')
pause

% SOLVE FOR INCIDENT HEAT FLUX
clc
wbar = waitbar(0,'Solving For Heat Flux Map...');
for j = 1:(Timesteps-1)                 %Cycle through time
    
    %Allocate space for temporary matrices used at this timestep
    tempflux = zeros(HeightPixels,WidthPixels); %incident heat flux at current timestep
    tempcondflux = tempflux;            %conduction flux at current timestep
    tempconvflux = tempflux;            %convection flux at current timestep
    tempradflux = tempflux;             %radiation flux at current timestep
    tempstoreflux = tempflux;           %storage eq. flux at current timestep
    temptimediff = tempflux;
    tempspacediff = tempflux;
    
    tempTemperature = Temperature(:,:,j);
    tempTemperatureFuture = Temperature(:,:,j+1);
    
    for i = 1:Pixels
        
        temptimediff(i) = tempTemperatureFuture(i)-tempTemperature(i);
        
        if i==1
            tempspacediff(i) = (tempTemperature(i+1)*DeltaX^2+tempTemperature(i+HeightPixels)*DeltaY^2-tempTemperature(i)*(DeltaX^2+DeltaY^2))/(DeltaX*DeltaY);
        elseif i==HeightPixels
            tempspacediff(i) = (tempTemperature(i-1)*DeltaX^2+tempTemperature(i+HeightPixels)*DeltaY^2-tempTemperature(i)*(DeltaX^2+DeltaY^2))/(DeltaX*DeltaY);
        elseif i==Pixels-HeightPixels+1
            tempspacediff(i) = (tempTemperature(i+1)*DeltaX^2+tempTemperature(i-HeightPixels)*DeltaY^2-tempTemperature(i)*(DeltaX^2+DeltaY^2))/(DeltaX*DeltaY);
        elseif i==Pixels
            tempspacediff(i) = (tempTemperature(i-1)*DeltaX^2+tempTemperature(i-HeightPixels)*DeltaY^2-tempTemperature(i)*(DeltaX^2+DeltaY^2))/(DeltaX*DeltaY);
        elseif i<HeightPixels
            tempspacediff(i) = (tempTemperature(i-1)*DeltaX^2+tempTemperature(i+1)*DeltaX^2+tempTemperature(i+HeightPixels)*DeltaY^2-tempTemperature(i)*(2*DeltaX^2+DeltaY^2))/(DeltaX*DeltaY);
        elseif i>Pixels-HeightPixels
            tempspacediff(i) = (tempTemperature(i-1)*DeltaX^2+tempTemperature(i+1)*DeltaX^2+tempTemperature(i-HeightPixels)*DeltaY^2-tempTemperature(i)*(2*DeltaX^2+DeltaY^2))/(DeltaX*DeltaY);
        elseif mod(i,HeightPixels)==0
            tempspacediff(i) = (tempTemperature(i-1)*DeltaX^2+tempTemperature(i+HeightPixels)*DeltaY^2+tempTemperature(i-HeightPixels)*DeltaY^2-tempTemperature(i)*(DeltaX^2+2*DeltaY^2))/(DeltaX*DeltaY);
        elseif mod(i,HeightPixels)==1
            tempspacediff(i) = (tempTemperature(i+1)*DeltaX^2+tempTemperature(i+HeightPixels)*DeltaY^2+tempTemperature(i-HeightPixels)*DeltaY^2-tempTemperature(i)*(DeltaX^2+2*DeltaY^2))/(DeltaX*DeltaY);
        else
            tempspacediff(i) = (tempTemperature(i-1)*DeltaX^2+tempTemperature(i+1)*DeltaX^2+tempTemperature(i-HeightPixels)*DeltaY^2+tempTemperature(i+HeightPixels)*DeltaY^2-tempTemperature(i)*(2*DeltaX^2+2*DeltaY^2))/(DeltaX*DeltaY);
        end
    end
    
    if sum(sum(temptimediff,1),2)==0
        temptimedifffilter = temptimediff;
    else
        temptimedifffilter = FilterTemperatureDifferential(temptimediff);
    end
    
    if sum(sum(tempspacediff,1),2)==0
        tempspacedifffilter = tempspacediff;
    else
        tempspacedifffilter = FilterTemperatureDifferential(tempspacediff);
    end
        
    for i = 1:Pixels
            %Calculate conduction flux at pixel
            tempcond = k*Thickness/(DeltaX*DeltaY)*tempspacedifffilter(i);
            %Caluclate convection flux at pixel
            %h = (0.0961*tempTemperature(i)-17.5797);
            tempconv = -h*(2)*(tempTemperature(i)-Tinf);
            %Calculate radiation flux at pixel
            temprad = -2*eps*sigma*(tempTemperature(i).^4-Tinf.^4);
            %Calculate energy storage term
            tempstored = ro*Thickness*cp*temptimedifffilter(i)/(TimeData(j+1)-TimeData(j));
            
            %Insert calculated pixel terms into global matrices for this
            %timestep
            tempflux(i) = (tempstored-tempcond-tempconv-temprad);
            tempcondflux(i) = tempcond;
            tempconvflux(i) = tempconv;
            tempradflux(i) = temprad;
            tempstoreflux(i) = tempstored;
    end %for i
    
    %Insert matrices for this timestep into full matrices containing all
    %fluxes for all times
    Flux(:,:,j) = tempflux;             %Incident Heat Flux (W/m^2)
    CondFlux(:,:,j) = tempcondflux;     %Conduction Flux between nodes (W/m^2)
    ConvFlux(:,:,j) = tempconvflux;     %Convective Losses on node (W/m^2)
    RadFlux(:,:,j) = tempradflux;       %Radiative Losses on node (W/m^2)
    StoreFlux(:,:,j) = tempstoreflux;   %Storage eq. flux
    TimeTDiff(:,:,j) = temptimediff;
    SpaceTDiff(:,:,j) = tempspacediff;
    waitbar(j/(Timesteps-1),wbar)
end %for j
close(wbar);