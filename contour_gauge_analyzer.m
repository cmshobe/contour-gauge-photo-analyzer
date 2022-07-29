%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is a MATLAB script used to analyze photos of contour gauges to measure
%rock surface roughness.

%Written 2013-2014 by Charlie Shobe while at the College of William and Mary

%Please cite (and read!) the following publication if you use this code:

%Shobe, C.M., Hancock, G.S., Eppes, M.C., and Small, E.E. (2017) Field evidence for 
%the influence of weathering on rock erodibility and channel form in bedrock rivers,
%Earth Surface Processes and Landforms, v. 42, 1997-2012, doi:10.1002/esp.4163.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Additional details about the code can be found in Charlie Shobe's BS thesis:
%https://publish.wm.edu/honorstheses/66.

%Input: name of pin gage image-- MUST be in MATLAB directory
%See Fig. 4 in Shobe et al (2017) for an example of an appropriate image
%Note that this script was written explicitly for use with the 15 cm Fowler contour gauge;
%other gauges may have different lengths and pin widths, thereby requiring code changes.

%Outputs: plot of roughness indices, .csv file of plot data

clear all
close all

%% 1) Photo import and edge detection
photoName=input('type photo name with extension','s'); %Get photo name
photo=imread(photoName); %Read in photo
photoBW=rgb2gray(photo); %Convert photo to black and white for edge detection
bwEdge=edge(photoBW,'canny'); %Detect edges using Canny algorithm
imshow(bwEdge); %Show image on screen

%% 2) Tilt point collection and slope calculation
disp('Click the upper left and upper right points on the gage to calculate tilt')
[lineClickX,lineClickY]=ginput(2); %Get two points, one from each end of the gage, to calculate tilt
upperLeftAnchorX=lineClickX(1,1); %Assign line anchor points
upperLeftAnchorY=lineClickY(1,1); 
upperRightAnchorX=lineClickX(2,1); 
upperRightAnchorY=lineClickY(2,1);
slope=(upperRightAnchorY-upperLeftAnchorY)/(upperRightAnchorX-upperLeftAnchorX); %Calculate line slope
angleInRad=atan(slope); %Calculate rotation angle in radians

%% 3) Get user to click left-most and right-most pixels on edge line
dotPinL=input('Ready to click leftmost pin? (y=1/n=2)'); %Starting pixel: User zooms to region of interest, then types "1"
if dotPinL==1
    [pixXL,pixYL,pixColL]=impixel; %Gets x-coordinate, y-coordinate, and color value of clicked pixel
else
end
dotPinR=input('Ready to click rightmost pin? (y=1/n=2)'); %Ending pixel: %User zooms to region of interest, then types "1"
if dotPinR==1
    [pixXR,pixYR,pixColR]=impixel; %Gets x-coordinate, y-coordinate, and color value of clicked pixel
else
end

%% 4) Follow edge (white) line from left to right selected pixel
searchSwitch=1; %Flag to keep thw while loop going until the right pixel is reached. Changes to 0 when reached
startPixX=pixXL; %Sets starting pixel X value to the X value of the left pixel
startPixY=pixYL; %Sets starting pixel Y value to the Y value of the left pixel
pixCoords=[pixXL pixYL]; %Initiates the matrix "pixCoords," which stores pixel x and y coordinates. Currently not pre-allocated... SLOW
figure
imshow(bwEdge)
hold on
plot(pixXL,pixYL,'b.'); %Puts blue dot on leftmost pixel
plot(pixXR,pixYR,'b.'); %Puts blue dot on rightmost pixel
while searchSwitch==1
    %Set coordinates of options for next pixel for code to move to
    upRightPixX=startPixX+1; %Set X value for pixel diagonally above and to the right of selected pixel
    upRightPixY=startPixY-1; %Set Y value for pixel diagonally above and to the right of selected pixel
    rightPixX=startPixX+1; %Set X value for pixel to the right of selected pixel
    rightPixY=startPixY; %Set Y value for pixel to the right of selected pixel
    downRightPixX=startPixX+1; %Set X value for pixel diagonally below and to the right of selected pixel
    downRightPixY=startPixY+1; %Set Y value for pixel diagonally below and to the right of selected pixel
    twoUpOneRightPixX=startPixX+1; %Set X value for pixel two up and one to the right of selected pixel
    twoUpOneRightPixY=startPixY-2; %Set Y value for pixel two up and one to the right of selected pixel
    twoDownOneRightPixX=startPixX+1; %Set X value for pixel two down and one to the right of selected pixel
    twoDownOneRightPixY=startPixY+2; %Set Y value for pixel two down and one to the right of selected pixel
    
    %Check color values to see which pixels are white. These if-blocks are
    %in an order such that the lower blocks will overwrite the upper ones
    %if more than one pixel is white, meaning that I gave the lower blocks
    %a higher priority based on the structure of the Canny edges.
    if impixel(bwEdge,rightPixX,rightPixY)==1 %If pixel is white...
        startPixX=rightPixX; %Change startPixX to new X value
        startPixY=rightPixY; %Change startPixY to new Y value
    elseif impixel(bwEdge,downRightPixX,downRightPixY)==1
        startPixX=downRightPixX;
        startPixY=downRightPixY;
    elseif impixel(bwEdge,upRightPixX,upRightPixY)==1
        startPixX=upRightPixX;
        startPixY=upRightPixY;
    elseif impixel(bwEdge,twoDownOneRightPixX,twoDownOneRightPixY)==1
        startPixX=twoDownOneRightPixX;
        startPixY=twoDownOneRightPixY;
    elseif impixel(bwEdge,twoUpOneRightPixX,twoUpOneRightPixY)==1
        startPixX=twoUpOneRightPixX;
        startPixY=twoUpOneRightPixY;
    else %If no surrounding pixel is white...
        dotNewPix=input('Want to dot a supplemental pixel? (y=1/n=2)'); %User zooms to area near red dot, types "1", clicks on next white pixel on edge.
        if dotNewPix==1
            [NewPixX,NewPixY,NewPixCol]=impixel; %Get color of the chosen pixel
            startPixX=NewPixX; %Set startPix to newly clicked value to enable loop to continue
            startPixY=NewPixY;
        else
        end
    end

    plot(startPixX,startPixY,'r.'); %Put a red dot on the last white pixel found
    drawnow
    
    newCoordSet=[startPixX startPixY]; %Make a row of coordinates for the chosen pixel
    pixCoords=[pixCoords; newCoordSet]; %Append that row onto the pixel coordinates matrix
    
    %Check to see if rightmost pixel has been reached. If so, change
    %the value of searchSwitch to 0 to quit loop.
    if startPixX==pixXR && startPixY==pixYR
        searchSwitch=0;
        disp('Reached rightmost pixel')
    else
    end
end

%% 5) Rotate pixel coordinates based on tilt line slope
newCoords=zeros(numel(pixCoords(:,1)),2); %Pre-allocate rotated coordinate matrix
for coordsIndex=1:numel(pixCoords(:,1)) %Perform coordinate system rotation
    newCoords(coordsIndex,1)=cos(angleInRad)*pixCoords(coordsIndex,1)+sin(angleInRad)*pixCoords(coordsIndex,2);
    newCoords(coordsIndex,2)=-sin(angleInRad)*pixCoords(coordsIndex,1)+cos(angleInRad)*pixCoords(coordsIndex,2);
end

%% 6) Calculate mm length of each pixel, normalize pixel x-values to zero
pixelLength=newCoords(numel(newCoords(:,1)))-newCoords(1,1); %Get X length of pixel matrix
mmLength=148.5; %Measured length of pin gage
mmPerPixel=mmLength/pixelLength; %mm per pixel
pixPerPin=pixelLength/185; %Number of pixels in each pin
normCoordsX=newCoords(:,1)-newCoords(1,1); %normalize x values by subtracting smallest X from whole array
normCoordsY=newCoords(:,2); %Don't change Y array
normCoords=[normCoordsX normCoordsY]; %Build new normalized coordinates matrix
%first clicked pixel is now x=0, which is center of first pin. pin diameter
%is .7874mm, so distance between two pin centers is .7874mm.

%% 7) Sample 185 pixels at the correct interval to hit each pin
sampledMatrix=zeros(185,2); %Pre-allocate matrix of sample coordinate values
for i=0:184
    sampleID=pixPerPin*i; %returns value of pixels corresponding to each pin
    [~,index]=min(abs(normCoords(:,1)-sampleID)); %Find value in coordinates array closest to sampleID
    pinXVal=normCoords(index,1); %Get X value of closest row
    pinYVal=normCoords(index,2); %Get Y value of closest row
    pinValrow=[pinXVal pinYVal]; %Put X and Y values into a row
    sampledMatrix(i+1,:)=pinValrow; %Add row to pre-allocated matrix
end
sampledMatrix(:,2)=-sampledMatrix(:,2);
figure
plot(sampledMatrix(:,1),sampledMatrix(:,2),'.')
set(gcf,'Position',[500 500 1000 100])

%% 8) Calculate roughness indices for 9 pin windows
pt9coords=sampledMatrix;
pt9diffmat=zeros(numel(pt9coords(:,1))-8,2);
pt9diffmat(:,1)=1:numel(pt9diffmat(:,2));
for pt9index=5:numel(pt9coords(:,1))-4 %Iterate through coordinate pairs
    y1=pt9coords(pt9index-4,2); %Set first y-value
    y2=pt9coords(pt9index-3,2); %Set second y-value
    y3=pt9coords(pt9index-2,2); %Set third y-value
    y4=pt9coords(pt9index-1,2); %Set fourth y-value--this is the one being tested in each run of the loop
    y5=pt9coords(pt9index,2); %Set fifth y-value
    y6=pt9coords(pt9index+1,2); %Set sixth y-value
    y7=pt9coords(pt9index+2,2); %Set seventh y-value
    y8=pt9coords(pt9index+3,2); %Set eighth y-value
    y9=pt9coords(pt9index+4,2); %Set ninth y-value
    ymat=[y1 y2 y3 y4 y5 y6 y7 y8 y9]; %matrix of the y-value of each pin in the window
    ymean=mean(ymat); %get the average y-value for the window of pins
    stdevrough=sqrt(((y1-ymean)^2+(y2-ymean)^2+(y3-ymean)^2+(y4-ymean)^2+(y5-ymean)^2+(y6-ymean)^2+(y7-ymean)^2+(y8-ymean)^2+(y9-ymean)^2)/8);
    pt9diffmat(pt9index-3,1)=pt9coords(pt9index,1); %Populate matrix with x-coordinate, mean, and max for window
    pt9diffmat(pt9index-3,2)=stdevrough;
end
pt9diffmat(1,:)=[]; %Kill first row of difference matrix
pt9diffmax=max(pt9diffmat(:,2)); %Maximum single-window standard deviation
pt9diffmean=mean(pt9diffmat(:,2)); %Mean of window standard deviations
clear pt9index y1 y2 y3 y4 y5 y6 y7 y8 y9 ymat ymax ymean pt9diff ydiff1 ydiff2 ydiff3 ydiff4 ydiff5 ydiff6 ydiff7 ydiff8 %Clear excess variables

%% 9) Convert results, plot on bar graph, and save as .csv file
windowMat=pt9diffmean; %This is the index I have chosen to use
windowMatmm=windowMat*mmPerPixel; %Change from pixels to mm
maxMat=pt9diffmax; %Maximum difference
maxMatmm=maxMat*mmPerPixel; %Change from pixels to mm
windowMatPlot=[maxMatmm windowMatmm]; %These will be the bars in the bar graph
printfig=figure; %Give the bar graph a handle
windowBarGraph=bar(windowMatPlot); %Make bar graph
set(gca,'XTickLabel',{'Max Window-Averaged STDEV','Average STDEV'}) %Label x axis
ylim([0 1.5]) %Set y axis limit
title('Moving Window Roughness Analysis Results') %Graph title
xlabel('Metric') %X label
ylabel('Avg. Difference Between Set of Pins (mm)') %Y label
grid on %Turn gridlines on
legcolors=cell(1,2); %Make legend colors array
legcolors{1}='9 Pin Windows'; legcolors{2}='Garbage'; %Set legend color labels
legend(windowBarGraph,legcolors) %Make legend
save([photoName '_coords.mat'],'normCoords','sampledMatrix')
xlswrite([photoName '_stdev.xls'],windowMatPlot) %Write xls file (always errors out to csv...)
%bar graph
finalfig=figure; %Do a final plot of the image with the detected edge overlain
imshow(photo)
hold on
plot(pixCoords(:,1),pixCoords(:,2),'r.')
hold on
plot(sampledMatrix(:,1),sampledMatrix(:,2),'b.')
clear  labels %mmPerPixel