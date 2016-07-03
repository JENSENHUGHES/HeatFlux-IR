function TemperatureFull = ReadTemperatureData_3_24_14(File, Frames, DataStart)

%Function will read a series of csv files that contain a matrix of
%temperatures from a FLIR camera.  Each file contains a 2D matrix of
%temperatures for a particular time.  A series of files corresponds to
%multiple times.  In this code ii is the number of included frames while jj
%is the total number of FLIR images.  If images are exported with skipped
%frames, ii and jj will be different (example:  if you have 1000 FLIR
%images but only export every other image, ii will end at 500 while jj
%reaches 1000.
h = waitbar(0,'Reading Temperature Data Files...');
%Determine 1st file name and existance
filename = sprintf('%s_%1.0f.csv',File,1);
existance = exist(filename);
jj = 1;

%Loop through all frames
for ii = 1:Frames
    status = sprintf('Reading Temperature Data Files...%2.0f%%',ii/Frames*100);
    waitbar(ii/Frames,h,status)
    %Determine existance of file
    kk = 0;
    while existance == 0;
        jj = jj + 1;
        filename = sprintf('%s_%1.0f.csv',File,jj);
        existance = exist(filename);
        if kk > 100
            warning('Over 100 frames skipped:  Check total number of frames')
            kk = 0;
        end
        kk = kk + 1;
    end
    %read files that exist
%    fprintf('Reading File:  %s_%1.0f.csv (%1.0f of %1.0f)\n',File,ii,ii,Frames)
    TemperatureFull(:,:,ii) = csvread(filename,DataStart(1)-1,DataStart(2));
    
    jj = jj + 1;
    filename = sprintf('%s_%1.0f.csv',File,jj);
    existance = exist(filename);
end
close(h)
%LEGACY CODE
% for i = 1:Frames-1
%     filename = sprintf('%s_%1.0f.csv',File,i*DropFrame+Beginning);
%     fprintf('Reading File:  %s_%1.0f.csv (%1.0f of %1.0f)\n',File,i*DropFrame+Beginning,i+1,Frames)
%     TemperatureFull(:,:,i+1) = csvread(filename,DataStart(1)-1,DataStart(2));
% end