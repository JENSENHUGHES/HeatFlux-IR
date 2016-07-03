function Temperature = CompressTemperatureData(TemperatureCropped,Compression, Frames)

%Function will compress the Temperature data by averaging temperatures of
%particular nodes within the compression distance within the matrix

h = waitbar(0,'Compressing Temperature Data...');

division = 2*Compression+1;
TSize = size(TemperatureCropped);
TempX = TSize(1)-mod(TSize(1),division);
TempY = TSize(2)-mod(TSize(2),division);

TemperatureAdjusted = TemperatureCropped(1:TempX,1:TempY,:);
TempSize = size(TemperatureAdjusted);
NewSize = size(TemperatureAdjusted)/division;

Temperature = zeros(NewSize(1),NewSize(2),Frames);

for i = 1:NewSize(2)
    %fprintf('Calculated compressed data for pixel row %1.0f of %1.0f\n',i,NewSize(2))
    for j = 1:NewSize(1)
        TempSet = TemperatureCropped((j-1)*division+1:j*division,(i-1)*division+1:i*division,:);
        for m = 1:Frames
            Temperature(j,i,m) = mean(mean(TempSet(:,:,m),1),2);
        end
    end
    waitbar(i/NewSize(2),h);
end
close(h);