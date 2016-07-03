function TemperatureFilter1 = FilterTemperatureData(TemperatureFull)

h = waitbar(0,'Filtering Temperature Data...');
TSize = size(TemperatureFull);
HeightPixels = TSize(1);
WidthPixels = TSize(2);
Frames = TSize(3);
f = 7;
Pixels = HeightPixels*WidthPixels;
TemperatureFilter1 = zeros(HeightPixels,WidthPixels,Frames);

for ii = 1:Frames
    waitbar(ii/Frames)
    if mod(ii,10) == 0
        %fprintf('Filtering Data for Timestep %3.0f of %3.0f\n',[ii Frames])
    end
    [~,noise] = wiener2(TemperatureFull(:,:,ii),[f f]);
    J = wiener2(TemperatureFull(:,:,ii),[f f],noise);
    [~,noise] = wiener2(J,[f f]);
    TemperatureFilter1(:,:,ii) = wiener2(J,[f f],noise);
end
close(h)
end