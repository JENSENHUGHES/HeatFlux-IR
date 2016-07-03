function TemperatureFilter1 = FilterFluxData(TemperatureFull)

h = waitbar(0,'Filtering Temperature Data...');
TSize = size(TemperatureFull);
HeightPixels = TSize(1);
WidthPixels = TSize(2);
Frames = TSize(3);

Pixels = HeightPixels*WidthPixels;
TemperatureFilter1 = zeros(HeightPixels,WidthPixels,Frames);

for ii = 1:Frames
    waitbar(ii/Frames)
    if mod(ii,10) == 0
        %fprintf('Filtering Data for Timestep %3.0f of %3.0f\n',[ii Frames])
    end
    [~,noise] = wiener2(TemperatureFull(:,:,ii),[2 2]);
    J = wiener2(TemperatureFull(:,:,ii),[2 2],noise);
    [~,noise] = wiener2(J,[2 2]);
    TemperatureFilter1(:,:,ii) = wiener2(J,[2 2],noise);
end
close(h)
end