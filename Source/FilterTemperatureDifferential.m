function TemperatureFilter1 = FilterTemperatureDifferential(TemperatureFull)

TSize = size(TemperatureFull);
HeightPixels = TSize(1);
WidthPixels = TSize(2);
f = 11;
Pixels = HeightPixels*WidthPixels;
TemperatureFilter1 = zeros(HeightPixels,WidthPixels);

    [~,noise] = wiener2(TemperatureFull,[f f]);
    J = wiener2(TemperatureFull,[f f],noise);
    [~,noise] = wiener2(J,[f f]);
    TemperatureFilter1 = wiener2(J,[f f],noise);
end