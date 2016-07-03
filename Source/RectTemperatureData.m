function TemperatureRect = RectTemperatureData(TemperatureFull,TopLeft,TopRight,BottomLeft,BottomRight)

%Function will turn any 4 sided ROI into rectangular matrix for processing
%using heat flux code.  Function works using linear interpolation between
%provided temperature matrix values.  WARNING:  Heat flux code will only
%work for ROIs with corner angles <180°

%Verify that provided corner coordinates are valid
if BottomLeft(1)-TopLeft(1)<=0 | BottomRight(1)-TopRight(1)<=0 | TopRight(2)-TopLeft(2)<=0 | BottomRight(2)-BottomLeft(2)<=0
    error('Provided corner coordinates not suitable for rectification. Ensure coordinate values correspond to correct locations. (i.e. top-right coordinate not further left than top-left coordinate)')
else
    %Continue with rectification
    
    %Calculate height and width of ROI
    HeightLeft = BottomLeft(1)-TopLeft(1);
    HeightRight = BottomRight(1)-TopRight(1);
    WidthTop = TopRight(2)-TopLeft(2);
    WidthBottom = BottomRight(2)-BottomLeft(2);
    
    TSize = size(TemperatureFull);
    TemperatureRect = zeros(TSize(1),TSize(2),TSize(3));
    
    h = waitbar(0,'Rectifying Temperature Data...');

    %loop through each value of rectified temperature matrix
    for mm = 1:TSize(3)
        %fprintf('Rectifying Frame %3.0f (%3.0f of %3.0f)\n',[mm mm TSize(3)])
        for ii = 1:TSize(1)
            fii = (ii-1)/(TSize(1)-1);
            for jj = 1:TSize(2)
                fjj = (jj-1)/(TSize(2)-1);
                %Calculate height and width fraction for particular pixel
                Height = fjj*(HeightRight-HeightLeft)+HeightLeft;
                Width = fii*(WidthBottom-WidthTop)+WidthTop;
                TopStart = fjj*(TopRight(1)-TopLeft(1))+TopLeft(1);
                LeftStart = fii*(BottomLeft(2)-TopLeft(2))+TopLeft(2);
                HeightPos = fii*Height+TopStart;
                WidthPos = fjj*Width+LeftStart;
                fh = mod(HeightPos,1);
                fw = mod(WidthPos,1);
                
                %Determine closest lower pixel from original temperature
                %matrix
                if mod(HeightPos,1) < 0.5
                    Top = int16(HeightPos);
                else
                    Top = int16(HeightPos)-1;
                end
                
                if mod(WidthPos,1) < 0.5
                    Left = int16(WidthPos);
                else
                    Left = int16(WidthPos)-1;
                end
                
                %Check for calculation of last pixel row and change to
                %access previous pixel grouping for interpolation
                if Top == TSize(1)
                    Top = Top-1;
                    fh = 1;
                end
                if Left == TSize(2)
                    Left = Left-1;
                    fw = 1;
                end
                
                %Calculate rectified pixel temperature using linear
                %interpolation
                Ttop = fw*(TemperatureFull(Top,Left+1,mm)-TemperatureFull(Top,Left,mm))+TemperatureFull(Top,Left,mm);
                Tbot = fw*(TemperatureFull(Top+1,Left+1,mm)-TemperatureFull(Top+1,Left,mm))+TemperatureFull(Top+1,Left,mm);
                TemperatureRect(ii,jj,mm) = fh*(Tbot-Ttop)+Ttop;
                
                %Display values for code validation
%                 if mm == 1 && ii == 479 && jj == 639
%                     fprintf('fii = %3.3f\n',fii)
%                     fprintf('fjj = %3.3f\n',fjj)
%                     fprintf('Height = %3.2f\n',Height)
%                     fprintf('Width = %3.2f\n',Width)
%                     fprintf('TopStart = %3.2f\n',TopStart)
%                     fprintf('LeftStart = %3.2f\n',LeftStart)
%                     fprintf('HeightPos = %3.2f\n',HeightPos)
%                     fprintf('WidthPos = %3.2f\n',WidthPos)
%                     fprintf('fh = %3.3f\n',fh)
%                     fprintf('fw = %3.3f\n',fw)
%                     fprintf('Top = %2.0f\n',Top)
%                     fprintf('Left = %2.0f\n',Left)
%                     fprintf('Ttop = %3.1f\n',Ttop)
%                     fprintf('Tbot = %3.1f\n',Tbot)
%                     fprintf('Trect = %3.1f\n',TemperatureRect(ii,jj,mm));
%                 end
            end %for jj
        end %for ii
        waitbar(mm/TSize(3),h);
    end %for mm
    close(h)
end %if