%this function adds lines infil
function [Supports_Points] = Add_Support_Infill(polyout, PrintParameters, BBox)


% define space between each line
line_spacing = PrintParameters.Support_Spacing;
segment_length = 0.2* line_spacing * sind(90-PrintParameters.Support_Angle);
segment_break = 0.7* line_spacing * sind(90-PrintParameters.Support_Angle);


% calculate distance along for each line
adjacent = abs(BBox.MaxY/tand(PrintParameters.Support_Angle));
if isinf(adjacent)
    adjacent = 0;
end

Start_list = [];
End_list = [];


%create infill lists of start and end points
if PrintParameters.Support_Angle == 0
    Start_list (:,2) = (BBox.MinY:PrintParameters.Spacing:BBox.MaxY)'; %y coordinate
    Start_list (:,1) = BBox.MinX;  %x coordinate
    End_list (:,2) = Start_list (:,2);
    End_list (:,1) = BBox.MaxX;
    
elseif (PrintParameters.Support_Angle > 0) && (PrintParameters.Support_Angle <= 90)
    count1 = 0;
    for start_ponit_x = BBox.MinX-adjacent:(line_spacing/sind(PrintParameters.Support_Angle)):BBox.MaxX
        Start_list_x1 = [];
        Start_list_y1 = [];

        Start_list_x2 = [];
        Start_list_y2 = [];
        %offset lines to avoid continuous gaps
        if ~mod(count1,2)
            Start_list_x1 = [Start_list_x1;(start_ponit_x:(segment_length+segment_break)*cosd(PrintParameters.Support_Angle):BBox.MaxX)'];        
            Start_list_y1 = tand(PrintParameters.Support_Angle)*(Start_list_x1-start_ponit_x);

            Start_list_x2 = Start_list_x1 + PrintParameters.Spacing*sind(PrintParameters.Support_Angle);
            Start_list_y2 = Start_list_y1 - PrintParameters.Spacing*cosd(PrintParameters.Support_Angle);
            Start_list = [Start_list; Start_list_x1, Start_list_y1; Start_list_x2, Start_list_y2];  
            count1 = count1 +1;
        else
            Start_list_x1 = [Start_list_x1;(start_ponit_x:(segment_length+segment_break)*cosd(PrintParameters.Support_Angle):BBox.MaxX)'];        
            Start_list_y1 =  tand(PrintParameters.Support_Angle)*(Start_list_x1-start_ponit_x) -sind(PrintParameters.Support_Angle) *segment_break; 
            
            Start_list_x2 = Start_list_x1 - PrintParameters.Spacing*sind(PrintParameters.Support_Angle);
            Start_list_y2 = Start_list_y1 + PrintParameters.Spacing*cosd(PrintParameters.Support_Angle);
            Start_list = [Start_list; Start_list_x1, Start_list_y1; Start_list_x2, Start_list_y2];  
            count1 = count1 +1;
        end 
    end  
    End_list = [Start_list(:,1) + segment_length*cosd(PrintParameters.Support_Angle), Start_list(:,2) + segment_length*sind(PrintParameters.Support_Angle)];        

elseif PrintParameters.Support_Angle >= 90
    Start_list (:,1) = (BBox.MinX:line_spacing:BBox.MaxX+adjacent)'; %x coordinate
    Start_list (:,2) = BBox.MinY; %y coordinate
    End_list (:,1) = Start_list (:,1) - adjacent;
    End_list (:,2) = BBox.MaxY;
end


Infill_lines1 = []; %initialize
check = 1; %switch for infill pattern

for k = 1:size(Start_list,1)  
    %find intersection points of infill pattern.
    [in,out] = intersect(polyout, [Start_list(k,:); End_list(k,:)]);
    out;
    if ~isempty(in)
        % flip direction to minimise air time
        if check == 1
            Infill_lines1 = [Infill_lines1; in; nan nan];
        else
            Infill_lines1 = [Infill_lines1; flip(in,1); nan, nan];
        end
        check = check * -1;
    end
end

%check for starting error
if isnan(Infill_lines1(1,1))
    Infill_lines1 = Infill_lines1(2:end, :);
end
Infill_lines1;
%% repeat at 90'
clear adjacent Start_list End_list start_ponit_x
% calculate distance along for each line
adjacent = abs(BBox.MaxY*tand(PrintParameters.Support_Angle));
segment_length = 0.4* line_spacing * sind(90-PrintParameters.Support_Angle);
segment_break = 0.5* line_spacing * sind(90-PrintParameters.Support_Angle);

if isinf(adjacent)
    adjacent = 0;
end

Start_list = [];
End_list = [];

%create infill lists of start and end points
if PrintParameters.Support_Angle == 0
    Start_list (:,2) = (BBox.MinY:PrintParameters.Spacing:BBox.MaxY)'; %y coordinate
    Start_list (:,1) = BBox.MinX;  %x coordinate
    End_list (:,2) = Start_list (:,2);
    End_list (:,1) = BBox.MaxX;
    
elseif (PrintParameters.Support_Angle > 0) && (PrintParameters.Support_Angle <= 90)
    count2 = 0;
    for start_ponit_x = BBox.MinX-adjacent:(line_spacing*cosd(PrintParameters.Support_Angle)):BBox.MaxX
        Start_list_x1 = [];
        Start_list_y1 = [];
        
        Start_list_x2 = [];
        Start_list_y2 = [];        
        %offset lines to avoid continuous gaps
        if ~mod(count2,2)
            Start_list_x1 = [Start_list_x1;(start_ponit_x:(segment_length+segment_break)*cosd(90-PrintParameters.Support_Angle):BBox.MaxX)'];        
            Start_list_y1 = tand(90-PrintParameters.Support_Angle)*(Start_list_x1-start_ponit_x); 

            Start_list_x2 = Start_list_x1 + PrintParameters.Spacing*cosd(PrintParameters.Support_Angle);
            Start_list_y2 = Start_list_y1 - PrintParameters.Spacing*sind(PrintParameters.Support_Angle);
            Start_list = [Start_list; Start_list_x1,Start_list_y1; Start_list_x2,Start_list_y2];  
            count2 = count2 +1;
        else
            Start_list_x1 = [Start_list_x1;(start_ponit_x:(segment_length+segment_break)*cosd(90-PrintParameters.Support_Angle):BBox.MaxX)'];        
            Start_list_y1 = cosd(PrintParameters.Support_Angle) *segment_break + tand(90-PrintParameters.Support_Angle)*(Start_list_x1-start_ponit_x); 

            Start_list_x2 = Start_list_x1 + PrintParameters.Spacing*cosd(PrintParameters.Support_Angle);
            Start_list_y2 = Start_list_y1 - PrintParameters.Spacing*sind(PrintParameters.Support_Angle);
            Start_list = [Start_list; Start_list_x1,Start_list_y1; Start_list_x2,Start_list_y2];  
            count2 = count2 +1;
        end
    end  
    End_list = [Start_list(:,1) + segment_length*cosd(90-PrintParameters.Support_Angle), Start_list(:,2) + segment_length*sind(90-PrintParameters.Support_Angle)];        
    
elseif PrintParameters.Support_Angle >= 90
%     Start_list (:,1) = (BBox.MinX:line_spacing:BBox.MaxX+adjacent)'; %x coordinate
%     Start_list (:,2) = BBox.MinY; %y coordinate
%     End_list (:,1) = Start_list (:,1) - adjacent;
%     End_list (:,2) = BBox.MaxY;
    for start_ponit_x = BBox.MinX:line_spacing:BBox.MaxX+adjacent
        Start_list = [Start_list;(start_ponit_x:segment_length+segment_break:BBox.MaxX+adjacent)',zeros(length((start_ponit_x:segment_length+segment_break:BBox.MaxX+adjacent)'),1)];
    end
        End_list = [(Start_list (:,1) + segment_length), zeros(length(Start_list),1)];        
        Start_list (:,2) = Start_list (:,1).*tand(PrintParameters.Support_Angle+90);
        End_list (:,2) = End_list (:,1).*tand(PrintParameters.Support_Angle+90);    
end

Infill_lines2 = []; %initialize
check = 1; %switch for infill pattern

for k = 1:size(Start_list,1)  
    %find intersection points of infill pattern.
    [in,~] = intersect(polyout, [Start_list(k,:); End_list(k,:)]);
    if ~isempty(in)
        % flip direction to minimise air time
        if check == 1
            Infill_lines2 = [Infill_lines2; in; nan nan];
        else
            Infill_lines2 = [Infill_lines2; flip(in,1); nan, nan];
        end
        check = check * -1;
    end
end

%check for starting error
if isnan(Infill_lines2(1,1))
    Infill_lines2 = Infill_lines2(2:end, :);
end
Infill_lines2;

%% combine lines

Supports_Points = [Infill_lines1; Infill_lines2];

end