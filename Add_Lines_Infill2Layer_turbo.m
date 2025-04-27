%this function adds lines infil
function [Layer_Infill] = Add_Lines_Infill2Layer_turbo(Pgon_Layer, PrintParameters, BBox, divisions)

% Define bounding box around internal part for infilling
Inner_temp_Wall = polybuffer(Pgon_Layer(PrintParameters.Number_of_Walls),-(PrintParameters.Spacing*0.499),'JointType','miter');

% define space between each line
if Inner_temp_Wall.NumRegions > 0
        line_spacing = PrintParameters.Spacing;


% calculate distance along for each line
adjacent = abs(BBox.MaxY/tand(PrintParameters.Infill_Angle));
if isinf(adjacent)
    adjacent = 0;
end

%create infill (for the whole bed) lists of start and end points
Start_list = [];
End_list = [];

if PrintParameters.Infill_Angle == 0
    %y coordinate
    Start_list (:,2) = repelem(flip((BBox.MinY-0.5*PrintParameters.Spacing:line_spacing:BBox.MaxY+0.5*PrintParameters.Spacing)'),divisions);
    End_list (:,2) = Start_list (:,2);
    %x coordinate
    Start_list (:,1) = repmat((BBox.MinX:(BBox.MaxX-BBox.MinX)/divisions:BBox.MaxX-((BBox.MaxX-BBox.MinX)/divisions))',numel(flip((BBox.MinY-0.5*PrintParameters.Spacing:line_spacing:BBox.MaxY+0.5*PrintParameters.Spacing)')),1);
    End_list (:,1) = repmat((BBox.MinX+((BBox.MaxX-BBox.MinX)/divisions):(BBox.MaxX-BBox.MinX)/divisions:BBox.MaxX)', numel(flip((BBox.MinY-0.5*PrintParameters.Spacing:line_spacing:BBox.MaxY+0.5*PrintParameters.Spacing)')),1);
  
 elseif (PrintParameters.Infill_Angle > 0) && (PrintParameters.Infill_Angle <= 90)
    %x coordinate
    Start_list (:,1) = flip((BBox.MinX-adjacent:line_spacing:BBox.MaxX)'); 
    End_list (:,1) = Start_list (:,1) + adjacent;
    %y coordinate
    Start_list (:,2) = BBox.MinY; 
    End_list (:,2) = BBox.MaxY;
  
elseif PrintParameters.Infill_Angle >= 90
    %x coordinate
    Start_list (:,1) = flip((BBox.MinX:line_spacing:BBox.MaxX+adjacent)');
    End_list (:,1) = Start_list (:,1) - adjacent;
    %y coordinate
    Start_list (:,2) = BBox.MinY; 
    End_list (:,2) = BBox.MaxY; 
end

Layer_Infill = []; %initialize


%     check = 1; %switch for infill pattern
for q = 1:size(Start_list,1)
    %find intersection points of infill pattern.
    [in,out] = intersect(Inner_temp_Wall, [Start_list(q,:); End_list(q,:)]);  % polyxpoly
    % out;

    if ~isempty(in)
        Layer_Infill = [Layer_Infill; in; nan nan];
    end
end
    
    %check for starting error
    if isnan(Layer_Infill(1,1))
        Layer_Infill = Layer_Infill(2:end, :);
    end

end

if Inner_temp_Wall.NumRegions == 0
    Layer_Infill = [nan nan];
end


end