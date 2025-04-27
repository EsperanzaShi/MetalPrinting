% This is the main slicing script to use the Aconity Metal 3D Printer
% Developed by Dr Connor Myant and Esperanza Shi, Imperial College London,
% Dyson School of Design Engineering, 2024.

% The script has not be computationally optimised nor can we garantee it is bug free. Please use, edit, add as you wish. Please share and communicate thoughts and development.

% STL files must be in mm!

%% INITIALISE
clear all
clc
close all

%% CREATE INPUTS, PRINT PARAMETER CONSTANTS
% this creates a large structure to store all the print parameters.
% you need to double check all parameters are correct for your machine etc.

% printer/build plate
PrintParameters.Bed_radius = 170; %size of printer build bed, mm
PrintParameters.Bed_X = PrintParameters.Bed_radius;
PrintParameters.Bed_Y = PrintParameters.Bed_radius;
PrintParameters.Bed_Z = 200; %mm
PrintParameters.BedLift = 3; %mm (3)

% General parameters
PrintParameters.SpotSize = 0.08; % laser spot diameter, mm (0.08)
PrintParameters.Layer_Thickness = 0.3; %mm (0.03)
PrintParameters.Scanningorder = ['sx';'vk';'vs';'kv']; % Aconity format - stands for supports, contours, infill, shape offsets (contours). The order you set here will be reflected 
PrintParameters.ParameterOrder = ['vk';'vs';'kv';'sx']; % IS THIS NEEDED?
PrintParameters.Spacing = 0.08; %mm (0.08)
PrintParameters.Rotate_Angle = 0; %rotation of scan strategy between layers

% single lines (sx)
PrintParameters.Line_LaserPower = 130; %Watts
PrintParameters.Line_ScanSpeed = 800; %Laser scan speed in mm/s
% vector support (st_sx)
PrintParameters.Vector_LaserPower = 135; %Watts
PrintParameters.Vector_ScanSpeed = 1000; %Laser scan speed in mm/s
% volume support (st_vs)
PrintParameters.Volume_LaserPower = 100; %Watts
PrintParameters.Volume_ScanSpeed = 1200; %Laser scan speed in mm/s
% contour (kv/vk)
PrintParameters.Contour_LaserPower = 130; %Watts
PrintParameters.Contour_ScanSpeed = 1000; %Laser scan speed in mm/s
% hatch (vs)
PrintParameters.Hatch_LaserPower = 260; %Watts (for CCM 140)
PrintParameters.Hatch_LaserPowerMin = 200; %Watts
PrintParameters.Hatch_LaserPowerMax = 300; %Watts
PrintParameters.Hatch_ScanSpeed = 400; %Laser scan speed in mm/s (for CCM 800)
PrintParameters.Hatch_ScanSpeedMin = 200; %Laser scan speed in mm/s (for CCM 800)
PrintParameters.Hatch_ScanSpeedMax = 400; %Laser scan speed in mm/s (for CCM 800)

PrintParameters.PointDistance = 200; %um
PrintParameters.ExposureTime = 400; %us

% toolpath parameters
PrintParameters.Number_of_Walls = 1;
PrintParameters.Infill_percentage = 1; %percentage out of 1
PrintParameters.Infill_Angle = 0;
PrintParameters.Support_Angle = 25; %degrees
PrintParameters.Support_Spacing = 1.2; %mm

% support parameters
PrintParameters.Minimum_Allowable_PrintedAngle = -60; % minimun printed angle of machine [(0 <Minimum_Allowable_PrintedAngle<90)recommneded degree is 40(the higher degree is, the more support material is required)]
PrintParameters.Support_LineSpace = 1;   % on/off switch of using diffrent air gap of infill in support

Part_ID = '001'; % required for CLT file output

%% on/off switches
% for compiler
CompileCLI_ILT = 1;
% for supports
Support_Generation = 1;
% make plots
Visualization = 1;

%% User selects file
% opens MATLAB GUI for user to select STL file for slicing. SLT MUST BE IN mm
[FileName, PathName] = uigetfile('*.stl','Select a File');
file_format = FileName(find(FileName=='.')+1:end);
TR = stlread([PathName, FileName]);

%% Translations and rotations of stl
% this insures the part is placed on the build plate and not hovering
% above, plus moves it to bed centre.

[BBox] = BoundingBoxXYZ (TR.Points); % Define bounding box around part
% insure part is on build platform
VerticesZ = TR.Points(:,3) - BBox.MinZ;

if Support_Generation == 1
    % lift part off build plate so it doesnt wield to build plate!
    VerticesZ = VerticesZ + PrintParameters.BedLift;
end

% move part to bed centre (0,0,0)
VerticesX = TR.Points(:,1) + abs(BBox.MinX);
VerticesY = TR.Points(:,2) + abs(BBox.MinY);

%Update
TR = triangulation(TR.ConnectivityList, [VerticesX, VerticesY, VerticesZ]);
[BBox] = BoundingBoxXYZ (TR.Points); % Define bounding box around part

%% Slice
% define Z heights for layer stack
Z_Height = 0:PrintParameters.Layer_Thickness:BBox.MaxZ;
Z_Height = round(Z_Height,9); %remove floating point errors
% Z_Height = linspace( 0, BBox.MaxZ,  BBox.MaxZ/PrintParameters.Layer_Thickness);

% set up progress bar
f = waitbar(0,'Slicing - forming raw layers');

% this function slices the STL into a layer stack of unsorted 2D
% coordinates. The coordinates are where the STL facets have intersected
% the layer/slice plane.
[Layer_Perimeters] = Find_Raw_Slices_Vectorised(TR,Z_Height);

%% Sort permiter coordinates, add walls, add infill. layer by layer
% now we take the unsorted coordinates and sort them into a CW/CCW order so
% that we can form polygons of each sliced layer. Once we have these we can
% define our infill patterns etc.



%loop through each layer and create the infill and contour (wall) scan
%strategies, ready for exporting into CLT files.
for Layer_Num=1:length(Layer_Perimeters)
   
    waitbar((1/size(Layer_Perimeters,1))*Layer_Num,f,'Sorting layers and adding infill') %update progress bar
    
    % start layer line counts
    Polylines_INTERNAL = [];
    Polylines_EXTERNAL = [];
    Infill = [];
    ForPlotting_INTERNAL = [];
    ForPlotting_EXTERNAL = [];
    ForPlotting_Infill =[];
    polyline_INT_count = 1;
    polyline_EXT_count = 1;
    Infill_count = 1;

    % connect the lines. Sort layer arrays into correct order
    if (~isempty(Layer_Perimeters{Layer_Num})) && (size(unique(single(Layer_Perimeters{Layer_Num}),'rows'),1)>1)

        %start by sorting the raw intersection points on each layer into
        %one closed polygon. the algorithm uses minimum 'flight time';
        %i.e. when a non-write ove is required it finds the nearest point.
        %Can easily be switched off inside the function (line 8)
        [Layer_Perimeters_Sorted{Layer_Num}] = Sort_LayerLines (Layer_Perimeters{Layer_Num}); %this will still contain points that are not required. i.e. additional points on a straight line. we only need turning points. This function is slow and should be vectorised but i couldnt work out how.

        % convert layers into polygons for simplification in MATLAB. this does the final sorting for perimeters (boundaries)
        Pgon_Layer(1) = polyshape(Layer_Perimeters_Sorted{Layer_Num},'SolidBoundaryOrientation','auto','Simplify',true,'KeepCollinearPoints',false);

        if ~isempty(Pgon_Layer(1).Vertices) % if layer is not empty create infills etc.

            %shrink polygon to account for walls
            if PrintParameters.Number_of_Walls > 1
                for k = 1:PrintParameters.Number_of_Walls-1
                    Pgon_Layer(k+1) = polybuffer(Pgon_Layer(1),-(k*PrintParameters.Spacing),'JointType','miter');
                end
            end

            % Add walls
            for q = 1:size(Pgon_Layer,2) % counts through the number of walls. each one is a reduced version of the real cross sectional perimenter
                if Pgon_Layer(q).NumRegions > 0 % check there are regions (vertices are not a singularity)
                    for k = 1:Pgon_Layer(q).NumRegions % check for multiple regions
                        Pgon_Layer_Regions = regions(Pgon_Layer(q));

                        %first identify internal walls (holes)
                        if Pgon_Layer_Regions(k).NumHoles > 0
                            Pgon_Region_Holes = holes(Pgon_Layer_Regions(k));
                            for i = 1:size(Pgon_Region_Holes,1) % each hole is given its own polyshape structure
                                [Xcord,Ycord] = boundary(Pgon_Region_Holes(i)); %extract boundary trace coordinates from polyshape object
                                if ~ispolycw(Xcord,Ycord) %check the polyline is clockwise for internal wall
                                    [Xcord,Ycord] = poly2cw(Xcord,Ycord);
                                end
                                Polylines_INTERNAL{polyline_INT_count} = reshape([Xcord,Ycord].',[],1);
                                ForPlotting_INTERNAL{polyline_INT_count} = [Xcord,Ycord]; % save for plotting is easier format
                                polyline_INT_count = polyline_INT_count + 1;
                            end
                        end

                        %Now add external walls (holes)
                        perimeter_ID = find(~ishole(Pgon_Layer_Regions(k))); % identify each external perimeter
                        for m = 1:size(perimeter_ID,1)
                            [Xcord,Ycord] = boundary(Pgon_Layer_Regions(k),perimeter_ID(m));
                            if ispolycw(Xcord,Ycord) %check the polyline is counter-clockwise for external wall
                                [Xcord,Ycord] = poly2ccw(Xcord,Ycord);
                            end
                            Polylines_EXTERNAL{polyline_EXT_count} = reshape([Xcord,Ycord].',[],1);
                            ForPlotting_EXTERNAL{polyline_EXT_count} = [Xcord,Ycord]; % save for plotting is easier format
                            polyline_EXT_count = polyline_EXT_count + 1;
                        end
                    end
                end
            end

            % Now add infill to the final polyshape (i.e. the most internal
            % wall). This function is a uni-directional scan strategy. Easy
            % enough to create other scan strategies and replace this
            % function as long as the output is a list of coordinates for
            % infill paths.
            divisions = 100; %defines how many the x axis tool path is to be divided into. 1 = no divisions. this makes changing the power or scan speeds along a path possible.
            [Layer_Infill] = Add_Lines_Infill2Layer_turbo(Pgon_Layer, PrintParameters, BBox, divisions);

            % add line types for exporting to CLT files
            Type = 4; %1 = polyline CW, 2= polyline CCW, 3= polyline open, 4= hatches

            % Now calculate laser power and speed for each scan
            Layer_Infill = Layer_Infill(:,1:2);
            % create min and max distance values of part from 0,0,0 coordinate
            rmin = sqrt(BBox.MinX^2 + BBox.MinY^2 + BBox.MinZ^2);
            rmax = sqrt(BBox.MaxX^2 + BBox.MaxY^2 + BBox.MaxZ^2);
            i = 1; %create count for loop
            r_average = []; % we use the average disctance of a scan, i.e. it's middle coordinate rather than start or end points to determine laser parameter value. easy to change.
            while i <= size(Layer_Infill,1)
                %catch lines with NAN on them
                if isnan(Layer_Infill(i))
                    Layer_Infill(i,3:4) = nan;
                    i = i+1;
                    if i> size(Layer_Infill,1)
                        break
                    end
                end

                %now create the corresponding power and speed values
                % change conditional statements based on a distance,
                % length, function etc as you want. i left a few examples
                % in.

                % a radial method from position 0,0,0. Power increases with
                % radial distance from position 0,0,0, and speed does the
                % opposite.
                r_average = sqrt( (((Layer_Infill(i+1,1) - Layer_Infill(i,1))/2) + min(Layer_Infill(i+1,1),Layer_Infill(i,1)))^2 + (((Layer_Infill(i+1,2) - Layer_Infill(i,2))/2) + min(Layer_Infill(i+1,2),Layer_Infill(i,2)))^2 + Z_Height(Layer_Num)^2);
                Layer_Infill(i,3) = PrintParameters.Hatch_LaserPowerMin + (((r_average-rmin)/(rmax-rmin)) * (PrintParameters.Hatch_LaserPowerMax-PrintParameters.Hatch_LaserPowerMin)); %d converts distance to power
                Layer_Infill(i,4) = PrintParameters.Hatch_ScanSpeedMax - (((r_average-rmin)/(rmax-rmin)) * (PrintParameters.Hatch_ScanSpeedMax-PrintParameters.Hatch_ScanSpeedMin)); %dconverts distance to speed

                % conditional statement on distance along x
                % distance_along_X = ((Layer_Infill(i+1,1) - Layer_Infill(i,1))/2) + min(Layer_Infill(i+1,1),Layer_Infill(i,1));
                % if distance_along_X >= 5
                %     Layer_Infill(i,3) = PrintParameters.Hatch_LaserPowerMax;
                % else
                %     Layer_Infill(i,3) = PrintParameters.Hatch_LaserPowerMin;
                % end

                %  % conditional statement on distance along y
                % distance_along_Y = ((Layer_Infill(i+1,2) - Layer_Infill(i,2))/2) + min(Layer_Infill(i+1,2),Layer_Infill(i,2));
                % if distance_along_Y >= 5
                %     Layer_Infill(i,3) = PrintParameters.Hatch_LaserPowerMax;
                %     Layer_Infill(i,4) = PrintParameters.Hatch_ScanSpeedMax;
                % else
                %     Layer_Infill(i,3) = PrintParameters.Hatch_LaserPowerMin;
                %     Layer_Infill(i,4) = PrintParameters.Hatch_ScanSpeedMin;
                % end

                % finish off table
                Layer_Infill(i+1,3) = Layer_Infill(i,3);
                i = i+2;
            end

            % sort infill hatches
            ForPlotting_Infill{Infill_count} = Layer_Infill; % save for plotting is easier format
            Layer_Infill(any(isnan(Layer_Infill), 2), :) = []; % remove nan lines for Aconity format
            Infill{Infill_count}.Points = reshape(Layer_Infill(:,1:2).',[],1);
            Infill{Infill_count}.Power = Layer_Infill(1:2:end,3);
            Infill{Infill_count}.Speed = Layer_Infill(1:2:end,4);
            Infill{Infill_count}.Type = Type;
            Infill_count = Infill_count + 1;
        end
    end

    % create structure for each layer
    Layer_toolpaths(Layer_Num).Walls_INTERNAL = Polylines_INTERNAL;
    Layer_toolpaths(Layer_Num).Walls_EXTERNAL = Polylines_EXTERNAL;
    Layer_toolpaths(Layer_Num).Infill = Infill;

    %for plotting
    Layer_plots(Layer_Num).Walls_INTERNAL = ForPlotting_INTERNAL;
    Layer_plots(Layer_Num).Walls_EXTERNAL = ForPlotting_EXTERNAL;
    Layer_plots(Layer_Num).Infill = ForPlotting_Infill;
end


%% Add support structures
if Support_Generation == 1

    waitbar(0,f,'creating support structures') %update progress bar

    % find all faces with a negative elivation > than min. allowable print angle
    FNorms = faceNormal(TR);
    [azimuth,elevation,r] = cart2sph(FNorms(:,1),FNorms(:,2),FNorms(:,3));
    ListF = find(elevation <= deg2rad(PrintParameters.Minimum_Allowable_PrintedAngle)); %list of face IDs that require supports

    %build sub structures of face groups needing supports
    New_Faces = TR.ConnectivityList(ListF,:);
    [Unique_Verts_ind_list,~,NewFaces_ind] = unique(New_Faces(:)); % Find all the vertex indices left in the new faces list.
    New_Faces = reshape(NewFaces_ind,[],3); % NewFaces_ind is the index corresponding each row in Unique_Verts_ind_list to New_Faces. so reshape it into a mxn array to get back to a connectivity list
    New_Verts = TR.Points(Unique_Verts_ind_list,:); % form new list of vertices
    New_TR = triangulation(New_Faces, New_Verts);

    %check for islands
    E = [New_Faces(:,1:2);New_Faces(:,2:3);New_Faces(:,[3,1])];
    E = unique([E;E(:,2:-1:1)],'rows');
    W = sparse(E(:,1),E(:,2),ones(1,size(E,1))); % adjacency matrix
    G = graph(W);
    [bin,binsize] = conncomp(G);
    %create supports for each island
    for i = 1:numel(binsize)
        sub_verts2keep_id = find(bin==i);
        [~,Indx_2Keep] = ismember(New_TR.ConnectivityList,sub_verts2keep_id);
        sub_New_Faces = New_TR.ConnectivityList(~any(Indx_2Keep == 0,2),:); %remove any row that contains a 0
        [Unique_Verts_ind_list,~,NewFaces_ind] = unique(sub_New_Faces(:)); % Find all the vertex indices left in the new faces list.
        sub_New_Faces = reshape(NewFaces_ind,[],3); % NewFaces_ind is the index corresponding each row in Unique_Verts_ind_list to New_Faces. so reshape it into a mxn array to get back to a connectivity list
        sub_New_Verts = New_TR.Points(Unique_Verts_ind_list,:); % form new list of vertices
        sub_New_TR = triangulation(sub_New_Faces, sub_New_Verts);
        flat_sub_New_TR = triangulation(sub_New_Faces, sub_New_Verts(:,1:2));

        Pgon_supports(i) = boundaryshape(flat_sub_New_TR);
        MaxZ(i) = max(New_Verts(bin==i,3));
    end

    %build  polyshapes of all layers
    clear Pgon_Layer
    for Layer_Num=1:length(Layer_Perimeters_Sorted)
        if ~isempty(Layer_Perimeters_Sorted{Layer_Num})
            Pgon_Layer (Layer_Num)= polyshape(Layer_Perimeters_Sorted{Layer_Num},'SolidBoundaryOrientation','auto','Simplify',true,'KeepCollinearPoints',false);
        end
    end

    for Layer_Num=1:length(Layer_Perimeters_Sorted)% start layer line counts

        log_MaxZ = MaxZ > Z_Height(Layer_Num);
        if any(log_MaxZ)  %any() determines if any array elements are nonzero
            Support_poly(Layer_Num) = union(Pgon_supports(log_MaxZ));
        else
            Support_poly(Layer_Num) = polyshape(0,0);
        end

        if ~isempty(Support_poly(Layer_Num).Vertices)
            % find the compound of all layers above this one
            Compound_Layer(Layer_Num) = union(Pgon_Layer(Layer_Num:end));
            %find the intersection with the compound poly
            polyout(Layer_Num) =  intersect(Support_poly(Layer_Num), Compound_Layer(Layer_Num));
            % subtract build layer
            polyout(Layer_Num) = subtract(polyout(Layer_Num),Pgon_Layer(Layer_Num));
            % shrink slightly to avoid full contact on sides
            polyout(Layer_Num) = polybuffer(polyout(Layer_Num),-PrintParameters.SpotSize*PrintParameters.Spacing,'JointType','miter');

            if polyout(Layer_Num).NumRegions > 0 && any(area(regions(polyout(Layer_Num))) >= 4) %need to add min area parameter to controls
                % add support lines
                [Supports_Points] = Add_Support_Infill(polyout(Layer_Num), PrintParameters, BBox);
                % add line types
                Type = 4; %1 = polyline CW, 2= polyline CCW, 3= polyline open, 4= hatches

                % sort infill hatches
                Supports_Points = Supports_Points(:,1:2);
                ForPlotting_Supports = Supports_Points; % save for plotting is easier format

                Supports_Points(any(isnan(Supports_Points), 2), :) = []; % remove nan lines for Aconity format
                Supports{1}.Points = reshape(Supports_Points.',[],1);
                Supports{1}.Type = Type;

                Layer_toolpaths(Layer_Num).Supports = Supports;
                Layer_plots(Layer_Num).Supports{1} = ForPlotting_Supports;
            end
        end
    end
else
    Supports=[];
    Layer_toolpaths(Layer_Num).Supports = [];
    Layer_plots(Layer_Num).Supports = Supports;
end

%% COMBINE ALL INTO A SINGLE TEXT FILE - create G code
if CompileCLI_ILT == 1

     waitbar(0.5,f,'compiling G Code') %update progress bar

    CLI_PLUS_ILT_Compiler_vk_2(PrintParameters, Layer_toolpaths, Part_ID) % use for Aconity versions <3. Results in larger output files but should do the same thing in the end.
    % CLI_PLUS_ILT_Compiler_vk_3(PrintParameters, Layer_toolpaths, Part_ID) % use if your Aconity version is 3.0 or better, this way you can use the POWERS function
end

%% VISUALISE PART AND TOOLPATHS
if Visualization == 1

    waitbar(0.75,f,'creating visualisation') %update progress bar

    figure
    % trisurf(TR,'FaceColor','y','EdgeColor','none', 'FaceAlpha',0.1);
    axis equal
    % axis tight
    axis vis3d
    grid on
    hold on
    for Layer_Num = 1:size(Layer_plots,2)
        % plot internal walls red
        for i = 1:size(Layer_plots(Layer_Num).Walls_INTERNAL,2)
            plot3(Layer_plots(Layer_Num).Walls_INTERNAL{i}(:,1),Layer_plots(Layer_Num).Walls_INTERNAL{i}(:,2),repmat(Z_Height(Layer_Num),size(Layer_plots(Layer_Num).Walls_INTERNAL{i},1),1),'color','r')
        end
        %         plot external walls blue
        for i = 1:size(Layer_plots(Layer_Num).Walls_EXTERNAL,2)
            plot3(Layer_plots(Layer_Num).Walls_EXTERNAL{i}(:,1),Layer_plots(Layer_Num).Walls_EXTERNAL{i}(:,2),repmat(Z_Height(Layer_Num),size(Layer_plots(Layer_Num).Walls_EXTERNAL{i},1),1),'color','b')
            scatter3(Layer_plots(Layer_Num).Walls_EXTERNAL{i}(1,1),Layer_plots(Layer_Num).Walls_EXTERNAL{i}(1,2),Z_Height(Layer_Num),'or')
        end
        %        plot internal hatches green
        for i = 1:size(Layer_plots(Layer_Num).Infill,2)
            plot3(Layer_plots(Layer_Num).Infill{i}(:,1),Layer_plots(Layer_Num).Infill{i}(:,2),repmat(Z_Height(Layer_Num),size(Layer_plots(Layer_Num).Infill{i},1),1),'color','g')
            scatter3(Layer_plots(Layer_Num).Infill{i}(1:3:end,1),Layer_plots(Layer_Num).Infill{i}(1:3:end,2),repmat(Z_Height(Layer_Num),length(Layer_plots(Layer_Num).Infill{i}(1:3:end,1)),1),20,Layer_plots(Layer_Num).Infill{i}(1:3:end,3))
            colormap(jet)
        end
        % plot support structure
        if isfield(Layer_plots(Layer_Num),'Supports')
            for i = 1:size(Layer_plots(Layer_Num).Supports,2)
                plot3(Layer_plots(Layer_Num).Supports{i}(:,1),Layer_plots(Layer_Num).Supports{i}(:,2),repmat(Z_Height(Layer_Num),size(Layer_plots(Layer_Num).Supports{i},1),1),'-*','color','c')
            end
        end
    end
end

close(f)

