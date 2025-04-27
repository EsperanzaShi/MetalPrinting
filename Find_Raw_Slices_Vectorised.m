% main slicing function
function [Layer_Lines] = Find_Raw_Slices_Vectorised (TR,Layer_Hieghts)
%plane-line intersect for all triangles with all slice planes

Slice_normal = [0,0,1]; % keep constant for now, will update later for non Z axis slicing
Spacing = repelem(Layer_Hieghts,size(TR.ConnectivityList,1)); %expand layer height list to repeat for the same number of triangles in TR structure on each layer.
V0 = [repelem(0,size(Spacing,2))',repelem(0,size(Spacing,2))',Spacing'];

%P0 and P1
Point1 = repmat(TR.Points(TR.ConnectivityList(:,1),:),size(Layer_Hieghts,2),1);
Point2 = repmat(TR.Points(TR.ConnectivityList(:,2),:),size(Layer_Hieghts,2),1);
Point3 = repmat(TR.Points(TR.ConnectivityList(:,3),:),size(Layer_Hieghts,2),1);

%for edge1
u = Point2-Point1;
w = Point1 - V0;
D = dot(repmat(Slice_normal,size(u,1),1),u,2);
N = -dot(repmat(Slice_normal,size(u,1),1),w,2);

sI = N ./ D;
I1 = Point1+ sI.*u;

% remove unwanted points
I1((abs(D) < 10^-7) & (N~=0), :) = nan;       % The segment is parallel to plane, but not or lies in the plane
I1(sI < 0 | sI > 1,:) = nan; % remove points that lie outside the segment, so there is no intersection

%for edge2
u = Point3-Point2;
w = Point2 - V0;
D = dot(repmat(Slice_normal,size(u,1),1),u,2);
N = -dot(repmat(Slice_normal,size(u,1),1),w,2);

sI = N ./ D;
I2 = Point2+ sI.*u;

% remove unwanted points
I2((abs(D) < 10^-7) & (N~=0), :) = nan;         % The segment is parallel to plane, but not on the plane
I2(sI < 0 | sI > 1,:) = nan; % remove points that lie outside the segment, so there is no intersection

%for edge1
u = Point1-Point3;
w = Point3 - V0;
D = dot(repmat(Slice_normal,size(u,1),1),u,2);
N = -dot(repmat(Slice_normal,size(u,1),1),w,2);

sI = N ./ D;
I3 = Point3+ sI.*u;

% remove unwanted points
I3((abs(D) < 10^-7) & (N~=0), :) = nan;        % The segment is parallel to plane, but not on the plane
I3(sI < 0 | sI > 1,:) = nan; % remove points that lie outside the segment, so there is no intersection

% we now have three lists of intersection points for all face edges from
% all faces in the TR. NOW we need to build the perimeter shapes (ordered)
% for each layer.
I = [I1,I2,I3]; % this now represents all intersections/faces from TR
I(all(isnan(I),2),:) = []; %remove all rows with only nan in them - no interestion from that face.
I(all(I(:,1:3)==I(:,4:6),2) | all(I(:,1:3)==I(:,7:9),2) | all(I(:,4:6)==I(:,7:9),2),:) = []; % remove rows where two points are the same, no line segment found, intersection occurred at face vertex.

% want to remove nan elemets and reduce to a n,6 matrix.
[~,X] = sort(isnan(I),2);
for k = 1:size(I,1)
    I(k,:) = I(k,X(k,:));
end
I = round(I(:,1:6),9); %remove floating point errors

% now put into cells
Layer_Lines = cell(size(Layer_Hieghts,2),1);
for i = 1:size(Layer_Hieghts,2)
    Layer_Lines{i} = I(I(:,6)==Layer_Hieghts(i),[1:2,4:5]);
    % Layer_Lines {i} = unique(Layer_Lines{i},'stable','rows'); %remove duplicates.  not sure why they are there
end
end