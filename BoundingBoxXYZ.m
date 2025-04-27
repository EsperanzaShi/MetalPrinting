%Function to find bounding box area/coordinates around array in cartisean space
function [BBox] = BoundingBoxXYZ (Points)
BBox.MaxX = max(Points(:,1));
BBox.MaxY = max(Points(:,2));
BBox.MinX = min(Points(:,1));
BBox.MinY = min(Points(:,2));

if size(Points,2)>2
    BBox.MaxZ = max(Points(:,3));
    BBox.MinZ = min(Points(:,3));
    
    BBox.Centre = [0.5*(BBox.MaxX-abs(BBox.MinX)),0.5*(BBox.MaxY-abs(BBox.MinY)),0.5*(BBox.MaxZ-abs(BBox.MinZ))];
else
    BBox.Centre = [0.5*(BBox.MaxX-abs(BBox.MinX)),0.5*(BBox.MaxY-abs(BBox.MinY))];
end

end