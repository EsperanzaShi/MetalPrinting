%this function takes the list of intersection points found in
%'Find_Raw_Slices' and re-organises them into a ordered list of points to
%trace out the layer perimeter

function [Layer_Lines_ordered] = Sort_LayerLines (Layer_Lines)

%find starting point - use nearest to 0,0,0
Idx_start = knnsearch([Layer_Lines(:,1:2);Layer_Lines(:,3:4)],[0, 0],'Distance','euclidean');
if Idx_start > size(Layer_Lines,1)
    Idx_start = Idx_start - size(Layer_Lines,1);
    Flag_start = 1;
    Last_Point = Layer_Lines(Idx_start,1:2);
else
    Flag_start = -1;
    Last_Point = Layer_Lines(Idx_start,3:4);
end

Layer_Lines_ordered = nan(size(Layer_Lines,1)+1,2); %intialise matrix

if Flag_start == -1   
    Layer_Lines_ordered(1:2,:) = [Layer_Lines(Idx_start,1:2); Layer_Lines(Idx_start,3:4)];   
else
    Layer_Lines_ordered(1:2,:) = [Layer_Lines(Idx_start,3:4); Layer_Lines(Idx_start,1:2)]; 
end

Layer_Lines(Idx_start,:) = NaN; %remove previous line from future searches

%now go through entire layer
count = 3;
for k = 2:size(Layer_Lines,1)
    index = 0; %index location of last move
    check = 0; %check counter for finding next move
    
    while index == 0
        if Flag_start == -1
            
            [~, index]=ismembertol(Layer_Lines_ordered(count-1,:),Layer_Lines(:,1:2),'ByRows',true);
            if index==0
                Flag_start = Flag_start*-1;
                check = check+1;
            else
                Layer_Lines_ordered(count,:) = Layer_Lines(index,3:4);
                count = count+1;
                Layer_Lines(index,:) = NaN;
            end
            
        else
            
            [~, index]=ismembertol(Layer_Lines_ordered(count-1,:),Layer_Lines(:,3:4),'ByRows',true);
            if index==0
                Flag_start = Flag_start*-1;
                check = check+1;
            else
                Layer_Lines_ordered (count,:) = Layer_Lines(index,1:2);
                count = count+1;
                Layer_Lines(index,:) = NaN;
            end
            
        end
        
        if check == 2  % no adjoining point can be found presumably because its a new section (internal window)           
            %find NEW starting point - use nearest to last point
            [Idx_start, ~] = knnsearch([Layer_Lines(:,1:2);Layer_Lines(:,3:4)],Layer_Lines_ordered(count-1,:),'Distance','euclidean');
            if Idx_start > size(Layer_Lines,1)
                Idx_start = Idx_start-size(Layer_Lines,1);
                Flag_start = 1;
            else
                 Flag_start = -1;
            end
            
            if Flag_start==-1
                Layer_Lines_ordered (count:count+2,:) = [Last_Point; [nan nan]; Layer_Lines(Idx_start,3:4)];
                count = count+3;
                Last_Point = Layer_Lines(Idx_start,3:4);
            else
                Layer_Lines_ordered (count:count+2,:) = [Last_Point; [nan nan]; Layer_Lines(Idx_start,1:2)];
                count = count+3;
                Last_Point = Layer_Lines(Idx_start,1:2);
            end
            
            Layer_Lines(Idx_start,:) = NaN;
            index = 1;
        end
        
    end
end
end