function [] = CLI_PLUS_ILT_Compiler_vk_3(PrintParameters, Layer_toolpaths, Part_ID)
% Updated version of CLI compilier to new CLI+ code
% compile CLI+ and ILT files for Aconity metal printer
% Connor Myant
% March 2024

% creates a command file for metal 3D printing on the Aconity Machine
% unlike other printers which use GCode format this use Common Layer Interface (CLI), these are then zipped into an ILT file
% this system is closer to DXF files and laser cutting


%% create CLI files
for Part_Num = 1:length(PrintParameters.ScanningOrder)
    % create a write enabled txt file
    if strcmp(PrintParameters.ScanningOrder(Part_Num,:),'sx')
        shellNum = 'st';
    else
        shellNum = 's1';
    end
    file_names{Part_Num} = ['modelsection_' Part_ID '_' shellNum '_' PrintParameters.ScanningOrder(Part_Num,:) '.cli'];
    fileID = fopen(file_names{Part_Num},'wt');
    
    % create header information
    fprintf(fileID,'%s\n', '$$HEADERSTART');
    fprintf(fileID,'%s\n', '$$ASCII' );
    fprintf(fileID,'%s\n', '$$UNITS/1.00000');
    fprintf(fileID,'%s\n', ['$$DATE/' datestr(now,'yyyymmdd')]);
    fprintf(fileID,'%s\n', ['$$LAYERS/' num2str(size(Layer_toolpaths,2))]);
    fprintf(fileID,'%s\n', '$$HEADEREND');
    fprintf(fileID,'%s\n', '$$GEOMETRYSTART');
    
    % WRITE WALLS AND INFILL TO LAYER TEXT FILE
    
    for Layer_Num = 1:size(Layer_toolpaths,2) %wall contours
        % layer header, not as layer number but layer height. height needs to be in the units as defined above
        fprintf(fileID,'%s\n', ['$$LAYER/' num2str(Layer_Num*PrintParameters.Layer_Thickness)]);
        
        if strcmp('sx',PrintParameters.ScanningOrder(Part_Num,:)) %supports
            fprintf(fileID,'%s\n', ['$$POWER/' num2str(PrintParameters.Vector_LaserPower)]);
            fprintf(fileID,'%s\n', ['$$SPEED/' num2str(PrintParameters.Vector_ScanSpeed)]);
            % Now infill hatches
            for i = 1:size(Layer_toolpaths(Layer_Num).Supports,2)
                if Layer_toolpaths(Layer_Num).Supports{i}.Type == 1 %polyline CW
                    fprintf(fileID,'%s\n', strcat(['$$POLYLINE/' num2str(str2num(Part_ID)) ',0,' num2str(size(Layer_toolpaths(Layer_Num).Supports{i}.Points,1)/2) ','], strjoin(string(Layer_toolpaths(Layer_Num).Supports{i}.Points.'),','))); % Polyline, ID, Direction, Number of points, Points (px1,py1,...pxn, pyn)
                elseif Layer_toolpaths(Layer_Num).Supports{i}.Type == 2 %polyline CCW
                    fprintf(fileID,'%s\n', strcat(['$$POLYLINE/' num2str(str2num(Part_ID)) ',1,' num2str(size(Layer_toolpaths(Layer_Num).Supports{i}.Points,1)/2) ','], strjoin(string(Layer_toolpaths(Layer_Num).Supports{i}.Points.'),','))); % Polyline, ID, Direction, Number of points, Points (px1,py1,...pxn, pyn)
                elseif Layer_toolpaths(Layer_Num).Supports{i}.Type == 3 %open polyline
                    fprintf(fileID,'%s\n', strcat(['$$POLYLINE/' num2str(str2num(Part_ID)) ',2,' num2str(size(Layer_toolpaths(Layer_Num).Supports{i}.Points,1)/2) ','], strjoin(string(Layer_toolpaths(Layer_Num).Supports{i}.Points.'),','))); % Polyline, ID, Direction, Number of points, Points (px1,py1,...pxn, pyn)
                elseif Layer_toolpaths(Layer_Num).Supports{i}.Type == 4 %hatches
                    fprintf(fileID,'%s\n', strcat(['$$HATCHES/' num2str(str2num(Part_ID)) ',' num2str(size(Layer_toolpaths(Layer_Num).Supports{i}.Points,1)/4) ','], strjoin(string(Layer_toolpaths(Layer_Num).Supports{i}.Points.'),','))); % Hatches, ID, number of hatches, points (H1startx, H1starty, H1endx, H1endy,...Hnstartx, Hnstarty, Hnendx, Hnendy)
                end
            end

        elseif strcmp('vk',PrintParameters.ScanningOrder(Part_Num,:))
            % Now any walls which are defined as polylines
            fprintf(fileID,'%s\n', ['$$POWER/' num2str(PrintParameters.Contour_LaserPower)]);
            fprintf(fileID,'%s\n', ['$$SPEED/' num2str(PrintParameters.Contour_ScanSpeed)]);
            % start with INTERNALS, CW direction = 0
            for i = 1:size(Layer_toolpaths(Layer_Num).Walls_INTERNAL,2)
                fprintf(fileID,'%s\n', strcat(['$$POLYLINE/' num2str(str2num(Part_ID)) ',0,' num2str(size(Layer_toolpaths(Layer_Num).Walls_INTERNAL{i},1)/2) ','], strjoin(string(Layer_toolpaths(Layer_Num).Walls_INTERNAL{i}.'),','))); % Polyline, ID, Direction, Number of points, Points (px1,py1,...pxn, pyn)
            end
            % next EXTERNALS, CCW direction = 1
            for i = 1:size(Layer_toolpaths(Layer_Num).Walls_EXTERNAL,2)
                fprintf(fileID,'%s\n', strcat(['$$POLYLINE/' num2str(str2num(Part_ID)) ',1,' num2str(size(Layer_toolpaths(Layer_Num).Walls_EXTERNAL{i},1)/2) ','], strjoin(string(Layer_toolpaths(Layer_Num).Walls_EXTERNAL{i}.'),','))); % Polyline, ID, Direction, Number of points, Points (px1,py1,...pxn, pyn)
            end            
        
        elseif strcmp('vs',PrintParameters.ScanningOrder(Part_Num,:)) %infills
            % Now infill hatches
            for i = 1:size(Layer_toolpaths(Layer_Num).Infill,2)
                if Layer_toolpaths(Layer_Num).Infill{i}.Type == 1 %polyline CW
                    fprintf(fileID,'%s\n', strcat(['$$POLYLINE/' num2str(str2num(Part_ID)) ',0,' num2str(size(Layer_toolpaths(Layer_Num).Infill{i}.Points,1)/2) ','], strjoin(string(Layer_toolpaths(Layer_Num).Infill{i}.Points.'),','))); % Polyline, ID, Direction, Number of points, Points (px1,py1,...pxn, pyn)
                elseif Layer_toolpaths(Layer_Num).Infill{i}.Type == 2 %polyline CCW
                    fprintf(fileID,'%s\n', strcat(['$$POLYLINE/' num2str(str2num(Part_ID)) ',1,' num2str(size(Layer_toolpaths(Layer_Num).Infill{i}.Points,1)/2) ','], strjoin(string(Layer_toolpaths(Layer_Num).Infill{i}.Points.'),','))); % Polyline, ID, Direction, Number of points, Points (px1,py1,...pxn, pyn)
                elseif Layer_toolpaths(Layer_Num).Infill{i}.Type == 3 %open polyline
                    fprintf(fileID,'%s\n', strcat(['$$POLYLINE/' num2str(str2num(Part_ID)) ',2,' num2str(size(Layer_toolpaths(Layer_Num).Infill{i}.Points,1)/2) ','], strjoin(string(Layer_toolpaths(Layer_Num).Infill{i}.Points.'),','))); % Polyline, ID, Direction, Number of points, Points (px1,py1,...pxn, pyn)
                elseif Layer_toolpaths(Layer_Num).Infill{i}.Type == 4 %hatches

                    fprintf(fileID,'%s\n', ['$$POWERS/' strjoin(string(reshape([(1:size(Layer_toolpaths(Layer_Num).Infill{i}.Points,1)/4)', round(Layer_toolpaths(Layer_Num).Infill{i}.Power)]',[],1) )',',') ] );
                    fprintf(fileID,'%s\n', ['$$SPEED/' num2str(PrintParameters.Hatch_ScanSpeed)]);
                    fprintf(fileID,'%s\n', strcat(['$$HATCHES/' num2str(str2num(Part_ID)) ',' num2str(size(Layer_toolpaths(Layer_Num).Infill{i}.Points,1)/4) ','], strjoin(string(Layer_toolpaths(Layer_Num).Infill{i}.Points.'),','))); % Hatches, ID, number of hatches, points (H1startx, H1starty, H1endx, H1endy,...Hnstartx, Hnstarty, Hnendx, Hnendy)
                end
            end
        elseif strcmp('kv',PrintParameters.ScanningOrder(Part_Num,:)) %shape offsets
        end
    end
    
    fprintf(fileID,'%s', '$$GEOMETRYEND');
    fclose(fileID);
end
%% create parameter-file
% create a write enabled txt file
file_names{Part_Num+1} = ['modelsection_' Part_ID '_param.txt'];
Parameter_file = fopen(file_names{Part_Num+1},'wt');

reOrder = [2; 3; 4; 1];
for i = 1:Part_Num   
    fprintf(Parameter_file,'%s\n', ['[' file_names{reOrder(i)} ']']);
    if strcmp('vk',PrintParameters.ParameterOrder(i,:)) %contour
        fprintf(Parameter_file,'%s\n', ['LaserSpeed = ' num2str(PrintParameters.Contour_ScanSpeed,'%.3f') ' mm/s']);
        fprintf(Parameter_file,'%s\n', ['LaserPower = ' num2str(PrintParameters.Contour_LaserPower,'%.3f') ' watt']);
        fprintf(Parameter_file,'%s\n', 'FocusShift = 0.00000 mm');
        fprintf(Parameter_file,'%s\n', ['PointDistance = ' num2str(PrintParameters.PointDistance,'%.3f') ' ' char(181) 'm']);
        fprintf(Parameter_file,'%s\r\n', ['ExposureTime = ' num2str(PrintParameters.ExposureTime,'%.3f') ' ' char(181) 's']);
    elseif strcmp('vs',PrintParameters.ParameterOrder(i,:)) %infills
        fprintf(Parameter_file,'%s\n', ['LaserSpeed = ' num2str(PrintParameters.Hatch_ScanSpeed,'%.3f') ' mm/s']);
        fprintf(Parameter_file,'%s\n', ['LaserPower = ' num2str(PrintParameters.Hatch_LaserPower,'%.3f') ' watt']);
        fprintf(Parameter_file,'%s\n', 'FocusShift = 0.00000 mm');
        fprintf(Parameter_file,'%s\n', ['PointDistance = ' num2str(PrintParameters.PointDistance,'%.3f') ' ' char(181) 'm']);
        fprintf(Parameter_file,'%s\r\n', ['ExposureTime = ' num2str(PrintParameters.ExposureTime,'%.3f') ' ' char(181) 's']);
    elseif strcmp('kv',PrintParameters.ParameterOrder(i,:)) %shape offsets
        fprintf(Parameter_file,'%s\n', ['LaserSpeed = ' num2str(PrintParameters.Contour_ScanSpeed,'%.3f') ' mm/s']);
        fprintf(Parameter_file,'%s\n', ['LaserPower = ' num2str(PrintParameters.Contour_LaserPower,'%.3f') ' watt']);
        fprintf(Parameter_file,'%s\n', 'FocusShift = 0.00000 mm');
        fprintf(Parameter_file,'%s\n', ['PointDistance = ' num2str(PrintParameters.PointDistance,'%.3f') ' ' char(181) 'm']);
        fprintf(Parameter_file,'%s\r\n', ['ExposureTime = ' num2str(PrintParameters.ExposureTime,'%.3f') ' ' char(181) 's']);
    elseif strcmp('sx',PrintParameters.ParameterOrder(i,:)) %supports
        fprintf(Parameter_file,'%s\n', ['LaserSpeed = ' num2str(PrintParameters.Line_ScanSpeed,'%.3f') ' mm/s']);
        fprintf(Parameter_file,'%s\n', ['LaserPower = ' num2str(PrintParameters.Line_LaserPower,'%.3f') ' watt']);
        fprintf(Parameter_file,'%s\n', 'FocusShift = 0.00000 mm');
        fprintf(Parameter_file,'%s\n', ['PointDistance = ' num2str(PrintParameters.PointDistance,'%.3f') ' ' char(181) 'm']);
        fprintf(Parameter_file,'%s\r\n', ['ExposureTime = ' num2str(PrintParameters.ExposureTime,'%.3f') ' ' char(181) 's']);
    end
end
fclose(Parameter_file);

% create folder to zip
%% create a txt file to output and write to
prompt = {'Enter the part name (NO NUMBERS)'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'Sample'};
PartName = inputdlg(prompt,dlg_title,num_lines,defaultans);


destination = ['ILT' PartName{1}];
movefile('modelsection*', destination);

end