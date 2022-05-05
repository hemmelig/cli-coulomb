function [x,y, z, model_elements, kode, poisson_ratio, youngs_mod, calculation_depth, coef_friction, regional_stress] = read_input_file(filename)
    % Parse the input information for the Okada deformation calculation from file.
    %
    % Args:
    %   filename : str
    %       Name of the input file to be parsed.
    %
    % Returns:
    %   x :
    %   y :
    %   z :
    %   el :
    %   kode :
    %   pois :
    %   young :
    %   cdepth : 
    %   fric :
    %   rstress :

    fid = fopen(filename, "r");
        % Parsing the header - first two lines are just meta information
        tmp = fgets(fid);
        tmp = fgets(fid);

        % Get number of faults - number after the '#fixed='
        line = textscan(fid, '%*s %*s %*s %*s %*s %4u16 %*s %*s', 1);
        num = int32([line{1}]);  % Max number of elements is 65,535

        % Get Poisson's ratio & calculation depth - numbers after 'PR1=' and 'DEPTH='
        line = textscan(fid, '%*s %15.3f32 %*s %*s %*s %15.3f32', 1);
        poisson_ratio = double([line{1}]);
        calculation_depth = double([line{2}]);

        % Get Young's modulus - number after 'E1='
        line = textscan(fid, '%*s %n %*s %*s', 1);
        youngs_mod = double([line{1}]);

        % Skip lines on symmetry
        line = textscan(fid, '%*s %*s %*s %*s', 1);

        % Get coefficient of friction - number after 'FRIC='
        line = textscan(fid, '%*s %15.3f32', 1);
        coef_friction = double([line{1}]);

        % Get regional stress field
        line = textscan(fid, '%*s %15.6f32 %*s %15.6f32 %*s %15.6f32 %*s %15.6f32', 3);
        regional_stress = [line{:}];

        % Skip next 2 header lines
        tmp = fgets(fid);
        tmp = fgets(fid);
        tmp = fgets(fid);
        tmp = fgets(fid);

        % Now on to reading in the deformation elements to be modelled!
        % Reads in the number of lines specified in the header + the buffer
        model = textscan(fid,...
            '%3u16 %10.9f32 %10.9f32 %10.9f32 %10.9f32 %3u16 %10.9f32 %10.9f32 %10.9f32 %10.9f32 %10.9f32',...
            num);

        id = uint16([model{1}]);
        if length(id) ~= num
            disp('************************************************************************');
            disp('**  Please change #fixed in the 3rd row in the input file afterward.  **');
            disp('************************************************************************');
        end

        num = length(id);

        % Parse out the important columns
        kode = uint16([model{6}]);
        % Model elements are - [x_start, y_start, x_fin, y_fin, lat_slip, dip_slip, dip, top, bottom]
        model_elements = [model{2:5} model{7:11}];
        model_elements = double(model_elements);

        % Parse grid parameters
        % Skip header line
        tmp = fgets(fid);
        tmp = fgets(fid);
        tmp = fgets(fid);
        grid = textscan(fid, '%*45c %16.7f32', 6);
        grid = double([grid{:}]);

        % Parse size parameters
        % Skip header line
        tmp = fgets(fid);
        tmp = fgets(fid);
        size = textscan(fid, '%*45c %16.7f32', 3);
        size = [size{:}];

        % Parse cross-section defaults
        % Skip header line
        tmp = fgets(fid);
        tmp = fgets(fid);
        tmp = fgets(fid);
        cross_section = textscan(fid, '%*45c %16.7f32', 7);
        cross_section = [cross_section{:}];
        if isempty(cross_section) == 1
            disp('   No info for a cross-section line is included in the input file.');
        end

        % Parse map information
        % Skip header line
        tmp = fgets(fid);
        tmp = fgets(fid);
        map_info = textscan(fid, '%*45c %16.7f32', 6);
        map_info = double([map_info{:}]);

        % Check if all parameters are properly read or not.
        if isempty(num) == 1
            errordlg('Number of modelling elements could not be read.');
            return;
        end
        if isempty(poisson_ratio) == 1
            errordlg('Poisson ratio could not be read.');
            return;
        end
        if isempty(youngs_mod) == 1
            errordlg('Youngs modulus could not be read.');
            return;
        end
        if isempty(coef_friction) == 1
            errordlg('Coefficient of friction could not be read.');
            return;
        end
        if isempty(regional_stress) == 1
            errordlg('Regional stress values could not be read.');
            return;
        end
        if isempty(grid) == 1
            errordlg('Grid info for study area could not be read properly.');
            return;
        end
    fclose (fid);

    x = grid(1):grid(5):grid(3);
    y = grid(2):grid(6):grid(4);
    z = calculation_depth;

    clear tmp;

    return