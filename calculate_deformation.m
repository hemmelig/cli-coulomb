function calculate_deformation(source_filename, output_filename, calculation_depth, mode)
    % Calculate the deformation induced for a number of model elements using
    % the equations of Okada, 1992, and write to file.
    %
    % Args:
    %   source_filename : str
    %       Name of the input file to be parsed.
    %   output_filename : str
    %       Name of output file containing calculated deformation.


    [x, y, ~, model_elements, kode, poisson_ratio, youngs_mod, ~, ~, ~] = read_input_file(source_filename);

    % Calculate in elastic half space of Okada (1992)
    [dc3d] = okada_elastic_halfspace(x, y, model_elements, youngs_mod, poisson_ratio, calculation_depth, kode, mode);

    fid = fopen(output_filename, 'wt');
        fprintf(fid, 'Input file selected: %s\n', source_filename);
        if mode == 'strain'
            fprintf(fid, 'x y z exx eyy ezz eyz exz exy dilatation\n');
            fprintf(fid, '(km) (km) (km) (-) (-) (-) (-) (-) (-) (-)\n');
        elseif mode == 'stress'
            fprintf(fid, 'x y z exx eyy ezz eyz exz exy dilatation\n');
            fprintf(fid, '(km) (km) (km) (-) (-) (-) (-) (-) (-) (-)\n');
        end
        for i = 1:size(dc3d, 1)
            fprintf(fid, '%18.12f', dc3d(i, 1:2));
            fprintf(fid, '%18.12f', calculation_depth);
            fprintf(fid, '%18.12f', dc3d(i, 9:14));
            fprintf(fid, '%18.12f', dc3d(i, 9) + dc3d(i, 10) + dc3d(i, 11));
            fprintf(fid, ' \n');
        end
    fclose(fid);
