function calculate_SHmax(source_filename, output_filename, calculation_depth)
    [x, y, ~, model_elements, kode, poisson_ratio, youngs_mod, ~, ~, ~] = read_input_file(source_filename);

    % Calculate in elastic half space of Okada (1992)
    [dc3d] = okada_elastic_halfspace(x, y, model_elements, youngs_mod, poisson_ratio, calculation_depth, kode, 'stress');

    % Stress tensor from Coulomb3
    ncells = length(dc3d(:, 1));
    A = zeros(ncells, 3, 3);
    A(:) = [ dc3d(:, 9), dc3d(:, 14), dc3d(:, 13);...
            dc3d(:, 14), dc3d(:, 10), dc3d(:, 12),;...
            dc3d(:, 13), dc3d(:, 12), dc3d(:, 11)];

    eigen_vals = zeros(ncells, 3);
    eigen_vecs = zeros(ncells, 3, 3);
    D = zeros(ncells, 3, 3);

    for i = 1:ncells
        eigen_vals(i, :) = eig(squeeze(A(i, :, :)));
        eigen_vals = -eigen_vals;
        [eigen_vecs(i,:,:), D(i,:,:)] = eig(squeeze(A(i, :, :)));
    end

    % From Lund and Townend, 2007
    % N=y, E=x
    s1n = eigen_vecs(:, 2, 1);
    s1e = eigen_vecs(:, 1, 1);
    s2n = eigen_vecs(:, 2, 2);
    s2e = eigen_vecs(:, 1, 2);

    R = (eigen_vals(:, 1) - eigen_vals(:, 2))./(eigen_vals(:, 1) - eigen_vals(:, 3));

    X = ((s1n.^2 - s1e.^2) + (1 - R).*(s2n.^2 - s2e.^2));
    Y = (2.*(s1n.*s1e + (1 - R).*s2n.*s2e));

    tan2alpha = Y./X;
    alpha = atan(tan2alpha)./2;

    % Use 2nd derivative to check min/max - negative -> max (wanted)
    dev2 = -2.0.*X.*cos(2.0.*alpha) - 2.0.*Y.*sin(2.0.*alpha);
    for i = 1:ncells
        if dev2(i) > 0
            % We found a minimum. Add 90 degrees to get the maximum.
            alpha(i) = alpha(i) + pi/2.0;
        end
    end

    % Calculate vectors in coordinate space for output
    xinc  = x(2) - x(1);
    yinc  = y(2) - y(1);
    ncell = length(x) * length(y);
    sl = sqrt(xinc*xinc + yinc*yinc)/4.0;
    
    angle = zeros(ncell, 1);
    angle(:, 1) = rad2deg(alpha);
    
    lg = 1.0;       % default
    shx = sl.*sin(deg2rad(angle)).*lg;
    shy = sl.*cos(deg2rad(angle)).*lg;

    xycoord = zeros(ncell, 2, 'double');
    for i = 1:length(x)
        xx = x(i);
        for j = 1:length(y)
            yy = y(j);
            nn = j + (i - 1) * length(y);
            xycoord(nn, 1) = xx;
            xycoord(nn, 2) = yy;
        end
    end

    xh1 = [rot90(xycoord(:, 1) - shx(:)); rot90(xycoord(:, 1) + shx(:))];
    yh1 = [rot90(xycoord(:, 2) - shy(:)); rot90(xycoord(:, 2) + shy(:))];

    % Save x1/y1 of SHmax vectors
    ars = zeros(1, length(xh1));
    for n = 1:length(ars)
        ars(n)='>';
    end
    xh2 = [xh1; ars];
    yh2 = [yh1; ars];
    xs = reshape(xh2, [], 1);
    ys = reshape(yh2, [], 1);
    SH = [xs ys];
    save(output_filename, 'SH', '-ascii', '-tabs');
