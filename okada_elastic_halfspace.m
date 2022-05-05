function [dc3d] = okada_elastic_halfspace(xgrid, ygrid, element, young, pois, calcdepth, kode, mode)
    % Calculate the resultant deformation of a given coordinate space using the formulation of Okada 1992.
    %
    %   Args:
    %       xgrid : array-like
    %           X-coordinates of volume of interest.
    %       ygrid : array-like
    %           Y-coordinates of volume of interest.
    %       element : array-like
    %           young:  young's modulus
    %           pois:   Poisson's ratio
    %           calcdepth: calculation depth (km, positive)
    %           kode:   number to specify the calc. type
    %
    %   CALL-OUT FUNCTION:
    %       ***** source fault coordinates *****
    %        Okada_DC3D: function to calculate rectangular source in Okada 1992
    %       Okada_DC3D0: function to calculate point source in Okada 1992
    
    alpha = 1.0/(2.0*(1.0 - pois));
    zd = calcdepth*(-1.0);
    
    if numel(xgrid) == 1 % one point calc.
        ncell = 1;
    else
        ncell = length(xgrid) * length(ygrid);
    end
    
    ux   = zeros(ncell, 1, 'double');
    uy   = zeros(ncell, 1, 'double');
    uz   = zeros(ncell, 1, 'double');
    uxx  = zeros(ncell, 1, 'double');
    uyx  = zeros(ncell, 1, 'double');
    uzx  = zeros(ncell, 1, 'double');
    uxy  = zeros(ncell, 1, 'double');
    uyy  = zeros(ncell, 1, 'double');
    uzy  = zeros(ncell, 1, 'double');
    uxz  = zeros(ncell, 1, 'double');
    uyz  = zeros(ncell, 1, 'double');
    uzz  = zeros(ncell, 1, 'double');
    iret = zeros(ncell, 1);
    xycoord = zeros(ncell, 2, 'double');
    
    dc3d  = zeros(ncell, 14, 'double');
    
    for i = 1:size(element, 1)
        depth = (element(i, 8) + element(i, 9))/2.0;  % depth should be positive
        for j = 1:length(xgrid)
            xx = xgrid(j);
            for k = 1:length(ygrid)
                yy = ygrid(k);
                nn = k + (j - 1)*length(ygrid);
                xycoord(nn, 1) = xx;
                xycoord(nn, 2) = yy;
            end
        end
        [c1, c2, c3, c4] = coord_conversion(xycoord(:, 1), xycoord(:, 2),...
                                            element(i, 1), element(i, 2),...
                                            element(i, 3), element(i, 4),...
                                            element(i, 8), element(i, 9),...
                                            element(i, 7));
        format long;
        aa = zeros(ncell, 1, 'double') + double(alpha);
        x  = zeros(ncell, 1, 'double') + double(c1);
        y  = zeros(ncell, 1, 'double') + double(c2);
        zz = zeros(ncell, 1, 'double') + double(zd);
        dp = zeros(ncell, 1, 'double') + double(depth);
        e7 = zeros(ncell, 1, 'double') + double(element(i, 7));
        al = zeros(ncell, 1, 'double') + double(c3);
        aw = zeros(ncell, 1, 'double') + double(c4);
        if kode(i) == 200
            e5 = zeros(ncell, 1, 'double') - double(element(i, 5));
            e6 = zeros(ncell, 1, 'double');
            zr = zeros(ncell, 1, 'double') + double(element(i, 6));
        elseif kode(i) == 300
            e5 = zeros(ncell, 1, 'double');
            e6 = zeros(ncell, 1, 'double') + double(element(i, 6));
            zr = zeros(ncell, 1, 'double') + double(element(i, 5));
        end
        if kode(i) == 400
            aw = zeros(ncell, 1, 'double') - double(element(i, 5));
            e5 = zeros(ncell, 1, 'double') + double(element(i, 6));
            e6 = zeros(ncell, 1, 'double');
            zr = zeros(ncell, 1, 'double');
        elseif kode(i) == 500
            aw = zeros(ncell, 1, 'double');
            e5 = zeros(ncell, 1, 'double');
            e6 = zeros(ncell, 1, 'double') + double(element(i, 5));
            zr = zeros(ncell, 1, 'double') + double(element(i, 6));
        end
    
        a = [aa x y zz dp e7 al al aw aw e5 e6 zr];

        if kode(i) == 100 | kode(i) == 200 | kode(i) == 300
            [ ux(:, 1),  uy(:, 1),  uz(:, 1),...
             uxx(:, 1), uyx(:, 1), uzx(:, 1),...
             uxy(:, 1), uyy(:, 1), uzy(:, 1),...
             uxz(:, 1), uyz(:, 1), uzz(:, 1),...
             iret(:, 1)] = Okada_DC3D( a(:, 1),  a(:, 2),  a(:, 3),...
                                       a(:, 4),  a(:, 5),  a(:, 6),...
                                       a(:, 7),  a(:, 8),  a(:, 9),...
                                      a(:, 10), a(:, 11), a(:, 12),...
                                      a(:, 13));
        elseif kode(i) == 400 | kode(i) == 500
            [ ux(:, 1),  uy(:, 1),  uz(:, 1),...
             uxx(:, 1), uyx(:, 1), uzx(:, 1),...
             uxy(:, 1), uyy(:, 1), uzy(:, 1),...
             uxz(:, 1), uyz(:, 1), uzz(:, 1),...
             iret(:, 1)] = Okada_DC3D0( a(:, 1),  a(:, 2),  a(:, 3),...
                                        a(:, 4),  a(:, 5),  a(:, 6),...
                                       a(:, 10), a(:, 11), a(:, 12),...
                                       a(:, 13));
        end
        
        x = a(:, 2);
        y = a(:, 3);
        z = a(:, 4);
        
        %-- Displacement Conversion from Okada's field to Given field
        sw = sqrt((element(i, 4) - element(i, 2))^2 + (element(i, 3) - element(i, 1))^2);
        sina = (element(i, 4) - element(i, 2))/double(sw);
        cosa = (element(i, 3) - element(i, 1))/double(sw);
        uxg = ux*cosa - uy*sina;
        uyg = ux*sina + uy*cosa;
        uzg = uz;

        if kode(i) == 400 | kode(i) == 500
            uxg = uxg./1000.0;
            uyg = uyg./1000.0;
            uzg = uzg./1000.0;
        end
    
        %-- Strain to Stress for the normal component
        % strain to stress
        sk = young/(1.0 + pois);
        gk = pois/(1.0 - 2.0*pois);
        vol = uxx + uyy + uzz;

        if mode == 'stress'
            sxx = sk*(gk*vol + uxx)*0.001;
            syy = sk*(gk*vol + uyy)*0.001;
            szz = sk*(gk*vol + uzz)*0.001;
            sxy = (young/(2.0*(1.0 + pois)))*(uxy + uyx)*0.001;
            sxz = (young/(2.0*(1.0 + pois)))*(uxz + uzx)*0.001;
            syz = (young/(2.0*(1.0 + pois)))*(uyz + uzy)*0.001;
        elseif mode == 'strain'
            sxx = uxx*0.001;
            syy = uyy*0.001;
            szz = uzz*0.001;
            sxy = uxy*0.001;
            sxz = uxz*0.001;
            syz = uyz*0.001;
        end

        ssxx = reshape(sxx, 1, ncell);
        ssyy = reshape(syy, 1, ncell);
        sszz = reshape(szz, 1, ncell);
        ssxy = reshape(sxy, 1, ncell);
        ssxz = reshape(sxz, 1, ncell);
        ssyz = reshape(syz, 1, ncell);
        
        s0 = [ssxx; ssyy; sszz; ssyz; ssxz; ssxy];

        %-- Strain Conversion from Okada's field to Given field
        s1 = tensor_trans(sina, cosa, s0, ncell);
    
        sxx = reshape(s1(1, :), ncell, 1);
        syy = reshape(s1(2, :), ncell, 1);
        szz = reshape(s1(3, :), ncell, 1);
        syz = reshape(s1(4, :), ncell, 1);
        sxz = reshape(s1(5, :), ncell, 1);
        sxy = reshape(s1(6, :), ncell, 1); 

        if i == 1
            dc3d0 = horzcat(xycoord, x, y, z, uxg, uyg, uzg, sxx, syy, szz, syz, sxz, sxy);
        else
            dc3d0 = horzcat(zeros(ncell, 1), zeros(ncell, 1),zeros(ncell, 1),...
                            zeros(ncell, 1), zeros(ncell, 1),...
                            uxg, uyg, uzg, sxx, syy, szz, syz, sxz, sxy);
        end
        dc3d = dc3d + dc3d0;
    end
end
