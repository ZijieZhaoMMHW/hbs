function adv = advection_gradient(lat, lon, T, u, v)
% ADVECTION_GRADIENT Calculate horizontal advection using gradient method
%
% This method assumes each grid cell is a cube and calculates advection
% using central differences of temperature gradients.
%
% Inputs:
%   lat - latitude vector or grid (degrees)
%   lon - longitude vector or grid (degrees)
%   T   - temperature field (same size as lat/lon grid)
%   u   - zonal velocity (eastward, m/s)
%   v   - meridional velocity (northward, m/s)
%
% Output:
%   adv - horizontal advection (°C/s or K/s)
%
% Formula: adv = -u*dT/dx - v*dT/dy
%
% Author: Based on user's original code
% Reference: Marin et al. (2022) Front. Clim. 4:788390

    % Create meshgrid if input is vectors
    if isvector(lat) && isvector(lon)
        [lon_grid, lat_grid] = meshgrid(lon, lat);
    else
        lon_grid = lon;
        lat_grid = lat;
    end
    
    % Calculate temperature gradients
    % cdtgradient handles spherical coordinates properly
    [dTdx, dTdy] = cdtgradient(lat_grid, lon_grid, T);
    
    % Calculate advection: -u*dT/dx - v*dT/dy
    % Negative sign because advection opposes the gradient
    adv = -u .* dTdx - v .* dTdy;
    
end
