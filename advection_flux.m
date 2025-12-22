function advH = advection_flux(lat, lon, T, u, v, dx, dy)
% ADVECTION_FLUX Calculate horizontal advection using flux formulation
%
% This method considers the four faces of each grid cell and calculates
% advection through mass fluxes, following Equation 2 from Marin et al. (2022).
% This approach removes dependence on zero-temperature reference.
%
% Inputs:
%   lat - latitude grid (degrees) [ny x nx]
%   lon - longitude grid (degrees) [ny x nx]
%   T   - temperature field (°C or K) [ny x nx]
%   u   - zonal velocity at cell faces (m/s) [ny x nx]
%   v   - meridional velocity at cell faces (m/s) [ny x nx]
%   dx  - zonal grid spacing (m) [ny x nx] or scalar
%   dy  - meridional grid spacing (m) [ny x nx] or scalar
%
% Output:
%   advH - horizontal advection (°C/s or K/s)
%
% Formula (horizontal part of Eq. 2):
%   AdvH = (1/V) * [∫∫_W u_W(T_W - T)dydz - ∫∫_E u_E(T_E - T)dydz +
%                   ∫∫_S v_S(T_S - T)dxdz - ∫∫_N v_N(T_N - T)dxdz]
%
% where W,E,S,N represent western, eastern, southern, northern faces
% and T is the volume-averaged temperature
%
% Reference: 
%   Marin et al. (2022) Front. Clim. 4:788390
%   Based on Lee et al. (2004) J. Phys. Oceanogr. 34:1936-1944

    [ny, nx] = size(T);
    
    % Initialize output
    advH = zeros(ny, nx);
    
    % Handle scalar dx, dy
    if isscalar(dx)
        dx = ones(ny, nx) * dx;
    end
    if isscalar(dy)
        dy = ones(ny, nx) * dy;
    end
    
    % Calculate cell volume (assuming unit depth for 2D)
    % For depth-integrated, multiply by depth (dz)
    dz = 1; % Unit depth for surface layer
    V = dx .* dy .* dz;
    
    % Interpolate temperature to cell faces
    % Western face (i-1/2, j)
    T_W = zeros(ny, nx);
    T_W(:, 2:end) = 0.5 * (T(:, 1:end-1) + T(:, 2:end));
    T_W(:, 1) = T(:, 1); % Boundary condition
    
    % Eastern face (i+1/2, j)
    T_E = zeros(ny, nx);
    T_E(:, 1:end-1) = 0.5 * (T(:, 1:end-1) + T(:, 2:end));
    T_E(:, end) = T(:, end); % Boundary condition
    
    % Southern face (i, j-1/2)
    T_S = zeros(ny, nx);
    T_S(2:end, :) = 0.5 * (T(1:end-1, :) + T(2:end, :));
    T_S(1, :) = T(1, :); % Boundary condition
    
    % Northern face (i, j+1/2)
    T_N = zeros(ny, nx);
    T_N(1:end-1, :) = 0.5 * (T(1:end-1, :) + T(2:end, :));
    T_N(end, :) = T(end, :); % Boundary condition
    
    % Calculate flux through each face
    % Note: u and v should already be on cell faces
    % If they're on cell centers, interpolate similarly
    
    % Flux through western face
    flux_W = zeros(ny, nx);
    flux_W(:, 2:end) = u(:, 2:end) .* (T_W(:, 2:end) - T(:, 2:end)) .* dy(:, 2:end) .* dz;
    
    % Flux through eastern face
    flux_E = zeros(ny, nx);
    flux_E(:, 1:end-1) = u(:, 1:end-1) .* (T_E(:, 1:end-1) - T(:, 1:end-1)) .* dy(:, 1:end-1) .* dz;
    
    % Flux through southern face
    flux_S = zeros(ny, nx);
    flux_S(2:end, :) = v(2:end, :) .* (T_S(2:end, :) - T(2:end, :)) .* dx(2:end, :) .* dz;
    
    % Flux through northern face
    flux_N = zeros(ny, nx);
    flux_N(1:end-1, :) = v(1:end-1, :) .* (T_N(1:end-1, :) - T(1:end-1, :)) .* dx(1:end-1, :) .* dz;
    
    % Calculate horizontal advection
    % Convention: positive flux = inflow
    advH = (1./V) .* (flux_W - flux_E + flux_S - flux_N);
    
end
