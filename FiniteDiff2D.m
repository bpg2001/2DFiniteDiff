% Parameters for Heat Map
Lx = 1.0; % Length of the box in x-direction
Ly = 1.0; % Length of the box in y-direction
nx = 50; % Number of grid points in x-direction
ny = 50; % Number of grid points in y-direction
dx = Lx / (nx - 1); % Grid spacing in x-direction
dy = Ly / (ny - 1); % Grid spacing in y-direction
dt = 0.0005; % Time step
nt = 1000; % Number of time steps

% Initial temperature distribution
T = zeros(ny, nx);

% Prompt user to pick boundary conditions for left, right, top and bottom.
% After picking Dirichlet or Neumann the user is prompted to pick
% temperature or heat flux value respectivley
disp('Choose boundary condition type for each boundary:');
disp('1. Dirichlet (fixed temperature)');
disp('2. Neumann (fixed heat flux)');

% Left boundary
left_bc = input('Left boundary condition type (1 or 2): ');
if left_bc == 1
    left_temp = input('Left boundary temperature: ');
else
    left_flux = input('Left boundary heat flux: ');
end

% Right boundary
right_bc = input('Right boundary condition type (1 or 2): ');
if right_bc == 1
    right_temp = input('Right boundary temperature: ');
else
    right_flux = input('Right boundary heat flux: ');
end

% Top boundary
top_bc = input('Top boundary condition type (1 or 2): ');
if top_bc == 1
    top_temp = input('Top boundary temperature: ');
else
    top_flux = input('Top boundary heat flux: ');
end

% Bottom boundary
bottom_bc = input('Bottom boundary condition type (1 or 2): ');
if bottom_bc == 1
    bottom_temp = input('Bottom boundary temperature: ');
else
    bottom_flux = input('Bottom boundary heat flux: ');
end

% Prompt user to pick the type of heat transfer to be tested
disp('Choose type of heat transfer to test:');
disp('1. Conduction');
disp('2. Convection');
disp('3. Radiation');
heat_transfer = input('Heat transfer type (1, 2, or 3): ');

% If conduction, user will be prompted to input thermal diffusivity
% If convection, user will be prompted to input velocity in x and y
% direction as well as heat transfer coefficient h
% If radiation, user will be prompted to input surrounding temperature
if heat_transfer == 1
    alpha = input('Enter thermal diffusivity: ');
elseif heat_transfer == 2
    vx = input('Enter convection velocity in x-direction: ');
    vy = input('Enter convection velocity in y-direction: ');
    h = input('Enter heat transfer coefficient: ');
elseif heat_transfer == 3
    sigma = 5.67e-8; % Stefan-Boltzmann constant
    T_surr = input('Enter surrounding temperature for radiation: ');
end

% Finite difference coefficients
% Coefficients for the diffusion terms (Dx, Dy) in the heat equation are 
% computed using the input thermal diffusivity and the grid spacings.
Dx = alpha * dt / dx^2;
Dy = alpha * dt / dy^2;

% Time-stepping loop. For each time step, a loop calculates the new temperature distribution
% Interior Update: The temperature at each internal grid point is updated using 
% the heat equation's finite difference approximation for conduction.
for n = 1:nt
    Tn = T;
    T(2:end-1, 2:end-1) = Tn(2:end-1, 2:end-1) + ...
        Dx * (Tn(2:end-1, 3:end) - 2 * Tn(2:end-1, 2:end-1) + Tn(2:end-1, 1:end-2)) + ...
        Dy * (Tn(3:end, 2:end-1) - 2 * Tn(2:end-1, 2:end-1) + Tn(1:end-2, 2:end-1));

    % Apply boundary conditions, temperature values at the boundaries are 
    % updated based on the user's input
    if left_bc == 1
        T(:, 1) = left_temp; % Left boundary (Dirichlet)
    else
        T(:, 1) = T(:, 2) - left_flux * dx / alpha; % Left boundary (Neumann)
    end
    
    if right_bc == 1
        T(:, end) = right_temp; % Right boundary (Dirichlet)
    else
        T(:, end) = T(:, end-1) + right_flux * dx / alpha; % Right boundary (Neumann)
    end
    
    if top_bc == 1
        T(1, :) = top_temp; % Top boundary (Dirichlet)
    else
        T(1, :) = T(2, :) - top_flux * dy / alpha; % Top boundary (Neumann)
    end
    
    if bottom_bc == 1
        T(end, :) = bottom_temp; % Bottom boundary (Dirichlet)
    else
        T(end, :) = T(end-1, :) + bottom_flux * dy / alpha; % Bottom boundary (Neumann)
    end
    
    % Convection: the temperature is further adjusted by the convection terms, 
    % calculated using a difference scheme that approximates convective derivatives.
    if heat_transfer == 2
        T = T + dt * (h * (vx * (circshift(T, [0 -1]) - T) / dx + ...
                              vy * (circshift(T, [-1 0]) - T) / dy));
    end
    
    % Radiation: the temperature is adjusted based on the radiation term, 
    % which involves the fourth power of the temperature difference between
    % the grid point and the environmental temperature.
    if heat_transfer == 3
        T = T - dt * sigma * (T.^4 - T_surr^4);
    end
end

% Plotting the final temperature distribution
figure;
contourf(T, 50, 'LineColor', 'none');
colormap hot;
colorbar;
title('Temperature Distribution');
xlabel('X');
ylabel('Y');
