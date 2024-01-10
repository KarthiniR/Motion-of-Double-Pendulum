function double_pendulum_final()
    % Constants
    g = 9.81;
    m1 = 1;
    m2 = 1;
    l1 = 1;
    l2 = 1;

    % Initial conditions
    theta1_0 = 0.2; %Initial angular positions
    theta2_0 = 0.4;
    theta1d_0 = 0.1;  %Initial Angular Velocities
    theta2d_0 = 0.3;

    % Time settings
    t0 = 0;
    tf = 10;
    td = 0.05; %Time step
    steps = floor((tf - t0) / td) + 1; %total no of time steps

    % Initialize arrays
    t = linspace(t0, tf, steps);
    theta1 = zeros(1, steps);
    theta2 = zeros(1, steps);
    theta1d = zeros(1, steps);
    theta2d = zeros(1, steps);

    % Set initial values for storing the angular positions and angular velocities
    theta1(1) = theta1_0;
    theta2(1) = theta2_0;
    theta1d(1) = theta1d_0;
    theta2d(1) = theta2d_0;

    % Initialize error arrays (GTE)
    gte_theta1_euler = zeros(1, steps);
    gte_theta2_euler = zeros(1, steps);
    gte_theta1_rk2 = zeros(1, steps);
    gte_theta2_rk2 = zeros(1, steps);
    gte_theta1_rk4 = zeros(1, steps);
    gte_theta2_rk4 = zeros(1, steps);
    gte_theta1d_euler = zeros(1, steps);
    gte_theta2d_euler = zeros(1, steps);
    gte_theta1d_rk2 = zeros(1, steps);
    gte_theta2d_rk2 = zeros(1, steps);
    gte_theta1d_rk4 = zeros(1, steps);
    gte_theta2d_rk4 = zeros(1, steps);


    % Reference solution
    steps_ref = 10 * steps; % You can adjust the factor for higher or lower resolution
    t_ref = linspace(t0, tf, steps_ref);

    % Solve ODEs for reference solution using a higher-order method (e.g., ode45)
    [~, Y_ref] = ode45(@(t, y) double_pendulum_ode(t, y, m1, m2, l1, l2, g), t_ref, [theta1_0, theta1d_0, theta2_0, theta2d_0]);
    theta1_ref = Y_ref(:, 1)';
    theta2_ref = Y_ref(:, 3)';

    % Solve ODEs using desired method (Euler, RK2 or RK4)
    for method = {'Euler', 'RK2', 'RK4'}
        % Reset initial values
        theta1(1) = theta1_0;
        theta2(1) = theta2_0;
        theta1d(1) = theta1d_0;
        theta2d(1) = theta2d_0;

        for i = 1:steps - 1
            switch method{1}
                case 'Euler'
                    %Calling Euler Method%
                    [d_theta1, d_theta2, d_theta1d, d_theta2d] = euler_method (theta1(i), theta2(i), theta1d(i), theta2d(i), td, g, m1, m2, l1, l2);
                case 'RK2'
                    %Calling heun/RK2 Method%
                    [d_theta1, d_theta2, d_theta1d, d_theta2d] = heun_method (theta1(i), theta2(i), theta1d(i), theta2d(i), td, g, m1, m2, l1, l2);
                case 'RK4'
                    %Calling RK4 Method%
                    [d_theta1, d_theta2, d_theta1d, d_theta2d] = rk4_method (theta1(i), theta2(i), theta1d(i), theta2d(i), td, g, m1, m2, l1, l2);
            end
            theta1(i+1) = theta1(i) + d_theta1;
            theta2(i+1) = theta2(i) + d_theta2;
            theta1d(i+1) = theta1d(i) + d_theta1d;
            theta2d(i+1) = theta2d(i) + d_theta2d;
        end

        % Resample the solution    %resampling the reference solution to match the time steps of the numerical solutions obtained using Euler, RK2, and RK4 methods.
        theta1_resampled = interp1(linspace(t0, tf, steps_ref), theta1_ref, t);
        theta2_resampled = interp1(linspace(t0, tf, steps_ref), theta2_ref, t);
        % Resample the solution for angular velocities
        theta1d_resampled = interp1(linspace(t0, tf, steps_ref), Y_ref(:, 2)', t);
        theta2d_resampled = interp1(linspace(t0, tf, steps_ref), Y_ref(:, 4)', t);

    % Calculate the global truncation error
        switch method{1}
            case 'Euler'
                gte_theta1_euler = calculate_global_truncation_error(theta1, theta1_resampled);
                gte_theta2_euler = calculate_global_truncation_error(theta2, theta2_resampled);
                gte_theta1d_euler = calculate_global_truncation_error(theta1d, theta1d_resampled);
                gte_theta2d_euler = calculate_global_truncation_error(theta2d, theta2d_resampled);
            case 'RK2'  
                gte_theta1_rk2 = calculate_global_truncation_error(theta1, theta1_resampled);
                gte_theta2_rk2 = calculate_global_truncation_error(theta2, theta2_resampled);
                gte_theta1d_rk2  = calculate_global_truncation_error(theta1d, theta1d_resampled );
                gte_theta2d_rk2  = calculate_global_truncation_error(theta2d, theta2d_resampled );
            case 'RK4'
                gte_theta1_rk4 = calculate_global_truncation_error(theta1, theta1_resampled);
                gte_theta2_rk4 = calculate_global_truncation_error(theta2, theta2_resampled);
                gte_theta1d_rk4 = calculate_global_truncation_error(theta1d, theta1d_resampled);
                gte_theta2d_rk4 = calculate_global_truncation_error(theta2d, theta2d_resampled);
        end
       
        % Plot results for the current method
        figure;
        sgtitle(sprintf('Results for %s method', method{1}));
        subplot(2,1,1);
        plot(t,theta1,'r',t,theta2,'b');
        title('Angular position');
        xlabel('Time (s)');
        ylabel('Angle (rad)');
        legend('Theta 1', 'Theta 2');

        subplot(2,1,2);
        plot(t,theta1d,'r',t,theta2d,'b');
        title('Angular velocity');
        xlabel('Time (s)');
        ylabel('Angular velocity (rad/s)');
        legend('Theta_d 1', 'Theta_d 2');
    end

    % Calculate stability using maximum allowable step size
    A = [     0, 0, 1, 0;                                              % The matrix A is 4x4 and represents the state space form of the linearized double pendulum system.
              0,0,0,1;
          -(m1 + m2)*g/l1,m2*g/l1,0,0;
            m1*g/l2,-(m1 + m2)*g/l2,0,0];
    [~, D] = eig(A);                                                   % the eigenvalues of the matrix A  are computed and stored them in the diagonal matrix D  
    eig_values = diag(D);
    max_eig_value = max(abs(eig_values));

    euler_stability = 2 / max_eig_value;
    rk2_stability = pi / max_eig_value;
    rk4_stability = 2 * sqrt(2) / max_eig_value;

    % Print the maximum allowable step size for each method
    fprintf('Euler method max step size: %.10f\n', euler_stability);
    fprintf('RK2 method max step size: %.10f\n', rk2_stability);
    fprintf('RK4 method max step size: %.10f\n', rk4_stability);


% Plot the global truncation error versus time for each method
figure;
subplot(3,1,1);
plot(t, gte_theta1_euler, 'r');
title('Global truncation error for theta1 (Euler)');
xlabel('Time (s)');
ylabel('Error');

subplot(3,1,2);
plot(t, gte_theta1_rk2, 'g');
title('Global truncation error for theta1 (RK2)');
xlabel('Time (s)');
ylabel('Error');

subplot(3,1,3);
plot(t, gte_theta1_rk4, 'b');
title('Global truncation error for theta1 (RK4)');
xlabel('Time (s)');
ylabel('Error');

figure;
subplot(3,1,1);
plot(t, gte_theta2_euler, 'r');
title('Global truncation error for theta2 (Euler)');
xlabel('Time (s)');
ylabel('Error');

subplot(3,1,2);
plot(t, gte_theta2_rk2, 'g');
title('Global truncation error for theta2 (RK2)');
xlabel('Time (s)');
ylabel('Error');

subplot(3,1,3);
plot(t, gte_theta2_rk4, 'b');
title('Global truncation error for theta2 (RK4)');
xlabel('Time (s)');
ylabel('Error');


% Plot the global truncation error versus time for each method - Angular velocities
figure;
subplot(3,1,1);
plot(t, gte_theta1d_euler, 'r');
title('Global truncation error for theta1d (Euler)');
xlabel('Time (s)');
ylabel('Error');

subplot(3,1,2);
plot(t, gte_theta1d_rk2, 'g');
title('Global truncation error for theta1d (RK2)');
xlabel('Time (s)');
ylabel('Error');

subplot(3,1,3);
plot(t, gte_theta1d_rk4, 'b');
title('Global truncation error for theta1d (RK4)');
xlabel('Time (s)');
ylabel('Error');

figure;
subplot(3,1,1);
plot(t, gte_theta2d_euler, 'r');
title('Global truncation error for theta2d (Euler)');
xlabel('Time (s)');
ylabel('Error');

subplot(3,1,2);
plot(t, gte_theta2d_rk2, 'g');
title('Global truncation error for theta2d (RK2)');
xlabel('Time (s)');
ylabel('Error');

subplot(3,1,3);
plot(t, gte_theta2d_rk4, 'b');
title('Global truncation error for theta2d (RK4)');
xlabel('Time (s)');
ylabel('Error');


end


% ODE45 solution
function dydt = double_pendulum_ode(t, y, m1, m2, l1, l2, g)
    theta1 = y(1);
    theta1d = y(2);
    theta2 = y(3);
    theta2d = y(4);

    delta_theta = theta2-theta1;

    M = [-(m1 + m2)*g*sin(theta1)-m2*l1*theta1d^2*sin(delta_theta),m2*l2*theta2d^2 * sin(delta_theta) + m2 * g * sin(theta2);
         
    
    m1 * g * sin(theta1) + m1 * l1 * theta1d^2 * sin(delta_theta), -(m1 + m2) * g * sin(theta2) - m1 * l2 * theta2d^2 * sin(delta_theta)];

    det_M = M(1, 1) * M(2, 2) - M(1, 2) * M(2, 1);
    inv_M = [M(2, 2), -M(1, 2); -M(2, 1), M(1, 1)] / det_M;

    C = [2*l1*l2*theta1d*theta2d*cos(delta_theta)*sin(delta_theta);
         -l1 * l2 * theta1d^2 * cos(delta_theta)*sin(delta_theta)];

    qdd = inv_M*(-C);

    dydt = [theta1d; qdd(1); theta2d; qdd(2)];
end

function gte = calculate_global_truncation_error(approx_solution, ref_solution)
    gte = abs(ref_solution - approx_solution);
end

function [d_theta1,d_theta2,d_theta1d,d_theta2d] = euler_method(theta1,theta2, theta1d, theta2d, dt, g, m1, m2, l1, l2)
% Calculate derivatives
[theta1dd,theta2dd] = double_pendulum_acceleration(theta1,theta2, theta1d, theta2d,g, m1, m2, l1, l2);
% Update using Euler's method
d_theta1 = dt * theta1d;
d_theta2 = dt * theta2d;                                                        %Angular positions calculations using Euler
d_theta1d = dt * theta1dd;
d_theta2d = dt * theta2dd;                                                      % Angular Velocity calculations using Euler
end

function [d_theta1,d_theta2,d_theta1d,d_theta2d] = heun_method(theta1,theta2,theta1d,theta2d, dt, g, m1, m2, l1, l2)
    % First step (Euler's method)
    [theta1dd_1, theta2dd_1] = double_pendulum_acceleration(theta1, theta2, theta1d, theta2d, g, m1, m2, l1, l2);
    k1_theta1 = dt * theta1d;
    k1_theta2 = dt * theta2d;
    k1_theta1d = dt * theta1dd_1;
    k1_theta2d = dt * theta2dd_1;                                 % Calculates angular acc, ang positions and ang. velocities same as Euler method

    % Second step
    [theta1dd_2, theta2dd_2] = double_pendulum_acceleration(theta1 + k1_theta1, theta2 + k1_theta2, theta1d + k1_theta1d, theta2d + k1_theta2d, g, m1, m2, l1, l2);
    k2_theta1 = dt*(theta1d + k1_theta1d);
    k2_theta2 = dt*(theta2d + k1_theta2d);
    k2_theta1d = dt*theta1dd_2;
    k2_theta2d = dt*theta2dd_2;                                             % Calculates again using updated values from the first step 

    % Update using Heun's method
    d_theta1 = 0.5*(k1_theta1 + k2_theta1);                                  % Performing weighted avg of the two estimates
    d_theta2 = 0.5*(k1_theta2 + k2_theta2);
    d_theta1d = 0.5*(k1_theta1d + k2_theta1d);
    d_theta2d = 0.5*(k1_theta2d + k2_theta2d);
end


function [d_theta1,d_theta2,d_theta1d,d_theta2d] = rk4_method(theta1,theta2,theta1d,theta2d, dt, g, m1, m2, l1, l2)
% First step
[theta1dd_1,theta2dd_1] = double_pendulum_acceleration(theta1,theta2,theta1d,theta2d, g, m1, m2, l1, l2);
k1_theta1 = dt*theta1d;
k1_theta2 = dt*theta2d;
k1_theta1d = dt* theta1dd_1;
k1_theta2d = dt * theta2dd_1;                                               % Same as Euler method
% Second step
[theta1dd_2,theta2dd_2] = double_pendulum_acceleration(theta1+0.5 *k1_theta1, theta2 + 0.5 * k1_theta2,theta1d + 0.5 * k1_theta1d, theta2d + 0.5 * k1_theta2d,g,m1,m2,l1,l2);
k2_theta1 = dt*(theta1d + 0.5*k1_theta1d);
k2_theta2 = dt*(theta2d + 0.5*k1_theta2d);
k2_theta1d = dt*theta1dd_2;
k2_theta2d = dt*theta2dd_2;
% Third step
[theta1dd_3,theta2dd_3] = double_pendulum_acceleration(theta1 + 0.5 * k2_theta1,theta2 + 0.5 * k2_theta2,theta1d + 0.5 * k2_theta1d,theta2d + 0.5 * k2_theta2d,g,m1,m2,l1,l2);
k3_theta1 = dt*(theta1d+0.5*k2_theta1d);
k3_theta2 = dt*(theta2d +0.5*k2_theta2d);
k3_theta1d = dt*theta1dd_3;
k3_theta2d = dt*theta2dd_3;

% Fourth step
[theta1dd_4,theta2dd_4] = double_pendulum_acceleration(theta1 + k3_theta1, theta2 + k3_theta2, theta1d + k3_theta1d, theta2d + k3_theta2d,g,m1,m2,l1,l2);
k4_theta1 = dt * (theta1d + k3_theta1d);
k4_theta2 = dt * (theta2d + k3_theta2d);
k4_theta1d = dt * theta1dd_4;
k4_theta2d = dt * theta2dd_4;

% Update using RK4 method
d_theta1 = (k1_theta1+2*k2_theta1+2*k3_theta1+k4_theta1)/6;                     % Weighted avf of the estimates
d_theta2 = (k1_theta2+2*k2_theta2 +2* k3_theta2 +k4_theta2)/6;
d_theta1d = (k1_theta1d+ 2 * k2_theta1d+ 2*k3_theta1d +k4_theta1d) / 6;
d_theta2d = (k1_theta2d + 2 * k2_theta2d + 2 * k3_theta2d +k4_theta2d) /6;
end

function [theta1dd,theta2dd]=double_pendulum_acceleration(theta1,theta2,  theta1d,theta2d,g,m1, m2,l1, l2)
delta =theta2-theta1;
M = m1+ m2;
num1 = -g*(2 * m1 + m2) *sin(theta1) - m2 * g * sin(theta1 - 2* theta2) - 2* sin(delta)* m2 *(theta2d^2 * l2 + theta1d^2*l1 *cos(delta));
den1 = l1 * (2 *m1 +m2 - m2*cos(2*delta));
theta1dd = num1 / den1;

num2 = 2 * sin(delta) * (theta1d^2 *l1 * M - g * (m1 + m2) * cos(theta1) - theta2d^2 * l2 * m2 * cos(delta));
den2 = l2 * (2 *m1 + m2 - m2 *cos(2 * delta));
theta2dd = num2/ den2;
end

