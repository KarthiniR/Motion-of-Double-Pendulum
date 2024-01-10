function double_pendulum_rk2()
    % Parameters
    g = 9.81; % acceleration due to gravity (m/s^2)
    m1 = 1; % mass of pendulum 1 (kg)
    m2 = 1; % mass of pendulum 2 (kg)
    L1 = 1; % length of rod 1 (m)
    L2 = 1; % length of rod 2 (m)

    % Time settings
    dt = 0.01; % time step (s)
    t_final = 10; % simulation time (s)
    t = 0:dt:t_final;

    % Initial conditions
    theta1_0 = pi/2; % initial angle of pendulum 1 (rad)
    theta2_0 = pi/4; % initial angle of pendulum 2 (rad)
    omega1_0 = 0; % initial angular velocity of pendulum 1 (rad/s)
    omega2_0 = 0; % initial angular velocity of pendulum 2 (rad/s)

    % RK2 integration
    state = [theta1_0; omega1_0; theta2_0; omega2_0];
    num_states = length(state);
    num_steps = length(t);
    states = zeros(num_states, num_steps);
    states(:, 1) = state;

    for i = 1:num_steps-1
        k1 = dt * pendulum_dynamics(state, m1, m2, L1, L2, g);
        k2 = dt * pendulum_dynamics(state + 0.5 * k1, m1, m2, L1, L2, g);
        state = state + k2;
        states(:, i+1) = state;
    end

    % Convert to Cartesian coordinates
    x1 = L1 * sin(states(1, :));
    y1 = -L1 * cos(states(1, :));
    x2 = x1 + L2 * sin(states(3, :));
    y2 = y1 - L2 * cos(states(3, :));

    % Plot results
    figure;
    plot(t, states(1, :), 'r', 'LineWidth', 1.5); hold on;
    plot(t, states(3, :), 'b', 'LineWidth', 1.5); hold off;
    xlabel('Time (s)');
    ylabel('Angle (rad)');
    legend('\theta_1', '\theta_2');
    title('Double Pendulum Angles vs Time');

    figure;
    plot(x1, y1, 'r', 'LineWidth', 1.5); hold on;
    plot(x2, y2, 'b', 'LineWidth', 1.5); hold off;
    xlabel('x (m)');
    ylabel('y (m)');
    legend('Pendulum 1', 'Pendulum 2');
    title('Double Pendulum Trajectories');
    axis equal;
end

function dstate = pendulum_dynamics(state, m1, m2, L1, L2, g)
    theta1 = state(1);
    omega1 = state(2);
    theta2 = state(3);
    omega2 = state(4);

    delta_theta = theta2 - theta1;
    M = m1 + m2;

    % Equations of motion (derived using Lagrange's equations)
    domega1 = (-(m2*L1*omega1^2*sin(delta_theta)*cos(delta_theta) + m2*g*sin(theta2)*cos(delta_theta) + m2*L2*omega2^2*sin(delta_theta) - M*g*sin(theta1)) / (L1*(M - m2*cos(delta_theta)^2)));
    domega2 = ((M*L1*omega1^2*sin(delta_theta)*cos(delta_theta) - m2*L1*omega1^2*sin(delta_theta) + M*g*sin(theta1)*cos(delta_theta) - M*L2*omega2^2*sin(delta_theta)*cos(delta_theta) + M*g*sin(theta2)) / (L2*(M - m2*cos(delta_theta)^2)));

    % State derivatives
    dtheta1 = omega1;
    dtheta2 = omega2;

    dstate = [dtheta1; domega1; dtheta2; domega2];
end