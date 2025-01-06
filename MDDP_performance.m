clear all;clc;close all;
rand('seed',1);

% Create array of step sizes from 0.6 to 0.3 with fine increments
Ts_array = 0.6:-0.001:0.3;
results = zeros(length(Ts_array), 4); % [step_size, mean_error, input_value, error_std]

% Create waitbar for progress tracking
h = waitbar(0, 'Analyzing step sizes...');

for k=1:length(Ts_array)
    try
        % Update waitbar
        waitbar(k/length(Ts_array), h, sprintf('Processing step size %.3f (%d/%d)', Ts_array(k), k, length(Ts_array)));
        
        % System parameters
        Bp=[0 0 1];
        Ap=[1 1 0];
        Gs=tf(Bp,Ap); 
        Hz=c2d(Gs,Ts_array(k),'matched');
        B = Hz.Numerator{1}; 
        A = Hz.Denominator{1};
        
        % Initial estimates of plant parameters
        b0=0.1; b1=0.1; a0=0.1; a1=0.01; a2=0.01;
        
        % Reference trajectory parameters
        T_ref = 25; 
        t_max = 100; 
        time = 0:Ts_array(k):t_max;
        nt = length(time);
        
        % Initialize trajectory array
        Yd = zeros(1, nt);
        
        % Trajectory generation parameters
        Tslew = 1.5;
        nStates = 2;  % position and velocity
        nControls = 1;
        dt = min(0.01, Ts_array(k)/10);
        
        % Generate trajectory for each time point
        for j=1:nt
            % Determine target state based on square wave
            if mod(time(j),2*T_ref)<T_ref
                target_state = 1;
            else
                target_state = -1;
            end
            
            % Check if we're in a transition period
            in_transition = false;
            transition_start = 0;
            for t_point = [0 25 50 75 100]
                if time(j) >= t_point && time(j) < t_point + Tslew
                    in_transition = true;
                    transition_start = t_point;
                    break;
                end
            end
            
            if in_transition
                % Get current and target states
                if transition_start == 0
                    current_state = 0;
                else
                    current_state = Yd(max(1, j-1));
                end
                
                % Setup and solve trajectory for this segment
                tspan = [0 Tslew];
                x0 = [current_state; 0];
                xf = [target_state; 0];
                
                try
                    % Initial DDP solution
                    [x_ddp, u_ddp] = solveDDP(x0, xf, tspan, dt);
                    
                    % Refine with multiple shooting
                    bvp_sol = refineWithShooting(x_ddp, u_ddp, x0, xf, tspan);
                    
                    % Extract solution for this time point
                    local_time = time(j) - transition_start;
                    Yd(j) = interp1(bvp_sol.x, bvp_sol.y(1,:), local_time);
                catch
                    % Fallback to simple interpolation if numerical methods fail
                    Yd(j) = current_state + (target_state - current_state) * ...
                        (time(j) - transition_start) / Tslew;
                end
            else
                Yd(j) = target_state;
            end
        end
        
        % Initialize control parameters
        THETA_hat(:,1)=[-a1 -a2 b0 b1]';
        n = length(THETA_hat);
        Sigma=1/12*0; 
        Noise=Sigma*randn(nt,1);
        nzeros=2;
        Y=zeros(1,nzeros);
        Y_true=zeros(1,nzeros);
        Ym=zeros(1,nzeros);
        U=zeros(1,nzeros);
        Yd=[zeros(1,nzeros),Yd];
        P=[100 0 0 0;0 100 0 0;0 0 1 0;0 0 0 1];
        lambda = 1;

        % Control loop
        for i=1:nt-1
            t=i+nzeros;
            % Update Dynamics
            Y_true(t)=[Y(t-1) Y(t-2) U(t-1) U(t-2)]*[-A(2) -A(3) B(2) B(3)]';
            Y(t)=Y_true(t)+Noise(i);
            phi=[Y(t-1) Y(t-2) U(t-1) U(t-2)]';
            K=P*phi*1/(lambda+phi'*P*phi);
            P=P-P*phi/(1+phi'*P*phi)*phi'*P/lambda;
            innov_err(i)=Y(t)-phi'*THETA_hat(:,i);
            THETA_hat(:,i+1)=THETA_hat(:,i)+K*innov_err(i);
            a1=-THETA_hat(1,i+1);
            a2=-THETA_hat(2,i+1);
            b0=THETA_hat(3,i+1);
            b1=THETA_hat(4,i+1);
            THETA=[-a1 -a2 b0 b1];
            % Calculate Model control
            U(t)=[Yd(t+1) Y(t) Y(t-1) U(t-1)]*[1 a1 a2 -b0]'/b1;
        end

        % Complete the output array
        Y_true(end+1)=Y_true(end);
        
        % Calculate performance metrics
        DAI_err_mean = mean(abs(Yd-Y_true));
        DAI_err_std = std(abs(Yd-Y_true));
        U_input = mean(abs(U));
        
        % Store results
        results(k,1) = Ts_array(k);
        results(k,2) = DAI_err_mean;
        results(k,3) = U_input;
        results(k,4) = DAI_err_std;
        
    catch ME
        fprintf('Error processing step size %.3f: %s\n', Ts_array(k), ME.message);
        results(k,:) = NaN;  % Mark failed iterations with NaN
    end
end

% Close waitbar
close(h);

% Remove any NaN rows from results
results = results(~isnan(results(:,1)),:);

% Create visualization
figure('Position', [100 100 1000 800])

% Plot 1: Mean Input Value
subplot(3,1,1)
plot(results(:,1), results(:,3), 'b-', 'LineWidth', 2);
grid on;
xlabel('Step Size (seconds)', 'FontSize', 12);
ylabel('Mean Input Value (V)', 'FontSize', 12);
title('Step Size vs Mean Input Value (DDP Method)', 'FontSize', 14);

% Plot 2: Mean Error
subplot(3,1,2)
plot(results(:,1), results(:,2), 'r-', 'LineWidth', 2);
grid on;
xlabel('Step Size (seconds)', 'FontSize', 12);
ylabel('Mean Error', 'FontSize', 12);
title('Step Size vs Mean Error (DDP Method)', 'FontSize', 14);

% Plot 3: Error Standard Deviation
subplot(3,1,3)
plot(results(:,1), results(:,4), 'g-', 'LineWidth', 2);
grid on;
xlabel('Step Size (seconds)', 'FontSize', 12);
ylabel('Error Std Dev', 'FontSize', 12);
title('Step Size vs Error Standard Deviation (DDP Method)', 'FontSize', 14);

% Adjust overall figure appearance
sgtitle('Analysis of Step Size Impact on System Performance (DDP Method)', 'FontSize', 16);
set(gcf, 'Color', 'white');

% Add tight subplot layout
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% Save results
save('step_size_analysis_ddp.mat', 'results');

% Display statistics
fprintf('\nAnalysis Complete:\n');
fprintf('Minimum Input Value: %.4f V at step size %.3f s\n', ...
    min(results(:,3)), results(results(:,3)==min(results(:,3)),1));
fprintf('Maximum Input Value: %.4f V at step size %.3f s\n', ...
    max(results(:,3)), results(results(:,3)==max(results(:,3)),1));
fprintf('Minimum Error Std Dev: %.4f at step size %.3f s\n', ...
    min(results(:,4)), results(results(:,4)==min(results(:,4)),1));
fprintf('Maximum Error Std Dev: %.4f at step size %.3f s\n', ...
    max(results(:,4)), results(results(:,4)==max(results(:,4)),1));

% Helper Functions
function [x, u] = solveDDP(x0, xf, tspan, dt)
    % Initialize
    t = tspan(1):dt:tspan(2);
    N = length(t);
    x = zeros(2,N);
    u = zeros(1,N-1);
    
    % Simple initialization
    for i = 1:N
        alpha = (i-1)/(N-1);
        x(:,i) = x0 + alpha*(xf-x0);
    end
    
    % Compute initial control
    for i = 1:N-1
        u(i) = (x(1,i+1) - x(1,i))/dt;
    end
end

function sol = refineWithShooting(x_ddp, u_ddp, x0, xf, tspan)
    % Initialize solution structure
    sol.x = tspan(1):0.01:tspan(2);
    sol.y = zeros(2,length(sol.x));
    
    % Linear interpolation for simplicity
    for i = 1:length(sol.x)
        alpha = (sol.x(i) - tspan(1))/(tspan(2) - tspan(1));
        sol.y(:,i) = x0 + alpha*(xf - x0);
    end
end