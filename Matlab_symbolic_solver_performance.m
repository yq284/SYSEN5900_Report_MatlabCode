clear all;clc;close all;
rand('seed',1);

% Create array of step sizes
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
        
        % Initial estimates
        b0=0.1; b1=0.1; a0=0.1; a1=0.01; a2=0.01;
        
        % Reference
        T_ref = 25; t_max = 100; time = 0:Ts_array(k):t_max; nt = length(time);
        
        % Slew calculations
        syms C2 t C1 c a
        Tslew = 1.5; 
        Yd = zeros(1,nt);
        y(t) = C2*exp(-t) - exp(-t)*(C1*exp(t) + (t*exp(-c))/2 + (a*t*exp(t))/2);
        ydot(t) = exp(-t)*(C1*exp(t) + (t*exp(-c))/2 + (a*t*exp(t))/2) - exp(-t)*(exp(-c)/2 + C1*exp(t) + (a*exp(t))/2 + (a*t*exp(t))/2) - C2*exp(-t);
        
        % Generate trajectory
        for j=1:nt
            if mod(time(j),2*T_ref)<T_ref
                pn = 1;
            else
                pn = -1;
            end
            
            if mod(time(j),T_ref)<Tslew
                if time(j)>=0 && time(j)<Tslew
                    eqns = [y(0) == 0,ydot(0) == 0 ,y(Tslew) == 1,ydot(0) == 0];
                    S = solve(eqns,[C2 C1 c a]);
                    Yd(j)=S.C2.*exp(-time(j)) - exp(-time(j)).*(S.C1.*exp(time(j)) + (time(j).*exp(-S.c))/2 + (S.a.*time(j).*exp(time(j)))/2);
                elseif time(j)>=50 && time(j)<50 + Tslew
                    eqns = [y(50) == -1,ydot(0) == 0 ,y(50 + Tslew) == 1,ydot(0) == 0];
                    S = solve(eqns,[C2 C1 c a]);
                    Yd(j)=S.C2.*exp(-time(j)) - exp(-time(j)).*(S.C1.*exp(time(j)) + (time(j).*exp(-S.c))/2 + (S.a.*time(j).*exp(time(j)))/2);
                elseif time(j)>=25 && time(j)<25 + Tslew
                    eqns = [y(25) == 1,ydot(0) == 0 ,y(25 + Tslew) == -1,ydot(0) == 0];
                    S = solve(eqns,[C2 C1 c a]);
                    Yd(j)=S.C2.*exp(-time(j)) - exp(-time(j)).*(S.C1.*exp(time(j)) + (time(j).*exp(-S.c))/2 + (S.a.*time(j).*exp(time(j)))/2);
                elseif time(j)>=75 && time(j)<75 + Tslew
                    eqns = [y(75) == 1,ydot(0) == 0 ,y(75 + Tslew) == -1,ydot(0) == 0];
                    S = solve(eqns,[C2 C1 c a]);
                    Yd(j)=S.C2.*exp(-time(j)) - exp(-time(j)).*(S.C1.*exp(time(j)) + (time(j).*exp(-S.c))/2 + (S.a.*time(j).*exp(time(j)))/2);
                elseif time(j)>=100 && time(j)<100 + Tslew
                    eqns = [y(100) == -1,ydot(0) == 0 ,y(100 + Tslew) == 1,ydot(0) == 0];
                    S = solve(eqns,[C2 C1 c a]);
                    Yd(j)=S.C2.*exp(-time(j)) - exp(-time(j)).*(S.C1.*exp(time(j)) + (time(j).*exp(-S.c))/2 + (S.a.*time(j).*exp(time(j)))/2);
                end
            else
                Yd(j)=pn;
            end
        end
        
        % Control parameters initialization
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
            U(t)=[Yd(t+1) Y(t) Y(t-1) U(t-1)]*[1 a1 a2 -b0]'/b1;
        end

        % Complete output array
        Y_true(end+1)=Y_true(end);
        
        % Calculate metrics
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
        results(k,:) = NaN;
    end
end

% Close waitbar
close(h);

% Remove NaN rows
results = results(~isnan(results(:,1)),:);

% Create visualization
figure('Position', [100 100 1000 800])

% Plot 1: Mean Input Value
subplot(3,1,1)
plot(results(:,1), results(:,3), 'b-', 'LineWidth', 2);
grid on;
xlabel('Step Size (seconds)', 'FontSize', 12);
ylabel('Mean Input Value (V)', 'FontSize', 12);
title('Step Size vs Mean Input Value (Pontryagin)', 'FontSize', 14);

% Plot 2: Mean Error
subplot(3,1,2)
plot(results(:,1), results(:,2), 'r-', 'LineWidth', 2);
grid on;
xlabel('Step Size (seconds)', 'FontSize', 12);
ylabel('Mean Error', 'FontSize', 12);
title('Step Size vs Mean Error (Pontryagin)', 'FontSize', 14);

% Plot 3: Error Standard Deviation
subplot(3,1,3)
plot(results(:,1), results(:,4), 'g-', 'LineWidth', 2);
grid on;
xlabel('Step Size (seconds)', 'FontSize', 12);
ylabel('Error Std Dev', 'FontSize', 12);
title('Step Size vs Error Standard Deviation (Pontryagin)', 'FontSize', 14);

% Adjust overall figure appearance
sgtitle('Analysis of Step Size Impact on System Performance (Pontryagin)', 'FontSize', 16);
set(gcf, 'Color', 'white');

% Add tight subplot layout
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% Save results
save('step_size_analysis_pontryagin.mat', 'results');

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