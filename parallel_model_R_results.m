%  parallel_model_R_results.m
%  Analysing battery behaviour with the cascading current apprpach and contact resistances added
%  Created on: 23 Jun 2024
%  Author(s): Nilsu Atlan, Dr. Ross Drummond
 
clc; close all; clear;


% This code run a simulation of n cells in parallel. It is based on the
% paper "esolving Kirchhoff's laws for state-estimator design of Li-ion
% battery packs connected in parallel" where the matrices are defined.

%% Load in the parameters for the batteries
battery_select = load('data.mat'); % load in a database containing battery model parametrs
bat_select_num = 9; % select what battery you want from the table
capacitance_Ah_all = [battery_select.data.capacity]; % Get the capacitances of all the batteries
R1_all = [battery_select.data.R1]; % Get the R1 values of all the batteries
C1_all = [battery_select.data.C1]; % Get the C1 values of all the batteries
r_all = [battery_select.data.R0]; % Get the series resistances of all the batteries

%% Open circuit voltage
ocv_coefs = battery_select.data(bat_select_num).OCV_coeff; % Load in the coefficients of the polynomial used to fit the OCV curves

ocv_points=battery_select.data(bat_select_num).OCV_points;

ocv_points= table2array(ocv_points);
N = 7;
ocv_coefs = polyfit(ocv_points(:,1),ocv_points(:,2),N);

nz = 1e2; % 
z_space= linspace(0,1,nz); %set of points for the soc to sample the ocv curve
for i = 1:nz
    OCV_plot(i) = ocv_eval(ocv_coefs,z_space(i)); % generate curves for the ocv
end

%% Model parameters
n_par = 5; % the number of cells in parallel

capacitance_Ah = capacitance_Ah_all(bat_select_num); % Define: capacitance in Ah.
capacitance = (3600*capacitance_Ah ); % Capacitance in As.

C_rate = 0.01; % Define the charging C-rate
I = n_par*C_rate*capacitance_Ah; % Define the charging current

R_nom = 1*1e-2;
%
r_nom = r_all(bat_select_num);  C_nom = C1_all(bat_select_num); % Set the nominal paraters for the 1st order circuit model of the battery
Q = capacitance*ones(n_par,1);
R = R_nom*ones(n_par,1); C = C_nom*ones(n_par,1); r =r_nom*ones(n_par,1); % Vectors of all the variables in the pack (as in stacked in parallel np times).

sd_vars = 5*1e-3; % standard deviation of the noise for the parameters
var_Q = sd_vars*capacitance; % Add noise to each of the parameters to simulate ageing.
var_R = sd_vars*R_nom*0;
var_F = sd_vars*R_nom;
var_r = sd_vars*r_nom;
var_C = sd_vars*C_nom;

for j = 1:n_par % add noise to each of the parametrs in the pack
    Q(j) = Q(j)+var_Q*randn;
    R(j) = R(j)+var_R*randn;
    r(j) = r(j)+var_r*randn;
    C(j) = C(j)+var_C*randn;
end

R(end) = 0;

R_nom = R1_all(bat_select_num);
F = R_nom*ones(n_par,1);
F_k = F;
% R = 1e-6*R_nom*ones(n_par,1); % makes error very high need to fix m
% R = 0e-5*R_nom*ones(n_par,1); % makes error NaN

tau = 1./(F_k.*C); % Define the time constant % make it for R2, make slight changes to r

% tau = R.*C; % Define the time constant

r(n_par) = r(n_par)+R_nom;
%% Compute functions
% [A11,A12,A21,A22,m] = compute_A11_A12_A21_A22(n_par,R,C,Q,tau); % Define the A11, A12, A21, A22 matrices from the Resolving ... paper.

[new_A11,new_A12,new_A21,new_A22,new_m] = new_compute_A11_A12_A21_A22(n_par,R,C,Q,tau,r); % Define the A11, A12, A21, A22 matrices from the parallel pack project paper.
% new_A11 = A11;
% new_A12 = A12;
% new_A21 = A21;
% new_A22 = A22;
% new_m = m;

%% Define the initial conditions of the model
x0 = zeros(2*n_par,1);
x0 = 0.2*ones(2*n_par,1);
% x0(1) = 0.05; x0(3) = 0.07; x0(5) = 0.1; % The initial state-of-charge of the model

for j = 1:n_par % the initial condition for v_rc
    x0(2*j) = 0;
end

%% Run the simulation
% nT_step = 80*3600;
t_f =  0.6 * 3600/(C_rate);
% tspan = [0,t_f]; % Define the time span of the simulation
% tspan = linspace(0,t_f,nT_step);
tspan = linspace(0, t_f);

global i_branch_store;
global i_branch_store_2;

% [t_ode,x_ode] = ode15s(@(t,x) parrallel_model(t,x,I,A11,A12,A21,A22,ocv_coefs), tspan, x0); % run the model
% [t_ode_no_inv,x_ode_no_inv] = ode15s(@(t,x) parrallel_model_no_mat_inverse(t,x,I,A11,A12,A21,A22,ocv_coefs,m), tspan, x0); % run the model with the analytic expression for the matrix inverse

tic
[t_ode_dae2, x_ode_dae2] = ode15s(@(t, x) parallel_model_dae6(t, x, I, new_A11, new_A12, new_A21, new_A22, ocv_coefs, R, r), tspan, x0);
time_dae2 = toc;

i_branch_values_analytic = i_branch_store;

% tic
% [t_ode_dae2, x_ode_dae2] = ode15s(@(t, x) parallel_model_dae5(t, x, I, new_A11, new_A12, new_A21, new_A22, ocv_coefs, R, r), tspan, x0);
% time_dae2 = toc;

i_branch_values_inverse = i_branch_store_2;

error_curr = norm(i_branch_values_analytic-i_branch_values_analytic)

x_ode = x_ode_dae2;
t_ode = t_ode_dae2;

size_t = size(x_ode,1); % get the number of time steps of the simulation
z_ode = zeros(size_t,n_par); w_ode = zeros(size_t,n_par);
for i = 1:n_par % extact the state-of-charge and relaxation voltage from the simulation
    z_ode(:,i) = x_ode(:,2*i-1);
    w_ode(:,i) = x_ode(:,2*i);
end

nT = max(size(t_ode));i_branch = zeros(n_par,nT);
phi = zeros(n_par,1); phi_init = zeros(n_par-1,1);

for j = 1:nT
    for i = 1:n_par
        OCV_ode(j,i) = ocv_eval(ocv_coefs,z_ode(j,i));
        volts(j,i) = OCV_ode(j,i)+ w_ode(j,i);
    end

    for i = 2:n_par
        phi_init(i-1) = volts(j,i)- volts(j,i-1);
    end
    phi =[I;phi_init];

    i_branch(:,j) = new_A22\phi; % Get the current going into each parallel branch
    current_sum(j) = sum(i_branch(:,j))/I;
end

for j = 1:nT
    for i = 1:n_par
        final_volts(j,i) = volts(j,i) + i_branch(i,j)*r(n_par);
    end
end
%% Plot the results
f_size = 24; 
f_size_gca = 24;
color_1 = 0.9 * [1, 1, 1];

results(z_space, OCV_plot, ocv_points, z_ode, t_ode, n_par, w_ode, i_branch, final_volts, f_size, f_size_gca, color_1);

fig1 = figure;
plot(z_space, OCV_plot,'k','LineWidth',2)
xlabel('State-of-charge ($z(t)$)','interpreter','latex','fontsize',f_size,'interpreter','latex')
ylabel('Open-circuit-voltage (OCV) [V]','interpreter','latex','fontsize',f_size,'interpreter','latex')
grid on
g = gca;
set(g, 'fontsize',f_size_gca);


fig2 = figure; hold on;
for j = 1:n_par
    plot(t_ode, z_ode(:,j),'-k','LineWidth',1.5,'color',color_1*j/n_par);
end
xlabel('Time (s)','interpreter','latex','fontsize',f_size,'interpreter','latex')
ylabel('State-of-charge ($z(t)$)','interpreter','latex','fontsize',f_size,'interpreter','latex')
grid on
leg = legend('Cell 1','Cell 2','location','northwest');
set(leg,'interpreter','latex','fontsize',f_size)
g = gca;
set(g, 'fontsize',f_size_gca); box;

fig3 = figure; hold on;
for j = 1:n_par
    plot(t_ode, w_ode(:,j),'-.k','LineWidth',1.5,'color',color_1*j/n_par);
end
xlabel('Time (s)','interpreter','latex','fontsize',f_size,'interpreter','latex')
ylabel('Relaxation voltage $v_{r}(t)$','interpreter','latex','fontsize',f_size,'interpreter','latex')
grid on
leg = legend('Cell 1','Cell 2','location','northeast');
set(leg,'interpreter','latex','fontsize',f_size)
g = gca;
set(g, 'fontsize',f_size_gca); box;


fig4 = figure; hold on;
for j = 1:n_par
    plot(t_ode/3600, i_branch(j,:),'-k','LineWidth',1.5,'color',color_1*j/n_par)
end
xlabel('Time (hrs)','interpreter','latex','fontsize',f_size,'interpreter','latex')
ylabel('Branch currents $i(t)$','interpreter','latex','fontsize',f_size,'interpreter','latex')
grid on
leg = legend('Cell 1','Cell 2','location','northwest');
set(leg,'interpreter','latex','fontsize',f_size)
g = gca;
set(g, 'fontsize',f_size_gca); box;

fig5 = figure; hold on;
for j = 1:n_par
    plot(t_ode/3600, current_sum,'-.k','LineWidth',1.5,'color',color_1*0/n_par)
end
xlabel('Time (hrs)','interpreter','latex','fontsize',f_size,'interpreter','latex')
ylabel('$\frac{1}{I(t)}\sum_{k = 1}^n i_k(t)$','interpreter','latex','fontsize',f_size,'interpreter','latex')
grid on
% leg = legend('Cell 1','Cell 2','location','northwest');
% set(leg,'interpreter','latex','fontsize',f_size)
g = gca;
set(g, 'fontsize',f_size_gca); box;
% axis([0 max(t_ode/3600) 0 2]);

%%
cell1_color = '#009083';
figure; hold on;
    h1 = plot(t_ode / 3600, final_volts(:,1), '-k', 'LineWidth', 3, 'color', cell1_color);
    for j = 2:n_par-1
        plot(t_ode / 3600, final_volts(:,j), '-k', 'LineWidth', 1.5, 'color', color_1*j/n_par)
    end
    h2 = plot(t_ode / 3600, final_volts(:,n_par), '-k', 'LineWidth', 1.5, 'color', color_1 * n_par / n_par);
    xlabel('Time (hrs)', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex');
    ylabel('Voltage at each branch $V(t)$', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex');
    title('Branch Voltage vs Time', 'interpreter', 'latex', 'fontsize', f_size);
    grid on;
    legend([h1 h2], {'Cell 1', ['Cell ', num2str(n_par)]}, 'Location', 'best');
    g = gca;
    set(g, 'fontsize', f_size_gca); 
    box;

%% Plot the results in subplots

function results(z_space, OCV_plot, ocv_points, z_ode, t_ode, n_par, w_ode, i_branch, final_volts, f_size, f_size_gca, color_1)
    
    cell1_color = '#009083';
    color1 = '#009083';
    color2 = '#1C5F71'; % Hexadecimal color code

    fig6 = figure; 

    % New subplot 1 (formerly subplot 2)
    subplot(2, 2, 1); hold on;
    h1 = plot(t_ode, z_ode(:, 1), '-k', 'LineWidth', 3, 'color', cell1_color); % Handle for first line
    for i = 2:n_par-1
        plot(t_ode, z_ode(:, i), '-k', 'LineWidth', 1.5, 'color', color_1 * i / n_par, 'HandleVisibility', 'off');
    end
    h2 = plot(t_ode, z_ode(:, n_par), '-k', 'LineWidth', 1.5, 'color', color_1 * n_par / n_par); % Handle for last line
    xlabel('Time (s)', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex');
    ylabel('State-of-charge ($z(t)$)', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex');
    title('SOC vs Time', 'interpreter', 'latex', 'fontsize', f_size);
    grid on;
    legend([h1 h2], {'Cell 1', ['Cell ', num2str(n_par)]}, 'Location', 'best');
    g = gca;
    set(g, 'fontsize', f_size_gca); 
    box;

    % New subplot 2 (formerly subplot 3)
    subplot(2, 2, 2); hold on;
    h1 = plot(t_ode, w_ode(:, 1), '-.k', 'LineWidth', 3, 'color', cell1_color); % Handle for first line
    for i = 2:n_par-1
        plot(t_ode, w_ode(:, i), '-.k', 'LineWidth', 1.5, 'color', color_1 * i / n_par, 'HandleVisibility', 'off');
    end
    h2 = plot(t_ode, w_ode(:, n_par), '-.k', 'LineWidth', 1.5, 'color', color_1 * n_par / n_par); % Handle for last line
    xlabel('Time (s)', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex');
    ylabel('Relaxation voltage $v_{r}(t)$', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex');
    title('Relaxation Voltage vs Time', 'interpreter', 'latex', 'fontsize', f_size);
    grid on;
    legend([h1 h2], {'Cell 1', ['Cell ', num2str(n_par)]}, 'Location', 'best');
    g = gca;
    set(g, 'fontsize', f_size_gca); 
    box;

    % Adding the plot of Branch Currents as new subplot 3
    subplot(2, 2, 3); hold on;
    h1 = plot(t_ode / 3600, i_branch(1, :), '-k', 'LineWidth', 3, 'color', cell1_color); % Handle for first line
    for i = 2:n_par-1
        plot(t_ode / 3600, i_branch(i, :), '-k', 'LineWidth', 1.5, 'color', color_1 * i / n_par, 'HandleVisibility', 'off');
    end
    h2 = plot(t_ode / 3600, i_branch(n_par, :), '-k', 'LineWidth', 1.5, 'color', color_1 * n_par / n_par); % Handle for last line
    xlabel('Time (hrs)', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex');
    ylabel('Branch currents $i(t)$', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex');
    title('Branch Currents vs Time', 'interpreter', 'latex', 'fontsize', f_size);
    grid on;
    legend([h1 h2], {'Cell 1', ['Cell ', num2str(n_par)]}, 'Location', 'best');
    g = gca;
    set(g, 'fontsize', f_size_gca); 
    box;

    % Adding the plot of OCV vs SOC as new subplot 4
    subplot(2, 2, 4); hold on;
    plot(z_space, OCV_plot, 'Color', color1, 'LineWidth', 4);
    plot(ocv_points(:, 1), ocv_points(:, 2), 'o', 'Color', color2, 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('State-of-charge ($z(t)$)', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex');
    ylabel('Open-circuit-voltage (OCV) [V]', 'interpreter', 'latex', 'fontsize', f_size, 'interpreter', 'latex');
    legend('OCV Curve', 'OCV Data Points', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', f_size);
    title('OCV vs SOC', 'interpreter', 'latex', 'fontsize', f_size);
    grid on;
    g = gca;
    set(g, 'fontsize', f_size_gca);
    box;

    
end
%% A function to evaluate the vector field of the parallel pack model
function xdot = parrallel_model(t,x,I,A11,A12,A21,A22,ocv_coefs)
n = max(size(x)/2);
Soc = zeros(n,1);  OCV_branch = zeros(n,1);  w = zeros(n); v_mod = zeros(n,1);

for i = 1:n
    Soc(i) = x(2*i-1)'; OCV_branch(i) = ocv_eval(ocv_coefs,Soc(i)); w(i) = x(2*i)'; % extract the soc and v_1
    v_mod(i) =  OCV_branch(i) + 1*w(i);
end

delta_v = zeros(n-1,1);
for i = 2:n
    delta_v(i-1) = v_mod(i)-v_mod(1); % compute the voltage difference with respect to cell one
end

phi = [-delta_v;-I];
i_branch = -A22\phi; % Compute the branch currents.
xdot = A11*x+A12*(-A22\phi); % the vector field
end

%% A function to evaluate the vector field of the parallel pack model
function xdot = parrallel_model_no_mat_inverse(t,x,I,A11,A12,A21,A22,ocv_coefs,m)
n = max(size(x)/2);
Soc = zeros(n,1);  OCV_branch = zeros(n,1);  w = zeros(n); v_mod = zeros(n,1);

for i = 1:n
    Soc(i) = x(2*i-1)'; OCV_branch(i) = ocv_eval(ocv_coefs,Soc(i)); w(i) = x(2*i)'; % extract the soc and v_1
    v_mod(i) =  OCV_branch(i) + 1*w(i);
end

delta_v = zeros(n-1,1);
for i = 2:n
    delta_v(i-1) = v_mod(i)-v_mod(1); % compute the voltage difference with respect to cell one
end

phi = [-delta_v;-I];
i_branch = -A22\phi; % Compute the branch currents.
xdot = A11*x+A12*(-m*phi); % the vector field. THIS IS THE MAIN DIFFERERNCE WITH RESPECT TO THE FUNCTION ABOVE. It uses m = inv(A) with m analytic in terms of the parameters.
end

%% A function which evaluates the polynomial vecotr field of the open circuit voltage
function ocv = ocv_eval(ocv_coefs,z)
ocv = polyval(ocv_coefs,z);
end


%%
function xdot = parallel_model_dae5(t, x, I, A11, A12, A21, A22, ocv_coefs, R, r)
global i_branch_store_2;
n = length(x) / 2;
Soc = zeros(n,1);  OCV_branch = zeros(n,1);  w = zeros(n,1); v_mod = zeros(n,1);

for i = 1:n
    Soc(i) = x(2*i-1)'; OCV_branch(i) = ocv_eval(ocv_coefs,Soc(i)); w(i) = x(2*i)'; % extract the soc and v_1
    v_mod(i) =  OCV_branch(i) + w(i);
end

i_branch = zeros(n, 1); theta = zeros(n, 1); rho = zeros(n, 1);

for i = 2:n
    theta(i) = (R(i)+r(i))/r(i-1);
    rho(i) = 1/r(i-1);
end

cn = 1;
for j = 2:n
    cn = cn+prod(theta(j:end));
end

f= zeros(n, 1);
for j = 1:n-2
    f(j) = rho(j+1)*(v_mod(j+1)-v_mod(j));
    for i = j+1:n-1
        f(j) = f(j)+prod(theta(i:end))*rho(i+1)*(v_mod(j+1)-v_mod(j));
    end
end
f(n-1) = rho(n)*(v_mod(n)-v_mod(n-1));

i_branch(n) = (I-sum(f))/cn;

for j = n-1:-1:1
    i_branch(j) = theta(j+1)* i_branch(j+1)+rho(j+1)*(v_mod(j+1)-v_mod(j));
end

i_branch_store = i_branch;
xdot = A11 * x + A12 * i_branch;
end

%%
function xdot = parallel_model_dae6(t, x, I, A11, A12, A21, A22, ocv_coefs, R, r)
global i_branch_store;
n = length(x) / 2;
Soc = zeros(n,1);  OCV_branch = zeros(n,1);  w = zeros(n,1); v_mod = zeros(n,1);

for i = 1:n
    Soc(i) = x(2*i-1)'; OCV_branch(i) = ocv_eval(ocv_coefs,Soc(i)); w(i) = x(2*i)'; % extract the soc and v_1
    v_mod(i) =  OCV_branch(i) + w(i);
end

i_branch = zeros(n, 1); theta = zeros(n, 1); rho = zeros(n, 1); omega = zeros(n, 1); alpha = zeros(n, 1); phi = zeros(n-1, 1);

for i = 2:n
    theta(i) = (r(i))/r(i-1);
    omega(i) =  (R(i))/r(i-1);
    rho(i) = (v_mod(i)-v_mod(i-1))/r(i-1);
    phi(i-1) =(v_mod(i)-v_mod(i-1));
    alpha(i) = 1+theta(i)+omega(i);
end

beta = zeros(n+2, 1);
beta(n+2) = 0;
beta(n+1) = 1;
beta(n) = alpha(n);
for j = n-1:-1:2
    beta(j) = alpha(j)*beta(j+1)-theta(j)*beta(j+2);
end

beta_hat = beta;
beta_hat(n+2) = 0;
beta_hat(n+1) = 1;
beta_hat(n) = 1;
for j = n-1:-1:2
    beta_hat(j) = alpha(j)*beta_hat(j+1)-theta(j)*beta_hat(j+2);
end

f = zeros(n+2, 1);
f(n) = rho(n); f(n+1) = 0; f(n+2) = 0;
for j = n-1:-1:2
    f(j) = beta_hat(j)*f(j+1)-theta(j)*f(j+2)+rho(j);
end

i_branch(n) = (I-f(2))/beta(2);

for j = n-1:-1:1
    i_branch(j) = theta(j+1)* i_branch(j+1)+omega(j+1)*sum(i_branch)+rho(j+1);
end

i_branch_store = i_branch;
%%
phi_mod  = [I; phi];
i_branch_invert = A22\phi_mod;
%%
check = sum(i_branch)/I;
check_invert = sum(i_branch_invert)/I;
%%

xdot = A11 * x + A12 * i_branch;
end

