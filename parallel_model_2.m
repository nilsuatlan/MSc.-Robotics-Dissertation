%  parallel_model_2.m
%  Analysing battery behaviour using 3 methodologies
%  Created on: 23 Jun 2024
%  Author(s): Nilsu Atlan, Dr. Ross Drummond

 
 clc; close all; clear;


% This code run a simulation of n cells in parallel. It is based on the
% paper "Resolving Kirchhoff's laws for state-estimator design of Li-ion
% battery packs connected in parallel" where the matrices are defined.

%% Load in the parameters for the batteries
battery_select = load('data.mat'); % load in a database containing battery model parametrs
bat_select_num = 8; % select what battery you want from the table
capacitance_Ah_all = [battery_select.data.capacity]; % Get the capacitances of all the batteries
R1_all = [battery_select.data.R1]; % Get the R1 values of all the batteries
C1_all = [battery_select.data.C1]; % Get the C1 values of all the batteries
r_all = [battery_select.data.R0]; % Get the series resistances of all the batteries

%% Open circuit voltage
ocv_coefs = battery_select.data(bat_select_num).OCV_coeff; % Load in the coefficients of the polynomial used to fit the OCV curves

ocv_pojnts=battery_select.data(bat_select_num).OCV_points

ocv_pojnts= table2array(ocv_pojnts)
N =7;
ocv_coefs = polyfit(ocv_pojnts(:,1),ocv_pojnts(:,2),N)

nz = 1e2; % 
z_space= linspace(0,1,nz); %set of points for the soc to sample the ocv curve
for i = 1:nz
    OCV_plot(i) = ocv_eval(ocv_coefs,z_space(i)); % generate curves for the ocv
end
%% Model parameters
n_par = 10; % the number of cells in parallel

capacitance_Ah = capacitance_Ah_all(bat_select_num); % Define: capacitance in Ah.
capacitance = (3600*capacitance_Ah ); % Capacitance in As.

C_rate = 0.01; % Define the charging C-rate
I = n_par*C_rate*capacitance_Ah; % Define the charging current

% R_nom = 1e-5; 
R_nom = 0*1e-5; 
r_nom = r_all(bat_select_num);  C_nom = C1_all(bat_select_num); % Set the nominal paraters for the 1st order circuit model of the battery
Q = capacitance*ones(n_par,1);
R = R_nom*ones(n_par,1); C = C_nom*ones(n_par,1); r =r_nom*ones(n_par,1); % Vectors of all the variables in the pack (as in stacked in parallel np times).


sd_vars = 1*1e-3; % standard deviation of the noise for the parameters
var_Q = sd_vars*capacitance; % Add noise to each of the parameters to simulate ageing.
var_R = sd_vars*R_nom;
var_F = sd_vars*R_nom;
var_r = sd_vars*r_nom;
var_C = sd_vars*C_nom;

for j = 1:n_par % add noise to each of the parametrs in the pack
    Q(j) = Q(j)+var_Q*randn;
    R(j) = R(j)+var_R*randn;
    r(j) = r(j)+var_r*randn;
    C(j) = C(j)+var_C*randn;
end

R_nom = R1_all(bat_select_num); 
F = R_nom*ones(n_par,1);
F_k = F;
% R = 1e-6*R_nom*ones(n_par,1); % makes error very high need to fix m
% R = 0e-5*R_nom*ones(n_par,1); % makes error NaN

tau = 1./(F_k.*C); % Define the time constant % make it for R2, make slight changes to r

% tau = R.*C; % Define the time constant

%% Compute functions
% [A11,A12,A21,A22,m] = compute_A11_A12_A21_A22(n_par,R,C,Q,tau); % Define the A11, A12, A21, A22 matrices from the Resolving ... paper. 

[new_A11,new_A12,new_A21,new_A22,new_m] = new_compute_A11_A12_A21_A22(n_par,R,C,Q,tau,r,F); % Define the A11, A12, A21, A22 matrices from the parallel pack project paper. 
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

% [t_ode,x_ode] = ode15s(@(t,x) parrallel_model(t,x,I,A11,A12,A21,A22,ocv_coefs), tspan, x0); % run the model
% [t_ode_no_inv,x_ode_no_inv] = ode15s(@(t,x) parrallel_model_no_mat_inverse(t,x,I,A11,A12,A21,A22,ocv_coefs,m), tspan, x0); % run the model with the analytic expression for the matrix inverse

tic
[t_ode_dae2, x_ode_dae2] = ode15s(@(t, x) parallel_model_dae5(t, x, I, new_A11, new_A12, new_A21, new_A22, ocv_coefs, R, r), tspan, x0);
time_dae2 = toc;

i_branch_values = i_branch_store;

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
        phi_init(i-1) = volts(j,i)- volts(j,1);
    end
    phi =[I;phi_init];
    % f = R(1)*I;
    % phi = [(R(1)*I - volts(j,1)); phi_init];
    i_branch(:,j) = new_A22\phi; % Get the current going into each parallel branch
    current_sum(j) = sum(i_branch(:,j))/I;
end

%% Plot the results
close all;
f_size = 16; f_size_gca = 13;
color_1 = 0.9*[1,1,1];

fig1 = figure; hold on
plot(z_space, OCV_plot,'k','LineWidth',2)
plot(ocv_pojnts(:,1),ocv_pojnts(:,2),'+','color',0.5*[1 1 1],'LineWidth',2)
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
axis([0 max(t_ode/3600) 0 2]);

%%
% print(fig1,'ocv','-depsc')
% print(fig2,'soc','-depsc')
% print(fig3,'w_relax','-depsc')
% print(fig4,'i_branch','-depsc')
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
%% A function to evaluate the vector field of the parallel pack model with a mass matrix
function xdot = parallel_model_dae1(t, x, I, A11, A12, A21, A22, ocv_coefs, R, r)
    global i_branch_store;
    n = length(x) / 2;
    Soc = zeros(n,1);  OCV_branch = zeros(n,1);  w = zeros(n,1); v_mod = zeros(n,1);

    for i = 1:n
        Soc(i) = x(2*i-1)'; OCV_branch(i) = ocv_eval(ocv_coefs,Soc(i)); w(i) = x(2*i)'; % extract the soc and v_1
        v_mod(i) =  OCV_branch(i) + w(i);
    end

    delta_v = zeros(n-1,1);

    for i = 2:n
        delta_v(i-1) = v_mod(i) - v_mod(i-1);
    end

    phi = [-delta_v;-I];
    i_branch = -A22\phi; % Get the current going into each parallel branch
    i_branch_store = i_branch;

    xdot = A11 * x + A12 * i_branch; 
end


%% A function to evaluate the vector field of the parallel pack model with a mass matrix
function xdot = parallel_model_dae2(t, x, I, A11, A12, A21, A22, ocv_coefs, R, r)
    global i_branch_store;
    n = length(x) / 2;
    Soc = zeros(n,1);  OCV_branch = zeros(n,1);  w = zeros(n,1); v_mod = zeros(n,1);

    for i = 1:n
        Soc(i) = x(2*i-1)'; OCV_branch(i) = ocv_eval(ocv_coefs,Soc(i)); w(i) = x(2*i)'; % extract the soc and v_1
        v_mod(i) =  OCV_branch(i) + w(i);
    end

    delta_v = zeros(n-1,1);

    for i = 2:n
        delta_v(i-1) = v_mod(i) - v_mod(i-1);
    end

    % phi = [-delta_v;-I];
    phi = [(R(1)*I - v_mod(1)); -delta_v];

    i_branch = -A22\phi; % Get the current going into each parallel branch
    i_branch_store = i_branch;

    xdot = A11 * x + A12 * i_branch; 
end

function xdot = parallel_model_dae3(t, x, I, A11, A12, A21, A22, ocv_coefs, R, r)
    global i_branch_store;
    n = length(x) / 2;
    Soc = zeros(n,1);  OCV_branch = zeros(n,1);  w = zeros(n,1); v_mod = zeros(n,1);

    for i = 1:n
        Soc(i) = x(2*i-1)'; OCV_branch(i) = ocv_eval(ocv_coefs,Soc(i)); w(i) = x(2*i)'; % extract the soc and v_1
        v_mod(i) =  OCV_branch(i) + w(i);
    end

    i_branch = zeros(n, 1);

    for i = n:-1:2
        if i == n
            % i_branch(i) = (v_mod(i) - v_mod(i-1)) / (r(i)); % ??
            i_branch(i) = I;
        % elseif i == (n-1)
        %     i_branch(i) = (v_mod(i+1) - v_mod(i) + ((R(i+1) + r(i+1)) * i_branch(i+1))) / (r(i));
        else
            i_branch(i) = (v_mod(i+1) - v_mod(i) + ((R(i+1) + r(i+1)) * i_branch(i+1)) + (R(i+1) * sum(i_branch(i+1:n))) ) / (r(i));
        end
    end

    % Calculate the first branch current to satisfy total current I
    i_branch(1) = (v_mod(2) - v_mod(1) + ((R(2) + r(2)) * i_branch(2)) + (R(2) * sum(i_branch(2:n)))) / r(1);

    i_branch_store = i_branch;

    xdot = A11 * x + A12 * i_branch; 
end

function xdot = parallel_model_dae4(t, x, I, A11, A12, A21, A22, ocv_coefs, R, r)
    global i_branch_store;
    n = length(x) / 2;
    Soc = zeros(n,1);  OCV_branch = zeros(n,1);  w = zeros(n,1); v_mod = zeros(n,1);

    for i = 1:n
        Soc(i) = x(2*i-1)'; OCV_branch(i) = ocv_eval(ocv_coefs,Soc(i)); w(i) = x(2*i)'; % extract the soc and v_1
        v_mod(i) =  OCV_branch(i) + w(i);
    end

    i_branch = zeros(n, 1);

    % Calculate the first branch current to satisfy total current I
    i_branch(1) = (R(1)*I - v_mod(i)) / r(1);
    % i_branch(1) = I*0.5

    for i = 2:n       
        i_branch(i) = (v_mod(i-1) - v_mod(i) + (r(i-1)* i_branch(i-1)) - (R(i) * sum(i_branch(i:n)))) / (r(i));
    end

    i_branch_store = i_branch;

    xdot = A11 * x + A12 * i_branch; 
end

%%
function xdot = parallel_model_dae5(t, x, I, A11, A12, A21, A22, ocv_coefs, R, r)
    global i_branch_store;
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


%% Final version that almost works --> current mismatch issues
% function xdot = parallel_model_dae2(t, x, I, A11, A12, A21, A22, ocv_coefs, R, r)
%     global i_branch_store;
%     n = length(x) / 2;
%     Soc = zeros(n,1);  OCV_branch = zeros(n,1);  w = zeros(n,1); v_mod = zeros(n,1);
% 
%     for i = 1:n
%         Soc(i) = x(2*i-1)'; OCV_branch(i) = ocv_eval(ocv_coefs,Soc(i)); w(i) = x(2*i)'; % extract the soc and v_1
%         v_mod(i) =  OCV_branch(i) + w(i);
%     end
% 
%     i_branch = zeros(n, 1);
% 
%     % for i = n:-1:2
%     %     if i == n
%     %         % i_branch(i) = (v_mod(i) - v_mod(i-1)) / (r(i)); % ??
%     %         i_branch(i) = I;
%     %     % elseif i == (n-1)
%     %     %     i_branch(i) = (v_mod(i+1) - v_mod(i) + ((R(i+1) + r(i+1)) * i_branch(i+1))) / (r(i));
%     %     else
%     %         i_branch(i) = (v_mod(i+1) - v_mod(i) + ((R(i+1) + r(i+1)) * i_branch(i+1)) + (R(i+1) * sum(i_branch(i+1:n))) ) / (r(i));
%     %     end
%     % end
% 
%     % % Calculate the first branch current to satisfy total current I
%     % i_branch(1) = (v_mod(2) - v_mod(1) + ((R(2) + r(2)) * i_branch(2)) + (R(2) * sum(i_branch(2:n)))) / r(1);
% 
%     % Calculate the first branch current to satisfy total current I
%     i_branch(1) = (R(1)*I - v_mod(i)) / r(1);
% 
%     for i = 2:n       
%         i_branch(i) = (v_mod(i-1) - v_mod(i) + (r(i-1)* i_branch(i-1)) - (R(i) * sum(i_branch(i:n)))) / (r(i));
%     end
% 
%     i_branch_store = i_branch;
% 
%     xdot = A11 * x + A12 * i_branch; 
% end

%% A function to evaluate the vector field of the parallel pack model with a mass matrix
% function xdot = parallel_model_dae2(t, x, I, A11, A12, A21, A22, ocv_coefs, tau, C, R)
%     % n = max(size(x)/2);
%     n = length(x) / 2;
%     Soc = zeros(n,1);  OCV_branch = zeros(n,1);  w = zeros(n,1); v_mod = zeros(n,1);
% 
%     for i = 1:n
%         Soc(i) = x(2*i-1)'; OCV_branch(i) = ocv_eval(ocv_coefs,Soc(i)); w(i) = x(2*i)'; % extract the soc and v_1
%         v_mod(i) =  OCV_branch(i) + w(i);
%     end
% 
%     i_branch = zeros(n, 1);
%     delta_v = zeros(n-1,1);
% 
%     for j = 2:n
%         delta_v(j-1) = v_mod(j)-v_mod(1); % compute the voltage difference with respect to cell one
%         sum_current = 0;
%         for m = j:n
%             sum_current = sum_current + i_branch(m);
%             % fprintf('Adding i_branch(%d) = %f to sum_current\n', m, i_branch(m));  % Debugging: print current sum
%         end
%         % fprintf('A22(%d,%d) = %f\n', j, j, A22(j, j));  % Check A22 before division
% 
%          % if A22(j, j) <= 0
%          %    % fprintf('Error: Non-positive value detected in A22(%d,%d) = %f\n', j, j, A22(j, j));
%          %    % i_branch(j) = 0;  % Set to zero or handle appropriately
%          %    i_branch(j) = -sum_current;
%          % else
%             % R1j(j) = 1 / (tau(j) * C(j));
%             % fprintf('R1j = %f, j = %d \n', R1j(j), j); 
%             % i_branch(j) = -sum_current + delta_v(j-1) / R1j(j);
%             i_branch(j) = -sum_current + delta_v(j-1) / R(j);
%          % end
%         % fprintf('i_branch(%d) = %f to \n', j, i_branch(j));  
% 
%     end
%     % l = length(i_branch);
%     % if l > 2
%         i_branch(1) = I - sum(i_branch(2:end));
%     % else if size(i_branch) == 2
%     % else
%     %     i_branch(1) = I - sum(i_branch(2));
%     % end
%     %     delta_v =zeros(n,1);
%     % for j =1:n-1
%     %     delta_v(n+1-j) = (v_mod(n+1-j)-v_mod(n-j))/R(n+1-j);
%     % end
%     % i_branch(n) =delta_v(n);
%     %   for j =1:n-1
%     %     i_branch(n+1-j) = -sum(i_branch) +delta_v(n-j);
%     % end
%     % i_branch(1) =I-sum(i_branch);
%     xdot = A11 * x + A12 * i_branch; 
%     % fprintf('xdot calculated as: %s\n', mat2str(xdot));
% end
% 

% %% A function to evaluate the vector field of the parallel pack model with a mass matrix
% function xdot = parrallel_model_dae(t,x,I,A11,A12,A21,A22,ocv_coefs,m)
% n = max(size(x)/2);
% Soc = zeros(n,1);  OCV_branch = zeros(n,1);  w = zeros(n); v_mod = zeros(n,1);
% 
% for i = 1:n
%     Soc(i) = x(2*i-1)'; OCV_branch(i) = ocv_eval(ocv_coefs,Soc(i)); w(i) = x(2*i)'; % extract the soc and v_1
%     v_mod(i) =  OCV_branch(i) + 1*w(i);
% end
% 
% delta_v = zeros(n-1,1);
% for i = 2:n
%     delta_v(i-1) = v_mod(i)-v_mod(1); % compute the voltage difference with respect to cell one
% end
% 
% phi = [-delta_v;-I];
% i_branch = -A22\phi; % Compute the branch currents.
% xdot = A11*x+A12*(-m*phi); % the vector field. THIS IS THE MAIN DIFFERERNCE WITH RESPECT TO THE FUNCTION ABOVE. It uses m = inv(A) with m analytic in terms of the parameters.
% 
% A = [A11 A12,A21,A22];
% 
% end



