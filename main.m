clear
close all


%% Initializing path variables and parsing input files
% add path to include directory
addpath(genpath([pwd,'/include']))

% parse inputs from input.ini file
inputs = parse_inputs('input.ini');

% read in temperature measurement inputs
T_input = read_input_data(inputs.inputs.T_input_path,inputs.inputs.downsample);
T_0 = T_input.T(1);

%% Preloading spatial and time integration parameters
% load quadrature weights and points into memory
legendre_general = quadrature_init('Legendre_roots.dat');

% Storing thickness of the TPS layer
L = inputs.(inputs.stackup.layers{1}).L;

% Calculating all volume/area/length numerical integration parameters
[L_coeff_vect,L_integrand_mat,L_weights_vect] = gauss_quadrature(L,0,inputs.params.N_x,legendre_general);

% Calculating length of integrand vector - this will be used may times
L_length = length(L_integrand_mat);

t_vect = T_input.norm_t';
t_vect_end_trunk = t_vect(2:end);
t_idx = 1:length(t_vect);

t_length = length(t_idx);

time_coeff_vect = zeros(length(t_idx));
[~,time_integrand_mat,time_weights_vect] = gauss_quadrature(t_vect(t_idx(2:end)),t_vect(t_idx(1:end-1)),inputs.params.N_time,legendre_general);

for i = 2:t_length
    [time_coeff_vect(i,1:i-1),~,~] = gauss_quadrature(t_vect(t_idx(2:i)),t_vect(t_idx(1:i-1)),inputs.params.N_time,legendre_general);
end

% Working vectorization
time_coeff_vect_trunk = time_coeff_vect(2:end,1:end-1);
time_coeff_vect_quick = time_coeff_vect_trunk(:,1);

t_mat_quick = repmat(t_vect(t_idx(2:end))',1,length(time_weights_vect));
time_integrand_mat_quick = repmat(time_integrand_mat(:,1)',t_length-1,1);

% 2nd Derivative operator - Good reference at https://dam.brown.edu/people/alcyew/handouts/numdiff.pdf
t_idx_truncate = t_idx(1:end-1);
L_coeff = 2./((t_vect(t_idx_truncate(3:end))-t_vect(t_idx_truncate(2:end-1)))'.*(t_vect(t_idx_truncate(2:end-1))-t_vect(t_idx_truncate(1:end-2)))'.*(t_vect(t_idx_truncate(3:end))-t_vect(t_idx_truncate(1:end-2)))');
L_op_2nd1 = [diag((t_vect(t_idx_truncate(2:end-1))-t_vect(t_idx_truncate(1:end-2)))'.*ones(t_length-3,1).*L_coeff,0), zeros(t_length-3,2)];
L_op_2nd2 = [zeros(t_length-3,1), -diag((t_vect(t_idx_truncate(3:end))-t_vect(t_idx_truncate(1:end-2)))'.*ones(t_length-3,1).*L_coeff,0), zeros(t_length-3,1)];
L_op_2nd3 = [zeros(t_length-3,2), diag((t_vect(t_idx_truncate(3:end))-t_vect(t_idx_truncate(2:end-1)))'.*ones(t_length-3,1).*L_coeff,0)];

L_op_2nd_forward = (t_vect(t_idx_truncate(2))-t_vect(t_idx_truncate(1)))^-3*[2, -5 ,4, -1, zeros(1,length(t_idx_truncate)-4)];
L_op_2nd_backward = fliplr((t_vect(t_idx_truncate(2))-t_vect(t_idx_truncate(1)))^-3*[2, -5 ,4, -1, zeros(1,length(t_idx_truncate)-4)]);
L_op_2nd_main = L_op_2nd1 + L_op_2nd2 + L_op_2nd3;
L_op_2nd = [L_op_2nd_forward;L_op_2nd_main;L_op_2nd_backward];

L_op_2nd_offset = [L_op_2nd, zeros(size(L_op_2nd,1),1);
    zeros(1,size(L_op_2nd,2)), 1];

%% Loading in material properties data
TPS_matl_str = inputs.(inputs.stackup.layers{1}).material;
matl_dir = [inputs.stackup.database,'/'];
air_dir = './materials/air/';

P_nom = inputs.env.P_nom;

% order of thermal conductivity polynomial fit
k_order = 5;

%cutoff temperature for calculating material curve fits
T_cutoff = 2000;


% Loading in kinetic parameters - these are static
kin_param = load_kinetics([matl_dir,TPS_matl_str,'/',TPS_matl_str,'_kinetics.xlsx']);

% Loading in pyrolysis gas properties for nominal pressure
% pyro_gas_options = sheetnames([matl_dir,matl_str,'_gas.xlsx']);

pyro_gas_h_raw = readmatrix([matl_dir,TPS_matl_str,'/',TPS_matl_str,'_gas.xlsx'],'Sheet','Gas_enthalpy');
pyro_gas_rho_raw = readmatrix([matl_dir,TPS_matl_str,'/',TPS_matl_str,'_gas.xlsx'],'Sheet','Gas_density');
pyro_gas_e_raw = readmatrix([matl_dir,TPS_matl_str,'/',TPS_matl_str,'_gas.xlsx'],'Sheet','Gas_energy');
pyro_gas_Cp_raw = readmatrix([matl_dir,TPS_matl_str,'/',TPS_matl_str,'_gas.xlsx'],'Sheet','Gas_Cp');

if size(pyro_gas_h_raw,2) == 2
    % If gas properties are only given for a single pressure, no need to
    % interpolate
    pyro_gas_h = pyro_gas_h_raw(2:end,2);
    pyro_gas_rho = pyro_gas_rho_raw(2:end,2);
    pyro_gas_e = pyro_gas_e_raw(2:end,2);
    pyro_gas_Cp = pyro_gas_Cp_raw(2:end,2);
else
    % Interpolate gas properties
    pyro_gas_h = interp1(pyro_gas_h_raw(1,2:end)',pyro_gas_h_raw(2:end,2:end)',P_nom,'linear','extrap')';
    pyro_gas_rho = interp1(pyro_gas_rho_raw(1,2:end)',pyro_gas_rho_raw(2:end,2:end)',P_nom,'linear','extrap')';
    pyro_gas_e = interp1(pyro_gas_e_raw(1,2:end)',pyro_gas_e_raw(2:end,2:end)',P_nom,'linear','extrap')';
    pyro_gas_Cp = interp1(pyro_gas_Cp_raw(1,2:end)',pyro_gas_Cp_raw(2:end,2:end)',P_nom,'linear','extrap')';
end

pyro_gas = table(pyro_gas_h_raw(2:end,1),pyro_gas_h,pyro_gas_rho,pyro_gas_e,pyro_gas_Cp);
pyro_gas.Properties.VariableNames = {'T','h','rho','e','Cp'};

% Loading in Air properties for nominal pressure
air_h_raw = readmatrix([air_dir,'air_properties.xlsx'],'Sheet','Air_enthalpy');
air_h = interp1(air_h_raw(1,2:end)',air_h_raw(2:end,2:end)',P_nom,'linear','extrap')';

air_rho_raw = readmatrix([air_dir,'air_properties.xlsx'],'Sheet','Air_density');
air_rho = interp1(air_rho_raw(1,2:end)',air_rho_raw(2:end,2:end)',P_nom,'linear','extrap')';

air_e_raw = readmatrix([air_dir,'air_properties.xlsx'],'Sheet','Air_energy');
air_e = interp1(air_e_raw(1,2:end)',air_e_raw(2:end,2:end)',P_nom,'linear','extrap')';

air_gas = table(air_h_raw(2:end,1),air_h,air_rho,air_e);
air_gas.Properties.VariableNames = {'T','h','rho','e'};

% Loading in virgin and char enthalpy tables -- solid is incompressible so
% these aren't pressure-dependent
h_virgin = readtable([matl_dir,TPS_matl_str,'/',TPS_matl_str,'_solid.xlsx'],'Sheet','Virgin_enthalpy');
h_virgin.Properties.VariableNames = {'T','h'};

h_char = readtable([matl_dir,TPS_matl_str,'/',TPS_matl_str,'_solid.xlsx'],'Sheet','Char_enthalpy');
h_char.Properties.VariableNames = {'T','h'};

% Calculating nominal thermal conductivity, heat capacity, and emissivity
% curve fits from polynomial coefficients.

%thermal conductivity could be pressure-dependent
% Calculating nominal curve fit for virgin conductivity
k_virgin_data_raw = readmatrix([matl_dir,TPS_matl_str,'/',TPS_matl_str,'_solid.xlsx'],'Sheet','Virgin_k');
if size(k_virgin_data_raw,2) == 2
    k_virgin_data = k_virgin_data_raw(2:end,:);
else
    k_virgin_data = [k_virgin_data_raw(2:end,1),interp1(k_virgin_data_raw(1,2:end)',k_virgin_data_raw(2:end,2:end)',P_nom,'linear','extrap')'];
end
k_virgin_data(k_virgin_data(:,1)>T_cutoff,:) = [];
k_virgin_data(isnan(k_virgin_data(:,1)),:) = [];
[k_virgin_nominal,k_virgin_coeffs] = matl_table_to_fun(k_virgin_data,k_order);

k_virgin_lin_data = k_virgin_data(k_virgin_data(:,1)<=1000,:);
[k_virgin_lin_nominal,k_virgin_T_lin_nominal] = matl_table_to_fun(k_virgin_lin_data,1);

% Calculating nominal curve fit for char conductivity
k_char_data_raw = readmatrix([matl_dir,TPS_matl_str,'/',TPS_matl_str,'_solid.xlsx'],'Sheet','Char_k');
if size(k_char_data_raw,2) == 2
    k_char_data = k_char_data_raw(2:end,:);
else
    k_char_data = [k_char_data_raw(2:end,1),interp1(k_char_data_raw(1,2:end)',k_char_data_raw(2:end,2:end)',P_nom,'linear','extrap')'];
end
k_char_data(k_char_data(:,1)>T_cutoff,:) = [];
k_char_data(isnan(k_char_data(:,1)),:) = [];
% [k_char_nominal,k_char_T_nominal] = matl_table_to_fun(k_char_data,1);
[k_char_nominal,k_char_coeffs] = matl_table_to_fun(k_char_data,k_order);

k_char_lin_data = k_char_data(k_char_data(:,1)<=1000,:);
[k_char_lin_nominal,k_char_T_lin_nominal] = matl_table_to_fun(k_char_lin_data,1);

% Calculating virgin and char heat capacity
Cp_virgin_data = readmatrix([matl_dir,TPS_matl_str,'/',TPS_matl_str,'_solid.xlsx'],'Sheet','Virgin_Cp');
Cp_virgin_data(isnan(Cp_virgin_data(:,1)),:) = [];
[Cp_virgin_nominal,~] = matl_table_to_fun(Cp_virgin_data,5);

Cp_char_data = readmatrix([matl_dir,TPS_matl_str,'/',TPS_matl_str,'_solid.xlsx'],'Sheet','Char_Cp');
Cp_char_data(isnan(Cp_char_data(:,1)),:) = [];
[Cp_char_nominal,~] = matl_table_to_fun(Cp_char_data,5);

% Calculating virgin and char emissivity
eps_virgin_data = readmatrix([matl_dir,TPS_matl_str,'/',TPS_matl_str,'_solid.xlsx'],'Sheet','Virgin_eps');
eps_virgin_data(isnan(eps_virgin_data(:,1)),:) = [];
[eps_virgin_nominal,~] = matl_table_to_fun(eps_virgin_data,2);

eps_char_data = readmatrix([matl_dir,TPS_matl_str,'/',TPS_matl_str,'_solid.xlsx'],'Sheet','Char_eps');
eps_char_data(isnan(eps_char_data(:,1)),:) = [];
[eps_char_nominal,~] = matl_table_to_fun(eps_char_data,1);

% Pulling virgin and char density
rho_virgin_nominal = readmatrix([matl_dir,TPS_matl_str,'/',TPS_matl_str,'_solid.xlsx'],'Sheet','Virgin_rho');
rho_char_nominal = readmatrix([matl_dir,TPS_matl_str,'/',TPS_matl_str,'_solid.xlsx'],'Sheet','Char_rho');

% Pullling properties from TPS stackup

% Load properties for the thermal resistance layer. In M2020 this is
% HT424
L_RTV = inputs.(inputs.stackup.layers{2}).L;
RTV_matl_str = inputs.(inputs.stackup.layers{2}).material;
rho_RTV = readmatrix([matl_dir,RTV_matl_str,'/',RTV_matl_str,'_solid.xlsx'],'Sheet','rho');
k_RTV = readmatrix([matl_dir,RTV_matl_str,'/',RTV_matl_str,'_solid.xlsx'],'Sheet','k');
Cp_RTV = readmatrix([matl_dir,RTV_matl_str,'/',RTV_matl_str,'_solid.xlsx'],'Sheet','Cp');
alpha_RTV = k_RTV/(rho_RTV*Cp_RTV);


% Load properties for "structure"- this would typically be the Al
% honeycomb
L_s = inputs.(inputs.stackup.layers{3}).L;
sub_matl_str = inputs.(inputs.stackup.layers{3}).material;
rho_structure = readmatrix([matl_dir,sub_matl_str,'/',sub_matl_str,'_solid.xlsx'],'Sheet','rho');
k_s = readmatrix([matl_dir,sub_matl_str,'/',sub_matl_str,'_solid.xlsx'],'Sheet','k');
Cp_s = readmatrix([matl_dir,sub_matl_str,'/',sub_matl_str,'_solid.xlsx'],'Sheet','Cp');
alpha_s = k_s/(rho_structure*Cp_s);


%% Calculating material properties and multi-parameter Cole-Hopf transformation
% placeholder for future integration of Monte-carlo analysis
scale_param = table(1,1,1,1,1,1,1,1,1,'VariableNames',{'k_v','k_c','Cp_v','Cp_c','rho_v','rho_c','eps_v','eps_c','TC_depth'});
m_iter = 1;

start_time = datetime('now');

% calculating fit functions for material properties
k_virgin = @(T_abs) scale_param.k_v(m_iter)*k_virgin_nominal(T_abs);

k_char = @(T_abs) scale_param.k_c(m_iter)*k_char_nominal(T_abs);

Cp_virgin = @(T_abs)  scale_param.Cp_v(m_iter)*Cp_virgin_nominal(T_abs);

Cp_char = @(T_abs)  scale_param.Cp_c(m_iter)*Cp_char_nominal(T_abs);

rho_virgin = scale_param.rho_v(m_iter)*rho_virgin_nominal;
rho_char = scale_param.rho_c(m_iter)*rho_char_nominal;

eps_virgin = @(T_abs) scale_param.eps_v(m_iter)*eps_virgin_nominal(T_abs);
eps_char = @(T_abs) scale_param.eps_c(m_iter)*eps_char_nominal(T_abs);

%---------- Everything below are calculated properties
char_mass_frac_fun = @(rho_var) (rho_char/(rho_virgin-rho_char))*((rho_virgin./rho_var)-1);
char_volume_frac_fun = @(rho_var) (rho_var-rho_virgin)/(rho_char-rho_virgin);

rho_inter_volume = @(char_frac) rho_virgin + char_frac*(rho_char-rho_virgin);

% Pyrolysis enthalpy lookup table
h_bar = @(T) ((rho_char*interp1(h_char.T,h_char.h,T,'linear','extrap'))-(rho_virgin*interp1(h_virgin.T,h_virgin.h,T,'linear','extrap')))/(rho_char-rho_virgin);

% Calculating reference material properties for Cole-Hopf transformation
% formulation
k_0 = k_virgin(T_0);
k_prime = @(T) k_virgin(T+T_0) - k_0;

% Expanding about virgin k_0 seems to produce the best result
k_char_0 = k_char(T_0);
k_char_prime = @(T) k_char(T+T_0) - k_char_0;

Cp_0 = Cp_virgin(T_0);
Cp_prime = @(T) Cp_virgin(T+T_0) - Cp_0;

Cp_char_0 = Cp_char(T_0);
Cp_char_prime = @(T) Cp_char(T+T_0) - Cp_char_0;

alpha_0 = k_0/(rho_virgin*Cp_0);
alpha_prime = @(T) (k_0 + k_prime(T))./(rho_virgin*(Cp_0+Cp_prime(T))) - alpha_0;

alpha_char_0 = k_char_0/(rho_char*Cp_char_0);
alpha_char_prime = @(T) (k_char_0 + k_char_prime(T))./(rho_char*(Cp_char_0+Cp_char_prime(T))) - alpha_char_0;

% Alpha using char volume fraction
alpha_pyro = @(char_frac,T_abs) ((1-char_frac).*k_virgin(T_abs)*rho_virgin./rho_inter_volume(char_frac) + char_frac.*k_char(T_abs)*rho_char./rho_inter_volume(char_frac))./...
    ((1-char_frac)*rho_virgin.*Cp_virgin(T_abs) + char_frac*rho_char.*Cp_char(T_abs));

alpha_prime_pyro = @(char_frac,T) alpha_pyro(char_frac,T+T_0) - alpha_0; % alpha prime is define relative to a reference alpha_0

% Switched to char mass fraction calculation of surface emissivity
eps_fun = @(T_abs,char_frac) (1-char_frac).*eps_virgin(T_abs).*rho_virgin./rho_inter_volume(char_frac) ...
    + char_frac.*eps_char(T_abs).*rho_char./rho_inter_volume(char_frac);

% Cole-Hopf Transformations
if contains(fieldnames(inputs.(inputs.stackup.layers{1})),'virgin_CH')
    switch fieldnames(inputs.(inputs.stackup.layers{1})).virgin_CH
        case 'true'
            virgin_CH = 1;
        case 'false'
            virgin_CH = 0;
        otherwise
            virgin_CH = 0;
    end
else
    virgin_CH = 0;
end

% Higher order transform for high temps
if k_order == 5
    sample_T_vect = linspace(0,1000,k_order+1)';

    % Keep in mind when using polyval that this is flipped vs. convention
    k_prime_sample_fit = [sample_T_vect,sample_T_vect.^2,sample_T_vect.^3,sample_T_vect.^4,sample_T_vect.^5];

    k_prime_coeffs = (k_prime_sample_fit'*k_prime_sample_fit)^-1*k_prime_sample_fit'*k_prime(sample_T_vect);
    
    % A is an abritrary constant - standardize definition to stay
    % consistent with other areas
    k_virgin_T = scale_param.k_v(m_iter)*k_virgin_T_lin_nominal(2);
    k_char_T = scale_param.k_c(m_iter)*k_char_T_lin_nominal(2);
    A = k_0^2/(2*k_virgin_T);

    % High order virgin transformation
    theta_T_transform_TPS = @(T) (k_0*T ...
        + k_prime_coeffs(1).*T.^2/2 ...
        + k_prime_coeffs(2).*T.^3/3 ...
        + k_prime_coeffs(3).*T.^4/4 ...
        + k_prime_coeffs(4).*T.^5/5 ...
        + k_prime_coeffs(5).*T.^6/6 ...
        )/A;

    % Linear direct transformation from theta to T
    T_theta_transform_TPS = @(theta) (k_0/k_virgin_T)*(sqrt(theta+1)-1);
else
    k_virgin_T = scale_param.k_v(m_iter)*k_virgin_T_lin_nominal(2);
    k_char_T = scale_param.k_c(m_iter)*k_char_T_lin_nominal(2);
    A = k_0^2/(2*k_virgin_T);
    theta_T_transform_TPS = @(T) T.*(k_0 + k_virgin_T.*T/2)/A;
    T_theta_transform_TPS = @(theta) (k_0/k_virgin_T)*(sqrt(theta+1)-1);
end

% Calculating drho/dT and dchi/dT -- note, here chi is defined as char MASS
% fraction
[drhodT,rho_local] = drhodT_fun(kin_param,t_vect(2:end),T_input.T(1) + (T_input.T(2:end)-T_input.T(1)));
dchidT_mass = -(rho_char/(rho_virgin-rho_char))*(rho_virgin./(rho_local.^2)).*drhodT;

% function for beta determination
beta_scale_fun = @(T,dchidT,char_frac_m) (char_frac_m.*(k_0 + k_prime(T))+ A*theta_T_transform_TPS(T).*dchidT)...
    ./(char_frac_m.*(k_char_0+k_char_prime(T)-(k_0+k_prime(T))));

% Calculating beta for the TC location
char_mass_frac_vect = char_mass_frac_fun(rho_local);
beta_indicator_vect = (beta_scale_fun((T_input.T(2:end)-T_input.T(1)),dchidT_mass',char_mass_frac_vect'));

% removing influence from beginning, where the surface is pure virgin. Beta
% only influences the char state.
start_char = find(abs(char_mass_frac_vect-0.02) == min(abs(char_mass_frac_vect-0.02)));
beta_cutoff = start_char;

% solving for beta that best approximates k_bar/k = 1 for the relevant time
% points
X_beta = beta_indicator_vect(beta_cutoff:end);

% Beta scaling parameter for TC probe
beta_scale_TC = ((X_beta'*X_beta)^-1)*X_beta'*ones(size(beta_indicator_vect(beta_cutoff:end)));

% Calculating multi-parameter C-H transform
theta_T_transform_TPS_general = @(char_frac,T,beta_scale) (1+beta_scale.*char_frac.*(rho_char./rho_inter_volume(char_frac))).*theta_T_transform_TPS(T);

% Gradient transform -- if temperature and char fraction are known
grad_T_theta_transform_TPS_shortcut = @(T_abs,char_frac,grad_theta) A./((1-char_frac).*(rho_virgin./rho_inter_volume(char_frac)).*k_virgin(T_abs) ...
    + char_frac.*(rho_char./rho_inter_volume(char_frac)).*(k_char(T_abs))).*grad_theta;

%% calculating Green's functions for stackup

% Solving eigenvalue problem for 3 part stackup
addpath('./Greens_functions/basis_functions/ins_stackup_full')
[V,gamma] = build_GF_stackup_full(inputs.params.N,...
    {k_0,k_virgin_T,alpha_0,L},...
    {k_s,alpha_s,L_s},...
    {k_RTV,alpha_RTV,L_RTV},1);

% Second, point to functions w/ eigenvectors/eigenvalues as inputs
addpath('./Greens_functions/build_stackup_full/')
GF_args = {'N',inputs.params.N,'V',V,'gamma',gamma,...
    'k_0',k_0,'k_T',k_virgin_T,'k_s',k_s,'k_RTV',k_RTV,...
    'A',A,'R_th',0,...
    'L_s',L_s,'L',L,'L_RTV',L_RTV};

%% Calculate Green's functions at integration points

if contains(fieldnames(inputs.params),'discretization_method')
    switch inputs.params.discretization_method
        case 'piecewise'
            q_interp_bool = 0;
        case 'interpolate'
            q_interp_bool = 1;
        otherwise
            q_interp_bool = 0;
    end
else
    q_interp_bool = 0;
end

G_q_virtual = cell(L_length,1);
G_q_del_virtual = cell(L_length,1);
G_q_del2_virtual = cell(L_length,1);

G_g_loc = cell(L_length,L_length); % rows for x, columns for x0
G_del_loc = cell(L_length,L_length);
G_del2_loc = cell(L_length,L_length);

zero_pad = zeros(1,t_length-2);

if q_interp_bool
    for j = 1:L_length

        % No interpolation for base Green's function
        G_q_virtual_kernel = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L,GF_args),2));
        G_q_virtual{j} = toeplitz(G_q_virtual_kernel,[G_q_virtual_kernel(1),zero_pad]);
        
        % G_q_virtual_kernel_back = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L,GF_args).*(t_vect(2)-time_integrand_mat_quick)/(t_vect(2)-t_vect(1)),2));  
        % G_q_virtual_kernel_central = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L,GF_args).*(1 + (time_integrand_mat_quick-t_vect(2))/(t_vect(2)-t_vect(1))),2));
        % G_q_virtual{j} = toeplitz([0;G_q_virtual_kernel_back(1:end-1)],[0,zero_pad]) ...
        %     + toeplitz(G_q_virtual_kernel_central,[G_q_virtual_kernel_central(1),zero_pad]);
        
        % No interpolation for Green's function 1st derivative
        G_q_del_virtual_kernel = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_del_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L,GF_args),2));
        G_q_del_virtual{j} = toeplitz(G_q_del_virtual_kernel,[G_q_del_virtual_kernel(1),zero_pad]);

        % G_q_del_virtual_kernel_back = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_del_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L,GF_args).*(t_vect(2)-time_integrand_mat_quick)/(t_vect(2)-t_vect(1)),2));  
        % G_q_del_virtual_kernel_central = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_del_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L,GF_args).*(1 + (time_integrand_mat_quick-t_vect(2))/(t_vect(2)-t_vect(1))),2));
        % G_q_del_virtual{j} = toeplitz([0;G_q_del_virtual_kernel_back(1:end-1)],[0,zero_pad]) ...
        %     + toeplitz(G_q_del_virtual_kernel_central,[G_q_del_virtual_kernel_central(1),zero_pad]);
    
        % 1st order interpolation for Green's function 2nd derivative
        % G_q_del2_virtual_kernel = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_del2_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L,GF_args),2));
        % G_q_del2_virtual{j} = toeplitz(G_q_del2_virtual_kernel,[G_q_del2_virtual_kernel(1),zero_pad]);
        
        G_q_del2_virtual_kernel_back = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_del2_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L,GF_args).*(t_vect(2)-time_integrand_mat_quick)/(t_vect(2)-t_vect(1)),2));  
        G_q_del2_virtual_kernel_central = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_del2_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L,GF_args).*(1 + (time_integrand_mat_quick-t_vect(2))/(t_vect(2)-t_vect(1))),2));
        G_q_del2_virtual{j} = toeplitz([0;G_q_del2_virtual_kernel_back(1:end-1)],[0,zero_pad]) ...
            + toeplitz(G_q_del2_virtual_kernel_central,[G_q_del2_virtual_kernel_central(1),zero_pad]);
    
        % Construct Green's functions for integration abscissas. No interpolation for base and 1st derivative Green's function. 1st order interpolation for 2nd derivative Greens' functions 
        for k = 1:L_length % Cycle through all r_0,z_0 in V1 ****BE CAREFUL NOT TO MIX UP J AND K INDICIES. THEY SWITCH IN THE ITERATIVE LOOP
    
            G_g_loc_kernel = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L_integrand_mat(k),GF_args),2));
            G_g_loc{j,k} = toeplitz(G_g_loc_kernel,[G_g_loc_kernel(1),zero_pad]);
            
            G_del_loc_kernel = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_del_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L_integrand_mat(k),GF_args),2));
            G_del_loc{j,k} = toeplitz(G_del_loc_kernel,[G_del_loc_kernel(1),zero_pad]);

            % G_del2_loc_kernel = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_del2_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L_integrand_mat(k),GF_args),2));
            % G_del2_loc{j,k} = toeplitz(G_del2_loc_kernel,[G_del2_loc_kernel(1),zero_pad]);

            G_del2_loc_kernel_back = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_del2_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L_integrand_mat(k),GF_args).*(t_vect(2)-time_integrand_mat_quick)/(t_vect(2)-t_vect(1)),2));  
            G_del2_loc_kernel_central = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_del2_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L_integrand_mat(k),GF_args).*(1 + (time_integrand_mat_quick-t_vect(2))/(t_vect(2)-t_vect(1))),2));
            G_del2_loc{j,k} = toeplitz([0;G_del2_loc_kernel_back(1:end-1)],[0,zero_pad]) ...
                + toeplitz(G_del2_loc_kernel_central,[G_del2_loc_kernel_central(1),zero_pad]);
        end
    end
else
    % If no linear interpolation selected
    % Calculating each toeplitz matrix kernel
    for j = 1:L_length

        G_q_virtual_kernel = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L,GF_args),2));        
        G_q_del_virtual_kernel = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_del_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L,GF_args),2));
        G_q_del2_virtual_kernel = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_del2_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L,GF_args),2));
    
        G_q_virtual{j} = toeplitz(G_q_virtual_kernel,[G_q_virtual_kernel(1),zero_pad]);
        G_q_del_virtual{j} = toeplitz(G_q_del_virtual_kernel,[G_q_del_virtual_kernel(1),zero_pad]);
        G_q_del2_virtual{j} = toeplitz(G_q_del2_virtual_kernel,[G_q_del2_virtual_kernel(1),zero_pad]);
    
        for k = 1:L_length % Cycle through all r_0,z_0 in V1 ****BE CAREFUL NOT TO MIX UP J AND K INDICIES. THEY SWITCH IN THE ITERATIVE LOOP    
            G_g_loc_kernel = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L_integrand_mat(k),GF_args),2));
            G_del_loc_kernel = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_del_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L_integrand_mat(k),GF_args),2));
            G_del2_loc_kernel = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_del2_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L_integrand_mat(j),L_integrand_mat(k),GF_args),2));
    
            G_g_loc{j,k} = toeplitz(G_g_loc_kernel,[G_g_loc_kernel(1),zero_pad]);
            G_del_loc{j,k} = toeplitz(G_del_loc_kernel,[G_del_loc_kernel(1),zero_pad]);
            G_del2_loc{j,k} = toeplitz(G_del2_loc_kernel,[G_del2_loc_kernel(1),zero_pad]);
           
        end    
    end
end

%% calculating TC location
if m_iter > 1
    error('TC_depth not scaled. pls fix')
end
TC_depth = L - inputs.inputs.TC_depth;

%% Calculating Green's functions at TC location
G_g_loc_TC = cell(1,L_length);

if q_interp_bool
    G_discrete_kernel_back = time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,TC_depth,L,GF_args).*(t_vect(2)-time_integrand_mat_quick)/(t_vect(2)-t_vect(1)),2);
    G_discrete_kernel_central = time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,TC_depth,L,GF_args).*(1 + (time_integrand_mat_quick-t_vect(2))/(t_vect(2)-t_vect(1))),2);
    G_discrete = toeplitz([0;G_discrete_kernel_back(1:end-1)],[0,zero_pad]) ...
        + toeplitz(G_discrete_kernel_central,[G_discrete_kernel_central(1),zero_pad]);
    
    for j = 1:L_length
        G_g_loc_TC_kernel_back = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,TC_depth,L_integrand_mat(j),GF_args).*(t_vect(2)-time_integrand_mat_quick)/(t_vect(2)-t_vect(1)),2));
        G_g_loc_TC_kernel_central = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,TC_depth,L_integrand_mat(j),GF_args).*(1 + (time_integrand_mat_quick-t_vect(2))/(t_vect(2)-t_vect(1))),2));
        G_g_loc_TC{j} = toeplitz([0;G_g_loc_TC_kernel_back(1:end-1)],[0,zero_pad]) ...
            + toeplitz(G_g_loc_TC_kernel_central,[G_g_loc_TC_kernel_central(1),zero_pad]);
    end
else
    G_discrete_kernel = time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,TC_depth,L,GF_args),2);
    G_discrete = toeplitz(G_discrete_kernel,[G_discrete_kernel(1),zero_pad]);
       
    for j = 1:L_length
        G_g_loc_TC_kernel = single(time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,TC_depth,L_integrand_mat(j),GF_args),2));
        G_g_loc_TC{j} = toeplitz(G_g_loc_TC_kernel,[G_g_loc_TC_kernel(1),zero_pad]);    
    end
end

% Transforming temperature input measurement using multi-parameter Cole Hopf transformation
char_frac_TC = char_volume_frac_fun(surf_pyro_r3(kin_param,t_vect(2:end),T_input.T(2:end),23))';
theta_TC = theta_T_transform_TPS_general(char_frac_TC,T_input.T(2:end)-T_input.T(1),beta_scale_TC);

L_op_2nd_double = L_op_2nd;

%% Iterative Solution Loop
% Initialize iteration loop
iter = 1;
iter_end = 20;

% Likely not needed- legacy
% addpath(function_path)

% Preallocate solution matrix
q_reconstruct_mat = zeros(t_length-1,iter_end);
theta_eff_mat = zeros(t_length-1,iter_end);

% Iterative solver parameters
converge_bool = 0;
last_step_bool = 0;
conv_crit = 1e8;
conv_crit_mat = [];

% Initialize virtual energy generation parameters
g_int = zeros(t_length-1,L_length);
g_del_int = g_int;
g_del2_int = g_int;
g_integrand = g_int;
g_int_TC = zeros(t_length-1,1);
char_frac_abis = g_int;

if contains(fieldnames(inputs.(inputs.stackup.layers{1})),'pyrolysis')
    switch fieldnames(inputs.(inputs.stackup.layers{1})).pyrolysis
        case 'on'
            pyro_bool = 1;
        case 'off'
            pyro_bool = 0;
        otherwise
            pyro_bool = 0;
    end
else
    pyro_bool = 1;
end

if contains(fieldnames(inputs.params),'NL_step_weight')
    NL_step_weight = inputs.params.NL_step_weight;
else
    NL_step_weight = 1;
end

% Nonlinear solution loop
while all([iter <= iter_end, converge_bool == 0])
   
  
    if iter > 1 %calculating non-linear correction g_prime
       
        % calculating temperature at absicias for temp-dependent material properties
        theta_abis = zeros(t_length-1,L_length);
        for j = 1:L_length
            theta_abis(:,j) = -(alpha_0./A).*(G_q_virtual{j}*q_reconstruct_finite)...
                + (alpha_0/A)*g_int(:,j);
        end

        % Calculating alpha_prime, either as a function of purely
        % temperature or as a function of both temperature and char
        % fraction
        if pyro_bool
           
            % on second iteration, initialize parameters for calculating
            % beta_sacle function
            if iter == 2
                T_prime_abis = T_theta_transform_TPS(theta_abis);
                drhodT_abis = zeros(size(theta_abis));
                rho_local_abis = zeros(size(theta_abis));
                start_abis_idx = 1;
            else
                start_abis_idx_raw = find(char_mass_frac_vect_abis(end,:) < 0.02);
                if length(start_abis_idx_raw) >= 3
                    start_abis_idx = start_abis_idx_raw(end-2);
                end
            end

            % Calculating beta_scale factor at integration abcissas
            
            for iii = start_abis_idx:L_length
                [drhodT_abis(:,iii),rho_local_abis(:,iii)] = drhodT_fun(kin_param,t_vect(2:end),T_0+T_prime_abis(:,iii));
            end

            dchidT_mass_abis = -(rho_char/(rho_virgin-rho_char))*(rho_virgin./(rho_local_abis.^2)).*drhodT_abis;
            
            % Using least squares to estimate beta instead, which better approximates
            % for the entire heating profile
            char_mass_frac_vect_abis = char_mass_frac_fun(rho_local_abis);
            beta_indicator_vect_abis = (beta_scale_fun(T_prime_abis,dchidT_mass_abis,char_mass_frac_vect_abis));
            
            beta_scale_abis = zeros(1,L_length);
            for iii = 1:L_length
                
                if char_mass_frac_vect_abis(end,iii) < 0.02
                    beta_scale_abis(:,iii) = 0;
                    continue
                end

                beta_cutoff_abis_local = find(abs(char_mass_frac_vect_abis(:,iii)-0.02) == min(abs(char_mass_frac_vect_abis(:,iii)-0.02)));
              
                % solving for beta that best approximates k_bar/k = 1 for the relevant time
                % points
                X_beta_abis_local = beta_indicator_vect_abis(beta_cutoff_abis_local:end,iii);
                
                if any(isinf(X_beta_abis_local))
                    beta_scale_abis(:,iii) = 0;
                else
                    beta_scale_abis(:,iii) = ((X_beta_abis_local'*X_beta_abis_local)^-1)*X_beta_abis_local'*ones(size(beta_indicator_vect_abis(beta_cutoff_abis_local:end,iii)));
                end
               
            end

            if all([~isscalar(beta_scale_abis(:,1)), size(beta_scale_abis,1) < 4],2)
                beta_scale_abis = [beta_scale_abis;zeros(4-size(beta_scale_abis,1),L_length)];
            end

            % If using multi-parameter Cole-Hopf transform, reverse transformation from theta-space to T-space becomes iterative.            
            theta_conv = 1;
            counter = 0;
            while all([theta_conv > 0.005, counter < 3])
                counter = counter + 1;

                % Calculate temperature rise at each integration point
                if virgin_CH
                    T_prime_abis = T_theta_transform_TPS(theta_abis);
                else
                    T_prime_abis = T_theta_transform_TPS_general(theta_abis,char_frac_abis,beta_scale_abis,theta_T_transform_TPS_general);
                end
            
                % Calculate density at each integration point
                rho_abis = zeros(size(theta_abis));
                for iii = 1:L_length
                    rho_abis(:,iii) = surf_pyro_r3(kin_param,t_vect(2:end),T_0 + T_prime_abis(:,iii),23)';
                end
                
                % Calculating char volume fraction at each integration point
                char_frac_abis = char_volume_frac_fun(rho_abis);
            
                % Depending on the transform used, either break the while-loop
                % or calculate convergence criteria
                if virgin_CH
                    theta_conv = 0;
                else
                    theta_conv = max(max(abs(theta_abis - theta_T_transform_TPS_general(char_frac_abis,T_prime_abis,beta_scale_abis))./theta_abis));
                end

                % keyboard
            
                if counter == 3
                    disp('Cole-Hopf Transform Might Not Have Converged --> check')
                end
            end

            % Calculating material properties at each integration point
            alpha_abis = alpha_prime_pyro(char_frac_abis,T_prime_abis);
            alpha_0_abis = alpha_0;

            % Calculating drho/dt for each integration point using 2-pt
            % central difference operator
            conv_kernel = [1 0 -1]./(t_vect(3)-t_vect(1));
            drhodt_abis = conv2(rho_abis',conv_kernel,'same')';
            drhodt_abis([1,end],:) = [drhodt_abis(2,:);drhodt_abis(end-1,:)];

            % Calculating energy storage term for both pyrolysis gas and
            % air (for the initially unchared state)
            pyro_storage_abis = kin_param.phi(1)*interp1(pyro_gas.T,pyro_gas.rho,T_0 + T_prime_abis,'linear','extrap').*interp1(pyro_gas.T,pyro_gas.e,T_0 + T_prime_abis,'linear','extrap');
            air_storage_abis = kin_param.phi(1)*interp1(air_gas.T,air_gas.rho,T_0 + T_prime_abis,'linear','extrap').*interp1(air_gas.T,air_gas.e,T_0 + T_prime_abis,'linear','extrap');
            
            % Calculating drhodt for pyro gas and air storage terms using
            % 2-pt central difference
            ddt_pyro_storage = conv2(pyro_storage_abis',conv_kernel,'same')';
            ddt_pyro_storage([1,end],:) = [ddt_pyro_storage(2,:);ddt_pyro_storage(end-1,:)];

            ddt_air_storage = conv2(air_storage_abis',conv_kernel,'same')';
            ddt_air_storage([1,end],:) = [ddt_air_storage(2,:);ddt_air_storage(end-1,:)];

            % Threshold to switch from pyrolysis to air
            pyro_threshold = 0.02;
            
            % Calculating delTheta at the abiscias for calculation of the
            % enthalpy gradient
            deltheta_abis = zeros(t_length-1,L_length);
            int_drhodt = zeros(t_length-1,L_length);
            for j = 1:L_length
                deltheta_abis(:,j) = -(alpha_0./A).*(G_q_del_virtual{j}*q_reconstruct_finite)...
                    + (alpha_0/A)*g_del_int(:,j);

                %Piggybacking on for loop to also calculate mass flux
                local_drhodt = drhodt_abis(:,1:j);
                int_drhodt(:,j) = trapz(L_integrand_mat(1:j)',local_drhodt,2);
                
            end

            % Calculating temperature gradient using shortcut
            grad_T_abis = grad_T_theta_transform_TPS_shortcut(T_0 + T_prime_abis,char_frac_abis,deltheta_abis);

            % Calculating enthalpy gradient at integration abscisas
            grad_h = interp1(pyro_gas.T,pyro_gas.Cp,T_0 + T_prime_abis,'linear','extrap').*grad_T_abis;

            % Calculating real volumetric energy generation term at each
            % integration point.
            % Production version- stable
            g_pyro = (interp1(pyro_gas.T,pyro_gas.h,T_0 + T_prime_abis,'linear','extrap')-h_bar(T_0 + T_prime_abis)).*drhodt_abis ... % Net enthalpy change due to pyrolysis and char formation
                - ((abs(char_frac_abis) >= pyro_threshold).*ddt_pyro_storage + (abs(char_frac_abis) < pyro_threshold).*ddt_air_storage) ... % Gas energy storage. Comment this line out of it is trying to solve for NaNs in the interpolated gas tables
                + grad_h.*int_drhodt ... % Pyrolysis gas convection term
                ;
         
        else
            % If no pyrolysis, just calculate change in material properties
            alpha_abis = alpha_prime_pyro(0,T_theta_transform_TPS(theta_abis));
            alpha_0_abis = alpha_0;
            char_frac_abis = 0;
            g_pyro = 0;
        end
        
        % Calculating del2 at the abiscias
        del2theta_abis = zeros(t_length-1,L_length);
        for j = 1:L_length
            del2theta_abis(:,j) = -(alpha_0./A).*(G_q_del2_virtual{j}*q_reconstruct_finite)...
                + (alpha_0/A)*g_del2_int(:,j);
        end
        
        % storing g_integrand from previous iteration step
        g_integrand_old = g_integrand;

        % Calculating new energy generation integrand
        g_integrand_new = A*(alpha_abis./alpha_0_abis).*del2theta_abis ...
            + ((alpha_0_abis + alpha_abis)./alpha_0_abis).*g_pyro; 
   
        % calculate g_integrand for cuirrent iteration step
        g_integrand = NL_step_weight*g_integrand_new + (1-NL_step_weight)*g_integrand_old;
                      
        % Calculating virtual energy generation terms for all integration
        % abscissas, as well as at the TC location
        g_spatial_integrand = cell(1,L_length);
        g_del_spatial_integrand = cell(1,L_length);
        g_del2_spatial_integrand = cell(1,L_length);

        g_spatial_integrand_TC = zeros(t_length-1,L_length);      

        for j = 1:L_length % Cycle through all r,z in v1
            g_spatial_integrand{j} = zeros(t_length-1,L_length);
            g_del_spatial_integrand{j} = zeros(t_length-1,L_length);
            g_del2_spatial_integrand{j} = zeros(t_length-1,L_length);

            g_spatial_integrand_TC(:,j) = G_g_loc_TC{j}*g_integrand(:,j);
            
            for k = 1:L_length % Cycle through all r_0,z_0 in V1
                g_spatial_integrand{j}(:,k) = G_g_loc{j,k}*g_integrand(:,k);
                g_del_spatial_integrand{j}(:,k) = G_del_loc{j,k}*g_integrand(:,k);
                g_del2_spatial_integrand{j}(:,k) = G_del2_loc{j,k}*g_integrand(:,k);
            end

            g_int(:,j) = L_coeff_vect*sum(L_weights_vect.*g_spatial_integrand{j},2);
            g_del_int(:,j) = L_coeff_vect*sum(L_weights_vect.*g_del_spatial_integrand{j},2);
            g_del2_int(:,j) = L_coeff_vect*sum(L_weights_vect.*g_del2_spatial_integrand{j},2);

        end % end for loop for L integration
        
        g_int_TC = L_coeff_vect*sum(L_weights_vect.*g_spatial_integrand_TC,2);
              
    end
 

    % Calculating effective transformed temperature for inversion
    theta_eff = theta_TC - (alpha_0/A).*g_int_TC;

    
    % Inverting the linear system to recover the hot-wall heat flux. These
    % are 4 different methods that were evaluated for speed

    % scaling_coeffs_finite = (G_discrete_mod'*G_discrete_mod + (lambda_reg_2nd_double.*L_op_2nd_double)'*(lambda_reg_2nd_double.*L_op_2nd_double))^(-1)*G_discrete_mod'*(-A./(alpha_0).*theta_eff);
    
    % General mldivide -> 2x speed
    % scaling_coeffs_finite = (G_discrete_mod'*G_discrete_mod + (lambda_reg_2nd_double.*L_op_2nd_double)'*(lambda_reg_2nd_double.*L_op_2nd_double))\(G_discrete_mod'*(-A./(alpha_0).*theta_eff));
    
    % direct LU factorization
    % [L_fact,U_fact,P_fact] = lu(G_discrete_mod'*G_discrete_mod + (lambda_reg_2nd_double.*L_op_2nd_double)'*(lambda_reg_2nd_double.*L_op_2nd_double));
    % LU_inter = L_fact\(P_fact*(G_discrete_mod'*(-A./(alpha_0).*theta_eff)));
    % scaling_coeffs_finite = U_fact\LU_inter;

    % QR factorization
    % [Q_fact,R_fact] = qr(G_discrete_mod'*G_discrete_mod + (lambda_reg_2nd_double.*L_op_2nd_double)'*(lambda_reg_2nd_double.*L_op_2nd_double));
    % scaling_coeffs_finite = R_fact^-1*Q_fact'*G_discrete_mod'*(-A./(alpha_0).*theta_eff);
    
    % toc
  
    % Recovering the Hot-wall heat flux
    switch inputs.params.regularization_method
        case 'Surrogate'
            scaling_coeffs_finite = (G_discrete_mod'*G_discrete_mod + (lambda_reg_2nd_double.*L_op_2nd_double)'*(lambda_reg_2nd_double.*L_op_2nd_double))\(G_discrete_mod'*(-A./(alpha_0).*theta_eff));
            q_reconstruct_finite = HW_scale*scaling_coeffs_finite;
        case 'Tikhonov'
            q_reconstruct_finite = (G_discrete'*G_discrete + inputs.params.lambda_reg_Tik^2*eye(size(G_discrete,2)))\(G_discrete'*(-A./(alpha_0).*theta_eff));
    end
    
    q_reconstruct_mat(:,iter) = q_reconstruct_finite;

    theta_eff_mat(:,iter) = theta_eff;

    % Logic to determine if the nonlinear solution has converged
    conv_crit_prev = conv_crit;
    if iter > 1
        conv_crit = abs(norm(q_reconstruct_finite-q_reconstruct_mat(:,iter-1))/mean(abs(q_reconstruct_finite)));
        conv_crit_mat = [conv_crit_mat,conv_crit];
        % keyboard
    end

    if conv_crit <= 0.02 % 2% change from iteration to iteration
        converge_bool = 1;
        disp(['Solution converged in ', num2str(iter),' iterations'])
    % Edited to prevent false divergence flags. >1 doesn't have much significance 
    elseif conv_crit-conv_crit_prev > 1
        if del2_reduce > 5
            converge_bool = 1;
            disp(['Solution began to diverge in ', num2str(iter),' iterations'])
        else
            del2_reduce = del2_reduce + 0.5;
            disp(['Solution began to diverge in ', num2str(iter),' iterations. Increase Laplacian Scale Factor'])
            q_reconstruct_finite = q_reconstruct_mat(:,iter-1);
        end
    else
        iter = iter + 1;
    end


end


% Calculating equivalent objective function at TC location for comparison
% Objective function = sum((TC_calc - TC_meas)^2)
theta_TC_residual = -(alpha_0/A).*G_discrete*q_reconstruct_finite + (alpha_0/A).*g_int_TC;

T_residual = T_theta_transform_TPS(theta_TC_residual) + T_0;

conv = 1;
iter = 0;
while all([conv > 0.02, iter <=5],2)
    iter = iter + 1;

    % Calculating estimated surface char fraction
    char_frac_residual = char_volume_frac_fun(surf_pyro_r3(kin_param,t_vect_end_trunk,T_residual)');
    
    % Storing previous surface temp calculation
    T_residual_prev = T_residual;

    if virgin_CH
        % End loop if only considereing virgin Cole-Hopf transform
        conv = 0;
    else
        % Calculating surface temperature using new char fraction estimation
        T_residual = T_theta_transform_TPS_general(theta_TC_residual,char_frac_residual,beta_scale_TC,theta_T_transform_TPS_general) + T_0;
    
        % While-loop convergence criteria
       conv = norm((T_residual-T_residual_prev)./T_residual_prev);
    end
end

obj_fn = sum((T_residual-T_input.T(2:end)).^2);
display(['Objective Function: ',num2str(obj_fn)]);

%% Calculating surface temperature

% Repeating all Green's function construction formulation for the surface
G_g_loc_surf = cell(1,length(L_integrand_mat));

if q_interp_bool
    G_discrete_surf_kernel_back = time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L,L,GF_args).*(t_vect(2)-time_integrand_mat_quick)/(t_vect(2)-t_vect(1)),2);
    G_discrete_surf_kernel_central = time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L,L,GF_args).*(1 + (time_integrand_mat_quick-t_vect(2))/(t_vect(2)-t_vect(1))),2);
    G_discrete_surf = toeplitz([0;G_discrete_surf_kernel_back(1:end-1)],[0,zero_pad]) ...
        + toeplitz(G_discrete_surf_kernel_central,[G_discrete_surf_kernel_central(1),zero_pad]);


    for j = 1:L_length
        G_g_loc_surf_kernel_back = time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L,L_integrand_mat(j),GF_args).*(t_vect(2)-time_integrand_mat_quick)/(t_vect(2)-t_vect(1)),2);
        G_g_loc_surf_kernel_central = time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L,L_integrand_mat(j),GF_args).*(1 + (time_integrand_mat_quick-t_vect(2))/(t_vect(2)-t_vect(1))),2);
        G_g_loc_surf{j} = toeplitz([0;G_g_loc_surf_kernel_back(1:end-1)],[0,zero_pad])...
            + toeplitz(G_g_loc_surf_kernel_central,[G_g_loc_surf_kernel_central(1),zero_pad]);
    end

else
    G_discrete_surf_kernel = time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L,L,GF_args),2);
    G_discrete_surf = toeplitz(G_discrete_surf_kernel,[G_discrete_surf_kernel(1),zero_pad]);
    
    for j = 1:L_length
        G_g_loc_surf_kernel = time_coeff_vect_quick.*sum(time_weights_vect.*G_gal_global_fun(alpha_0,t_mat_quick,time_integrand_mat_quick,L,L_integrand_mat(j),GF_args),2);
        G_g_loc_surf{j} = toeplitz(G_g_loc_surf_kernel,[G_g_loc_surf_kernel(1),zero_pad]);
    end
end

g_spatial_integrand_surf = zeros(length(theta_TC),length(L_integrand_mat));

for j = 1:L_length % cycle through all r_0,z_0 in V2
    g_spatial_integrand_surf(:,j) = G_g_loc_surf{j}*g_integrand(:,j);
end
g_int_surf = L_coeff_vect*sum(L_weights_vect.*g_spatial_integrand_surf,2);

% Calculating transformed temperature at the surface
theta_surf = - (alpha_0/A)*G_discrete_surf*q_reconstruct_finite + (alpha_0/A)*g_int_surf;

% Simultaneously calculating surface temperature and surface char fraction

% Initial guess using virgin Cole-Hopf transform
T_surf = T_theta_transform_TPS(theta_surf) + T_0;

% Setting up loop to calculate surface temperature from theta-inversion
surf_conv = 1;
iter_surf = 0;
while all([surf_conv > 0.02, iter_surf <=5],2)
    iter_surf = iter_surf + 1;

    [drhodT_surf,rho_local_surf] = drhodT_fun(kin_param,t_vect(2:end),T_surf);

    %-------------------------------------------------------- in-situ
    %calculation of CH transformation at the surface
    
    dchidT_mass_surf = -(rho_char/(rho_virgin-rho_char))*(rho_virgin./(rho_local_surf.^2)).*drhodT_surf;
    char_mass_frac_vect_surf = char_mass_frac_fun(rho_local_surf);

    beta_indicator_vect_surf = (beta_scale_fun(T_surf-T_0,dchidT_mass_surf',char_mass_frac_vect_surf'));   
   
    start_char_surf = find(abs(char_mass_frac_vect_surf-0.02) == min(abs(char_mass_frac_vect_surf-0.02)));
    beta_cutoff_surf = start_char_surf;

    % solving for beta that best approximates k_bar/k = 1 for the relevant time
    % points
    X_beta_surf = beta_indicator_vect_surf(beta_cutoff_surf:end);

    if any(isinf(X_beta_surf))
        beta_scale_surf = 0;
    else
        anchor_soln = ones(size(beta_indicator_vect_surf(beta_cutoff_surf:end)));
        beta_scale_surf = ((X_beta_surf'*X_beta_surf)^-1)*X_beta_surf'*anchor_soln;
    end

    % if all([~isscalar(beta_scale_surf), size(beta_scale_surf,1) < 4],2)
    %     beta_scale_surf = [beta_scale_surf;zeros(4-size(beta_scale_surf,1),1)];
    % end

    % Calculating estimated surface char fraction
    char_frac_surf = char_volume_frac_fun(rho_local_surf');
    % Storing previous surface temp calculation
    T_surf_prev = T_surf;

    if virgin_CH
        % End loop if only considereing virgin Cole-Hopf transform
        surf_conv = 0;
    else
        % Calculating surface temperature using new char fraction estimation
        T_surf = T_theta_transform_TPS_general(theta_surf,char_frac_surf,beta_scale_surf,theta_T_transform_TPS_general) + T_0;
    
        % While-loop convergence criteria
        surf_conv = norm((T_surf-T_surf_prev)./T_surf_prev);
    end
end

% Calculating surface emissivity
eps_surf = eps_fun(T_surf,char_frac_surf);
sigma = 5.67E-8;

% T_infinity for radiative emission calculation
T_rad_inf = inputs.env.T_rad_inf;

% Calculating surface radiation emission
q_emit = -eps_surf.*sigma.*(T_surf.^4-T_rad_inf.^4);

% Calculating incident heat flux
q_incident = q_reconstruct_finite + q_emit;

% total_runtime = duration(end_time-start_time);
end_time = datetime('now');
disp('Total Runtime')
disp([num2str(milliseconds(duration(end_time-start_time))/1000),' seconds'])
%% Outputting results
mkdir('results')

% convert q to W/cm2
q_output = table(t_vect(1:end-1)', -q_reconstruct_finite/10000, -q_incident/10000, T_surf,...
    'VariableNames',{'t (s)','q_cond (W/cm2)','q_inc (W/cm2)','T_s (K)'});

writetable(q_output,'./results/solution.dat')
writematrix(L_integrand_mat,'./results/integration_points.dat')
writematrix(T_0 + T_prime_abis,'./results/T_global.dat')
writematrix(char_frac_abis,'./results/char_frac_global.dat')


