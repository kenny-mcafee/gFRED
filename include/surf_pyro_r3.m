function rho_tot = surf_pyro_r3(kin_param,t_vect,T,varargin)
%{
Revision: 1
Revision Date: 9/6/2024
Shifted gamma (resin volume fraction) and phi (porosity) upstream to the
kinetic tables.
Generalized calculation of virgin and char densities.
Added interpolation to time and temperature vectors for more accurate ODE
solver --> reverted this change

Revision: 2
Revision Date: 9/9/2024
Switched to more robust ODE45 solver for decomposition reaction

Revision: 3
Revision Date: 9/13/2024
BIG REVISION - NOW OUTPUTS DENSITY

Revision: 3.1
Revision Date: 9/22/24
- Switch to lower order ODE23 method - runs 20% faster

Revision: 4
Revision Date: 7/17/25
- added option to select ode solver
%}
    ode_opt = 23;
    if ~isempty(varargin)
        ode_opt = varargin{1};
    end

    tol_opt = odeset('RelTol',1e-4,'AbsTol',1e-7);

    % Interpolating between measurements for better ODE accuracy (maybe?)
    % need to change inputs to t_vect_raw and T_raw

    % keyboard
    % d_rho = @(rho_0,rho_r,A,psi,T_act,T_min,T,rho) (rho>rho_r)*(T>=T_min)*(-A*rho_0*exp(-T_act/T)*((rho-rho_r)/rho_0)^psi);

    rho_indi = zeros(height(kin_param),length(t_vect));
    
    for i = 1:height(kin_param)
        temp_fun = @(T,rho) d_rho(kin_param.rho_0(i),kin_param.rho_r(i),kin_param.A(i),kin_param.psi(i)...
            ,kin_param.T_act(i),kin_param.T_min(i),T,rho);
        
        rho_indi(i,1) = kin_param.rho_0(i);

        if ode_opt == 23
            [~,rho_indi_lcl] = ode23(@(t,rho) subfn(t,t_vect,kin_param.rho_0(i),kin_param.rho_r(i),kin_param.A(i),kin_param.psi(i)...
                ,kin_param.T_act(i),kin_param.T_min(i),T,rho),t_vect,kin_param.rho_0(i),tol_opt);
        elseif ode_opt == 45
            [~,rho_indi_lcl] = ode45(@(t,rho) subfn(t,t_vect,kin_param.rho_0(i),kin_param.rho_r(i),kin_param.A(i),kin_param.psi(i)...
                ,kin_param.T_act(i),kin_param.T_min(i),T,rho),t_vect,kin_param.rho_0(i),tol_opt);
        elseif ode_opt == 89
            [~,rho_indi_lcl] = ode89(@(t,rho) subfn(t,t_vect,kin_param.rho_0(i),kin_param.rho_r(i),kin_param.A(i),kin_param.psi(i)...
                ,kin_param.T_act(i),kin_param.T_min(i),T,rho),t_vect,kin_param.rho_0(i),tol_opt);
        end

        rho_indi(i,:) = rho_indi_lcl;
        rho_indi(i,rho_indi(i,:)<kin_param.rho_r(i)) = kin_param.rho_r(i);
    end

    rho_tot = (1-kin_param.phi(1))*(kin_param.gamma_virgin(1)*(rho_indi(1,:)+rho_indi(2,:)) + (1-kin_param.gamma_virgin(1))*rho_indi(3,:));

    % rho_virgin = rho_tot(1);
    % rho_char = (1-kin_param.phi(1))*(kin_param.gamma_virgin(1)*(kin_param.rho_r(1)+kin_param.rho_r(2)) + (1-kin_param.gamma_virgin(1))*kin_param.rho_r(3));
    % 
    % char_frac = 1-((rho_virgin./(rho_virgin-rho_char)).*(1-rho_char./rho_tot));

function drhodt = subfn(t,t_vect,rho_0,rho_r,A,psi,T_act,T_min,T,rho)
    T_lcl = interp1(t_vect,T,t);
    drhodt = (rho>rho_r).*(T_lcl>=T_min).*(-A*rho_0*exp(-T_act./T_lcl).*((rho-rho_r)/rho_0).^psi);
end
end