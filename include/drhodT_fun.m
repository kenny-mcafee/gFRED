function [drhodT_tot,rho_tot] = drhodT_fun(kin_param,t_vect,T)
%{
Revision: 1
Revision Date: 2/23/2024

%}

    % Interpolating between measurements for better ODE accuracy (maybe?)
    % need to change inputs to t_vect_raw and T_raw
    rho_indi = zeros(height(kin_param),length(t_vect));
    drhodt_indi = rho_indi;
    drhodT_indi = rho_indi;

    %ODE23 tolerance options
    tol_opt = odeset('RelTol',[],'AbsTol',[]);
    
    for i = 1:height(kin_param)
        temp_fun = @(T,rho) d_rho(kin_param.rho_0(i),kin_param.rho_r(i),kin_param.A(i),kin_param.psi(i)...
            ,kin_param.T_act(i),kin_param.T_min(i),T,rho);
        
        rho_indi(i,1) = kin_param.rho_0(i);
        % for j = 2:length(t_vect)
        %     rho_indi(i,j) = rho_indi(i,j-1)+(t_vect(j)-t_vect(j-1))*temp_fun(T(j-1),rho_indi(i,j-1));
        % end
        [~,rho_indi_lcl] = ode23(@(t,rho) subfn1(t,t_vect,kin_param.rho_0(i),kin_param.rho_r(i),kin_param.A(i),kin_param.psi(i)...
            ,kin_param.T_act(i),kin_param.T_min(i),T,rho),t_vect,kin_param.rho_0(i),tol_opt);

        rho_indi(i,:) = rho_indi_lcl;
        rho_indi(i,rho_indi(i,:)<kin_param.rho_r(i)) = kin_param.rho_r(i);

        drhodt_indi(i,:) = subfn1(t_vect,t_vect,kin_param.rho_0(i),kin_param.rho_r(i),kin_param.A(i),kin_param.psi(i)...
            ,kin_param.T_act(i),kin_param.T_min(i),T,rho_indi(i,:));

        [~,drhodT_indi_lcl] = ode23(@(t,drhodT) subfn2(t,t_vect,kin_param.rho_0(i),kin_param.rho_r(i),kin_param.psi(i)...
            ,kin_param.T_act(i),T,rho_indi(i,:)',drhodt_indi(i,:)',drhodT),t_vect,0,tol_opt);
        drhodT_indi(i,:) = drhodT_indi_lcl;
    end

    rho_tot = (1-kin_param.phi(1))*(kin_param.gamma_virgin(1)*(rho_indi(1,:)+rho_indi(2,:)) + (1-kin_param.gamma_virgin(1))*rho_indi(3,:));

    drhodT_tot = (1-kin_param.phi(1))*(kin_param.gamma_virgin(1)*(drhodT_indi(1,:)+drhodT_indi(2,:)) + (1-kin_param.gamma_virgin(1))*drhodT_indi(3,:));


function drhodt = subfn1(t,t_vect,rho_0,rho_r,A,psi,T_act,T_min,T,rho)
    T_lcl = interp1(t_vect,T,t);
    drhodt = (rho>rho_r).*(T_lcl>=T_min).*(-A*rho_0*exp(-T_act./T_lcl).*((rho-rho_r)/rho_0).^psi);
end

function drhodT_dt = subfn2(t,t_vect,rho_0,rho_r,psi,T_act,T,rho,drhodt,drhodT)
    T_lcl = interp1(t_vect,T,t);
    rho_lcl = interp1(t_vect,rho,t);
    drhodt_lcl = interp1(t_vect,drhodt,t);
    drhodT_dt = drhodt_lcl.*(T_act./(T_lcl.^2) + (psi/rho_0).*(((rho_lcl-rho_r)/rho_0).^-1).*drhodT);
    if isnan(drhodT_dt)
        if drhodt_lcl == 0
            drhodT_dt = 0;
        end
    end
end
end