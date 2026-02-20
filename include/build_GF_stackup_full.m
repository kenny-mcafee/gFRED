function [V,gamma] = build_GF_stackup_full(N,TPS_props,Struct_props,RTV_props,theta_scale) 

    [k_0,k_T,alpha_0,L] = deal(TPS_props{1},TPS_props{2},TPS_props{3},TPS_props{4});
    [k_s,alpha_s,L_s] = deal(Struct_props{1},Struct_props{2},Struct_props{3});
    [k_RTV,alpha_RTV,L_RTV] = deal(RTV_props{1},RTV_props{2},RTV_props{3});

    % Load in legendre polynomial abscissa's and quadrature weights
    % sheet_names = sheetnames("Legendre_roots.xlsx");
    % legendre_general = cellfun(@(x) readmatrix('Legendre_roots.xlsx','Sheet',x),sheet_names,'UniformOutput',false);
    legendre_general = quadrature_init('Legendre_roots.dat');

    % Calculating all volume/area/length numerical integration parameters
    [L_coeff_vect_eig,L_integrand_mat_eig,L_weights_vect_eig] = gauss_quadrature(L,0,4*(N-1),legendre_general);
    [L_s_coeff_vect_eig,L_s_integrand_mat_eig,L_s_weights_vect_eig] = gauss_quadrature(-L_RTV,-L_s-L_RTV,4*(N-1),legendre_general);
    [L_RTV_coeff_vect_eig,L_RTV_integrand_mat_eig,L_RTV_weights_vect_eig] = gauss_quadrature(0,-L_RTV,4*(N-1),legendre_general);

    f1_vect = f1((L_s_integrand_mat_eig'),N,L_s,L_RTV);
    grad_f1_vect = f1_del((L_s_integrand_mat_eig'),N,L_s,L_RTV);

    f15_vect = f15((L_RTV_integrand_mat_eig'),N,L_s,L_RTV,k_s,k_RTV);
    grad_f15_vect = f15_del((L_RTV_integrand_mat_eig'),N,L_s,L_RTV,k_s,k_RTV);

    f2_vect = f2((L_integrand_mat_eig'),N,L_s,L_RTV,k_0,k_s,k_RTV);
    grad_f2_vect = f2_del((L_integrand_mat_eig'),N,L_s,L_RTV,k_0,k_s,k_RTV);

    a1 = -(k_s)*L_s_coeff_vect_eig*((L_s_weights_vect_eig.*grad_f1_vect)*grad_f1_vect');
    a15 = -k_RTV*L_RTV_coeff_vect_eig*((L_RTV_weights_vect_eig.*grad_f15_vect)*grad_f15_vect');
    a2 = -(k_0)*L_coeff_vect_eig*((L_weights_vect_eig.*grad_f2_vect)*grad_f2_vect');

    a = a1 + a15 + a2;

    b1 = (k_s/alpha_s)*L_s_coeff_vect_eig*((L_s_weights_vect_eig.*f1_vect)*f1_vect');
    b15 = (k_RTV/alpha_RTV)*L_RTV_coeff_vect_eig*((L_RTV_weights_vect_eig.*f15_vect)*f15_vect');
    b2 = (k_0/alpha_0)*L_coeff_vect_eig*((L_weights_vect_eig.*f2_vect)*f2_vect');

    b = b1 + b15 + b2;
    
    a = (a+a')/2;
    b = (b+b')/2;
    
    [R,flag] = chol(b);
    if flag ~= 0
        error('b matrix was not symmetric positive definite.')
    end

    %Computing eigenvalue problem
    [V_bar,D] = eig(-(R'^-1*a*(R)^-1));

    %note: outputs have to be floats
    V = double((R)^-1*V_bar);
    gamma = double(diag(D));
    

    % %-----------diagnostics
    % ortho_test = double(V'*b*V)
    % ortho_test = V*(V')*b
    % test1 = (b^-1*a) - (b^-1*a)';
    % test2 = theta_scale*V'*(((theta_scale*V*b)')^-1) - V'*(((V*b)')^-1)

    % psi1 = V'*f1_vect;
    % psi15 = V'*f15_vect;
    % psi2 = V'*f2_vect;
    % 
    % grad_psi1 = V'*grad_f1_vect;
    % grad_psi15 = V'*grad_f15_vect;
    % grad_psi2 = V'*grad_f2_vect;

    % figure(10001)
    % clf(10001)
    % tiledlayout('flow')
    % 
    % for i = 1:size(psi1,1)
    %     nexttile
    %     hold off
    %     plot(L_s_integrand_mat_eig',psi1(i,:))
    %     hold on
    %     plot(L_RTV_integrand_mat_eig',psi15(i,:))
    %     plot(L_integrand_mat_eig',psi2(i,:))
    % end
    % 
    % figure(10002)
    % clf(10002)
    % for i = 1:size(psi1,1)
    %     nexttile
    %     hold off
    %     plot(L_s_integrand_mat_eig',(k_s/k_0)*grad_psi1(i,:))
    %     hold on
    %     plot(L_RTV_integrand_mat_eig',(k_RTV/k_0)*grad_psi15(i,:))
    %     plot(L_integrand_mat_eig',grad_psi2(i,:))
    % end

end

