function G_gal_global = G_gal_global_fun(alpha_scale,t,tau,x,x_0,varargin)
  
    varargin = varargin{:};
    while ~isempty(varargin)
        switch varargin{1}
            case 'N'
                N = varargin{2};
            case 'V'
                V = varargin{2};
            case 'gamma'
                gamma = varargin{2};
            case 'k_0'
                k_0 = varargin{2};
            case 'k_T'
                k_T = varargin{2};
            case 'k_s'
                k_s = varargin{2};
            case 'k_RTV'
                k_RTV = varargin{2};
            case 'A'
                A = varargin{2};
            case 'R_th'
                R_th = varargin{2};
            case 'L_s'
                L_s = varargin{2};
            case 'L'
                L = varargin{2};
            case 'L_RTV'
                L_RTV = varargin{2};
        end
        varargin(1:2) = [];
    end
    
    f = f2(x,N,L_s,L_RTV,k_0,k_s,k_RTV);
    f_0 = f2(x_0,N,L_s,L_RTV,k_0,k_s,k_RTV);

    psi = V'*f;
    psi_0 = V'*f_0;
  
    G_gal_global_flat = (k_0/alpha_scale)*(psi.*psi_0)'*exp(-gamma*(reshape(t,1,[])-reshape(tau,1,[])));
    G_gal_global = reshape(G_gal_global_flat,size(t,1),[]);
  
end

