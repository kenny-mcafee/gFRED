function out = big_gamma_23(N,L_RTV,L_s,k_0,k_s,k_RTV)
    out = (L_s+L_RTV)*(((2*(1:N)'-2)./((1:N)')).*(L_s+L_RTV).^(2*(1:N)'-3) ...
        - big_gamma_12(N,L_s,k_RTV,k_s)*L_s^-1.*(1+L_RTV/L_s).^((1:N)'-1)).*(1-k_RTV/k_0);

end