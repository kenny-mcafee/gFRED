function f_del_out = f2_del(x,N,L_s,L_RTV,k_0,k_s,k_RTV)
    f_del_out = f15_del(x,N,L_s,L_RTV,k_s,k_RTV) ...
        - (L_s+L_RTV)^-1*big_gamma_23(N,L_RTV,L_s,k_0,k_s,k_RTV).*((1:N)').*(1+x/(L_s+L_RTV)).^((1:N)'-1);
end
