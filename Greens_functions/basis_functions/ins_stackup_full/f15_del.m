function f_out = f15_del(x,N,L_s,L_RTV,k_s,k_RTV)
    f_out = f1_del(x,N,L_s,L_RTV) ...
        - L_s^-1*((1:N)').*big_gamma_12(N,L_s,k_RTV,k_s).*((1+L_RTV/L_s+x/L_s).^((1:N)'-1));
end