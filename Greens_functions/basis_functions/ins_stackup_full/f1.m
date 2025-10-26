function f_out = f1(x,N,L_s,L_RTV)
    f_out = (x+L_s+L_RTV).^(2*(1:N)'-2);
end