function [matl_fun,coeffs_out] = matl_table_to_fun(data,n)
%{
Revision: 1
Revision Date: 10/16/24

- Added check to reduce n if not enough data
size(data,1) <= n

Revision: 2
Revision Date: 6/16/25
- Changed K_T_out to coeffs_out
%}
    if size(data,1) <= n
        n = size(data,1) - 1;
    end
    
    X = data(:,1).^(0:n);
    Y = data(:,2);
    
    coeffs = (X'*X)^-1*X'*Y;
    
    fun_str = '@(T_abs) ';
    for i = 1:n+1
        if i == 1
            fun_str = [fun_str,num2str(coeffs(i))];
        else
        fun_str = [fun_str,' + ',num2str(coeffs(i)),'*T_abs.^',num2str(i-1)];
        end
    end

    matl_fun = str2func(fun_str);
    if n >= 1
        coeffs_out = coeffs;
    else
        coeffs_out = [];
    end
end