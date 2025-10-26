function T = T_theta_transform_TPS_general(theta_input,char_frac,beta_scale,theta_transform)
    T = zeros(size(theta_input));
    for j = 1:size(theta_input,2)
        for i = 1:size(theta_input,1)
    
            if i == 1
                T(i,j)= fzero(@(T) theta_transform(char_frac(i,j),T,beta_scale(:,j))-theta_input(i,j),0);
            else
  
                T(i,j)= fzero(@(T) theta_transform(char_frac(i,j),T,beta_scale(:,j))-theta_input(i,j),T(i-1,j));
                % if any([T(i,j) <= -20, T(i,j) > 400])
                %     T(i,j) = T(i-1,j);
                % end
            end
        end
        % keyboard
    end
end
