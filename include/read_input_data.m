function T_input = read_input_data(path,downsample)

    % if strcmp(path(1),'.')
    %     path = ['.', path];
    % end
    T_input = readtable(path);
    T_input.norm_t = T_input.t - T_input.t(1);
    T_input = T_input(1:downsample:end,:);
end