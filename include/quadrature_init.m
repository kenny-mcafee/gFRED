function legendre_general = quadrature_init(fname)
    f = fopen(fname,'r');
    legendre_general = {};
    i=1;
    while ~feof(f)
        header = str2double(fgetl(f));
        if i~=header
            continue
        end
        local_points = str2double(split(fgetl(f)))';
        local_points(isnan(local_points)) = [];

        local_weights = str2double(split(fgetl(f)))';
        local_weights(isnan(local_weights)) = [];
        legendre_general{i} = [local_points;local_weights];
        i=i+1;
    end
    fclose(f);
end