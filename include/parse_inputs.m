function inputs = parse_inputs(fname)
    inputs = [];
    f = fopen(fname,'r');

    while ~feof(f)
        s = strip(fgetl(f));
        if isempty(s)
            continue
        elseif strcmp(s(1),'#')
            continue
        end

        switch s(1)
            case '['
               local_field = genvarname(s(2:end-1));
               inputs.(local_field) = [];
            otherwise
                [field,val] = strtok(s,'=');
                field = strip(field);
                val = strip(val);

                % Remove '=' delimiter
                if strcmp(val(1),'=')
                    val(1) = [];
                    val = strip(val);
                end
                
                % Check if input is array
                if strcmp(val(1),'[')
                    val = val(2:end-1);
                    val = split(val,',');
                end

                if ~isnan(str2double(val))
                    val = str2double(val);
                end

                inputs.(local_field).(genvarname(field)) = val;
        end

        
end