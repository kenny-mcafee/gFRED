function kin_param = load_kinetics(fname)
%{
Revision: 1
Revision Date: 9/6/2024

Included definition of solid porosity and resin volume fraction
%}
    kin_param_imperial = readtable(fname);

    if length(kin_param_imperial.Properties.VariableNames) == 7
        kin_param_imperial.Properties.VariableNames = {'mat','rho_0','rho_r','A','psi','T_act','T_min'};
    elseif length(kin_param_imperial.Properties.VariableNames) == 9
        kin_param_imperial.Properties.VariableNames = {'mat','rho_0','rho_r','A','psi','T_act','T_min','gamma_virgin','phi'};
    end

    kin_param = table();
    kin_param.mat = kin_param_imperial.mat;
    kin_param.rho_0 = kin_param_imperial.rho_0*16.0185;
    kin_param.rho_r = kin_param_imperial.rho_r*16.0185;
    kin_param.A = kin_param_imperial.A;
    kin_param.psi = kin_param_imperial.psi;
    kin_param.T_act = kin_param_imperial.T_act/1.8;
    kin_param.T_min = kin_param_imperial.T_min/1.8;

    if length(kin_param_imperial.Properties.VariableNames) == 7 %defaults
        kin_param.phi = 0*ones(height(kin_param),1);
        kin_param.gamma_virgin = 0.5*ones(height(kin_param),1);
    elseif length(kin_param_imperial.Properties.VariableNames) == 9
        kin_param.phi = kin_param_imperial.phi;
        kin_param.gamma_virgin = kin_param_imperial.gamma_virgin;
    end

end