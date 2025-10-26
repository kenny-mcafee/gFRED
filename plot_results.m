results_data = readtable('./results/solution.dat');
results_data.Properties.VariableNames = {'t','q_cond','q_inc','T_s'};

figure(1)
clf(1)
plot(results_data.t,results_data.q_cond,'Color',"#0072BD",'LineWidth',1,'DisplayName','q_{cond}')
hold on
plot(results_data.t,results_data.q_inc,'Color',"#D95319",'LineWidth',1,'DisplayName','q_{inc}')
xlabel('Time (s)')
ylabel('Heat Flux (W/cm^2)')
set(gca,'FontSize',16,'FontName','serif')

figure(2)
clf(2)
plot(results_data.t,results_data.T_s,'Color','Black','LineWidth',1,'DisplayName','T_s')
xlabel('Time (s)')
ylabel('Surface Temperature (K)')
set(gca,'FontSize',16,'FontName','serif')