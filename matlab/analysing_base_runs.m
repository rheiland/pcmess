clear all
close all

%% load in base runs
baserun_1 = load("baserun_1.mat");
baserun_2 = load("baserun_2.mat");
baserun_3 = load("baserun_3.mat");
baserun_4 = load("baserun_4.mat");
baserun_5 = load("baserun_5.mat");
baserun_6 = load("baserun_6.mat");
baserun_7 = load("baserun_7.mat");

fig1 = figure
%% Plotting activated cells over time 

mat_activated = [baserun_1.ActiveT_all(1:50),...
    baserun_2.ActiveT_all(1:50),...
    baserun_3.ActiveT_all(1:50),...
    baserun_4.ActiveT_all(1:50),...
    baserun_5.ActiveT_all(1:50),...
    baserun_6.ActiveT_all(1:50),...
    baserun_7.ActiveT_all(1:50)];

figure(fig1)
options.handle = fig1;
options.alpha = 0.5;
options.color_area = [0.8 0.2 0.1];
options.color_line = [0.8 0.2 0.1];
options.line_width = 2;
options.error = 'std';
hold on
plot_areaerrorbar(mat_activated',options)
xlabel('time')
ylabel('cells')
set(gca,'FontSize',18)


%% Plotting naive cells over time 

mat_naive = [baserun_1.InactiveT_all(1:50),...
    baserun_2.InactiveT_all(1:50),...
    baserun_3.InactiveT_all(1:50),...
    baserun_4.InactiveT_all(1:50),...
    baserun_5.InactiveT_all(1:50),...
    baserun_6.InactiveT_all(1:50),...
    baserun_7.InactiveT_all(1:50)];

figure(fig1)
options.handle = fig1;
options.color_area = [0.4 0.2 0.8];
options.color_line = [0.4 0.2 0.8];
options.alpha = 0.5;
options.line_width = 2;
options.error = 'std';
hold on
plot_areaerrorbar(mat_naive',options)
xlabel('time')
ylabel('cells')
set(gca,'FontSize',18)

%% Total IL-2 over time

fig2 = figure

mat_IL2 = [baserun_1.IL2_mass(1:50),...
    baserun_2.IL2_mass(1:50),...
    baserun_3.IL2_mass(1:50),...
    baserun_4.IL2_mass(1:50),...
    baserun_5.IL2_mass(1:50),...
    baserun_6.IL2_mass(1:50),...
    baserun_7.IL2_mass(1:50)];

figure(fig2)
options.handle = fig2;
options.color_area = [0.1 0.8 0.8];
options.color_line = [0.1 0.8 0.8];
options.alpha = 0.5;
options.line_width = 2;
options.error = 'std';
hold on
plot_areaerrorbar(mat_IL2',options)
xlabel('time')
ylabel('IL-2 mass')
set(gca,'FontSize',18)

%% Plotting contact time distribution at different time points for one simulation

tcount_1day = 10; % NEED TO WORK OUT THE CORRECT VALUES FOR THESE, NOW JUST USING DUMMY VALUES
tcount_2day = 20; % NEED TO WORK OUT THE CORRECT VALUES FOR THESE, NOW JUST USING DUMMY VALUES
tcount_3day = 30; % NEED TO WORK OUT THE CORRECT VALUES FOR THESE, NOW JUST USING DUMMY VALUES
tcount_4day = 40; % NEED TO WORK OUT THE CORRECT VALUES FOR THESE, NOW JUST USING DUMMY VALUES

xlim_max = max(baserun_1.contact_time_distribution{end});

figure
subplot(1,4,1)
histogram(baserun_1.contact_time_distribution{tcount_1day})
xlim([0 xlim_max])
xlabel('Contact time (mins)')
ylabel('Population frequency')
set(gca,'FontSize',18)
title('1 day')

subplot(1,4,2)
histogram(baserun_1.contact_time_distribution{tcount_2day})
xlim([0 xlim_max])
xlabel('Contact time (mins)')
ylabel('Population frequency')
set(gca,'FontSize',18)
title('2 day')

subplot(1,4,3)
histogram(baserun_1.contact_time_distribution{tcount_3day})
xlim([0 xlim_max])
xlabel('Contact time (mins)')
ylabel('Population frequency')
set(gca,'FontSize',18)
title('3 day')

subplot(1,4,4)
histogram(baserun_1.contact_time_distribution{tcount_4day})
xlim([0 xlim_max])
title('4 day')
xlabel('Contact time (mins)')
ylabel('Population frequency')
set(gca,'FontSize',18)




%% plotting Mason's metric at different time points 

tcount_1day = 10; % NEED TO WORK OUT THE CORRECT VALUES FOR THESE, NOW JUST USING DUMMY VALUES
tcount_2day = 20; % NEED TO WORK OUT THE CORRECT VALUES FOR THESE, NOW JUST USING DUMMY VALUES
tcount_3day = 30; % NEED TO WORK OUT THE CORRECT VALUES FOR THESE, NOW JUST USING DUMMY VALUES
tcount_4day = 40; % NEED TO WORK OUT THE CORRECT VALUES FOR THESE, NOW JUST USING DUMMY VALUES

[cluster_cells_1, g_cc] = clustering_metric(baserun_1);
[cluster_cells_2, g_cc] = clustering_metric(baserun_2);
[cluster_cells_3, g_cc] = clustering_metric(baserun_3);
[cluster_cells_4, g_cc] = clustering_metric(baserun_4);
[cluster_cells_5, g_cc] = clustering_metric(baserun_5);
[cluster_cells_6, g_cc] = clustering_metric(baserun_6);
[cluster_cells_7, g_cc] = clustering_metric(baserun_7);

%% plotting the clustering metric

cluster_cells_mat = [cluster_cells_1;...
    cluster_cells_2;...
    cluster_cells_3;...
    cluster_cells_4;...
    cluster_cells_5;...
    cluster_cells_6;...
    cluster_cells_7];

figure 
plot_areaerrorbar(cluster_cells_mat)
ylabel('Clustering metric (Mason)')
xlabel('Time')
title('average clustering metric over time')
set(gca,'FontSize',18)

figure
plot(g_cc_store{10})
hold on 
plot(g_cc_store{20})
plot(g_cc_store{30})
plot(g_cc_store{40})
legend('time step 10', 'time step 20', 'time step 30', 'time step 40')
xlabel('radius')
ylabel('g_cc')
title('Pair clustering function at different time points')
set(gca,'FontSize',18)

