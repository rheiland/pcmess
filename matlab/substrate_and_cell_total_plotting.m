%% Getting Data from Output

% NOTE: change this to your directory
base_dir = 'C:\Users\mason\OneDrive - Queensland University of Technology\PhD Notes\Publications\PhysiCell paper\PhysiCell\output'; 

plot_IL2_2D = false; % true to make animation of IL-2, false to create population plots faster
measure_cluster = false; % true to measure clustering
index = 1;

%automating the loading of the outputs
A1 = 'output0000000';
A2 = 'output000000';
A3 = 'output00000';
A4 = 'output0000';
B = '.xml';

timetotal = 0; % initialise total number of output files
runDir = base_dir; % get the directory of the first output (requires that at least one simulation was run)
while true % while the output files continue to exist
    if timetotal < 10
        K = [A1 num2str(timetotal,'%d') B]; % get the file name for the next output (assuming it exists)
    elseif timetotal < 100
        K = [A2 num2str(timetotal,'%d') B];
    elseif timetotal < 1000
        K = [A3 num2str(timetotal,'%d') B];
    else
        K = [A4 num2str(timetotal,'%d') B];
    end

    if isfile(fullfile(runDir, K)) % the current file is a file that exists in the output folder
        timetotal = timetotal + 1; % increment the total number of outputs
    else % the file does not exist
        break; % exit while loop
    end
end


cell_type_0 = zeros(timetotal,1);
cell_type_1 = zeros(timetotal,1);
InactiveT_all = zeros(timetotal,1);
ActiveT_all = zeros(timetotal,1);
IL2_mass = zeros(timetotal,1);
t_contact_time = 0;

avgdist_cells = zeros(timetotal,1); % average distance to closest neighbour
avgdist_rods = zeros(timetotal,1);
avgdist_both = zeros(timetotal,1);

% opening a video
if plot_IL2_2D
    v = VideoWriter('video1.avi');
    v.FrameRate = 10;
    open(v)
    figure('Position', [50 50 900 700])
    hold on 
end
for tcount = 1:timetotal

    % set up output name
    if tcount<11
        K = [A1 num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    elseif tcount<1001
        K = [A3 num2str(tcount-1,'%d') B];
    else
        K = [A4 num2str(tcount-1,'%d') B];
    end

    % reading output
    MCDS = read_MultiCellDS_xml(K,base_dir);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 

    cell_type_0(tcount) = length(find(MCDS.discrete_cells.metadata.type==0));
    cell_type_1(tcount) = length(find(MCDS.discrete_cells.metadata.type==1));
 	 

    %active v inactive v exhausted cells 

    types = MCDS.discrete_cells.metadata.type;
    T_c = (types == 0);
    state = MCDS.discrete_cells.custom.state;
    t_states = state(T_c);
    
    %ExhaustedT_all(tcount)   = sum(t_states > 1.5);
    InactiveT_all(tcount) = sum(t_states < 0.5);
    ActiveT_all(tcount) = sum(t_states == 1); 
    totalTcells(tcount) = length(T_c);
    
    IL2_conc = MCDS.continuum_variables.data;
    IL2_mass(tcount) = sum(IL2_conc(:) .* [MCDS.mesh.voxels.volume]'); % sum(sum(MCDS.continuum_variables(1).data(:,:,k)))*(MCDS.mesh.X_coordinates(2)-MCDS.mesh.X_coordinates(1))*(MCDS.mesh.Y_coordinates(2)-MCDS.mesh.Y_coordinates(1));

    % Contact time 

    contact_time = MCDS.discrete_cells.custom.attached_time;
    t_contact_time = contact_time(T_c);


    % Measure clustering
    if measure_cluster
        agent_pos = MCDS.discrete_cells.state.position(:,1:2);
        for i = 1:length(types)
            closest_cells = Inf;
            closest_rods = Inf;
            closest_both = Inf;
            for j = [1:i-1 i+1:length(types)] % for each pair of agents (except agent i)
                curr_dist = norm(agent_pos(i,:)-agent_pos(j,:)); % current distance between both agents
                if types(i) == 0
                    if types(j) == 0 % both T cells
                        if curr_dist < closest_cells % current distance is closer than the prev closest
                            closest_cells = curr_dist;
                        end
                    else % i is T cell, j is rod
                        if curr_dist < closest_both % current distance is closer than the prev closest
                            closest_both = curr_dist;
                        end
                    end
                elseif types(j) == 1 % both rods
                    if curr_dist < closest_rods % current distance is closer than the prev closest
                        closest_rods = curr_dist;
                    end
                end
            end
            if types(i) == 0
                avgdist_cells(tcount) = avgdist_cells(tcount) + closest_cells; % add closest distance to current timestep
                avgdist_both(tcount) = avgdist_both(tcount) + closest_both;
            else
                avgdist_rods(tcount) = avgdist_rods(tcount) + closest_rods;
            end
        end
        avgdist_cells(tcount) = avgdist_cells(tcount)/sum(types==0); % average the results
        avgdist_rods(tcount) = avgdist_rods(tcount)/sum(types==1);
        avgdist_both(tcount) = avgdist_both(tcount)/sum(types==0);
    end


    %plotting contour plot of substrate
    if plot_IL2_2D
        contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
            MCDS.continuum_variables(index).data(:,:,k), 20 );
        axis image;
        colorbar; 
        xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
        ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
        title( sprintf('%s (%s) at t = %3.2f days', MCDS.continuum_variables.name , ...
            MCDS.continuum_variables.units , ...
            MCDS.metadata.current_time/24/60 ) ); 
        frame = getframe(gcf);
        writeVideo(v,frame);
    end

    fprintf('%f%% completed.', round(tcount/timetotal*100,2))

end


simulation_time = MCDS.metadata.current_time/24/60; % number of days in simulation



% Plotting

% total cells 
%figure
%plot(linspace(0,simulation_time,timetotal),cell_type_0,'LineWidth',2)
%xlabel('Time (days)')
%ylabel('Cell type 0')
%set(gca,'FontSize',18)

% total rod 
%figure
%plot(linspace(0,simulation_time,timetotal),cell_type_1,'LineWidth',2)
%xlabel('Time (days)')
%ylabel('Cell type 1')
%set(gca,'FontSize',18)

% Active v inactive v exhuasted cells 
figure 
plot(linspace(0,simulation_time,timetotal),InactiveT_all, 'LineWidth',2)
hold on
plot(linspace(0,simulation_time,timetotal),ActiveT_all, 'LineWidth',2)
%hold on 
%plot(linspace(0,simulation_time,timetotal),ExhaustedT_all, '.', 'Linewidth',2)
%xticks([1 25 49 73 97 121 145 169 193 217 241 265 289 313 336])
%xticklabels({'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14'})
%title('Activated v Naive T Cells', 'a_{max} = 1, s = 1')
xlabel('Time (days)')
ylabel('T cell count')
legend('Naive', 'Activated', 'Location','best')
set(gca,'FontSize',18)

figure 
plot(linspace(0,simulation_time,timetotal),IL2_mass, 'LineWidth',2)
yline(3.47e-12*MCDS.mesh.X(end)*MCDS.mesh.Y(end), '--')
%title('Activated v Naive T Cells', 'a_{max} = 1, s = 1')
xlabel('Time (days)')
ylabel('IL-2 mass (ng)')
set(gca,'FontSize',18)

% Contact time histogram 
figure, 
histogram(t_contact_time, 20)
hold on 
xline(720)
xlabel('Contact time')
ylabel('Number of T cells')
title('Total contact time of T cells with micro-rods')

% Clustering plots
if measure_cluster
    figure 
    plot(linspace(0,simulation_time,timetotal),avgdist_cells, 'LineWidth',2)
    hold on
    plot(linspace(0,simulation_time,timetotal),avgdist_rods, 'LineWidth',2)
    hold on
    plot(linspace(0,simulation_time,timetotal),avgdist_both, 'LineWidth',2)
    xlabel('Time (days)')
    ylabel('Average distance to nearest neighbour')
    legend('T cells', 'Rods', 'Both', 'Location','best')
    set(gca,'FontSize',18)
end





%% testing rod position and rotation

num_cells = 500;
num_rods = 100; 

types = zeros(num_cells+num_rods,1); types(1:num_rods) = 1;

cell_mean = 300; cell_var = 1;
rod_mean = 100; rod_var = 100;

cell_pos = cell_mean + cell_var*randn(num_cells,2);
rod_pos = rod_mean + rod_var*randn(num_rods,2);
agent_pos = [rod_pos; cell_pos];

%agent_pos = MCDS.discrete_cells.state.position(:,1:2);
%cell_pos = agent_pos(types==0,:);
%rod_pos = agent_pos(types==1,:);
%rod_rot = MCDS.discrete_cells.state.orientation(types==1,:);

avgdist_cells = 0;
avgdist_rods = 0;
avgdist_both = 0;
for i = 1:length(types)
    closest_cells = Inf;
    closest_rods = Inf;
    closest_both = Inf;
    for j = [1:i-1 i+1:length(types)] % for each pair of agents (except agent i)
        curr_dist = norm(agent_pos(i,:)-agent_pos(j,:)); % current distance between both agents
        if types(i) == 0
            if types(j) == 0 % both T cells
                if curr_dist < closest_cells % current distance is closer than the prev closest
                    closest_cells = curr_dist;
                end
            else % i is T cell, j is rod
                if curr_dist < closest_both % current distance is closer than the prev closest
                    closest_both = curr_dist;
                end
            end
        elseif types(j) == 1 % both rods
            if curr_dist < closest_rods % current distance is closer than the prev closest
                closest_rods = curr_dist;
            end
        end
    end
    if types(i) == 0
        avgdist_cells = avgdist_cells + closest_cells; % add closest distance to current timestep
        avgdist_both = avgdist_both + closest_both;
    else
        avgdist_rods = avgdist_rods + closest_rods;
    end
end
avgdist_cells = avgdist_cells/sum(types==0); % average the results
avgdist_rods = avgdist_rods/sum(types==1);
avgdist_both = avgdist_both/sum(types==0);
