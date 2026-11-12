

base_dir = 'C:/Users/adria/Documents/PhysiMeSS project data/';
specific_folder = 'att_12hrs_b08_fs1_cca04';

runDir = [base_dir,specific_folder];

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
    MCDS = read_MultiCellDS_xml(K,runDir);
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
    
    % Contact time 

    contact_time = MCDS.discrete_cells.custom.attached_time;
    t_contact_time = contact_time(T_c);

    %if tcount = 1*60*24%,3, 5
    contact_time_distribution{tcount} = t_contact_time;

    Tcellagent_position_matrix{tcount} = MCDS.discrete_cells.state.position(T_c,:);
    R_c = (types == 1);
    Rodagent_pos_orien_length_matrix{tcount} = [MCDS.discrete_cells.state.position(R_c,:),...
        MCDS.discrete_cells.state.orientation(R_c,:),MCDS.discrete_cells.custom.f_length(R_c)'];
    %end
end

save(specific_folder, 'Tcellagent_position_matrix','Rodagent_pos_orien_length_matrix',...
    'totalTcells', 'ActiveT_all', 'InactiveT_all', ...
    'contact_time_distribution');


%% plotting to check

figure
scatter(Tcellagent_position_matrix{tcount}(:,1),Tcellagent_position_matrix{tcount}(:,2))
hold on 
for i = 1:length(MCDS.discrete_cells.custom.f_length(R_c))
x1 = Rodagent_pos_orien_length_matrix{tcount}(i,1)-Rodagent_pos_orien_length_matrix{tcount}(i,7)/2*Rodagent_pos_orien_length_matrix{tcount}(i,4);
y1 = Rodagent_pos_orien_length_matrix{tcount}(i,2)-Rodagent_pos_orien_length_matrix{tcount}(i,7)/2*Rodagent_pos_orien_length_matrix{tcount}(i,5);
x2 = Rodagent_pos_orien_length_matrix{tcount}(i,1)+Rodagent_pos_orien_length_matrix{tcount}(i,7)/2*Rodagent_pos_orien_length_matrix{tcount}(i,4);
y2 = Rodagent_pos_orien_length_matrix{tcount}(i,2)+Rodagent_pos_orien_length_matrix{tcount}(i,7)/2*Rodagent_pos_orien_length_matrix{tcount}(i,5);

plot([x1 x2],[y1 y2],'b')
end


