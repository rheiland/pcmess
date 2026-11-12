/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"



void create_cell_types( void )
{
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 

	cell_defaults.functions.instantiate_cell = instantiate_physimess_cell;	
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = physimess_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
    
    cell_defaults.functions.plot_agent_SVG = fibre_agent_SVG;
	cell_defaults.functions.plot_agent_legend = fibre_agent_legend;
	
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	cell_defaults.functions.cell_division_function = tcell_division_function; 
	cell_defaults.functions.calculate_distance_to_membrane = custom_distance_to_domain; 
    
	for (auto* pCD: *getFibreCellDefinitions()){
		pCD->functions.instantiate_cell = instantiate_physimess_fibre;
		pCD->functions.plot_agent_SVG = fibre_agent_SVG;
		pCD->functions.plot_agent_legend = fibre_agent_legend;
        pCD->functions.update_phenotype = fibre_time_secretion_function;   // <-- new
	
	}

	for (auto* pCD: cell_definitions_by_index){
		if (!isFibre(pCD) && pCD->custom_data.find_variable_index("fibre_custom_degradation") > 0){	
			if (pCD->custom_data["fibre_custom_degradation"] > 0.5)
				pCD->functions.instantiate_cell = instantiate_physimess_cell_custom_degrade;	
		}
	}
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 

    // code to compute x_max, y_max, and number_of_fibres from Brynn's digitised scaffold 
    microenvironment.mesh.bounding_box[0] = 0; microenvironment.mesh.bounding_box[3] = default_microenvironment_options.X_range[1]; // reassign x_max and y_max in here since these default to 0.5 for some reason
    microenvironment.mesh.bounding_box[1] = 0; microenvironment.mesh.bounding_box[4] = default_microenvironment_options.Y_range[1];
    double Xmax = microenvironment.mesh.bounding_box[3]; double Ymax = microenvironment.mesh.bounding_box[4]; 
    std::cout<< "Maximum x and y: " << Xmax <<std::endl;

    // load in scaffold spreadsheet and add up the mass of rods (code stolen from PhysiCell_geometry -> load_cells_csv_v1)
    pugi::xml_node node = physicell_config_root.child( "initial_conditions" ).child( "cell_positions" ); // find where the file name is defined
    std::string folder = xml_get_string_value( node, "folder" ); 
	std::string filename = xml_get_string_value( node, "filename" ); 
	std::string input_filename = folder + "/" + filename; // define file name/directory
    std::ifstream file( input_filename, std::ios::in ); // open file
    std::string line;
    double rod_mass = 0;
	while (std::getline(file, line)) // for each row in the file
	{
		std::vector<double> data;
		csv_to_vector( line.c_str() , data ); // extract the row and convert to array

        parameters.ints("number_of_fibres") = parameters.ints("number_of_fibres") + 1; // increment the number of rods
        rod_mass = rod_mass + 6.28318530718*parameters.doubles("rod_radius")*(parameters.doubles("rod_radius")+data[4])/parameters.doubles("rod_surface_area_ratio"); // increment total mass of rods by dividing the average surface area by the ratio of surface area to mass, using length of individual rod
	}
	file.close();
    std::cout<< "Number of rods: " << parameters.ints("number_of_fibres") <<std::endl;
    std::cout<< "Total mass of rods: " << rod_mass <<std::endl;

    //parameters.ints("number_of_fibres") = 411; // define number of rods from Brynn's digitised scaffold
    //parameters.doubles("rod_mass") = 6.28318530718*parameters.doubles("rod_radius")*(parameters.doubles("rod_radius")+parameters.doubles("rod_length"))/parameters.doubles("rod_surface_area_ratio"); // compute mass of individual rod by dividing the average surface area by the ratio of surface area to mass

    parameters.doubles("scaling_factor") = parameters.doubles("tot_rod_mass")/rod_mass; // compute scaling factor by comparing total mass of rods in our simulation to the total mass of rods used in experiments
    std::cout<< "Scaling factor: " << parameters.doubles("scaling_factor") <<std::endl;

    //parameters.doubles("rod_mass") = parameters.doubles("rod_mass")/parameters.ints("number_of_fibres"); // replace rod_mass with the average mass of rods in domain (for use in secretion calculation)
    //std::cout<< "Average rod mass: " << parameters.doubles("rod_mass") <<std::endl;

    microenvironment.mesh.bounding_box[2] = 0; microenvironment.mesh.bounding_box[5] = parameters.doubles("exp_vol")/(Xmax*Ymax*parameters.doubles("scaling_factor")); // compute scaling in z direction by comparing simulation domain to scaling factor
    double Zmax = microenvironment.mesh.bounding_box[5]; 
    default_microenvironment_options.dz = Zmax; 
    std::cout<< "Domain height (z direction): " << Zmax <<std::endl;

    //parameters.doubles("scaling_factor") = parameters.doubles("exp_vol")/(Xmax*Ymax*Zmax); // compute scaling factor as the number of our simulation domains to represent the full experimental volume
    parameters.ints("number_of_cells") = parameters.doubles("exp_cells")/parameters.doubles("scaling_factor")+0.5; // scale the experimental cell count using scaling factor, rounds to nearest int (truncates the decimals, so +0.5 to round up or down)
    std::cout<< "Number of cells (before rounding): " << parameters.doubles("exp_cells")/parameters.doubles("scaling_factor") <<std::endl;
    if (parameters.ints("number_of_cells") == 0){ parameters.ints("number_of_cells") = 1; } // make sure there is at least one cell


    //parameters.ints("number_of_cells") = 5; // hardcode number of initial cells - for testing only
    //exit(-1); // stops running here - only for testing parameter values
	

	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin;

    load_cells_from_pugixml();

    // new fibre related parameters and bools
    bool isFibreFromFile = false;
    
    for( int i=0; i < (*all_cells).size(); i++ ){

        if (isFibre((*all_cells)[i]))
        {
            /* fibre positions are given by csv
               assign fibre orientation and test whether out of bounds */
            isFibreFromFile = true;
			static_cast<PhysiMeSS_Fibre*>((*all_cells)[i])->assign_fibre_orientation();    
            //static_cast<PhysiMeSS_Fibre*>((*all_cells)[i])->check_out_of_bounds((*all_cells)[i]->position);
			
        } 
    }

    /* agents have not been added from the file but do want them
       create some of each agent type */

    if(!isFibreFromFile){
        Cell* pC;
        std::vector<double> position = {0, 0, 0};

        for( int k=0; k < cell_definitions_by_index.size() ; k++ ) {

            Cell_Definition *pCD = cell_definitions_by_index[k];
            // std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
            
            if (!isFibre(pCD))
            {
                for (int n = 0; n < parameters.ints("number_of_cells"); n++) {

                    position[0] = Xmin + UniformRandom() * Xrange;
                    position[1] = Ymin + UniformRandom() * Yrange;
                    position[2] = Zmin + UniformRandom() * Zrange;

                    pC = create_cell(*pCD);
                                        
                    pC->assign_position(position);
                }
            } 
            
            else 
            {
                for ( int nf = 0 ; nf < parameters.ints("number_of_fibres") ; nf++ ) {

                    position[0] = Xmin + UniformRandom() * Xrange;
                    position[1] = Ymin + UniformRandom() * Yrange;
                    position[2] = Zmin + UniformRandom() * Zrange;

                    pC = create_cell(*pCD);

                    static_cast<PhysiMeSS_Fibre*>(pC)->assign_fibre_orientation();
                    static_cast<PhysiMeSS_Fibre*>(pC)->check_out_of_bounds(position);

                    pC->assign_position(position);
                }
            }
        }
    }
    else { // ML: if fibres are loaded in by file, add in the T cells (not in the files)
        Cell* pC;
        std::vector<double> position = {0, 0, 0};
        Cell_Definition *pCD = cell_definitions_by_index[0];
        for (int n = 0; n < parameters.ints("number_of_cells"); n++) {
            position[0] = Xmin + UniformRandom() * Xrange;
            position[1] = Ymin + UniformRandom() * Yrange;
            position[2] = Zmin + UniformRandom() * Zrange;

            pC = create_cell(*pCD);    
            pC->assign_position(position);
        }
    }

    remove_physimess_out_of_bounds_fibres();
    
}


std::vector<std::string> paint_by_cell_type_and_state( Cell* pCell )
{
    
	std::vector< std::string > output( 0);
	
    int k_state = pCell->custom_data.find_variable_index("state");


    if(pCell->custom_data[k_state]<0.5)
    { return std::vector<std::string>(4, "lightsteelblue");}
    //else
    //{ return std::vector<std::string>(4, "black");}
   
    int k_time  = pCell->custom_data.find_variable_index("attached_time");

    double time_attached = pCell->custom_data[k_time];

    // Get values
    double activation_time_thresh = parameters.doubles("activation_time_thresh");
    double max_time = 3000;

    int color = (int) round( ((time_attached)/max_time) * 255 );
    //std::cout<<"time attached: "<<time_attached<<" col: "<<color<<std::endl;
	if(color > 255){
		color = 255;
	}
	char szTempString [128];
	sprintf( szTempString , "rgb(%u,0,%u)", color, 255 - color);
	output.push_back( szTempString );
	output.push_back( szTempString );
	output.push_back( szTempString );
	output.push_back( szTempString );

    return output;
		

}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ 
	return paint_by_cell_type_and_state(pCell); 
}
std::string my_coloring_function_for_substrate( double concentration, double max_conc, double min_conc )
{ return paint_by_density_percentage( concentration,  max_conc,  min_conc); }

void my_cellcount_function(char* string){
	int nb_fibres = 0;
	for (Cell* cell : *all_cells) {
		if (isFibre(cell))
			nb_fibres++;
	}

	sprintf( string , "%lu cells, %u fibres" , all_cells->size()-nb_fibres, nb_fibres ); 
}

void custom_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ 

    // adding code to bump cells back into the domain if they go outside
    check_cell_position_bump_back(pCell);

    static int rod_type = get_cell_definition("ecm").type;

    //std::cout<<pCell->type<<std::endl;

    static int t_type   = get_cell_definition("cell").type; 
	int k_state = pCell->custom_data.find_variable_index("state");
	double THRESH_MIN = parameters.doubles("activation_time_thresh");
 
    if( pCell->type == t_type && pCell->custom_data[k_state] > 0.5) // if cell is already active
    {
		cell_proliferation_based_on_IL2(pCell, phenotype, dt);
		//check_cell_contact(pCell,phenotype,dt);
	}
	else if(pCell->type == t_type && pCell->custom_data["attached_time"]>THRESH_MIN )
    {
        pCell->custom_data[k_state] = 1.0;

		//#pragma omp critical
		//std::cout << "[Tcell ACTIVATED]" << " id =" << pCell->ID << " total_time =" << pCell->custom_data["attached_time"] << " t =" << PhysiCell_globals.current_time << std::endl;
    }
    //if(pCell->type == t_type) // cell is a T cell but not active
	//{
	//	check_cell_contact(pCell,phenotype,dt);
	//	check_for_activation(pCell,phenotype,dt);
	//}
		
	return; 
}

void check_cell_position_bump_back(Cell* pCell)
{
    
    double Xmin = microenvironment.mesh.bounding_box[0]; 
    double Xmax = microenvironment.mesh.bounding_box[3]; 
    double Ymin = microenvironment.mesh.bounding_box[1]; 
    double Ymax = microenvironment.mesh.bounding_box[4];

    /*
    std::cout<<pCell->position<<std::endl;

    if(pCell->position[0]<Xmin || pCell->position[0]>Xmax
     || pCell->position[1]>Ymax || pCell->position[1]<Ymin)
    {  
        std::cout<<"cell outside domain: "<<pCell->position<<std::endl;
        
        if(pCell->position[0]<Xmin)
        {pCell->position[0] = Xmin+1;}
        if(pCell->position[0]>Xmax)
        {pCell->position[0] = Xmax-1;}
        if(pCell->position[1]<Ymin)
        {pCell->position[1] = Ymin+1;}
        if(pCell->position[1]>Ymax)
        {pCell->position[1] = Ymax-1;}

    }   
    */

    return;
}


void cell_proliferation_based_on_IL2( Cell* pCell , Phenotype& phenotype, double dt )
{

	static int cycle_start_index = live.find_phase_index(PhysiCell_constants::live);
	static int cycle_end_index   = live.find_phase_index(PhysiCell_constants::live);

	static int IL2_index = microenvironment.find_density_index("nutrient");
	double IL2 = pCell->nearest_density_vector()[IL2_index];

	double rPmax = parameters.doubles("rPmax");
	double IP    = parameters.doubles("IP")/(microenvironment.mesh.bounding_box[5]-microenvironment.mesh.bounding_box[2]); // divide by width in z direction to get volumetric density

	phenotype.cycle.data.transition_rate( cycle_start_index, cycle_end_index ) = rPmax*IL2/(IP+IL2);

    //std::cout<<"IL2 conc: "<<IL2<<std::endl;
    //std::cout<<"I_P: "<<IP<<std::endl;
    //std::cout<<"rpmax scale: "<<IL2/(IP+IL2)<<std::endl;
	
	return;
}


void check_cell_contact( Cell* pCell , Phenotype& phenotype, double dt )
{

    static int rod_type = get_cell_definition("ecm").type;

    int k_touch = pCell->custom_data.find_variable_index("touching_rod");
    int k_time  = pCell->custom_data.find_variable_index("attached_time");
    int k_state = pCell->custom_data.find_variable_index("state");

    if( k_touch < 0 || k_time < 0 || k_state < 0 ) return;

    pCell->custom_data[k_touch] = 0.0;

    auto nearby = pCell->nearby_interacting_cells();

 
    //set migration speed at 1 to start
    
    // AJ ADDED 
    //pCell->phenotype.motility.migration_speed = 1;

    for( Cell* other : nearby )
    {
        if( other->type != rod_type ) continue;

        double dx = pCell->position[0] - other->position[0];
        double dy = pCell->position[1] - other->position[1];
        double dz = pCell->position[2] - other->position[2];
        double d  = std::sqrt(dx*dx + dy*dy + dz*dz);

        double contact_dist = pCell->phenotype.geometry.radius + other->custom_data["fibre_radius"];

        if( d <= contact_dist )
        {
            pCell->custom_data[k_touch] = 1.0;
        	pCell->custom_data[k_time] += dt;

            // add that if the cells are in contact, the speed of the T cells slows to almost nothing
            //pCell->phenotype.motility.migration_speed = 0.05;

            //if (pCell->custom_data[k_state]>0.5) {  // AJ ADDED - if cell is activated and touching other cells, then slow it down
            //    pCell->phenotype.motility.migration_speed = 0.1;
            //}

            break;
        }
      
    }
	
	return;
}

void check_for_activation( Cell* pCell , Phenotype& phenotype, double dt )
{
	

    int k_time  = pCell->custom_data.find_variable_index("attached_time");
    int k_state = pCell->custom_data.find_variable_index("state");
	double THRESH_MIN = parameters.doubles("activation_time_thresh");

	if( pCell->custom_data[k_time] >= THRESH_MIN )
	{
		pCell->custom_data[k_state] = 1.0;

		#pragma omp critical
		std::cout
		<< "[Tcell ACTIVATED]"
		<< " id =" << pCell->ID
		<< " total_time =" << pCell->custom_data[k_time]
		<< " t =" << PhysiCell_globals.current_time
		<< std::endl;
	}
    
	return;
}
void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

void fibre_time_secretion_function( Cell* pCell, Phenotype& phenotype, double dt )
{
    static int substrate_index = microenvironment.find_density_index( "nutrient" ); // your substrate name

    double t = PhysiCell_globals.current_time;

    double v   = parameters.doubles("v"); 
    double q   = parameters.doubles("q");

    static double rod_length = pCell->custom_data["f_length"]; 
    static double rod_radius = pCell->custom_data["fibre_radius"];
    static double rod_mass = pCell->custom_data["rod_mass"]; 
    double rho = rod_mass/(rod_length*3.14159265358979*rod_radius*rod_radius); // mass density of a single rod (ug/um^3) 

    phenotype.secretion.secretion_rates[substrate_index] = v*q*rho*exp(-v*t);
	if (phenotype.secretion.secretion_rates[substrate_index]<0)
	{std::cout<<"error in secretion rate"<<std::endl;}


    return;
}

void tcell_division_function( Cell* pParent, Cell* pDaughter )
{


    static int t_type   = get_cell_definition("cell").type; 
	
    if( pParent == nullptr || pDaughter == nullptr ) return;

    if( pParent->type != t_type ) return;

    int k_state = pParent->custom_data.find_variable_index("state");
    int k_time  = pParent->custom_data.find_variable_index("attached_time");
    int k_touch = pParent->custom_data.find_variable_index("touching_rod");

    if( k_state < 0 )
	{
		std::cout<<"k_state index zero?"<<std::endl;
		return;
	}

    if( pParent->custom_data[k_state] > 0.5 )
    {
        pDaughter->custom_data[k_state] = 0.0;

        if( pDaughter->custom_data[k_time]  >= 0 ) pDaughter->custom_data[k_time]  = 0.0;
        if( pDaughter->custom_data[k_touch] >= 0 ) pDaughter->custom_data[k_touch] = 0.0;

        //#pragma omp critical
        //std::cerr << "[DIV] parent active -> daughter naive | parent id=" << pParent->ID << " daughter id=" << pDaughter->ID << " t=" << PhysiCell_globals.current_time << "\n";
    }

   // clamp_cell_to_domain(pParent);
   // clamp_cell_to_domain(pDaughter);
	return;
}

Cell* instantiate_physimess_cell() { return new PhysiMeSS_Cell; }
Cell* instantiate_physimess_fibre() { return new PhysiMeSS_Fibre; }
Cell* instantiate_physimess_cell_custom_degrade() { return new PhysiMeSS_Cell_Custom_Degrade; }


void PhysiMeSS_Cell_Custom_Degrade::degrade_fibre(PhysiMeSS_Fibre* pFibre)
{
	// Here this version of the degrade function takes cell pressure into account in the degradation rate
    double distance = 0.0;
    pFibre->nearest_point_on_fibre(position, displacement);
    for (int index = 0; index < 3; index++) {
        distance += displacement[index] * displacement[index];
    }
    distance = std::max(sqrt(distance), 0.00001);
    
    
        // Fibre degradation by cell - switched on by flag fibre_degradation
        double stuck_threshold = this->custom_data["fibre_stuck_time"];
        double pressure_threshold = this->custom_data["fibre_pressure_threshold"];
        if (this->custom_data["fibre_degradation"] > 0.5 && (stuck_counter >= stuck_threshold
                                                        || state.simple_pressure > pressure_threshold)) {
            // if (stuck_counter >= stuck_threshold){
            //     std::cout << "Cell " << ID << " is stuck at time " << PhysiCell::PhysiCell_globals.current_time
            //                 << " near fibre " << pFibre->ID  << std::endl;;
            // }
            // if (state.simple_pressure > pressure_threshold){
            //     std::cout << "Cell " << ID << " is under pressure of " << state.simple_pressure << " at "
            //                 << PhysiCell::PhysiCell_globals.current_time << " near fibre " << pFibre->ID  << std::endl;;
            // }
            displacement *= -1.0/distance;
            double dotproduct = dot_product(displacement, phenotype.motility.motility_vector);
            if (dotproduct >= 0) {
                double rand_degradation = PhysiCell::UniformRandom();
                double prob_degradation = this->custom_data["fibre_degradation_rate"];
                if (state.simple_pressure > pressure_threshold){
                    prob_degradation *= state.simple_pressure;
                }
                if (rand_degradation <= prob_degradation) {
                    //std::cout << " --------> fibre " << (*other_agent).ID << " is flagged for degradation " << std::endl;
                    // (*other_agent).parameters.degradation_flag = true;
                    pFibre->flag_for_removal();
                    // std::cout << "Degrading fibre agent " << pFibre->ID << " using flag for removal !!" << std::endl;
                    stuck_counter = 0;
                }
            }
        }
    // }
}

// just copying the default "distance_to_domain_edge" standardfunction, but you can modify it to do whatever you want.
double custom_distance_to_domain(Cell* pCell, Phenotype& phenotype, double dummy)
{
	static double tolerance = 1e-7;
	static double one_over_sqrt_2 = 0.70710678118;
	static double one_over_sqrt_3 = 0.57735026919; 
	
		
	double min_distance = 9e99; 
	int nearest_boundary = -1; 
	
	// check against xL and xU
	double temp_distance = pCell->position[0] - microenvironment.mesh.bounding_box[0]; 
	if( temp_distance < min_distance )
	{
		min_distance = temp_distance; 
		nearest_boundary = 0; 
	}
	temp_distance = microenvironment.mesh.bounding_box[3] - pCell->position[0]; 
	if( temp_distance < min_distance )
	{
		min_distance = temp_distance; 
		nearest_boundary = 1; 
	}
	
	// check against yL and yU
	temp_distance = pCell->position[1] - microenvironment.mesh.bounding_box[1]; 
	if( temp_distance < min_distance )
	{
		min_distance = temp_distance; 
		nearest_boundary = 2; 
	}
	temp_distance = microenvironment.mesh.bounding_box[4] - pCell->position[1]; 
	if( temp_distance < min_distance )
	{
		min_distance = temp_distance; 
		nearest_boundary = 3; 
	}	
	
	if( default_microenvironment_options.simulate_2D == false )
	{
		// if in 3D, check against zL and zU
		temp_distance = pCell->position[2] - microenvironment.mesh.bounding_box[2]; 
		if( temp_distance < min_distance )
		{
			min_distance = temp_distance; 
			nearest_boundary = 4; 
		}
		temp_distance = microenvironment.mesh.bounding_box[5] - pCell->position[2]; 
		if( temp_distance < min_distance )
		{
			min_distance = temp_distance; 
			nearest_boundary = 5; 
		}			
		
		// check for 3D exceptions 
		
		// lines 
		if( fabs( (pCell->position[0]) - (pCell->position[1]) ) < tolerance && 
			fabs( (pCell->position[1]) - (pCell->position[2]) ) < tolerance && 
			fabs( (pCell->position[0]) - (pCell->position[2]) ) < tolerance )
		{
			if( pCell->position[0] > 0 )
			{
				if( pCell->position[0] > 0 && pCell->position[1] > 0 )
				{ pCell->displacement = { -one_over_sqrt_3 , -one_over_sqrt_3 , -one_over_sqrt_3 }; }
				if( pCell->position[0] < 0 && pCell->position[1] > 0 )
				{ pCell->displacement = { one_over_sqrt_3 , -one_over_sqrt_3 , -one_over_sqrt_3 }; }
				
				if( pCell->position[0] > 0 && pCell->position[1] < 0 )
				{ pCell->displacement = { -one_over_sqrt_3 , one_over_sqrt_3 , -one_over_sqrt_3 }; }
				if( pCell->position[0] < 0 && pCell->position[1] < 0 )
				{ pCell->displacement = { one_over_sqrt_3 , one_over_sqrt_3 , -one_over_sqrt_3 }; }
			} 
			else
			{
				if( pCell->position[0] > 0 && pCell->position[1] > 0 )
				{ pCell->displacement = { -one_over_sqrt_3 , -one_over_sqrt_3 , one_over_sqrt_3 }; }
				if( pCell->position[0] < 0 && pCell->position[1] > 0 )
				{ pCell->displacement = { one_over_sqrt_3 , -one_over_sqrt_3 , one_over_sqrt_3 }; }
				
				if( pCell->position[0] > 0 && pCell->position[1] < 0 )
				{ pCell->displacement = { -one_over_sqrt_3 , one_over_sqrt_3 , one_over_sqrt_3 }; }
				if( pCell->position[0] < 0 && pCell->position[1] < 0 )
				{ pCell->displacement = { one_over_sqrt_3 , one_over_sqrt_3 , one_over_sqrt_3 }; }				
			}
			return min_distance; 
		}
		
		// planes - let's not worry for today 
		
	}
	else
	{
		// check for 2D  exceptions 
		
		if( fabs( (pCell->position[0]) - (pCell->position[1]) ) < tolerance )
		{
			if( pCell->position[0] > 0 && pCell->position[1] > 0 )
			{ pCell->displacement = { -one_over_sqrt_2 , -one_over_sqrt_2 , 0 }; }
			if( pCell->position[0] < 0 && pCell->position[1] > 0 )
			{ pCell->displacement = { one_over_sqrt_2 , -one_over_sqrt_2 , 0 }; }
			
			if( pCell->position[0] > 0 && pCell->position[1] < 0 )
			{ pCell->displacement = { -one_over_sqrt_2 , one_over_sqrt_2 , 0 }; }
			if( pCell->position[0] < 0 && pCell->position[1] < 0 )
			{ pCell->displacement = { one_over_sqrt_2 , one_over_sqrt_2 , 0 }; }
			return min_distance; 
		}
	}
	
	//std::cout << "nearest_boundary = " << nearest_boundary << std::endl;	 // dummy debug print

	// no exceptions 
	switch(nearest_boundary)
	{
		case 0:
			pCell->displacement = {1,0,0}; 
			return min_distance; 
		case 1:
			pCell->displacement = {-1,0,0}; 
			return min_distance;
		case 2:
			pCell->displacement = {0,1,0}; 
			return min_distance; 
		case 3: 
			pCell->displacement = {0,-1,0}; 
			return min_distance; 
		case 4: 
			pCell->displacement = {0,0,1}; 
			return min_distance; 
		case 5: 
			pCell->displacement = {0,0,-1}; 
			return min_distance; 
		default:
			pCell->displacement = {0,0,0};
			return 9e99; 
	}
	
	pCell->displacement = {0,0,0};
	return 9e99; 
}	
