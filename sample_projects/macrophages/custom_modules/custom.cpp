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
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
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
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <string> 
#include <unistd.h>
#include <vector>
#include <sstream>

#include "./custom.h"

// declare cell definitions here 

Cell_Definition macro0; 
Cell_Definition macro1;
Cell_Definition macro2;
Cell_Definition neutro;
//Cell_Definition MSC;

using namespace std;

vector<vector<float>> TNFamech;
vector<vector<float>> IL10mech;

bool mech_do=false;

int voxelRegion;

double ecf=1e12; //exponential correction factor due to small concentration values

double xmin=-1000;
double xmax=1000;
double ymin=-1000;
double ymax=1000;
double zmin=-1000;
double zmax=1000;
	

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 

	// housekeeping 
	
	initialize_default_cell_definition();
	//create_standard_cell_death_models();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

	// Name the default cell type 
	
	cell_defaults.type = 0; 
	cell_defaults.name = "cell"; 

	int debris_index = microenvironment.find_density_index( "debris" );
	int TNF_index = microenvironment.find_density_index( "TNFa" );
	int TGF_index = microenvironment.find_density_index( "TGFb" );
	int IL10_index = microenvironment.find_density_index( "IL10" );
	int IFN_index = microenvironment.find_density_index( "IFNg" );

	// set default cell cycle model 
	
	cell_defaults.functions.cycle_model = live; 
	cell_defaults.phenotype.cycle.data.transition_rate(0,0) = 0.0; //parameters.doubles("macro_proliferation_ratio"); //proliferation ratio
	cell_defaults.phenotype.volume.fluid_change_rate=0.0;
	//cell_defaults.phenotype.molecular.internalized_total_substrates[debris_index]=1.0;

	// set default_cell_functions; 

	//cell_defaults.functions.update_phenotype = default_function; 
	
	//cell_defaults.functions.update_phenotype = macrophage0_function;
	//cell_defaults.functions.update_migration_bias = macrophage_chemotaxis; 

	// needed for a 2-D simulation: 
	
	// grab code from heterogeneity 
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = false; 

	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment );
	//cell_defaults.phenotype.sync_to_functions(macro0.functions ); 
	//cell_deafults.phenotype.secretion.secretion_rates[...] =0;

	// set the rate terms in the default phenotype 

	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int napoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Neutrophil Apoptosis" );

	//int G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
	//int S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );

	// initially no death 
	cell_defaults.phenotype.death.rates[apoptosis_model_index] = 0.0; 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 
	cell_defaults.phenotype.death.rates[napoptosis_model_index] = 0.0; 

	// initially no birth 
	//cell_defaults.phenotype.cycle.data.transition_rate(G0G1_index, S_index ) = 0.0 ; 

	// not motile 
	cell_defaults.phenotype.motility.is_motile = true;

	// add custom data here, if any 
	
	// release virus at death or when eaten 
	cell_defaults.phenotype.molecular.fraction_released_at_death[debris_index]= 1.0; //100% of debris released into microenvironment

	// Now, let's define another cell type. 
	// It's best to just copy the default and modify it.
	
	//M0 PROPERTIES
	macro0 = cell_defaults; 
	macro0.type = 1; 
	macro0.name = "macrophage0";

	macro0.functions.update_phenotype = macrophage0_function; 
	macro0.phenotype.sync_to_functions( cell_defaults.functions );
	double mvol=parameters.doubles("macro_volume"); //um^3
	macro0.phenotype.volume.multiply_by_ratio(mvol/macro0.phenotype.volume.total); 
	//macro0.parameters.pReference_live_phenotype = &( macro0.phenotype ); 

	macro0.functions.cycle_model = live; 
	macro0.phenotype.cycle.data.transition_rate(0,0) = parameters.doubles("k_P0"); // min^-1
	macro0.phenotype.volume.fluid_change_rate=0.0;
	macro0.phenotype.molecular.fraction_released_at_death[debris_index]= 1.0; //100% of debris released into microenvironment

	macro0.functions.update_migration_bias = macrophage_chemotaxis; 

	// enable random motility 
	macro0.phenotype.motility.is_motile = true; 
	macro0.phenotype.motility.persistence_time = 1.0; //parameters.doubles( "macrophage_persistence_time" ); // min
	macro0.phenotype.motility.migration_speed = parameters.doubles( "m_ms" ); // um min^-1
	macro0.phenotype.motility.migration_bias = parameters.doubles( "m_mb" ); // -
	
	macro0.phenotype.death.rates[apoptosis_model_index] = parameters.doubles("k_A0"); // min^-1

	// Set cell-cell adhesion relative to other cells 
	//macrophage.phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "macrophage_relative_adhesion" ); 
	
	// no birth 
	//macrophage.phenotype.cycle.data.transition_rate(G0G1_index, S_index ) = 0.05;
	
	// SECRETION + UPTAKE
	//macro0.phenotype.secretion.uptake_rates[debris_index] = 0.0;
	macro0.phenotype.secretion.secretion_rates[TNF_index] = parameters.doubles("k_T0"); // ng min^-1
	macro0.phenotype.secretion.saturation_densities[TNF_index] = 1.00; // ng um^-3
	macro0.phenotype.secretion.secretion_rates[TGF_index] = parameters.doubles("k_TG0"); // ng min^-1
	macro0.phenotype.secretion.saturation_densities[TGF_index] = 1.00; // ng um^-3
	macro0.phenotype.secretion.secretion_rates[IL10_index] = parameters.doubles("k_I0"); // ng min^-1
	macro0.phenotype.secretion.saturation_densities[IL10_index] = 1.00; // ng um^-3

	macro0.phenotype.secretion.secretion_rates[IFN_index] = parameters.doubles("k_IF0"); // ng min^-1
	macro0.phenotype.secretion.saturation_densities[IFN_index] = 1.00; // ng um^-3

	// MACROPHAGES M1

	macro1 = macro0; 
	macro1.type = 2; 
	macro1.name = "macrophage1";

	macro1.functions.update_phenotype = macrophage1_function; 
	macro1.phenotype.sync_to_functions( cell_defaults.functions ); 
	
	macro1.parameters.pReference_live_phenotype = &( macro1.phenotype ); 
	macro1.phenotype.secretion.secretion_rates[TNF_index] = parameters.doubles("k_T1"); // ng min^-1
	macro1.phenotype.secretion.secretion_rates[TGF_index] = parameters.doubles("k_TG1"); // ng min^-1
	macro1.phenotype.secretion.secretion_rates[IL10_index] = parameters.doubles("k_I1"); // ng min^-1
	macro1.phenotype.secretion.secretion_rates[IFN_index] = parameters.doubles("k_IF1"); // ng min^-1

	// MACROPHAGES M2

	macro2 = macro0; 
	macro2.type = 3; 
	macro2.name = "macrophage2";

	macro2.functions.update_phenotype = macrophage2_function; 
	macro2.phenotype.sync_to_functions( cell_defaults.functions ); 
	
	macro2.parameters.pReference_live_phenotype = &( macro2.phenotype ); 
	macro2.phenotype.secretion.secretion_rates[TNF_index] = parameters.doubles("k_T2"); // ng min^-1
	macro2.phenotype.secretion.secretion_rates[TGF_index] = parameters.doubles("k_TG2"); // ng min^-1
	macro2.phenotype.secretion.secretion_rates[IL10_index] = parameters.doubles("k_I2"); // ng min^-1

	//cell mechanics

 	//macro0.phenotype.mechanics.cell_cell_adhesion_strength*=0.0;
	//macro1.phenotype.mechanics.cell_cell_adhesion_strength*=0.0;
	//macro2.phenotype.mechanics.cell_cell_adhesion_strength*=0.0; 

 	//macro0.phenotype.mechanics.cell_cell_repulsion_strength*=3.0;
	//macro1.phenotype.mechanics.cell_cell_repulsion_strength*=3.0;
	//macro2.phenotype.mechanics.cell_cell_repulsion_strength*=3.0; 
	

	// NEUTROPHILS PMN
	neutro = cell_defaults; 
	neutro.type = 4; 
	neutro.name = "neutrophils";

	neutro.functions.update_phenotype = neutrophils_function;
	neutro.phenotype.sync_to_functions( cell_defaults.functions );
	double nvol=parameters.doubles("neutro_volume"); // um^3
	neutro.phenotype.volume.multiply_by_ratio(nvol/neutro.phenotype.volume.total); 

	neutro.functions.cycle_model = live; 
	neutro.phenotype.cycle.data.transition_rate(0,0) = 0.0; // min^-1
	neutro.phenotype.volume.fluid_change_rate=0.0;
	neutro.phenotype.molecular.internalized_total_substrates[debris_index]=1.6e-3;
	neutro.phenotype.molecular.fraction_released_at_death[debris_index]= 1.0; //100% of debris released into microenvironment

	neutro.functions.update_migration_bias = neutrophils_chemotaxis;

	// enable random motility 
	neutro.phenotype.motility.is_motile = true; 
	neutro.phenotype.motility.persistence_time = 1.0; // min;
	neutro.phenotype.motility.migration_speed = parameters.doubles( "n_ms" ); // um min^-1
	neutro.phenotype.motility.migration_bias = parameters.doubles( "n_mb" ); // -
	
	neutro.phenotype.death.rates[apoptosis_model_index] = 0.0;
	neutro.phenotype.death.rates[napoptosis_model_index] = parameters.doubles("k_AN"); 

	// Set cell-cell adhesion relative to other cells 
	//macrophage.phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "macrophage_relative_adhesion" ); 
	
	// no birth 
	//macrophage.phenotype.cycle.data.transition_rate(G0G1_index, S_index ) = 0.05;
	
	// SECRETION + UPTAKE
	neutro.phenotype.secretion.secretion_rates[TNF_index] = 1.00*parameters.doubles("k_TN"); // ng um^-3 min^-1
	neutro.phenotype.secretion.saturation_densities[TNF_index] = 1.00; // ng um^-3
	neutro.phenotype.secretion.secretion_rates[IFN_index] = 1.00*parameters.doubles("k_IFN"); // ng um^-3 min^-1
	neutro.phenotype.secretion.saturation_densities[IFN_index] = 1.00; // ng um^-3

	build_cell_definitions_maps(); 

	// MECHANICS 

	if (mech_do){
	string fname="TNFa_mech.csv";
 
	vector<float> row;
	string line, word;
 
	fstream file;
	file.open(fname);
	if(file.is_open())
	{
		while(getline(file, line))
		{
			//std::cout << "GO" << std::endl; 
			row.clear();
 
			stringstream str(line);
 
			while(getline(str, word, ',')){
				//std::cout << word << std::endl; 
				row.push_back(stod(word.c_str()));
			}
			TNFamech.push_back(row);
		}
	}
	else
		std::cout<<"Could not open the file\n";

	fname="IL10_mech.csv";

	fstream file2;
	file2.open(fname);
	if(file2.is_open())
	{
		while(getline(file2, line))
		{
			//std::cout << "GO" << std::endl; 
			row.clear();
 
			stringstream str(line);
 
			while(getline(str, word, ',')){
				//std::cout << word << std::endl; 
				row.push_back(stod(word.c_str()));
			}
			IL10mech.push_back(row);
		}
	}
	else
		std::cout<<"Could not open the file\n";
	
	std::cout<<"---- MECHANOREGULATION INCLUDED ----"<<std::endl;	
	}

}

void setup_microenvironment( void )
{

	pugi::xml_node node = xml_find_node( physicell_config_root , "geometry" );

	double dA=xml_get_double_value(node, "dA");
	double dB=xml_get_double_value(node, "dB");
	double dC=xml_get_double_value(node, "dC");
	double dD=xml_get_double_value(node, "dD");
	double dE=xml_get_double_value(node, "dE");
	double dF=xml_get_double_value(node, "dF");

	node = node.parent();
	
	node = xml_find_node( physicell_config_root , "domain" );
	
	double dx = xml_get_double_value( node, "dx" ); 
	double dy = xml_get_double_value( node, "dy" ); 
	double dz = xml_get_double_value( node, "dz" ); 

	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == true )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = false; 
	}

	default_microenvironment_options.calculate_gradients = true; 
	
	//set Dirichlet conditions 

	default_microenvironment_options.outer_Dirichlet_conditions = true;
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector( 5 , 0.0); // ng um^-3
	std::vector<bool> bc_activ(5, true); // USER
	std::vector<bool> bc_inactiv(5, false); // USER
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
	
	// set initial conditions 
	default_microenvironment_options.initial_condition_vector = bc_vector; 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment();

	//microenvironment.name = "tissue"; 
	
	int debris_ID = microenvironment.find_density_index( "debris" ); 

	// set a new initial condition for debris

	double xDebris;
	double yDebris;
	double zDebris;

	double debrisconc = 1.00e-4; // ng um^-3

	for( int n=0; n < microenvironment.number_of_voxels(); n++ )
	{
		// x coordinate of the nth voxel's center
		xDebris = microenvironment.mesh.voxels[n].center[0];
		yDebris = microenvironment.mesh.voxels[n].center[1];
		zDebris = microenvironment.mesh.voxels[n].center[2];

		if (((xDebris<dB)&&(xDebris>-dB)&&(dC>sqrt((yDebris*yDebris)+(zDebris*zDebris))))||((xDebris<dD)&&(xDebris>-dD)&&((dC+dE)>sqrt((yDebris*yDebris)+(zDebris*zDebris))))||(((dC+dE)<sqrt((yDebris*yDebris)+(zDebris*zDebris)))&&((dA-(dA-dE-dC)*(((xDebris)*(xDebris))/((dD+dF)*(dD+dF))))>sqrt((yDebris*yDebris)+(zDebris*zDebris)))))
		{
			microenvironment(n)[debris_ID] = debrisconc*dx*dy*dz; // ng
			voxelRegion++;
		}
	}

	// display the microenvironment again 
	
	microenvironment.display_information( std::cout ); 
	
	return; 
}

void release_new_macrophages(void)
{

	pugi::xml_node node = xml_find_node( physicell_config_root , "geometry" );

	double dA=xml_get_double_value(node, "dA");
	double dB=xml_get_double_value(node, "dB");
	double dC=xml_get_double_value(node, "dC");
	double dD=xml_get_double_value(node, "dD");
	double dE=xml_get_double_value(node, "dE");
	double dF=xml_get_double_value(node, "dF");

	//std::cout << voxelRegion << std::endl;

	node = node.parent();
	
	node = xml_find_node( physicell_config_root , "domain" );
	double dx = xml_get_double_value( node, "dx" ); 
	double dy = xml_get_double_value( node, "dy" ); 
	double dz = xml_get_double_value( node, "dz" ); 

	std::cout << "--> Macrophages releasing" << std::endl;
	
	Cell* pC;

	double volume_total= voxelRegion*dx*dy*dz; // um^3 
	//std::cout << " volume total: " << volume_total << " um3"<<std::endl;

	int cell_total=0;
 	for(int i=1; i<(*all_cells).size();i++)
	{
		if(((*all_cells)[i]->type==1)||((*all_cells)[i]->type==2)||((*all_cells)[i]->type==3))
		{
			cell_total++; //cell
		}
	}
	double macro_concentration=cell_total/volume_total; // cell um^-3
	
	int debris_index = microenvironment.find_density_index( "debris" );
	double debris_particles=0.00; // ng
	for( int n=0; n < microenvironment.number_of_voxels(); n++ )
	{
		debris_particles+=microenvironment(n)[debris_index]; // ng
	}
	double debris_density=debris_particles/volume_total; // ng um^-3 

	double Mmax=parameters.doubles("M0_max"); // cell um^-3
	double kmax=parameters.doubles("k_Rmax"); // min^-1  
	double RM=kmax*(1-(macro_concentration/Mmax))*debris_density*volume_total; // [cell min^-1]

	double fraction_marrow=parameters.doubles("frac_marrow"); // -

	//std::cout<< "RMacro = " << RM << std::endl;

	int perimRelease_tissues;
	int perimRelease_marrow;
	int xPosition;
	int yPosition;
	int zPosition;

	int released_M0=ceil(RM); // [cell min^-1]
	released_M0*=60.0; // [cell h^-1]

	for(int z=0;z<released_M0;z++)
	{
		perimRelease_tissues=dD+dF-1.0; //um
		perimRelease_marrow=dC-1.0; //um
		int zone=rand()%8;
		int turn=rand()%2;
		int side=rand()%100;
		int ytop, ztop;
		if(side < (1-fraction_marrow)*100)
		{
			pC=create_cell(macro0);
			xPosition=rand()%perimRelease_tissues+1;
			if(turn==0){
				ytop=(dA-1.0)-(dA-dC-dE)*((xPosition)*(xPosition)/((dD+dF-1.0)*(dD+dF-1.0)));
				yPosition=rand()%ytop+1;
				zPosition=sqrt((ytop*ytop)-(yPosition*yPosition));
			}
			else{
				ztop=(dA-1.0)-(dA-dC-dE)*((xPosition)*(xPosition)/((dD+dF-1.0)*(dD+dF-1.0)));
				zPosition=rand()%ztop+1;
				yPosition=sqrt((ztop*ztop)-(zPosition*zPosition));
			}
			switch(zone){
			case 0:
				pC->assign_position(xPosition,yPosition,zPosition);
				break;
			case 1:
				pC->assign_position(-xPosition,yPosition,zPosition);
				break;
			case 2:
				pC->assign_position(-xPosition,-yPosition,zPosition);
				break;
			case 3:
				pC->assign_position(xPosition,-yPosition,zPosition);
				break;
			case 4:
				pC->assign_position(xPosition,yPosition,-zPosition);
				break;
			case 5:
				pC->assign_position(-xPosition,yPosition,-zPosition);
				break;
			case 6:
				pC->assign_position(-xPosition,-yPosition,-zPosition);
				break;
			case 7:
				pC->assign_position(xPosition,-yPosition,-zPosition);
				break;
			}
		}
		else 
		{
			pC=create_cell(macro0);
			xPosition=dB-1.0;
			yPosition=ymax+1;
			zPosition=zmax+1;
			while(perimRelease_marrow<sqrt((yPosition*yPosition)+(zPosition*zPosition))){
				yPosition=rand()%perimRelease_marrow;
				zPosition=rand()%perimRelease_marrow;
			}
			int zone=rand()%8;
			switch(zone){
			case 0:
				pC->assign_position(xPosition,yPosition,zPosition);
				break;
			case 1:
				pC->assign_position(-xPosition,yPosition,zPosition);
				break;
			case 2:
				pC->assign_position(-xPosition,-yPosition,zPosition);
				break;
			case 3:
				pC->assign_position(xPosition,-yPosition,zPosition);
				break;
			case 4:
				pC->assign_position(xPosition,yPosition,-zPosition);
				break;
			case 5:
				pC->assign_position(-xPosition,yPosition,-zPosition);
				break;
			case 6:
				pC->assign_position(-xPosition,-yPosition,-zPosition);
				break;
			case 7:
				pC->assign_position(xPosition,-yPosition,-zPosition);
				break;
			}
		}
	}

	std::cout << " DONE" << std::endl;
	std::cout << released_M0 << " macrophages released! <--" << std::endl;
	
	return; 
}

void release_new_neutrophils( void )
{
	std::cout << "--> Neutrophils releasing";

	pugi::xml_node node = xml_find_node( physicell_config_root , "geometry" );

	double dA=xml_get_double_value(node, "dA");
	double dB=xml_get_double_value(node, "dB");
	double dC=xml_get_double_value(node, "dC");
	double dD=xml_get_double_value(node, "dD");
	double dE=xml_get_double_value(node, "dE");
	double dF=xml_get_double_value(node, "dF");

	node = node.parent();
	
	node = xml_find_node( physicell_config_root , "domain" );
	double dx = xml_get_double_value( node, "dx" );
	double dy = xml_get_double_value( node, "dy" );
	double dz = xml_get_double_value( node, "dz" );

	Cell* pC;

	double volume_total=voxelRegion*dx*dy*dz; //um^3 

	int cell_total=0;
 	for(int i=1; i<(*all_cells).size();i++)
	{
		if(((*all_cells)[i]->type==4))
		{
			cell_total++; // cell
		}
	}
	double neutro_concentration=cell_total/volume_total; // cell um^-3
	
	int debris_index = microenvironment.find_density_index( "debris" );
	double debris_particles=0;
	for( int n=0; n < microenvironment.number_of_voxels(); n++ )
	{
		debris_particles+=microenvironment(n)[debris_index]; // ng
	}
	double debris_density=debris_particles/volume_total; // ng um^-3 

	double PMNmax=parameters.doubles("PMN_max"); // cell um^-3
	double kmax=parameters.doubles("k_RNmax"); // min^-1 
	double RM=kmax*(1-(neutro_concentration/PMNmax))*debris_density*volume_total; // [cell min^-1]

	int perimRelease_marrow;
	int perimPosition;

	int xPosition;
	int yPosition;
	int zPosition;

	int released_PMN=ceil(RM); // [cell min^-1]
	released_PMN*=60.0; // [cell h^-1]

	perimRelease_marrow=dC-1.0; // um

	for(int z=0;z<released_PMN;z++)
	{
		pC=create_cell(neutro);
		xPosition=dB-1.0;
		yPosition=ymax+1;
		zPosition=zmax+1;
		while(perimRelease_marrow<sqrt((yPosition*yPosition)+(zPosition*zPosition))){
			yPosition=rand()%perimRelease_marrow;
			zPosition=rand()%perimRelease_marrow;
		}
		int zone=rand()%8;
		switch(zone){
		case 0:
			pC->assign_position(xPosition,yPosition,zPosition);
			break;
		case 1:
			pC->assign_position(-xPosition,yPosition,zPosition);
			break;
		case 2:
			pC->assign_position(-xPosition,-yPosition,zPosition);
			break;
		case 3:
			pC->assign_position(xPosition,-yPosition,zPosition);
			break;
		case 4:
			pC->assign_position(xPosition,yPosition,-zPosition);
			break;
		case 5:
			pC->assign_position(-xPosition,yPosition,-zPosition);
			break;
		case 6:
			pC->assign_position(-xPosition,-yPosition,-zPosition);
			break;
		case 7:
			pC->assign_position(xPosition,-yPosition,-zPosition);
			break;
		}
	}

	std::cout << " DONE" << std::endl;
	std::cout << released_PMN << " neutrophils released! <--" << std::endl;

	return; 
}

void setup_tissue( void )
{

	std::cout << "--> Tissue Setup: " << std::endl;
	int debris_index = microenvironment.find_density_index( "debris" ); 
	// create some cells near the origin

	pugi::xml_node node = xml_find_node( physicell_config_root , "geometry" );

	double dA=xml_get_double_value(node, "dA");
	double dB=xml_get_double_value(node, "dB");
	double dC=xml_get_double_value(node, "dC");
	double dD=xml_get_double_value(node, "dD");
	double dE=xml_get_double_value(node, "dE");
	double dF=xml_get_double_value(node, "dF");

	node = node.parent();
	
	node = xml_find_node( physicell_config_root , "domain" );
	double dx = xml_get_double_value( node, "dx" );
	double dy = xml_get_double_value( node, "dy" );
	double dz = xml_get_double_value( node, "dz" );

	double volume_total=voxelRegion*dx*dy*dz; //um^3
	
	double length_x = microenvironment.mesh.bounding_box[3] - 
		microenvironment.mesh.bounding_box[0]; 
		
	double length_y = microenvironment.mesh.bounding_box[4] - 
		microenvironment.mesh.bounding_box[1]; 

	double length_z = microenvironment.mesh.bounding_box[5] - 
		microenvironment.mesh.bounding_box[2]; 


	// create some cells near the origin
	
	Cell* pC;
	int initial_M0=parameters.doubles("initial_M0"); // cell mm^-3
	int initial_PMN=parameters.doubles("initial_PMN"); // cell mm^-3

	int num_res_M0=ceil(initial_M0 * volume_total * 1.00e-9); // cell
	int num_PMN=ceil(initial_PMN * volume_total * 1.00e-9); // cell

	std::cout << "-> Macrophages initialization";
	int macrox;
	int macroy;
	int macroz;
	int lsize;

	for(int z=0;z<num_res_M0;z++)
	{
		macrox=(rand()%int(length_x))+int(microenvironment.mesh.bounding_box[0]);
		macroy=(rand()%int(length_y))+int(microenvironment.mesh.bounding_box[1]);
		macroz=(rand()%int(length_z))+int(microenvironment.mesh.bounding_box[2]);

		if (((macrox<dB)&&(macrox>-dB)&&(dC>sqrt((macroy*macroy)+(macroz*macroz))))||((macrox<dD)&&(macrox>-dD)&&((dC+dE)>sqrt((macroy*macroy)+(macroz*macroz))))||(((dC+dE)<sqrt((macroy*macroy)+(macroz*macroz)))&&((dA-(dA-dE-dC)*(((macrox)*(macrox))/((dD+dF)*(dD+dF))))>sqrt((macroy*macroy)+(macroz*macroz)))))
		{
			pC=create_cell(macro0);
			pC->assign_position(macrox,macroy,macroz);
		}
		else
		{
			z--;
		} 
	}

	std::cout << num_res_M0 << " macrophages released! <-" << std::endl;

	std::cout << "-> Neutrophils initialization";
	int neux;
	int neuy;
	int neuz;
	//int lsize;

	for(int z=0;z<num_PMN;z++)
	{
		neux=(rand()%int(length_x))+int(microenvironment.mesh.bounding_box[0]);
		neuy=(rand()%int(length_y))+int(microenvironment.mesh.bounding_box[1]);
		neuz=(rand()%int(length_z))+int(microenvironment.mesh.bounding_box[2]);

		if (((neux<dB)&&(neux>-dB)&&(dC>sqrt((neuy*neuy)+(neuz*neuz))))||((neux<dD)&&(neux>-dD)&&((dC+dE)>sqrt((neuy*neuy)+(neuz*neuz))))||(((dC+dE)<sqrt((neuy*neuy)+(neuz*neuz)))&&((dA-(dA-dE-dC)*(((neux)*(neux))/((dD+dF)*(dD+dF))))>sqrt((neuy*neuy)+(neuz*neuz)))))
		{
			pC=create_cell(neutro);
			pC->assign_position(neux,neuy,neuz);
		}
		else
		{
			z--;
		} 
	}

	std::cout << num_PMN << " neutrophils released! <-" << std::endl;
	
	std::cout << " Tissue Setup COMPLETE <--" << std::endl;
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	
	std::vector<std::string> output =  false_cell_coloring_live_dead(pCell); 

	if(pCell->type == 2)
	{
		 output[0] = "red";  
	}

	if(pCell->type == 3)
	{
		 output[0] = "blue";  
	}

	output[1]=output[0];
		
	if(pCell->type == 4) 
	{
		 output[0] = "yellow"; 
		 output[1] = "black"; 
	}
	
	return output; 
}

std::vector<std::string> immunefluorescence_function( Cell* pCell )
{
	
	std::vector<std::string> output =  false_cell_coloring_live_dead(pCell); 

	output[2]="rgb(49,89,242)";
	output[3]=output[2];

	if(pCell->type == 1)
	{
		 output[0] = "rgb(102,255,0)";  
	}

	if(pCell->type == 2)
	{
		 output[0] = "rgb(235,236,240)";  
	}

	if(pCell->type == 3)
	{
		 output[0] = "rgb(255,22,12)";  
	}
		
	if(pCell->type == 4)
	{
		 output[0] = "rgb(49,89,242)"; 
	}
	
	output[1]=output[0];
	return output; 
}


void macrophage0_function( Cell* pCell, Phenotype& phenotype, double dt )
{

	pugi::xml_node node = xml_find_node( physicell_config_root , "domain" );
	
	double dx = xml_get_double_value( node, "dx" ); 
	double dy = xml_get_double_value( node, "dy" ); 
	double dz = xml_get_double_value( node, "dz" ); 

	int debris_index = microenvironment.find_density_index("debris");	
	int TNF_index = microenvironment.find_density_index("TNFa");	
	int IL10_index = microenvironment.find_density_index("IL10");
	int TGF_index = microenvironment.find_density_index("TGFb");

	pCell->phenotype.molecular.internalized_total_substrates[debris_index]=1.0e-3*pCell->phenotype.volume.total; // ng
	
	Cell* pTestCell = NULL; 
	std::vector<Cell*> neighbors = get_possible_neighbors(pCell);

	// SURVIVAL

	int nsize=0;

	for( int n=0; n < neighbors.size() ; n++ )
	{
		pTestCell = neighbors[n]; 
		std::vector<double> displacement = pTestCell->position;
		displacement -= pCell->position;
		double distance = norm( displacement ); 
		
		double max_distance = pCell->phenotype.geometry.radius + pTestCell->phenotype.geometry.radius; 
		max_distance *= 1.1; 

		if(distance < max_distance)
			{
				nsize++; 
			} 
	}
	
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
 	
	pCell->phenotype.death.rates[apoptosis_model_index] = (nsize+1) * parameters.doubles("k_A0"); // min^-1
	pCell->phenotype.cycle.data.transition_rate(0,0)=parameters.doubles("k_P0")/(nsize+1);  // min^-1

	std::vector<double> CellPosition = pCell->position;

	bool FoundPosition=false;
	int PositionFinal=0;

	double xVoxel;
	double yVoxel;
	double zVoxel;

	for(int n=0; n< microenvironment.number_of_voxels() && FoundPosition==false; n++ )
	{
		xVoxel=microenvironment.mesh.voxels[n].center[0];
		yVoxel=microenvironment.mesh.voxels[n].center[1];
		zVoxel=microenvironment.mesh.voxels[n].center[2];
		
		if(abs(CellPosition[0]-xVoxel)<=microenvironment.mesh.dx/2 && abs(CellPosition[1]-yVoxel)<=microenvironment.mesh.dy/2 && abs(CellPosition[2]-zVoxel)<=microenvironment.mesh.dz/2)
		{
			FoundPosition=true;
			PositionFinal=n;
		}
	}

	//UPDATE DEBRIS UPTAKE

	double ke1=parameters.doubles("k_e1"); // min^-1 
	double aed=parameters.doubles("a_ed"); // ng um^-3
	double RD=(1.00+exp(1.00)); //-
	RD/=(1.00+exp(((microenvironment(PositionFinal)[debris_index]/(dx*dy*dz))+aed)/(microenvironment(PositionFinal)[debris_index]/(dx*dy*dz)))); // - 

	pCell->phenotype.secretion.uptake_rates[debris_index] = ke1*RD; // min^-1

	// PROLIFERATION

	int IFN_index = microenvironment.find_density_index( "IFNg" );
	pCell->phenotype.cycle.data.transition_rate(0,0)*=(parameters.doubles("k_0")+((1-parameters.doubles("k_0"))/(1.00+exp((microenvironment(PositionFinal)[IFN_index]/(dx*dy*dz))*log(3)/parameters.doubles("a_0")))));; //min^-1

	// CYTOKINE SECRETION

	float TNFfact=1.0;
	float IL10fact=1.0;

	if (mech_do){
		int indexmin=0;
		float distmin=100.0;
		float dist;
		float xCell=CellPosition[1]/1000.0;
		float yCell=CellPosition[2]/1000.0;
		float zCell=CellPosition[0]/1000.0;

		if ((abs(xCell)<0.90)&&(abs(yCell)<0.5)&&(abs(zCell)<0.5)){
			for(int n=0; n<TNFamech.size(); n++ ){
				dist=sqrt((TNFamech[n][1]-xCell)*(TNFamech[n][1]-xCell)+(TNFamech[n][2]-yCell)*(TNFamech[n][2]-yCell)+(TNFamech[n][3]-zCell)*(TNFamech[n][3]-zCell));
				if (dist<distmin){
					indexmin=n;
					distmin=dist;
				}
			} 
			TNFfact=1.00*TNFamech[indexmin][5]/112.81;
			IL10fact=1.00*IL10mech[indexmin][5]/100.0;
		}
	}

	// SECRETION RATIOS

	// TNFa
	pCell->phenotype.secretion.secretion_rates[TNF_index]=parameters.doubles("k_T0")*(TNFfact); // min^-1

	pCell->phenotype.secretion.secretion_rates[TNF_index]*=parameters.doubles("k_TIL")*exp(-parameters.doubles("a_TIL")*(ecf*microenvironment(PositionFinal)[IL10_index]/(dx*dy*dz)))+parameters.doubles("b_TIL"); // min^-1 // TNFa affected by IL10
	pCell->phenotype.secretion.secretion_rates[TNF_index]*=parameters.doubles("k_TTGF")*exp(-parameters.doubles("a_TTGF")*(ecf*microenvironment(PositionFinal)[TGF_index]/(dx*dy*dz)))+parameters.doubles("b_TTGF"); // min^-1 // TNFa affected by TGFb
	
	// CHECK // pCell->phenotype.secretion.secretion_rates[TNF_index]*=(1.00+(parameters.doubles("k_TI")/(1.00+exp(ecf*(parameters.doubles("a_TI")-microenvironment(PositionFinal)[IFN_index])))));
	
	// IFNg //check
	pCell->phenotype.secretion.secretion_rates[IFN_index]=parameters.doubles("k_IF0")*exp(-parameters.doubles("a_IF0")*microenvironment(PositionFinal)[TNF_index]/(dx*dy*dz));
	
	// IL10
	pCell->phenotype.secretion.secretion_rates[IL10_index]=parameters.doubles("k_I0")*(IL10fact); 
	
	// DIFFERENTIATE INTO M1 OR M2

 	bool polarized=false;

	double k01=parameters.doubles("k_01"); // min^-1
	double a01=parameters.doubles("a_01"); // ng um^-3
	double g1=k01*microenvironment(PositionFinal)[TNF_index]/(dx*dy*dz);
	g1/=(a01+(microenvironment(PositionFinal)[TNF_index]/(dx*dy*dz))); // min^-1

	double gtest1=(double)rand() / RAND_MAX;

	if(gtest1<=g1*dt) 
	{
		pCell->type=2;
		polarized=true;
		pCell->functions.update_phenotype = macrophage1_function; 
	}

	double k02=parameters.doubles("k_02"); // min^-1
	double a02=parameters.doubles("a_02"); // ng um^-3
	double g2=k02*microenvironment(PositionFinal)[IL10_index]/(dx*dy*dz);
	g2/=(a02+(microenvironment(PositionFinal)[IL10_index]/(dx*dy*dz))); // min^-1

	double gtest2=(double)rand() / RAND_MAX;

	if(gtest2<=g2*dt && polarized==false) 
	{
		pCell->type=3;
		pCell->functions.update_phenotype = macrophage2_function; 
	} 

	boundaries(pCell);
	return; 
}

void boundaries(Cell* pCell)
{	
	Cell* pC = NULL; 

	std::vector<double> CellPosition = pCell->position;

	int debris_index = microenvironment.find_density_index("debris");

	pugi::xml_node node = xml_find_node( physicell_config_root , "geometry" );

	double dA=xml_get_double_value(node, "dA");
	double dB=xml_get_double_value(node, "dB");
	double dC=xml_get_double_value(node, "dC");
	double dD=xml_get_double_value(node, "dD");
	double dE=xml_get_double_value(node, "dE");
	double dF=xml_get_double_value(node, "dF");
	
	return; 
}

void macrophage1_function( Cell* pCell, Phenotype& phenotype, double dt )
{

	pugi::xml_node node = xml_find_node( physicell_config_root , "domain" );
	
	double dx = xml_get_double_value( node, "dx" ); 
	double dy = xml_get_double_value( node, "dy" ); 
	double dz = xml_get_double_value( node, "dz" ); 

	int debris_index = microenvironment.find_density_index("debris"); 
	pCell->phenotype.molecular.internalized_total_substrates[debris_index]=1.0e-3*pCell->phenotype.volume.total; // ng
	
	Cell* pTestCell = NULL; 
	std::vector<Cell*> neighbors = get_possible_neighbors(pCell);

	// SURVIVAL

	int nsize=0;

	for( int n=0; n < neighbors.size() ; n++ )
	{
		pTestCell = neighbors[n]; 
		std::vector<double> displacement = pTestCell->position;
		displacement -= pCell->position;
		double distance = norm( displacement ); 
		
		double max_distance = pCell->phenotype.geometry.radius + pTestCell->phenotype.geometry.radius; 
		max_distance *= 1.1; 

		if(distance < max_distance)
			{
				nsize++; 
			} 
	}
	
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	pCell->phenotype.death.rates[apoptosis_model_index] = (nsize+1) * parameters.doubles("k_A1"); // min^-1
	pCell->phenotype.cycle.data.transition_rate(0,0)=parameters.doubles("k_P1")/(nsize+1); // min^-1
	
	std::vector<double> CellPosition = pCell->position;

	bool FoundPosition=false;
	int PositionFinal=0;

	double xVoxel;
	double yVoxel;
	double zVoxel;

	for(int n=0; n< microenvironment.number_of_voxels() && FoundPosition==false; n++ )
	{
		xVoxel=microenvironment.mesh.voxels[n].center[0];
		yVoxel=microenvironment.mesh.voxels[n].center[1];
		zVoxel=microenvironment.mesh.voxels[n].center[2];
		
		if(abs(CellPosition[0]-xVoxel)<=microenvironment.mesh.dx/2 && abs(CellPosition[1]-yVoxel)<=microenvironment.mesh.dx/2 && abs(CellPosition[2]-zVoxel)<=microenvironment.mesh.dx/2)
		{
			FoundPosition=true;
			PositionFinal=n;
		}
	}

	//UPDATE DEBRIS UPTAKE

	double ke1=parameters.doubles("k_e1"); // min^-1
	double aed=parameters.doubles("a_ed"); // ng um^-3
	double RD=(1.00+exp(1.00)); //-
	RD/=(1.00+exp(((microenvironment(PositionFinal)[debris_index]/(dx*dy*dz))+aed)/(microenvironment(PositionFinal)[debris_index]/(dx*dy*dz)))); // - 

	pCell->phenotype.secretion.uptake_rates[debris_index] = ke1*RD;

	// Cytokine secretion influence 
	
	int TNF_index = microenvironment.find_density_index("TNFa");
	int TGF_index = microenvironment.find_density_index( "TGFb" );
	int IL10_index = microenvironment.find_density_index( "IL10" ); 
	int IFN_index = microenvironment.find_density_index( "IFNg" );


	//DYNAMICS TO CHECK

	float TNFfact=1.0;
	float IL10fact=1.0;

	if (mech_do){
		int indexmin=0;
		float distmin=100.0;
		float dist;
		float xCell=CellPosition[1]/1000.0;
		float yCell=CellPosition[2]/1000.0;
		float zCell=CellPosition[0]/1000.0;

		if ((abs(xCell)<0.90)&&(abs(yCell)<0.5)&&(abs(zCell)<0.5)){
			for(int n=0; n<TNFamech.size(); n++ ){
				dist=sqrt((TNFamech[n][1]-xCell)*(TNFamech[n][1]-xCell)+(TNFamech[n][2]-yCell)*(TNFamech[n][2]-yCell)+(TNFamech[n][3]-zCell)*(TNFamech[n][3]-zCell));
				if (dist<distmin){
					indexmin=n;
					distmin=dist;
				}
			} 
			TNFfact=1.00*TNFamech[indexmin][5]/112.81;
			IL10fact=1.00*IL10mech[indexmin][5]/100.0;
		}
	}

	// SECRETION RATIOS
	
	// TNFa
	pCell->phenotype.secretion.secretion_rates[TNF_index]=parameters.doubles("k_T0")*(TNFfact); // min^-1

	pCell->phenotype.secretion.secretion_rates[TNF_index]*=parameters.doubles("k_TIL")*exp(-parameters.doubles("a_TIL")*(ecf*microenvironment(PositionFinal)[IL10_index]/(dx*dy*dz)))+parameters.doubles("b_TIL"); // min^-1 // TNFa affected by IL10
	pCell->phenotype.secretion.secretion_rates[TNF_index]*=parameters.doubles("k_TTGF")*exp(-parameters.doubles("a_TTGF")*(ecf*microenvironment(PositionFinal)[TGF_index]/(dx*dy*dz)))+parameters.doubles("b_TTGF"); // min^-1 // TNFa affected by TGFb
	// CHECK // pCell->phenotype.secretion.secretion_rates[TNF_index]*=(1.00+(parameters.doubles("k_TI")/(1.00+exp(ecf*(parameters.doubles("a_TI")-microenvironment(PositionFinal)[IFN_index])))));
	
	// IFNg // check
	pCell->phenotype.secretion.secretion_rates[IFN_index]=parameters.doubles("k_IF0")*exp(-parameters.doubles("a_IF0")*microenvironment(PositionFinal)[TNF_index]/(dx*dy*dz));
	
	// IL10
	pCell->phenotype.secretion.secretion_rates[IL10_index]=parameters.doubles("k_I0")*(IL10fact); 

	//PROLIFERATION
	pCell->phenotype.cycle.data.transition_rate(0,0)*=(parameters.doubles("k_1")+((1-parameters.doubles("k_1"))/(1.00+exp((microenvironment(PositionFinal)[IFN_index]/(dx*dy*dz))*log(3)/parameters.doubles("a_1")))));; //min^-1

	// DIFFERENTIATE INTO M2

	double k12=parameters.doubles("k_12"); // min^-1
	double a12=parameters.doubles("a_12"); // ng um^-3
	double g2=k12*microenvironment(PositionFinal)[IL10_index]/(dx*dy*dz);
	g2/=(a12+(microenvironment(PositionFinal)[IL10_index]/(dx*dy*dz))); // min^-1

	double gtest2=(double)rand() / RAND_MAX;

	if(gtest2<=g2*dt)
	{
		pCell->type=3;
		pCell->functions.update_phenotype = macrophage2_function; 
	}

	boundaries(pCell);
	return;
}

void macrophage2_function( Cell* pCell, Phenotype& phenotype, double dt )
{

	pugi::xml_node node = xml_find_node( physicell_config_root , "domain" );
	
	double dx = xml_get_double_value( node, "dx" ); 
	double dy = xml_get_double_value( node, "dy" ); 
	double dz = xml_get_double_value( node, "dz" ); 

	int debris_index = microenvironment.find_density_index("debris");	
	pCell->phenotype.molecular.internalized_total_substrates[debris_index]=1.0e-3*pCell->phenotype.volume.total; // ng
	
	Cell* pTestCell = NULL; 
	std::vector<Cell*> neighbors = get_possible_neighbors(pCell);

	// SURVIVAL

	int nsize=0;

	for( int n=0; n < neighbors.size() ; n++ )
	{
		pTestCell = neighbors[n]; 
		std::vector<double> displacement = pTestCell->position;
		displacement -= pCell->position;
		double distance = norm( displacement ); 
		
		double max_distance = pCell->phenotype.geometry.radius + pTestCell->phenotype.geometry.radius; 
		max_distance *= 1.1; 

		if(distance < max_distance)
			{
				nsize++; 
			} 
	}
	
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	pCell->phenotype.death.rates[apoptosis_model_index] = (nsize+1) * parameters.doubles("k_A2"); // min^-1
	pCell->phenotype.cycle.data.transition_rate(0,0) = parameters.doubles("k_P2")/(nsize+1); // min^-1

	std::vector<double> CellPosition = pCell->position;

	bool FoundPosition=false;
	int PositionFinal=0;

	double xVoxel;
	double yVoxel;
	double zVoxel;

	for(int n=0; n< microenvironment.number_of_voxels() && FoundPosition==false; n++ )
	{
		xVoxel=microenvironment.mesh.voxels[n].center[0];
		yVoxel=microenvironment.mesh.voxels[n].center[1];
		zVoxel=microenvironment.mesh.voxels[n].center[2];
		
		if(abs(CellPosition[0]-xVoxel)<=microenvironment.mesh.dx/2 && abs(CellPosition[1]-yVoxel)<=microenvironment.mesh.dx/2 && abs(CellPosition[2]-zVoxel)<=microenvironment.mesh.dx/2)
		{
			FoundPosition=true;
			PositionFinal=n;
		}
	}

	//UPDATE DEBRIS UPTAKE

	double ke2=parameters.doubles("k_e2"); // min^-1
	double aed=parameters.doubles("a_ed"); // ng min^-1
	double RD=(1.00+exp(1.00)); //-
	RD/=(1.00+exp(((microenvironment(PositionFinal)[debris_index]/(dx*dy*dz))+aed)/(microenvironment(PositionFinal)[debris_index]/(dx*dy*dz)))); // - 

	pCell->phenotype.secretion.uptake_rates[debris_index] = ke2*RD; // min^-1

	//Cytokine secretion
	//DYNAMICS TO CHECK
	int TNF_index = microenvironment.find_density_index("TNFa");
	int TGF_index = microenvironment.find_density_index( "TGFb" );
	int IL10_index = microenvironment.find_density_index( "IL10" );
	int IFN_index = microenvironment.find_density_index( "IFNg" );

	float TNFfact=1.0;
	float IL10fact=1.0;

	if (mech_do){
		int indexmin=0;
		float distmin=100.0;
		float dist;
		float xCell=CellPosition[1]/1000.0;
		float yCell=CellPosition[2]/1000.0;
		float zCell=CellPosition[0]/1000.0;

		if ((abs(xCell)<0.90)&&(abs(yCell)<0.5)&&(abs(zCell)<0.5)){
			for(int n=0; n<TNFamech.size(); n++ ){
				dist=sqrt((TNFamech[n][1]-xCell)*(TNFamech[n][1]-xCell)+(TNFamech[n][2]-yCell)*(TNFamech[n][2]-yCell)+(TNFamech[n][3]-zCell)*(TNFamech[n][3]-zCell));
				if (dist<distmin){
					indexmin=n;
					distmin=dist;
				}
			} 
			TNFfact=1.00*TNFamech[indexmin][5]/112.81;
			IL10fact=1.00*IL10mech[indexmin][5]/100.0;
		}
	}

	// SECRETION TO CHECK

	// TNFa
	pCell->phenotype.secretion.secretion_rates[TNF_index]=parameters.doubles("k_T2")*(TNFfact); // min^-1
	pCell->phenotype.secretion.secretion_rates[TNF_index]*=(1.00+(parameters.doubles("k_TI")/(1.00+1.00*exp(ecf*(parameters.doubles("a_TI")-microenvironment(PositionFinal)[IFN_index]/(dx*dy*dz)))))); //min^-1
	
	pCell->phenotype.secretion.secretion_rates[TGF_index] =1.00*parameters.doubles("k_TG2"); // min^-1
	
	pCell->phenotype.secretion.secretion_rates[IL10_index] = 1.00*parameters.doubles("k_I2")*IL10fact; // min^-1
	
	pCell->phenotype.secretion.secretion_rates[IFN_index] = 0.00;

	// PROLIFERATION
	pCell->phenotype.cycle.data.transition_rate(0,0)*=(parameters.doubles("k_2")+((1-parameters.doubles("k_2"))/(1.00+exp((microenvironment(PositionFinal)[IFN_index]/(dx*dy*dz))*log(3)/parameters.doubles("a_2")))));; //min^-1

	//DIFFERENTIATE INTO M1

	double k21=parameters.doubles("k_21"); // min^-1
	double a21=parameters.doubles("a_21"); // ng um^-3
	double g1=k21*microenvironment(PositionFinal)[TNF_index]/(dx*dy*dz);
	g1/=(a21+(microenvironment(PositionFinal)[TNF_index]/(dx*dy*dz))); // min^-1

	double gtest2=(double)rand() / RAND_MAX;

	if(gtest2<=g1*dt)  
	{
		pCell->type=2;
		pCell->functions.update_phenotype = macrophage1_function; 
	}

	boundaries(pCell);
	return;
}

void neutrophils_function( Cell* pCell, Phenotype& phenotype, double dt )
{

	pugi::xml_node node = xml_find_node( physicell_config_root , "domain" );
	
	double dx = xml_get_double_value( node, "dx" ); 
	double dy = xml_get_double_value( node, "dy" ); 
	double dz = xml_get_double_value( node, "dz" ); 

	int debris_index = microenvironment.find_density_index("debris");	
	int TNF_index = microenvironment.find_density_index("TNFa");	
	int TGF_index = microenvironment.find_density_index( "TGFb" );
	int IL10_index = microenvironment.find_density_index( "IL10" );
	int IFN_index = microenvironment.find_density_index( "IFNg" );
	
	Cell* pTestCell = NULL; 
	std::vector<Cell*> neighbors = get_possible_neighbors(pCell);

	// More debris it eats, faster it dies	

	int napoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Neutrophil Apoptosis" );
	pCell->phenotype.death.rates[napoptosis_model_index] = parameters.doubles("k_AN"); 
	
	double neutrodeb=pCell->phenotype.molecular.internalized_total_substrates[debris_index]/pCell->phenotype.volume.total;

	if (neutrodeb-1.6e-3>0)
	{
		pCell->phenotype.death.rates[napoptosis_model_index] *= (1.00+(neutrodeb-1.6e-3));
	}

	std::vector<double> CellPosition = pCell->position;

	bool FoundPosition=false;
	int PositionFinal=0;

	double xVoxel;
	double yVoxel;
	double zVoxel;

	for(int n=0; n< microenvironment.number_of_voxels() && FoundPosition==false; n++ )
	{
		xVoxel=microenvironment.mesh.voxels[n].center[0];
		yVoxel=microenvironment.mesh.voxels[n].center[1];
		zVoxel=microenvironment.mesh.voxels[n].center[2];
		
		if(abs(CellPosition[0]-xVoxel)<=microenvironment.mesh.dx/2 && abs(CellPosition[1]-yVoxel)<=microenvironment.mesh.dx/2 && abs(CellPosition[2]-zVoxel)<=microenvironment.mesh.dx/2)
		{
			FoundPosition=true;
			PositionFinal=n;
		}
	}

	//UPDATE DEBRIS UPTAKE // work on this

	double ken=parameters.doubles("k_en");
	double aedn=parameters.doubles("a_edn");
	double RD=(1.00+exp(1.00)); //-
	RD/=(1.00+exp(((microenvironment(PositionFinal)[debris_index]/(dx*dy*dz))+aedn)/(microenvironment(PositionFinal)[debris_index]/(dx*dy*dz)))); // - 

	pCell->phenotype.secretion.uptake_rates[debris_index] = ken*RD;
	
	// CYTOKINE SECRETION  

	pCell->phenotype.secretion.secretion_rates[TNF_index]=1.00*parameters.doubles("k_TN");
	pCell->phenotype.secretion.secretion_rates[IL10_index]=0.00;
	pCell->phenotype.secretion.secretion_rates[TGF_index]=0.00;
	pCell->phenotype.secretion.secretion_rates[IFN_index]=1.00*parameters.doubles("k_IFN"); //  /(1+exp(2.5e-12-microenvironment(PositionFinal)[IFN_index]));

	boundaries(pCell);
	return; 
}

void macrophage_chemotaxis( Cell* pCell, Phenotype& phenotype, double dt )
{
	static double bias = parameters.doubles("m_mb");
	static int debris_index = microenvironment.find_density_index( "debris" ); 
	
	phenotype.motility.migration_bias = bias; 
	
	phenotype.motility.migration_bias_direction = pCell->nearest_gradient( debris_index ); 
	double denominator =  norm( phenotype.motility.migration_bias_direction ) + 1e-17; 
	
	phenotype.motility.migration_bias_direction /= denominator; 
	
	return; 
}

void neutrophils_chemotaxis( Cell* pCell, Phenotype& phenotype, double dt )
{
	static double bias = parameters.doubles("n_mb");
	static int debris_index = microenvironment.find_density_index( "debris" ); 
	
	phenotype.motility.migration_bias = bias; 
	
	phenotype.motility.migration_bias_direction = pCell->nearest_gradient( debris_index ); 
	double denominator =  norm( phenotype.motility.migration_bias_direction ) + 1e-17; 
	
	phenotype.motility.migration_bias_direction /= denominator; 
	
	return; 
}

std::vector<double> integrate_total_substrates( void )
{
	// start with 0 vector 
	std::vector<double> out( microenvironment.number_of_densities() , 0.0 ); 

	// integrate extracellular substrates 
	for( unsigned int n = 0; n < microenvironment.number_of_voxels() ; n++ )
	{
		// out = out + microenvironment(n) * dV(n) 
		axpy( &out , microenvironment.mesh.voxels[n].volume , microenvironment(n) ); 
	}

	// inte
	for( unsigned int n=0; n < (*all_cells).size(); n++ )
	{
		Cell* pC = (*all_cells)[n];
		out += pC->phenotype.molecular.internalized_total_substrates;
	}
	
	return out; 
}

std::vector<Cell*> get_possible_neighbors( Cell* pCell )
{
	std::vector<Cell*> neighbors = {}; 

	// First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end =
		pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for( neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{ neighbors.push_back( *neighbor ); }

	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end; 
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{ neighbors.push_back( *neighbor ); }
	}
	
	return neighbors; 
}