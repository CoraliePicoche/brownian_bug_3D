//FB 2023-11-04 Testing the initial move function after Coralie's code. We just make it move around in 3D. 

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include "basic_particle.h"
#include <random>
#include <vector>
#include <array>
#include <fstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <time.h>
#include <iomanip>

gsl_rng *rgslbis2 = gsl_rng_alloc(gsl_rng_mt19937);

extern const int num_simu=1; //Simulation identifier

extern const double pi=3.14159265;
extern const int print_distrib=1; //If 1, we output the position of each particle at the end of the simulation. This is not recommended for huge populations.
extern const int tmax=10; //length of the simulation. Tmax is negative if we only want the initial distribution

//All variables are defined as a function of the duration of \tau (or U?)
extern const double tau=0.00028; //in day
extern const double Utot=0.5; //advection, corresponds to U\tau/3
//extern const double Utot=0.0; //advection, corresponds to U\tau/3

//Diatoms
extern const double radius=25*pow(10,-6);
extern const double growth_rate=1; //in day^-1
extern const double Lmax=pow(1000,1.0/3.0); //size of the grid: adapted here to always use size_pop=10 000, but different concentrations


//Community definition
extern const int nb_species=1; //not even used

//Environment
extern const double volume=Lmax*Lmax*Lmax; 
extern const double k=2*pi; //could be 2pi/Lmax, but then scaling leads to another flow which does not have the same properties 


using namespace std;

// Main program
int main()
{
	int i,j,t,s;
	double a_x,a_y,a_z,phi,theta,psi,a_n,tmp_pop;
	std::vector<basic_particle> Part_table;
	std::ofstream f_space;
	f_space << std::fixed << setprecision(10) << endl;

	gsl_rng_set(rgslbis2, num_simu); // addendum FB 02/11/2024. Default seed is zero, for previous simulations. 


	//Open the file in which we will have the x, y, parent of each particle
	f_space.open("Spatial_distribution_"+std::to_string(num_simu)+".txt");

	//Initialize
	j=0;
	a_x=gsl_rng_uniform(rgslbis2)*Lmax;
	a_y=gsl_rng_uniform(rgslbis2)*Lmax;
	a_z=gsl_rng_uniform(rgslbis2)*Lmax;
	s=1; // species ID
	i=1; //individual ID
	Part_table.push_back(basic_particle(a_x,a_y,a_z,a_y,i,s));
	
	//Run the simulation
	for(t=0;t<=tmax;t++)
	{
	
		//For each particle, move
		for(j=0; j<Part_table.size(); j++) // normally should only select j=1
		{
			Part_table[j].initial_move(2*radius, Lmax);
			
					
		//Printing
		f_space<< j <<";";
		f_space<< Part_table[j].get_x()  << ';';
		f_space<< Part_table[j].get_y()  << ';';
		f_space<< Part_table[j].get_z()  << ';';
		f_space<< Part_table[j].get_yfirst()  << ";";
		f_space<< Part_table[j].get_firstparent() << ";";
		f_space<< 2*radius << ";";
		f_space<< Part_table[j].get_species()  << std::endl;
		
		} //end j, i.e. the particle mvt
			
		
	} //end t, i.e the whole simulation
	std::cout<<"End simulation"<<std::endl;
	f_space.close();
	return 0;
}

int main();
