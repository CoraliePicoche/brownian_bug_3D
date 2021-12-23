//24/03/21 CP: added species in the model by Young et al. 2001

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

gsl_rng *rgslbis2 = gsl_rng_alloc(gsl_rng_mt19937);

extern const int num_simu=100;

//Define constant for simulation
extern const char type_simul='P'; // P for Poisson distribution, T for Thomas distribution, B for Brownian Bug Model
extern const char type_init='P'; // Initial distribution of particles : uniform (Poisson) or already aggregated (Thomas)
extern const double pi=3.14159265;
extern const int print_distrib=0; //If 1, we output the position of each particle at the end of the simulation. This is not recommended for huge populations.
extern const int tmax=-10; //length of the simulation. Tmax is negative if we only want the initial distribution

//All variables are defined as a function of the duration of \tau (or U?)
extern const double tau=0.0002; //in day
//extern const double Utot=0.5; //advection, corresponds to U\tau/2
extern const double Utot=0.0; //advection, corresponds to U\tau/2

//Define variables to compute diffusivity
double R=8.314, T=293, Na=6.0225*pow(10,23), eta=pow(10,-3);
double factor=pow(R*T/(Na*3*pi*eta),0.5);

//Diatoms
extern const double radius=25*pow(10,-6);
extern const double growth_rate=1; //in day^-1

//Nanophytoplankton
//extern const double radius=1.5*pow(10,-6);
//extern const double growth_rate=2.5; //in day^-1

extern const double Delta=factor*pow(tau/radius*(3600*24),0.5)*pow(10,2); //diffusion. The factor 10^2 is here because the length unit is cm and the 3600*24 is the conversion from day to second for tau
extern const double proba_death=growth_rate*tau; //Death and birth probability
extern const double proba_repro=growth_rate*tau; //Death and birth probability

//Community definition
extern const int nb_species=3;
//extern const std::vector<double> size_pop={55000,43000,41000,18000,6500,6300,2400,2000,1500,600,400}; 
extern const std::vector<double> size_pop={10000,10000,10000}; 
extern const int N_parent_init=200;
extern const int N_children_init=0;
extern const double sigma=0.01;
//extern const std::vector<double> size_pop={N_children_init*N_parent_init,N_children_init*N_parent_init,N_parent_init*N_children_init};

//Environment
extern const double Lmax=pow(1000,1.0/3.0); //size of the grid
extern const double volume=Lmax*Lmax*Lmax; 
extern const double k=2*pi; //could be 2pi/Lmax, but then scaling leads to another flow which does not have the same properties 

//PCF computation
extern const double pow_max=-1;
extern const double pow_min=-4;
extern const int nb_r_pcf=5000; //Number of values for r when computing pcf
extern const int delta_spatstat=0; //delta is the bandwidth for the computation. If the boolean is 1, we used delta=0.26/lambda^(1/3). If not, we use a fixed delta
extern const double delta_fixed=pow(10,-5); //Only used if delta_spatstat==0


using namespace std;

//This function stores the values used in this specific simulation to avoid mistakes
void write_parameters(std::ofstream& f_param)
{
	int i=0;

	f_param<<"type_simu="<<type_simul<<std::endl;
	f_param<<"type_init="<<type_init<<std::endl;
	f_param<<"tau="<<tau<<std::endl;
	f_param<<"Utot="<<Utot<<std::endl;
	f_param<<"radius="<<radius<<std::endl;
	f_param<<"growth_rate="<<growth_rate<<std::endl;
	f_param<<"Delta="<<Delta<<std::endl;
	f_param<<"proba_repro="<<proba_repro<<std::endl;
	f_param<<"proba_death="<<proba_death<<std::endl;
	f_param<<"nb_species="<<nb_species<<std::endl;
	for (i=0; i<nb_species; i++){
	f_param<<"init_size "<<i<<"="<<size_pop[i]<<std::endl;
	}
	f_param<<"N_parent="<<N_parent_init<<std::endl;
	f_param<<"N_children="<<N_children_init<<std::endl;
	f_param<<"sigma="<<sigma<<std::endl;
	f_param<<"L="<<Lmax<<std::endl;
	f_param<<"volume="<<volume<<std::endl;
	f_param<<"L="<<Lmax<<std::endl;
	f_param<<"tmax="<<tmax<<std::endl;
}

std::vector<basic_particle> initialize_thomas(std::vector<basic_particle> Part_table,int s)
{
        int i,j,N_children_tmp,N_parent,init_size;
        double a_x,a_y,a_z;

        N_parent=N_parent_init;
        init_size=Part_table.size();
        //Initialize the parent distribution
        for(i=init_size-1;i>=(init_size-N_parent);i--)
        {
                if(Part_table.at(i).get_species()==s){
                N_children_tmp=gsl_ran_poisson(rgslbis2,N_children_init);
                for(j=0;j<N_children_tmp;j++)
        {
                a_x=Part_table.at(i).get_x()+gsl_ran_gaussian(rgslbis2,sigma);
                a_y=Part_table.at(i).get_y()+gsl_ran_gaussian(rgslbis2,sigma);
                a_z=Part_table.at(i).get_y()+gsl_ran_gaussian(rgslbis2,sigma);
                Part_table.push_back(basic_particle(a_x,a_y,a_z,a_y,i,s)); //Children are put after their parents
                Part_table.at(Part_table.size()-1).check_boundaries(Lmax);
        } //end n_children
        } //end sp(i)==s

        }
        return Part_table;
}


//This function is mostly a copy-paste of pcf3est in the spatstat package
void K_PCF_kernel_spatstat(std::vector<basic_particle> Part_table, int nb_indiv[], double pcf[nb_species][nb_species][nb_r_pcf], double dominance[nb_species][nb_r_pcf],std::ofstream& f_pcf, std::ofstream& f_param)
{
double vx,vy,vz,dx,dy,dz,dist,lmin,lmax,invweight,frac,kernel,coef,delta,tval,dt,rondel,bias;
basic_particle temp, current;
int l,p1=0,p2=0,s1=0,s2=0;
double Concentration, h, Concentration_square;
double lambda_K[nb_species][nb_species][nb_r_pcf],K[nb_species][nb_species][nb_r_pcf];
clock_t t1, t2;
double temps, num, denom;

//Debug variables
int p2_min=0,p1_min=0;

dt=(pow(10,pow_max)-pow(10,pow_min))/(nb_r_pcf-1);

//Initialize
for(p1=0;p1<nb_species;p1++){
        for(l=0;l<nb_r_pcf;l++){
                dominance[p1][l]=0.0;
                for(p2=0;p2<nb_species;p2++){
                        pcf[p1][p2][l]=0.0;
                	lambda_K[p1][p2][l]=0.0;
                	K[p1][p2][l]=0.0;
                }
        }
}
	t1=clock();
    // double loop (or loop on all pairs)
    for(p1=0;p1<Part_table.size();p1++)
{
    current=Part_table.at(p1);
    s1=Part_table.at(p1).get_species();
    if(p1==0.25*Part_table.size()){
	std::cout<<"Reach 25%";
        t2 = clock();
     	temps = (float)(t2-t1)/CLOCKS_PER_SEC;
        std::cout<<"BBM PCF LOOP TEMPS NUMBER ONE = "<< temps<<std::endl;
    }
    if(p1==0.5*Part_table.size()){
	std::cout<<"Reach 50%";
        t2 = clock();
     	temps = (float)(t2-t1)/CLOCKS_PER_SEC;
        std::cout<<"BBM PCF LOOP TEMPS NUMBER ONE = "<< temps<<std::endl;
    }
    if(p1==0.75*Part_table.size()){
	std::cout<<"Reach 75%";
        t2 = clock();
     	temps = (float)(t2-t1)/CLOCKS_PER_SEC;
        std::cout<<"BBM PCF LOOP TEMPS NUMBER ONE = "<< temps<<std::endl;
    }
    if(p1==Part_table.size()){
	std::cout<<"Reach 100%";
        t2 = clock();
     	temps = (float)(t2-t1)/CLOCKS_PER_SEC;
        std::cout<<"BBM PCF LOOP TEMPS NUMBER ONE = "<< temps<<std::endl;
    }
    for (p2=0;p2<Part_table.size();p2++)
{

        if(p1!=p2)
{
                temp=Part_table.at(p2);
                s2=Part_table.at(p2).get_species();

                Concentration_square=(nb_indiv[s1]/volume)*(nb_indiv[s2]/volume);
		if (delta_spatstat==1){
        		delta=0.26*pow(pow(Concentration_square,0.5),-1.0/3.0); //Used in pcf3est spatstat
			if(p1==0 & p2==1){
				f_param<<"delta_spatstat=TRUE"<<std::endl;
				f_param<<"delta="<<delta<<std::endl;
			}
		}else{
  			delta=delta_fixed;
			if(p1==0 & p2==1){
				f_param<<"delta_spatstat=FALSE"<<std::endl;
				f_param<<"delta="<<delta<<std::endl;
			}
		}
		
		coef = (3.0/(4.0 * delta)) * 1/(Concentration_square);
		dx = temp.get_x() - current.get_x();
		dy = temp.get_y() - current.get_y();
		dz = temp.get_z() - current.get_z();
		dist = pow(dx * dx + dy * dy + dz * dz,0.5);
		//Using periodic conditions as in Illian et al. (2008) p. 184: JUST A TEST, NOT RECOMMENDED
		//dist=pow(min(abs(dx),Lmax-abs(dx)),2)+pow(min(abs(dy),Lmax-abs(dy)),2)+pow(min(abs(dz),Lmax-abs(dz)),2);
		//dist=pow(dist, 0.5);
		//Only corresponding to num 102 and 103 in simulations

	        lmin = std::ceil( ((dist - delta) - pow(10,pow_min))/dt);
	        lmax = std::floor( ((dist + delta) - pow(10,pow_min))/dt);   
		
		//First, we need to compute dominance outside of the loop with lmin and lmax which depends on delta; that should not be the case for dominance
		for(l = 0; l <nb_r_pcf; l++) {
			      tval = pow(10,pow_min) + l * dt;
		if(lmax >= 0 && lmin < nb_r_pcf) {
	  /* kernel centred at 'dist' has nonempty intersection 
	     with specified range of t values */
	  /* compute intersection */
	 		 if(lmin < 0)
	   			 lmin = 0;
			 if(lmax >= nb_r_pcf)
	    			lmax = nb_r_pcf - 1;
	  		/* compute (inverse) edge correction weight */
			  vx = Lmax - (dx > 0 ? dx : -dx);
			  vy = Lmax - (dy > 0 ? dy : -dy);
			  vz = Lmax - (dz > 0 ? dz : -dz);
			  invweight = vx * vy * vz * 4*pi * dist * dist;
			  if(invweight > 0.0) {
		////	    for(l = lmin; l <nb_r_pcf; l++) {
		////	      tval = pow(10,pow_min) + l * dt;
				if(dist<tval){
                                	lambda_K[s1][s2][l]+=1.0/(vx*vy*vz);
                                	//lambda_K[s1][s2][l]+=1.0; // Only 102 and 103
                                	K[s1][s2][l]+=1.0/(Concentration_square*vx*vy*vz);
                                	//K[s1][s2][l]+=1.0/(Concentration_square); //Only 102 and 103
				}//test on dist
				if(l>=lmin){ ////
	      /* unnormalised Epanechnikov kernel with halfwidth delta */
			      frac = (dist - tval)/delta;
			      kernel = (1 - frac * frac);
			      if(kernel > 0){
				rondel=tval/delta;
				if (rondel<=1){ //This correction is present in https://github.com/spatstat/spatstat.core/blob/master/R/fgk3.R
					bias=(3.0/4.0)*(rondel + 2.0/3.0 - (1.0/3.0)*pow(rondel,3));
				}else{
					bias=1;
				}
				pcf[s1][s2][l] = pcf[s1][s2][l] + coef*kernel / invweight*1/bias;
				//pcf[s1][s2][l] = pcf[s1][s2][l] + coef*kernel /bias * 1/(4*pi * dist * dist); //Only 102 and 103
	   		 	} //end kernel>0
				} /// end test on lmin
	  		} //end loop on l
			} //end invweight > 0
             }// end checks on lmax and lmin
	} //check p1\neq p2
    	} //end loop on p2
  	 }//end loop on p1
        for(s1=0;s1<nb_species;s1++){
                for(l=0;l<nb_r_pcf;l++){
			tval=pow(10,pow_min) + l * dt;
			num=lambda_K[s1][s1][l];
                        denom=0;
			for(s2=0;s2<nb_species;s2++){
				denom=denom+lambda_K[s1][s2][l];
			}
                        dominance[s1][l]=num/denom;
			
                        for(s2=0;s2<nb_species;s2++){
                                if(s2==s1){
                                        f_pcf << tval<<";"<< s1 <<";"<< s2 <<";"<<pcf[s1][s2][l]<<";"<<dominance[s1][l]<<";"<<lambda_K[s1][s2][l]<<";"<<K[s1][s2][l]<<std::endl;
                                }else{
                                        f_pcf << tval<<";"<< s1 <<";"<< s2 <<";"<<pcf[s1][s2][l]<<";NA;"<<lambda_K[s1][s2][l]<<";"<<K[s1][s2][l]<<std::endl;
                                }
                        }
                }
        }
}

//This function ouptuts a distribution of distances between pairs of particles (mostly for debugging purposes)
void distrib_distance(std::vector<basic_particle> Part_table, int repart[])
{
double d2,dt2,dtt2;
int ki,kj,i,kh;
basic_particle temp, current;
int iter=0;
int p1=0,p2=0;
int pow_dist,id_pow;
std::ofstream f3;


//f3.open("distance_table_Poisson.txt");

    for(p1=0;p1<Part_table.size();p1++)
{
    current=Part_table.at(p1);
    p2=p1;
    // double loop (or loop on all pairs)
    while (p2<(Part_table.size()-1))
{
	p2++;
                temp=Part_table.at(p2);
                dtt2=volume;
                for (ki=-1;ki<=1;ki++)
            {
                    for (kj=-1;kj<=1;kj++)
                    {
                    	for (kh=-1;kh<=1;kh++)
                    	{
                        	dt2 = pow((temp.get_x() - current.get_x()+ki*Lmax),2) + pow((temp.get_y() - current.get_y()+kj*Lmax),2) + pow((temp.get_z() - current.get_z()+kh*Lmax),2) ;
                        	d2 = std::min(dtt2,dt2);
                        	dtt2 = d2;
			}
                      }
             }

                        pow_dist=int(std::max(int(-10),int(round(log10(pow(d2,0.5))))));
                        id_pow=-1*pow_dist;
                        repart[id_pow]=repart[id_pow]+1;
//                        f3<<p1<<";"<<p2<<";"<<pow(d2,0.5)<<std::endl;

}
}
//f3.close();
}


//Births and deaths in a community. The current vector of particles is updated with the new particles (at the end of the table), then we remove the dead ones
void branching_process(std::vector<basic_particle> &part_1,double proba_repro, double proba_death)
{
	double proba;	
	std::vector <int> id_to_remove;
	int size=part_1.size();
	basic_particle tmp_part;

	for(int j=0;j<size;j++)
	{
		proba=gsl_rng_uniform(rgslbis2); //Computes a number between 0 and 1, then compare it with the different probabilites of death, birth, or just staying alive
		if (proba<proba_repro) //new particle at the same place, and the parent is still there
		{
			tmp_part=basic_particle(part_1.at(j).get_x(),part_1.at(j).get_y(),part_1.at(j).get_z(),part_1.at(j).get_yfirst(),part_1.at(j).get_firstparent(),part_1.at(j).get_species());
			part_1.push_back(tmp_part);

               } else if( (proba_repro<=proba) and (proba<(proba_repro+proba_death)) )
		{
			id_to_remove.push_back(j); //Here, we just do a list of the indices we will need to remove at the end. We cannot remove them on the fly because it would modify the sequence of rows
		}
	}
	for(int j=id_to_remove.size()-1;j>=0;j--) //We remove dead individuals starting at the end of the table. Indeed, if we were to start by the beginning, removing row number 2 would shift all following rows (element 4 from the old table would be element 3, for example)
	{
			part_1.erase(part_1.begin()+id_to_remove[j]);
	}
} 

// Main program
int main()
{
	int i,j,t,s,s1,s2;
	double a_x,a_y,a_z,phi,theta,psi,a_n,tmp_pop;
	std::vector<basic_particle> Part_table;
	std::ofstream f_pcf,f_param,f_space,f_end_simu;
	int nb_indiv[nb_species];
	int repart[11];
	double pcf[nb_species][nb_species][nb_r_pcf];
	double dominance[nb_species][nb_r_pcf];

	if(type_simul != 'B' & tmax>0){
		std::cout<<"Conflicting parameters type_simul and tmax"<<std::endl;
		return 0;
	}
	if(type_simul == 'T' & N_parent_init<=0){
		std::cout<<"Conflicting parameters type_simul and N_parent_init"<<std::endl;
		return 0;
	}

	//Open the file in which we will have the x, y, parent of each particle
	f_space.open("Spatial_distribution_"+std::to_string(num_simu)+".txt");
	f_end_simu.open("nb_indiv_"+std::to_string(num_simu)+".txt");
	f_pcf.open("lambda_K_"+std::to_string(num_simu)+".txt");
	f_param.open("param_"+std::to_string(num_simu)+".txt");

	write_parameters(f_param);

	//Initialize
	j=0;
	for (s=0; s< nb_species ; s++)
	{
		nb_indiv[s]=0;

	if (type_init=='P'){ //Poisson Distribution
		for(i=0; i < size_pop.at(s); i++)
		{
			a_x=gsl_rng_uniform(rgslbis2)*Lmax;
			a_y=gsl_rng_uniform(rgslbis2)*Lmax;
			a_z=gsl_rng_uniform(rgslbis2)*Lmax;
			Part_table.push_back(basic_particle(a_x,a_y,a_z,a_y,i,s));
	                nb_indiv[Part_table.at(j).get_species()]=nb_indiv[Part_table.at(j).get_species()]+1;
			j++;
		}
	} else if (type_init=='T'){ //Thomas distribution
	        for(i=0; i <N_parent_init; i++)
        	{
                	a_x=gsl_rng_uniform(rgslbis2)*Lmax;
	                a_y=gsl_rng_uniform(rgslbis2)*Lmax;
        	        a_z=gsl_rng_uniform(rgslbis2)*Lmax;
                	Part_table.push_back(basic_particle(a_x,a_y,a_z,a_y,i,s));
	                nb_indiv[s]=nb_indiv[s]+1;
                	j++;
	        }
                tmp_pop=Part_table.size();
                Part_table=initialize_thomas(Part_table,s);
                nb_indiv[s]=nb_indiv[s]+Part_table.size()-tmp_pop;
                std::cout<<"For sp="<<s<<" End parent="<<nb_indiv[s]<<std::endl;

	} //End test on initial distribution
		std::cout<<"SPECIES="<<s<<" INDIV="<<nb_indiv[s]<<std::endl;
	} //End loop on species


	//Run the simulation
	for(t=0;t<=tmax;t++)
	{
	//Compute the phase in x and y for the turbulent flow from Pierrehumbert. These phases are common to each particle as they correspond to a unique flow
	phi=gsl_rng_uniform(rgslbis2)*2*pi;
	theta=gsl_rng_uniform(rgslbis2)*2*pi;
	psi=gsl_rng_uniform(rgslbis2)*2*pi;

	//First step: birth and death	
	branching_process(Part_table,proba_repro,proba_death);

	//For each particle, diffusion, then advection
	for(j=0; j<Part_table.size(); j++)
	{
		Part_table[j].diffusion(Delta, Lmax);
		Part_table[j].pierrehumbert_flow(Utot, k, phi,theta,psi, Lmax);

	} //end j, i.e. the particle mvt
	} //end t, i.e the whole simulation
	std::cout<<"End simulation"<<std::endl;

	//Reinitialize nb species to use it when computing PCF at the end
	for (s=0; s< nb_species ; s++)
	{
		nb_indiv[s]=0;
	}
	
	//This is not optimized, but makes the code a bit cleaner by separating the simulation part and the output part
	for(j=0; j<Part_table.size(); j++)
        {
                nb_indiv[Part_table.at(j).get_species()]=nb_indiv[Part_table.at(j).get_species()]+1;
                if(print_distrib==1)
                {
        	        f_space<< j <<";";
                	f_space<< Part_table[j].get_x()  << ';';
	                f_space<< Part_table[j].get_y()  << ';';
        	        f_space<< Part_table[j].get_z()  << ';';
                	f_space<< Part_table[j].get_yfirst()  << ";";
	                f_space<< Part_table[j].get_firstparent() << ";";
        	        f_space<< Part_table[j].get_species()  << std::endl;
                }
        }

	//End of the simulation

//Distribution of distance between pairs of particles, for debugging purposes
 //       for(i=0;i<11;i++){
 //               repart[i]=0;
 //       }
//	distrib_distance(Part_table, repart);
  /*      for(i=0;i<11;i++){
              f0<<Utot<<";"<<i<<";"<<repart[i]<<std::endl;
       } */

        for (s1=0;s1<nb_species;s1++)
        {
        	f_end_simu<<s1<<";"<<nb_indiv[s1]<<std::endl;
	}

        K_PCF_kernel_spatstat(Part_table, nb_indiv, pcf, dominance,f_pcf,f_param);

	f_space.close();
	f_end_simu.close();
	f_pcf.close();
	f_param.close();
	return 0;
}

int main();
