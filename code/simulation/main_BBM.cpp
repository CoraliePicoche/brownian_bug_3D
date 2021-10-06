//24/03/21 CP: added species in the model by Young et al. 2001

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
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

//Define constant for simulation
extern const double pi=3.14159265;

//Diatoms
extern const double Delta=7*pow(10,-5); //diffusion
extern const double pow_min=-3;
extern const double pow_max=-0.6;
extern const double proba_death=2*pow(10,-4); //Death and birth probability
extern const double proba_repro=2*pow(10,-4);

//Nanophytoplankton
//extern const double Delta=3*pow(10,-4); //diffusion
//extern const double pow_min=-3.5;
//extern const double pow_max=-2;
//extern const double proba_death=5*pow(10,-4); //Death and birth probability
//extern const double proba_repro=5*pow(10,-4);

//Community definition
extern const int nb_species=2;
extern const std::vector<double> size_pop={100000,100000}; 

//Environment
extern const double Lmax=pow(10.0,1.0/3.0); //size of the grid
extern const double volume=Lmax*Lmax*Lmax; 
extern const double Utot=0.5; //advection, corresponds to U\tau/2 in Young et al. 2001
extern const double k=2*pi; //could be 2pi/Lmax, but then scaling leads to another flow which does not have the same properties 

//Simulation duration
extern const int tmax=1000; //length of the simulation

//Number of points for pcf computation
extern const int nb_r_pcf=1000; //Number of values for r when computing pcf

using namespace std;

//This function is mostly a copy-paste of pcf3est in the spatstat package
void PCF_kernel_spatstat(std::vector<basic_particle> Part_table, int nb_indiv[], double pcf[nb_species][nb_species][nb_r_pcf], double dominance[nb_species][nb_r_pcf],std::ofstream& f0)
{
double vx,vy,vz,dx,dy,dz,dist,lmin,lmax,invweight,frac,kernel,coef,delta,tval,dt,rondel,bias;
basic_particle temp, current;
int l,p1=0,p2=0,s1=0,s2=0;
double Concentration, h, Concentration_square;
double mingling[nb_species][2][nb_r_pcf];

//Debug variables
int p2_min=0,p1_min=0;

dt=(pow(10,pow_max)-pow(10,pow_min))/(nb_r_pcf-1);

//Initialize
for(p1=0;p1<nb_species;p1++){
        for(l=0;l<nb_r_pcf;l++){
                mingling[p1][0][l]=0.0;
                mingling[p1][1][l]=0.0;
                dominance[p1][l]=0.0;
                for(p2=0;p2<nb_species;p2++){
                        pcf[p1][p2][l]=0.0;
                }
        }
}

    // double loop (or loop on all pairs)
    for(p1=0;p1<Part_table.size();p1++)
{
    current=Part_table.at(p1);
    s1=Part_table.at(p1).get_species();
    for (p2=0;p2<Part_table.size();p2++)
{

        if(p1!=p2)
{
                temp=Part_table.at(p2);
                s2=Part_table.at(p2).get_species();

                Concentration_square=(nb_indiv[s1]/volume)*(nb_indiv[s2]/volume);
        	//delta=0.26*pow(pow(Concentration_square,0.5),-1.0/3.0); //Used in pcf3est spatstat
  		delta=0.0001;
		coef = (3.0/(4.0 * delta)) * 1/(Concentration_square);
		dx = temp.get_x() - current.get_x();
		dy = temp.get_y() - current.get_y();
		dz = temp.get_z() - current.get_z();
		dist = pow(dx * dx + dy * dy + dz * dz,0.5);

	        lmin = std::ceil( ((dist - delta) - pow(10,pow_min))/dt);
	        lmax = std::floor( ((dist + delta) - pow(10,pow_min))/dt);    
                
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
			    for(l = lmin; l <nb_r_pcf; l++) {
			      tval = pow(10,pow_min) + l * dt;
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
	   		 	} //end kernel>0
				
                        	if(dist<tval){
					if(s1==s2){ //Mingling index _ii
                                		mingling[s1][0][l]+=1.0;
                        		}else{
                                		mingling[s1][1][l]+=1.0;
                        		} //test on s1==s2
				}//test on dist
	  		} //end loop on l
			} //end invweight > 0
             }// end checks on lmax and lmin
	} //check p1\neq p2
    	} //end loop on p2
  	 }//end loop on p1
        for(s1=0;s1<nb_species;s1++){
                for(l=0;l<nb_r_pcf;l++){
			tval=pow(10,pow_min) + l * dt;
                        dominance[s1][l]=mingling[s1][0][l]/(mingling[s1][0][l]+mingling[s1][1][l]);
			
                        for(s2=0;s2<nb_species;s2++){
                                if(s2==s1){
                                        f0 << tval<<";"<< s1 <<";"<< s2 <<";"<<pcf[s1][s2][l]<<";"<<dominance[s1][l]<<std::endl;
                                }else{
                                        f0 << tval<<";"<< s1 <<";"<< s2 <<";"<<pcf[s1][s2][l]<<";NA"<<std::endl;
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
	double a_x,a_y,a_z,phi,theta,psi,a_n;
	std::vector<basic_particle> Part_table;
	std::ofstream f0,f1,f2;
	int nb_indiv[nb_species];
	clock_t t1, t2;
	double temps;
	int repart[11];
	double pcf[nb_species][nb_species][nb_r_pcf];
	double dominance[nb_species][nb_r_pcf];


	//Open the file in which we will have the x, y, parent of each particle
	//f2.open("Spatial_distribution_BBM_kernel_2sp_realistic_values_10000_U0p5.txt");
	f1.open("nb_indiv_BBM_kernel_2sp_realistic_values_20000_U0p5_test_pcf.txt");
	f0.open("pcf_BBM_kernel_2sp_realistic_values_20000_U0p5_test_pcf.txt");

	//Initialize
	j=0;
	for (s=0; s< nb_species ; s++)
	{
		nb_indiv[s]=0;
	for(i=0; i < size_pop.at(s); i++)
	{
		a_x=gsl_rng_uniform(rgslbis2)*Lmax;
		a_y=gsl_rng_uniform(rgslbis2)*Lmax;
		a_z=gsl_rng_uniform(rgslbis2)*Lmax;
		Part_table.push_back(basic_particle(a_x,a_y,a_z,a_y,i,s));
                nb_indiv[Part_table.at(j).get_species()]=nb_indiv[Part_table.at(j).get_species()]+1;
                //f2<< "0" << ";"<< Part_table[j].get_x()  << ";" << Part_table[j].get_y()  << ";" << Part_table[j].get_z() << ';' << Part_table[j].get_yfirst()  << ";" << Part_table[j].get_firstparent() << ";" << Part_table[j].get_species()  << std::endl;
		j++;

	}
		std::cout<<"SPECIES="<<s<<" INDIV="<<nb_indiv[s]<<std::endl;
	}


	//Reinitialize nb species to use it when computing PCF at the end
	for (s=0; s< nb_species ; s++)
	{
		nb_indiv[s]=0;
	}
	
	//Run the simulation
	for(t=0;t<=tmax;t++)
	{
	std::cout<<"TIME="<<t<<std::endl;
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
		if(t==tmax){
			nb_indiv[Part_table.at(j).get_species()]=nb_indiv[Part_table.at(j).get_species()]+1;
		}
		/*if(t==tmax) 
                {
                f2<< t <<";";
                f2<< Part_table[j].get_x()  << ';';
                f2<< Part_table[j].get_y()  << ';';
                f2<< Part_table[j].get_z()  << ';';
                f2<< Part_table[j].get_yfirst()  << ";";
                f2<< Part_table[j].get_firstparent() << ";";
                f2<< Part_table[j].get_species()  << std::endl;
                }*/


	} //end j, i.e. the particle mvt
	} //end t, i.e the whole simulation

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
        	f1<<s1<<";"<<nb_indiv[s1]<<std::endl;
	}

	t1=clock();
        PCF_kernel_spatstat(Part_table, nb_indiv, pcf, dominance,f0);
        t2 = clock();
     	temps = (float)(t2-t1)/CLOCKS_PER_SEC;
        std::cout<<"BBM PCF LOOP TEMPS NUMBER ONE = "<< temps<<std::endl;

//	f2.close();
	f0.close();
	f1.close();
	return 0;
}

int main();
