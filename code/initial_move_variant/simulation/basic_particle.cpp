#include <stdio.h>
#include <stdlib.h>
#include "basic_particle.h"
#include <random>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>

gsl_rng *rgslbis = gsl_rng_alloc(gsl_rng_mt19937);

using namespace std;

//Creator and destructor
basic_particle::basic_particle ()
{
  x=0;
  y=0;
  z=0;
  y_first=y;
  first_parent=-1; //default value is wrong
}


basic_particle::~basic_particle()
{
}

basic_particle::basic_particle (double pos_x , double pos_y, double pos_z, double pos_yfirst, int first, int sp)
{
  x=pos_x;
  y=pos_y;
  z=pos_z;
  y_first=pos_yfirst;
  first_parent=first;
  species=sp;
}

//Use private variables
void basic_particle::set_x (double pos_x)
{
  x=pos_x;
}

void basic_particle::set_y (double pos_y)
{
  y=pos_y;
}

void basic_particle::set_z (double pos_z)
{
  z=pos_z;
}

void basic_particle::set_yfirst (double pos_yfirst)
{
  y_first=pos_yfirst;
}

void basic_particle::set_firstparent (int first)
{
  first_parent=first;
}

double basic_particle::get_x()
{
  return x;
}

double basic_particle::get_y()
{
  return y;
}

double basic_particle::get_z()
{
  return z;
}

double basic_particle::get_yfirst()
{
  return y_first;
}

int basic_particle::get_firstparent()
{
  return first_parent;
}

int basic_particle::get_species()
{
  return species;
}


//Movement of each particle
//Diffusion
void basic_particle::diffusion(double Delta, double Lmax)
{
        double d_x,d_y,d_z;
	d_x=gsl_ran_gaussian(rgslbis,Delta);
	d_y=gsl_ran_gaussian(rgslbis,Delta);
	d_z=gsl_ran_gaussian(rgslbis,Delta);
        x=x+d_x;
        y=y+d_y;
        z=z+d_z;
	check_boundaries(Lmax);
}

//Periodic conditions
void basic_particle::check_boundaries(double Lmax)
{
	x = fmod(x, Lmax);
	y = fmod(y, Lmax);
	z = fmod(z, Lmax);
	if (x>Lmax)
	{
                x=x-Lmax;
	}
        else if (x<0)
	{
                x=x+Lmax;
	}

        if (y>Lmax)
	{
                y=y-Lmax;
	}
        else if (y<0)
	{
                y=y+Lmax;
	}

        if (z>Lmax)
        {
                z=z-Lmax;
        }
        else if (z<0)
        {
                z=z+Lmax;
        }

}

//Advective flow, with turbulence described by the Pierrehumbert map (1991)
void basic_particle::pierrehumbert_flow(double U_tot, double k, double phi, double theta, double psi, double Lmax)
{
        x=x+U_tot*cos(k*y+phi);
        y=y+U_tot*cos(k*z+theta);
	z=z+U_tot*cos(k*x+psi);
        check_boundaries(Lmax);
}

// Initial move of fixed distance in a random direction in 3D
void basic_particle::initial_move(double diameter, double Lmax)
{
        double phi, theta, pi=3.14159265;
        phi=gsl_rng_uniform(rgslbis)*2*pi; // standard notation for angle, unrelated to Pierrehumbert algo
        theta=gsl_rng_uniform(rgslbis)*pi; // between 0 and Pi, standard notation for spherical coordinates
       //here we use rgslbis not rgslbis2

        x=x+diameter*cos(phi)*sin(theta); //physicist notation
        y=y+diameter*sin(phi)*sin(theta); //theta=azimuthal, phi = polar angle
        z=z+diameter*cos(theta);
        check_boundaries(Lmax);
}
