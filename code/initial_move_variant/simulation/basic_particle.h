#include <string>

class basic_particle
{
  double x;
  double y;
  double z;
  double y_first; //y_first is the initial y position of the first parent, tracked for plotting reasons
  int first_parent;
  int species; //id of the species

 public:
  //Creator, destructor
  basic_particle();
  basic_particle (double pos_x, double pos_y, double pos_z, double pos_yfirst, int first, int sp);
  ~basic_particle();

  //Use private variables
  void set_x (double pos_x);
  void set_y (double pos_y);
  void set_z (double pos_z);
  void set_yfirst (double pos_yfirst);
  void set_firstparent (int first);
  double get_x();
  double get_y();
  double get_z();
  double get_yfirst();
  int get_firstparent();
  int get_species();

  //Movement of each particle
  void diffusion(double Delta,double Lmax); //new values of x and y after particle random walk
  void check_boundaries(double Lmax); //checking that particles remain in a Lmax x Lmax grid
  void pierrehumbert_flow(double U_tot, double k, double phi, double theta, double psi, double Lmax); //new values of x and y after advection
  void initial_move(double diameter, double Lmax); //new values of x,y and z after initial movement
};
