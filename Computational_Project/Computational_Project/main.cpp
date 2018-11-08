// This code is writen by : Maedeh Seyedolmohadesin
// Computational Physics Project Fall 2018, Northeastern University, Boston, MA
// This code will simulate a monoatomic Lennard jones system consists of N atoms in an (N,V,E) ensemble.
// BCC Lattice is used for the initial positions of the atoms
// Reduced unit is used here
// RDF stands for Radial Distribution Function

// **************************************************** Include Libraries ****************************************************

#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <fstream>

using namespace std;

//********************************************************* Settings *********************************************************

//system settings

int N0 =7 ; // number of one type atom in each direction (BCC lattic)
int N =  2 * pow (N0 , 3); // total number of atoms
double rho=1.2; // density of the system
double sigma1=1.0; // Lennard jones parameter
double mass1=1.0;
double sigma2=1.0;
double mass2=1.0;
double K=1.0; //Boltzman const.
double epsilon=1.0; // Lennard jones parameter
double d=pow((N/(rho)),(1.00/3.00)); // length of box in each dimension

//MD settings

int nsteps=1000; // number of steps
double dt=0.001; // time element
int nstxout=100; // number of steps that elapse between writing coordinates to output trajectory file
double force_cut=3.5 ; // after this radius the force won't be calulated
double T=1.2; // initial temperature

// RDF setting

int nstrdf = 1000; // number of steps that elapse between writing the RDF outputs
int rdfpoints= 1000; // Number of data points on RDF plot

// Energy  variables

double E=0 ; // total energy of the system
double KE=0 ; // total kinetic energy of the system
double U=0 ; // total potential energy of the system

//****************************************************** Global Variables ******************************************************

// RDF Variables

double binsize = (d / 2) / rdfpoints ;
double *bin = new double[rdfpoints] ;

// Anderson Thermostat

double temp=1.2; // Temperature
int nsttemp=1; // number of step

// ******************************************************* class Atom **********************************************************

class Atom
{
public:
    
    double x;
    double y;        // positions
    double z;
    double vx;
    double vy;       // velocities
    double vz;
    double fx;
    double fy;       // forces
    double fz;
    double u;  // potential energy of atom
    double ke; // kinetic energy of atom
    
};

//****************************************************** List of Functions ******************************************************

Atom* gen_atom ();                                        // Generates objects for class Atom

void init_pos( Atom *p );                                 // Initializes the positions

double normal_dist (double sigma);                        // generates numbers from normal dist. with variance sigma^2 and mean 0

void init_vel( Atom *p );                                 // Initializes the velocities

Atom* initialize ();                                      // Initializes the system

void calculate_forces (Atom *p);                          // Calculates forces on atoms

void MD_run (Atom* p);                                    // Runs the simulations and does the calculations

void write_trj (ofstream &file, Atom*p, int n);           // Writes the trajectories of atoms Atom*p into 'file'

void check_boundry (Atom* p, int i);                      // Apllies the periodic boundry condition

//******************************************************* main() function *******************************************************

int main()
{
    srand(time(NULL));
    
    Atom *atom_pointer; // defining pointer to store class object's addresses
    
    atom_pointer=initialize();
    
    MD_run(atom_pointer);
    
    // ******** test *******
    
    return 0;
}


