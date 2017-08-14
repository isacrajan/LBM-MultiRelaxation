/*
	Class Header File for LBM-MRT
	Single Sided Lid Driven Cavity
	Author: Isac Rajan
	Start Date:	09 August 2017
	Email:	isacrajan@gmail.com
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <ctime>
#include <cstdio>
#include <string>
#include <sstream>

#define NX 201 // No of grids in X direction
#define NY 201 // No of grids in Y direction

using std::vector;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::ostringstream;

typedef double mydoub;
typedef int myint;
typedef vector<vector<vector<double> > > vdouble3;
typedef vector<vector<double> > vdouble2;
typedef vector<double> vdouble1;

class LBM_MRT
{
	private:
		//Using the D2Q9 Lattice Model
		myint cx[9];
		myint cy[9];
		mydoub w[9]; //weights
		mydoub rhoo;
		mydoub uo;
		mydoub Re; //Reynolds Number
		mydoub nu; //Kinematic Viscosity
		mydoub omega;
		
		vdouble3 meq, m, f;
		vdouble2 u, v, rho, u_old, v_old, S;
		
	public:
		myint t;
		myint TSTEPS; //Max number of timesteps
		LBM_MRT(); //constructor
		void collision();
		void stream();
		void boundaryCondition();
		void computeField();
};

LBM_MRT::LBM_MRT(void)
{
	//Initializing the variables
	TSTEPS = 150000; //Max number of timesteps
	rhoo = 1.0; 
	Re = 100; //Reynolds Number
	uo = 0.01;
	
	//Initialization of lattice weights and lattice velocities
	cx[0] = 0;
	cx[1] = 1;
	cx[2] = 0;
	cx[3] = -1;
	cx[4] = 0;
	cx[5] = 1;
	cx[6] = -1;
	cx[7] = -1;
	cx[8] = 1;
	cy[0] = 0;
	cy[1] = 0;
	cy[2] = 1;
	cy[3] = 0;
	cy[4] = -1;
	cy[5] = 1;
	cy[6] = 1;
	cy[7] = -1;
	cy[8] = -1;
	w[0] = 4./9;
	w[1] = w[2] = w[3] = w[4] = 1./9;
	w[5] = w[6] = w[7] = w[8] = 1./36;

	nu = (mydoub)uo*NY/Re; //Kinematic Viscosity
	omega = 1.0/(3.0*nu + 0.5); // inverse of tau, as in formulation
	
	//Initializing rho, u, v vectors
	rho.resize(NX);
	u.resize(NX); 
	u_old.resize(NX);
	v.resize(NX); 
	v_old.resize(NX);
	for(int i=0; i<NX; i++)
	{
		rho[i].resize(NY);
		u[i].resize(NY); 
		u_old[i].resize(NY);
		v[i].resize(NY); 
		v_old[i].resize(NY);
	}
	for(int i=0; i<NX; i++)
		for(int j=0; j<NY; j++)
			{
				rho[i][j] = rhoo;
				u[i][j] = 0.0; 
				u_old[i][j] = 0.0; 
				v[i][j] = 0.0; 
				v_old[i][j] = 0.0;
			}
	
	//Initializing the velocity of the lid
	for (int i=0; i<NX; i++)
	{
		u[i][NY-1]=uo;
		v[i][NY-1]=0.0;
	}
	
	//Initialization of m, meq
	m.resize(9);
	meq.resize(9);
	f.resize(9);
	for(int k=0; k<9; k++)
	{
		m[k].resize(NX);
		meq[k].resize(NX);
		f[k].resize(NX);
	}
	for(int k=0; k<9; k++)
		for(int i=0; i<NX; i++)
		{
			m[k][i].resize(NY);
			meq[k][i].resize(NY);
			f[k][i].resize(NY);
		}
	for(int k=0; k<9; k++)
		for(int i=0; i<NX; i++)
			for(int j=0; j<NY; j++)
			{
				switch(k)
				{
					case 0:
						meq[k][i][j] = rho[i][j];
						break;
					case 1:
						meq[k][i][j] = -2*rho[i][j] + 3(u[i][j]*u[i][j]+v[i][j]*v[i][j]);
						break;
					case 2:
						meq[k][i][j] = rho[i][j] - 3(u[i][j]*u[i][j]+v[i][j]*v[i][j]);
						break;
					case 3:
						meq[k][i][j] = rho[i][j]*u[i][j];
						break;
					case 4:
						meq[k][i][j] = -u[i][j];
						break;
					case 5:
						meq[k][i][j] = rho[i][j]*v[i][j];
						break;
					case 6:
						meq[k][i][j] = -v[i][j];
						break;
					case 7:
						meq[k][i][j] = u[i][j]*u[i][j] - v[i][j]*v[i][j];
						break;
					case 8:
						meq[k][i][j] = u[i][j]*v[i][j];
						break;
				}
				m[k][i][j] = meq[k][i][j];
			}
	
	//Intialization of M & S
	vdouble2 M {{1, 1, 1, 1, 1, 1, 1, 1, 1},
							{-4, -1, -1, -1, -1, 2, 2, 2, 2},
							{4, -2, -2, -2, -2, 1, 1, 1, 1},
							{0, 1, 0, -1, 0, 1, -1, -1, 1},
							{0, -2, 0, 2, 0, 1, -1, -1, 1},
							{0, 0, 1, 0, -1, 1, 1, -1, -1},
							{0, 0, -2, 0, 2, 1, 1, -1, -1},
							{0, 1, -1, 1, -1, 0, 0, 0, 0},
							{0, 0, 0, 0, 0, 1, -1, 1, -1}
						 };
	S.resize(9);
	for(int i=0; i<9; i++)
		S[i].resize(9);
	for(int i=0; i<9; i++)
		for(int j=0; j<9; j++)
		{
			
		}
	
	
	
	cout << "All the variables & arrays are initialized! We are good to go!" << endl;
}










