#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include "meep.hpp"
#include "ctl-math.h"
#include <mpb.h>
#include "ctlgeom.h"
#include "meepgeom.hpp"
#include <math.h>
#ifndef DATADIR
#define DATADIR "./"
#endif

using namespace meep;
using namespace std;





typedef std::complex<double> cdouble;

double epsr=0.9999;
//SI Conversion Factors
double a0=1e-3;//1mm
double c0=2.99792458e8;//Speed of Light (m/s)
double f0=c0/a0;//300GHz
double t0=1/f0;//0.33e-11 (s)
double mu0=4*pi*(1e-7);// (H/m)
double eps0=8.854187817e-12;// (F/m)
double I0=1; //(A)
double E0=I0/(a0*eps0*c0);//Electric Field
double D0=I0/(a0*c0);//Electric Displacement Field
double B0=I0/(a0*eps0*c0*c0);//Magnetic Field
double H0=I0/(a0);//Magnetizing Field
double sigmaD0=(epsr*eps0*c0)/a0;//Electric Conductivity
double J0=I0/(a0*a0);//Electric Current Density
double u0=(I0*I0)/(eps0*c0*c0*a0*a0);//Energy Density
double S0=(I0*I0)/(eps0*c0*a0*a0);//Poynting Vector
double Sc0=1/(c0);//Courant Factor

double xcen=0.0, ycen=0.0, zcen=0.0;
double dxmin=0.0, dxmax=0.25, dymin=0.0, dymax=0.25, dzmin=0.0, dzmax=1.0;
double wcore=2.0*dymax;
double winding_thickness_p=0.25, insulation_thickness_p=0.25,pml_thickness=1.0;  
//double sigma_Cu=100000000;
double mu_core=10000;
int Np=1;//must be odd
int Ns=1;//must be odd
double margin=0.2;   
double amplitude=1.0;
double divisions=5;

int numcoord=0;
int corecoord=0;
double k=2*pi*Np/(2*dzmin);
double theta=0.0;
double rp=(wcore/2.0)+(winding_thickness_p/2.0)+insulation_thickness_p;
double * xpcoord=new double [1000000];
double * ypcoord=new double [1000000];
double * zpcoord=new double [1000000];
double * xscoord=new double [1000000];
double * yscoord=new double [1000000];
double * zscoord=new double [1000000];
double * xccoord=new double [1000000];
double * yccoord=new double [1000000];
double * zccoord=new double [1000000];


//Copper
//metal_range = mp.FreqRange(min=um_scale/12.398, max=um_scale/.20664)
double um_scale = 1000.0;//1000um
double eV_um_scale = 1.0/1.23984193;
double Cu_plasma_frq = 10.83*eV_um_scale;
double Cu_f0 = 0.575;
double Cu_frq0 = 1e-10;
double Cu_gam0 = 0.030*eV_um_scale;
double Cu_sig0 = Cu_f0*(Cu_plasma_frq*Cu_plasma_frq)/(Cu_frq0*Cu_frq0);
double Cu_f1 = 0.061;
double Cu_frq1 = 0.291*eV_um_scale;
double Cu_gam1 = 0.378*eV_um_scale;
double Cu_sig1 = Cu_f1*(Cu_plasma_frq*Cu_plasma_frq)/(Cu_frq1*Cu_frq1);
double Cu_f2 = 0.104;
double Cu_frq2 = 2.957*eV_um_scale;      
double Cu_gam2 = 1.056*eV_um_scale;
double Cu_sig2 = Cu_f2*(Cu_plasma_frq*Cu_plasma_frq)/(Cu_frq2*Cu_frq2);
double Cu_f3 = 0.723;
double Cu_frq3 = 5.300*eV_um_scale;      
double Cu_gam3 = 3.213*eV_um_scale;
double Cu_sig3 = Cu_f3*(Cu_plasma_frq*Cu_plasma_frq)/(Cu_frq3*Cu_frq3);
double Cu_f4 = 0.638;
double Cu_frq4 = 11.18*eV_um_scale;
double Cu_gam4 = 4.305*eV_um_scale;
double Cu_sig4 = Cu_f4*(Cu_plasma_frq*Cu_plasma_frq)/(Cu_frq4*Cu_frq4);
double sigma_Cu=Cu_sig0;


//NiFe
//NiFe_range = mp.FreqRange(min=um_scale/0.83, max=um_scale/0.25)
double NiFe_frq = 1/(0.0838297450980392*um_scale);
double NiFe_gam = 1/(0.259381156903766*um_scale);
double NiFe_sig = 1;



double annulus(const vec & v, double xcen, double ycen, double zcen, double dxmin, double dxmax, double dymin, double dymax, double dzmin, double dzmax, double control, double special, double std = 0.0)
{
  double mu=control;

  double dx=v.x() - xcen;
  double dy=v.y() - ycen;
  double dz=v.z() - zcen;

  if ( (abs(dx)<=dxmax) && (abs(dy)<=dymax) && (abs(dz)<=dzmax) )
  {
    {mu= gaussian_random(special,std);}
  }

  if ( (abs(dx)<=dxmin) && (abs(dy)<=dymin) && (abs(dz)<=dzmin) )
  {
    mu= control;
  }

    return mu;  
}

bool inside_box(const vec & v, double xcen, double ycen, double zcen, double dxmax, double dymax, double dzmax)
{
  bool ans=false;

  double dx=v.x() - xcen;
  double dy=v.y() - ycen;
  double dz=v.z() - zcen;


  if ( (abs(dx)<=dxmax) && (abs(dy)<=dymax) && (abs(dz)<=dzmax) )
  {
    {ans=true;}
  }

  return ans;
}

  
cdouble line_integral_x(fields & f, component C, double dx, double xmin, double xmax, double y, double z)
{
  cdouble sum(0.0,0.0);
  cdouble deltax(dx,0.0);
  for (double x=xmin; x<=xmax; x=x+dx)
  {
    monitor_point p;
    f.get_point(&p, vec(x,y,z));
    cdouble dF = p.get_component(C);
    sum += dF*deltax;
  }
  return sum;
}

cdouble line_integral_y(fields & f, component C, double dy, double ymin, double ymax, double x, double z)
{
  cdouble sum(0.0,0.0);
  cdouble deltay(dy,0.0);
  for (double y=ymin; y<=ymax; y=y+dy)
  {
    monitor_point p;
    f.get_point(&p, vec(x,y,z));
    cdouble dF = p.get_component(C);
    sum += dF*deltay;
  }
  return sum;
}

cdouble line_integral_z(fields & f, component C, double dz, double zmin, double zmax, double x, double y)
{
  cdouble sum(0.0,0.0);
  cdouble deltaz(dz,0.0);
  for (double z=zmin; z<=zmax; z=z+dz)
  {
    monitor_point p;
    f.get_point(&p, vec(x,y,z));
    cdouble dF = p.get_component(C);
    //cout<<dF.real()<<" , "<<dF.imag()<<endl;
    sum += dF*deltaz;
    //cout<<dF.real()<<" , "<<dF.imag()<<endl;
  }
  return sum;
}

cdouble compute_Im(fields & f, double z)
{
  cdouble Iyf=line_integral_y(f,Ey,0.001,ycen-dymax,ycen+dymax,xcen+dxmax,z);
    cdouble Iyb=line_integral_y(f,Ey,0.001,ycen-dymax,ycen+dymax,xcen-dxmax,z);
    cdouble Ixt=line_integral_x(f,Ex,0.001,xcen-dxmax,xcen+dxmax,ycen-dymax,z);
    cdouble Ixb=line_integral_x(f,Ex,0.001,xcen-dxmax,xcen+dxmax,ycen+dymax,z);
    cdouble Im=Ixt+Iyf-Ixb-Iyb;
    return Im;
}

cdouble compute_Ie(fields & f, double z)
{
  cdouble Iyf=line_integral_y(f,Hy,0.001,ycen-dymax,ycen+dymax,xcen+dxmax,z);
    cdouble Iyb=line_integral_y(f,Hy,0.001,ycen-dymax,ycen+dymax,xcen-dxmax,z);
    cdouble Ixt=line_integral_x(f,Hx,0.001,xcen-dxmax,xcen+dxmax,ycen-dymax,z);
    cdouble Ixb=line_integral_x(f,Hx,0.001,xcen-dxmax,xcen+dxmax,ycen+dymax,z);
    cdouble Ie=Ixt+Iyf-Ixb-Iyb;
  return Ie;
}

cdouble compute_Vm(fields & f, double z)
{
  //cdouble Vy=line_integral_y(f,Hy,0.001,ycen-dymin,y,xcen,zcen-dzmin);
  cdouble Vz=line_integral_z(f,Hz,0.001,zcen-dzmax,z,xcen,ycen);
  //cdouble Vm=Vy+Vz;
  return Vz;
}

cdouble compute_Ve(fields & f, double z)
{
  //cdouble Vy=line_integral_y(f,Ey,0.001,ycen-dymin,y,xcen,zcen-dzmin);
  cdouble Vz=line_integral_z(f,Ez,0.001,zcen-dzmax,z,xcen,ycen);
  //cdouble Ve=Vy+Vz;
  return Vz;
}

double mu(const vec &v) 
{
  //return annulus(v,xcen,ycen,zcen,dxmin,dxmax,dymin,dymax,dzmin,dzmax,1.0,mu_core,0.0);
  return 1.0;
}

double eps(const vec &v)
{
  //return annulus(v,xcen,ycen,zcen,dxmin,dxmax,dymin,dymax,dzmin,dzmax,1.0,3.5,0.0);
  return 1.0;
}

double conductivity(const vec &v) 
{ 
}

class core_material : public material_function {
  public:
};

class winding_material : public material_function {
  public:
};


typedef struct my_material_func_data {
  double rxInner, ryInner, rOuter;
  bool with_susceptibility;
} my_material_func_data;

void my_material_func(vector3 p, void *user_data, meep_geom::medium_struct *m) {
  my_material_func_data *data = (my_material_func_data *)user_data;
  
    bool in_middle=false;

    double dx=p.x - xcen;
    double dy=p.y - ycen;
    double dz=p.z - zcen;

	for (int i=0; i<corecoord; i++) {
	double dxp=p.x - xccoord[i];
    double dyp=p.y - yccoord[i];
    double dzp=p.z - zccoord[i];  	
  	double drp=sqrt(dxp*dxp+dyp*dyp+dzp*dzp);
  	if(drp<=(wcore/2.0)){
  		in_middle=true;
  	}
  	}

  // set permittivity and permeability
  double nn = in_middle ? sqrt(mu_core) : 1.0;
  double mm = in_middle ? sqrt(10.0) : 1.0;
  m->epsilon_diag.x = m->epsilon_diag.y = m->epsilon_diag.z = mm * mm;
  //m->epsilon_offdiag.x.re = m->epsilon_offdiag.x.im = epsilon_offdiag.y.re = m->epsilon_offdiag.y.im = epsilon_offdiag.z.re = m->epsilon_offdiag.z.im = nn * nn;
  m->mu_diag.x = m->mu_diag.y = m->mu_diag.z = nn*nn;

  //m->mu_offdiag.x.re = m->mu_offdiag.x.im = mu_offdiag.y.re = m->mu_offdiag.y.im = mu_offdiag.z.re = m->mu_offdiag.z.im = nn * nn;
  m->E_chi2_diag.x = m->E_chi2_diag.y = m->E_chi2_diag.z = 1.0;
  m->E_chi3_diag.x = m->E_chi3_diag.y = m->E_chi3_diag.z = 1.0;
  m->H_chi2_diag.x = m->H_chi2_diag.y = m->H_chi2_diag.z = 1.0;
  m->H_chi3_diag.x = m->H_chi3_diag.y = m->H_chi3_diag.z = 1.0;
  //m->D_conductivity_diag.x = m->D_conductivity_diag.y = m->D_conductivity_diag.z = nn * nn;
  //m->B_conductivity_diag.x = m->B_conductivity_diag.y = m->B_conductivity_diag.z = 0.0;

  if (in_middle)
  {
    m->H_susceptibilities.num_items = 1;
    m->H_susceptibilities.items = new meep_geom::susceptibility[1];

    m->H_susceptibilities.items[0].sigma_offdiag.x = 0.0;
    m->H_susceptibilities.items[0].sigma_offdiag.y = 0.0;
    m->H_susceptibilities.items[0].sigma_offdiag.z = 0.0;
    m->H_susceptibilities.items[0].sigma_diag.x = gaussian_random(1.0,1.0);//NiFe_sig;
    m->H_susceptibilities.items[0].sigma_diag.y = gaussian_random(1.0,1.0);//NiFe_sig;
    m->H_susceptibilities.items[0].sigma_diag.z = gaussian_random(1.0,1.0);//NiFe_sig;
    m->H_susceptibilities.items[0].bias.x = gaussian_random(1.0,1.0);
    m->H_susceptibilities.items[0].bias.y = gaussian_random(1.0,1.0);
    m->H_susceptibilities.items[0].bias.z = gaussian_random(1.0,1.0);
    m->H_susceptibilities.items[0].frequency = (0.2e6)/f0;//NiFe_frq;
    m->H_susceptibilities.items[0].gamma = 100.0*((0.2e6)/f0);//NiFe_gam;
    m->H_susceptibilities.items[0].alpha = gaussian_random(1.0,1.0);//NiFe_alpha;
    m->H_susceptibilities.items[0].noise_amp = 0.01;
    m->H_susceptibilities.items[0].drude = false;
    m->H_susceptibilities.items[0].saturated_gyrotropy = true;
    m->H_susceptibilities.items[0].is_file = false;
  }

  for (int i=0; i<numcoord; i++) {
	double dxp=p.x - xpcoord[i];
    double dyp=p.y - ypcoord[i];
    double dzp=p.z - zpcoord[i];  	
  	double drp=sqrt(dxp*dxp+dyp*dyp+dzp*dzp);
  	if(drp<=(winding_thickness_p/2.0)){
  		m->D_conductivity_diag.x = m->D_conductivity_diag.y = m->D_conductivity_diag.z = 5.8e7/sigmaD0;
  	}
  	double dxs=p.x - xscoord[i];
    double dys=p.y - yscoord[i];
    double dzs=p.z - zscoord[i];
    double drs=sqrt(dxs*dxs+dys*dys+dzs*dzs);
  	if(drs<=(winding_thickness_p/2.0)){
  		m->D_conductivity_diag.x = m->D_conductivity_diag.y = m->D_conductivity_diag.z = 5.8e7/sigmaD0;
  	}
  }

  if ( (p.z<(zcen-dzmax))||(p.z>(zcen+dzmax)) ) {
  	//m->D_conductivity_diag.x = m->D_conductivity_diag.y = m->D_conductivity_diag.z = 1e20;
  	//m->B_conductivity_diag.x = m->B_conductivity_diag.y = m->B_conductivity_diag.z = 1e20;
  	m->epsilon_diag.x = m->epsilon_diag.y = m->epsilon_diag.z = 20.0;
	  //m->mu_diag.x = m->mu_diag.y = m->mu_diag.z = 1e20;
  }

}



int main(int argc, char *argv[]) {
  
  initialize mpi(argc,argv);
  const char *mydirname = "MMTL-out";
    std::ofstream Time;
    std::ofstream Space;
    std::ofstream FieldsIn;
    std::ofstream FieldsOut;
    std::ofstream Fluxes;
    std::ofstream Skin;
    Time.open ("TimeEvolution.txt");
    Space.open ("SpaceEvolution.txt");
    FieldsIn.open ("FieldEvolutionIn.txt");
    FieldsOut.open ("FieldEvolutionOut.txt");
    Fluxes.open ("Flux.txt");
    Skin.open ("Skin.txt");
    //trash_output_directory(mydirname);
    double xsize=5;
    double ysize=10;
    double zsize=10;
    //double xsize=6, ysize=6, zsize=6;
  	
  	k=2*pi*Np/(0.50);

    for (double z = -0.5; z <= 0.5; z=z+0.01)
    {
    	theta=k*(z+0.50);
    	xpcoord[numcoord]=-rp*sin(theta);
    	ypcoord[numcoord]=rp*cos(theta);
    	zpcoord[numcoord]=z;
    	xscoord[numcoord]=-rp*sin(theta);
    	yscoord[numcoord]=rp*cos(theta);;
    	zscoord[numcoord]=z+0.5;
    	numcoord++;
    }

	for(double y=-1.25;y<=1.25;y=y+1.25)
	{for (double z=-0.75;z<=0.75;z=z+0.01){
		xccoord[corecoord]=xcen;
		yccoord[corecoord]=y;
		zccoord[corecoord]=z;
		corecoord++;
	}}

	for(double z=-0.75;z<=0.75;z=z+1.5)
	{for (double y=-1.25;y<=1.25;y=y+0.01){
		xccoord[corecoord]=xcen;
		yccoord[corecoord]=y;
		zccoord[corecoord]=z;
		corecoord++;
	}}


    grid_volume gv = vol3d(xsize, ysize, zsize, divisions);
    //grid_volume vol3d(double xsize, double ysize, double zsize, double a);
    gv.center_origin();
    //void center_origin(void) { shift_origin(-icenter()); }
  

    structure transformer(gv, eps, pml(pml_thickness));
    //transformer.set_epsilon(eps,true);
    //transformer.set_mu(mu,true);
    //transformer.add_susceptibility(core_mat, H_stuff, lorentzian_susceptibility(2*pi*NiFe_frq, NiFe_gam, true));
    //transformer.add_susceptibility(core_mat, H_stuff, gyrotropic_susceptibility(vec(0.0,0.0,0.0),2*pi*NiFe_frq, NiFe_gam, true));
    //gyrotropic_susceptibility(const vec &bias, double omega_0, double gamma, double alpha = 0.0,gyrotropy_model model = GYROTROPIC_LORENTZIAN);
    transformer.set_output_directory(mydirname);
    //transformer->set_chi2(mu);
    //transformer->set_chi3(mu);
    
    meep_geom::medium_struct my_medium_struct;
    my_medium_struct.epsilon_diag.x = 1.0;
    my_medium_struct.epsilon_diag.y = 1.0;
    my_medium_struct.epsilon_diag.z = 1.0;
    my_medium_struct.mu_diag.x=1.0;
    my_medium_struct.mu_diag.y=1.0;
    my_medium_struct.mu_diag.z=1.0;
    /*my_medium_struct.mu_offdiag.x=mu_core;
    my_medium_struct.mu_offdiag.y=mu_core;
    my_medium_struct.mu_offdiag.z=mu_core;
    */my_medium_struct.H_chi2_diag.x=1.0;
    my_medium_struct.H_chi2_diag.y=1.0;
    my_medium_struct.H_chi2_diag.z=1.0;
    my_medium_struct.H_chi3_diag.x=1.0;
    my_medium_struct.H_chi3_diag.y=1.0;
    my_medium_struct.H_chi3_diag.z=1.0;

    my_medium_struct.E_chi2_diag.x=1.0;
    my_medium_struct.E_chi2_diag.y=1.0;
    my_medium_struct.E_chi2_diag.z=1.0;
    my_medium_struct.E_chi3_diag.x=1.0;
    my_medium_struct.E_chi3_diag.y=1.0;
    my_medium_struct.E_chi3_diag.z=1.0;
    //m->epsilon_offdiag.x.re = m->epsilon_offdiag.x.im = epsilon_offdiag.y.re = m->epsilon_offdiag.y.im = epsilon_offdiag.z.re = m->epsilon_offdiag.z.im = nn * nn;
    //m->mu_offdiag.x.re = m->mu_offdiag.x.im = mu_offdiag.y.re = m->mu_offdiag.y.im = mu_offdiag.z.re = m->mu_offdiag.z.im = nn * nn;
    //m->E_chi2_diag.x = m->E_chi2_diag.y = m->E_chi2_diag.z = mu_core;
    //m->E_chi3_diag.x = m->E_chi3_diag.y = m->E_chi3_diag.z = mu_core;
    //m->H_chi2_diag.x = m->H_chi2_diag.y = m->H_chi2_diag.z = mu_core;
    //m->H_chi3_diag.x = m->H_chi3_diag.y = m->H_chi3_diag.z = mu_core;
    //m->D_conductivity_diag.x = m->D_conductivity_diag.y = m->D_conductivity_diag.z = nn * nn;
    //m->B_conductivity_diag.x = m->B_conductivity_diag.y = m->B_conductivity_diag.z = 0.0;

    my_medium_struct.H_susceptibilities.num_items = 1;
    my_medium_struct.H_susceptibilities.items = new meep_geom::susceptibility[1];

    my_medium_struct.H_susceptibilities.items[0].sigma_offdiag.x = 0.0;
    my_medium_struct.H_susceptibilities.items[0].sigma_offdiag.y = 0.0;
    my_medium_struct.H_susceptibilities.items[0].sigma_offdiag.z = 0.0;
    my_medium_struct.H_susceptibilities.items[0].sigma_diag.x = gaussian_random(1.0,1.0);
    my_medium_struct.H_susceptibilities.items[0].sigma_diag.y = gaussian_random(1.0,1.0);
    my_medium_struct.H_susceptibilities.items[0].sigma_diag.z = gaussian_random(1.0,1.0);
    my_medium_struct.H_susceptibilities.items[0].bias.x = gaussian_random(1.0,1.0);
    my_medium_struct.H_susceptibilities.items[0].bias.y = gaussian_random(1.0,1.0);
    my_medium_struct.H_susceptibilities.items[0].bias.z = gaussian_random(1.0,1.0);
    my_medium_struct.H_susceptibilities.items[0].frequency = (0.2e6/f0);
    my_medium_struct.H_susceptibilities.items[0].gamma = 100.0*(0.2e6/f0);
    my_medium_struct.H_susceptibilities.items[0].alpha = gaussian_random(1.0,1.0);
    my_medium_struct.H_susceptibilities.items[0].noise_amp = 0.01;
    my_medium_struct.H_susceptibilities.items[0].drude = false;
    my_medium_struct.H_susceptibilities.items[0].saturated_gyrotropy = true;
    my_medium_struct.H_susceptibilities.items[0].is_file = false;
    
    my_material_func_data data;
    data.with_susceptibility = true;


  meep_geom::material_type default_material =meep_geom::make_dielectric(1.0);    
  //default_material->which_subclass = material_data::MEDIUM;
  default_material->medium=my_medium_struct;
  default_material->user_func = my_material_func;
  default_material->user_data = (void *) &data;
  default_material->do_averaging = false;  

    meep_geom::material_type my_user_material =meep_geom::make_user_material(my_material_func, (void *)&data, false);
    //vector3 center = {0, 0, 0};
    //geometric_object go = ctlgeom::geometric_object(my_material,center);
    
    geometric_object objects[5];
  vector3 center = {0.0, 0.0, 0.0};  
  vector3 center1 = {0.0, -1.25, 0.0};
  vector3 center2 = {0.0, 0.0, 0.0};
  vector3 center3 = {0.0, 1.25, 0.0};
  vector3 center4 = {0.0, 0.0, 0.75};
  vector3 center5 = {0.0, 0.0, -0.75};
  double radius = 3.0;
  double height = 1.0e20;
  vector3 xhat1 = {1.0, 0.0, 0.0};
  vector3 yhat1 = {0.0, 1.0, 0.0};
  vector3 zhat1 = {0.0, 0.0, 1.0};
  vector3 size1 = {wcore, wcore, 2.0};
  vector3 size2 = {wcore, wcore, 2.0};
  vector3 size3 = {wcore, wcore, 2.0};
  vector3 size4 = {wcore, 3.0, wcore};
  vector3 size5 = {wcore, 3.0, wcore};
  //objects[0] = make_block(my_material, center, radius, height, zhat);
  objects[0] = make_block(my_user_material, center1, xhat1, yhat1, zhat1, size1);
  objects[1] = make_block(my_user_material, center2, xhat1, yhat1, zhat1, size2);
  objects[2] = make_block(my_user_material, center3, xhat1, yhat1, zhat1, size3);
  objects[3] = make_block(my_user_material, center4, xhat1, yhat1, zhat1, size4);
  objects[4] = make_block(my_user_material, center5, xhat1, yhat1, zhat1, size5);
  geometric_object_list g = {5, objects};
  


    //geometric_object_list g = {0,0};
    //g.num_items=1;
    meep_geom::absorber_list al= new meep_geom::absorber_list_type;

    meep_geom::material_type_list mtl= meep_geom::material_type_list();    
    mtl.num_items=1;
    mtl.items=new meep_geom::material_type;
    mtl.items[0]=my_user_material;
    cout<<"MAIN "<<mtl.num_items<<endl;
    bool use_anisotropic_averaging = false;
    bool ensure_periodicity = true;
    set_materials_from_geometry(&transformer, g, center, use_anisotropic_averaging,DEFAULT_SUBPIXEL_TOL, DEFAULT_SUBPIXEL_MAXEVAL,ensure_periodicity, default_material,al,mtl);

    /*    
    winding_material winding_mat;
    structure winding(gv, winding_mat);
    //winding.set_epsilon(eps);
    winding.set_conductivity(Ex,conductivity);
    winding.set_conductivity(Ey,conductivity);
    winding.set_conductivity(Ez,conductivity);    
    winding.add_susceptibility(winding_mat, E_stuff, lorentzian_susceptibility(2*pi*Cu_frq0, Cu_gam0, true));
    winding.add_susceptibility(winding_mat, E_stuff, lorentzian_susceptibility(2*pi*Cu_frq1, Cu_gam1));
    winding.add_susceptibility(winding_mat, E_stuff, lorentzian_susceptibility(2*pi*Cu_frq2, Cu_gam2));
    winding.add_susceptibility(winding_mat, E_stuff, lorentzian_susceptibility(2*pi*Cu_frq3, Cu_gam3));
    winding.add_susceptibility(winding_mat, E_stuff, lorentzian_susceptibility(2*pi*Cu_frq4, Cu_gam4));
    //lorentzian_susceptibility(double omega_0, double gamma, bool no_omega_0_denominator = false): omega_0(omega_0), gamma(gamma), no_omega_0_denominator(no_omega_0_denominator)
    */      


    fields f(& transformer);
    //fields(structure *, double m = 0, double beta = 0, bool zero_fields_near_cylorigin = true);
  

 for (double fp=0.0;fp<=3.0*(1e9/f0);fp=fp+0.5*(1e9/f0))
    { cout<<fp<<endl;
      for (double yp = -2.0; yp <= 2.0; yp=yp+1.0)
      {
        cout<<f.get_mu(vec(0.0,yp,(dzmin+dzmax)/2.0),fp)-10000<<" ";
      }
      cout<<endl;
    }
    

    double xcenp=xcen; 
    double ycenp=ycen-dymin-(0.5*wcore); 
    double zcenp=zcen;
    double dxminp=dxmax+insulation_thickness_p;
    double dxmaxp=dxminp+winding_thickness_p;
    double dyminp=(0.5*wcore)+insulation_thickness_p; 
    double dymaxp=dyminp+winding_thickness_p;
    double dzminp=0; 
    double dzmaxp=0.5*winding_thickness_p;    
    zcenp=zcen-(0.5*double(Np-1)*insulation_thickness_p)-(0.5*double(Np-1)*winding_thickness_p);

    

    /*for (int i=0;i<Np;i++)
    {
      const volume vsrc1 =volume(vec(xcenp+dxmaxp,ycenp-dymaxp,zcenp+dzmaxp), vec(xcenp+dxminp,ycenp+dymaxp,zcenp-dzmaxp));
      const volume vsrc2 =volume(vec(xcenp+dxmaxp,ycenp+dymaxp,zcenp+dzmaxp), vec(xcenp-dxmaxp,ycenp+dyminp,zcenp-dzmaxp));
      const volume vsrc3 =volume(vec(xcenp-dxmaxp,ycenp+dymaxp,zcenp+dzmaxp), vec(xcenp-dxminp,ycenp-dymaxp,zcenp-dzmaxp));
      const volume vsrc4 =volume(vec(xcenp-dxmaxp,ycenp-dyminp,zcenp+dzmaxp), vec(xcenp+dxmaxp,ycenp-dymaxp,zcenp-dzmaxp));
      f.add_volume_source(Ey, src, vsrc1, cdouble(amplitude,0));
      f.add_volume_source(Ex, src, vsrc2, cdouble(-amplitude,0));
      f.add_volume_source(Ey, src, vsrc3, cdouble(-amplitude,0));
      f.add_volume_source(Ex, src, vsrc4, cdouble(amplitude,0));

      zcenp=zcenp+insulation_thickness_p+winding_thickness_p;     
    }*/

    //f_range;1e-5,1e-2
    double fcen = (1e10)/f0; // ; pulse center frequency
    double df = 0.999999*((1e10)/f0);    // ; df
    //continuous_src_time src(cdouble(fcen,0));
    //gaussian_src_time src(fcen,df);

    for (int i=0;i<numcoord;i++)
    {
    	theta=k*(zpcoord[i]+0.5);
      const volume vsrc1 =volume(vec(xpcoord[i],ypcoord[i],zpcoord[i]), vec(xpcoord[i],ypcoord[i],zpcoord[i]));
      //for (double fr = 100000000.0; fr <= 1000000000.0; fr=fr+100000000.0)
      {
        //continuous_src_time src(cdouble(fr/f0,0.0));
        gaussian_src_time src(fcen,df);
        f.add_volume_source(Ex, src, vsrc1, cdouble(amplitude*cos(theta),0));
        f.add_volume_source(Ey, src, vsrc1, cdouble(-amplitude*sin(theta),0));
      }
    
    }
    //void add_point_source(component c, double freq, double width, double peaktime, double cutoff, const vec &, std::complex<double> amp = 1.0, int is_continuous = 0);
    //void add_volume_source(component c, const src_time &src, const volume &, std::complex<double> amp = 1.0);


    xcenp=xcen; 
    ycenp=ycen-dymin-(0.5*wcore); 
    zcenp=zcen;
    volume box1( vec(dxmax,-dymax,-0.5), vec(-dxmax,dymax,0.0) );
    ycenp=ycen+dymin+(0.5*wcore); 
    volume box2( vec(dxmax,-dymax,0.0), vec(-dxmax,dymax,0.5) );

    fcen = (1e10)/f0; // ; pulse center frequency
    df = 0.001*(fcen/f0);    // ; df
  
    double fmin = 0, fmax = (1e11)/f0;
    int Nfreq = 1000;
    dft_flux flux1 = f.add_dft_flux_box(box1, fmin, fmax, Nfreq);
    dft_flux flux2 = f.add_dft_flux_box(box2, fmin, fmax, Nfreq);
    double init_energy = f.field_energy_in_box(box1);
   


  /*vec lb(vec()), rb(vec());
  vec lt(vec()), rt(vec());
  flux_vol *left = f.add_flux_plane(lb, lt);
  flux_vol *right = f.add_flux_plane(rb, rt);
  flux_vol *bottom = f.add_flux_plane(lb, rb);
  flux_vol *top = f.add_flux_plane(lt, rt);
  long double fluxL = 0;
  */
  //integral of flux = change in energy of box
	f.step();
  f.step();
  f.step();
  f.step();
  f.step();
    volume vxy=volume(vec(-xsize+1,-ysize+1,0),vec(xsize-1,ysize-1,0));
    volume vxz=volume(vec(-xsize+1,0,-zsize+1),vec(xsize-1,0,zsize-1));
    volume vyz=volume(vec(0,-ysize+1,-zsize+1),vec(0,ysize-1,zsize-1));
    
    /*h5file * fEx=f.open_h5file("fEx",h5file::WRITE,0,false);    
    h5file * fEy=f.open_h5file("fEy",h5file::WRITE,0,false);
    h5file * fEz=f.open_h5file("fEz",h5file::WRITE,0,false);
    h5file * fDx=f.open_h5file("fDx",h5file::WRITE,0,false);    
    h5file * fDy=f.open_h5file("fDy",h5file::WRITE,0,false);
    h5file * fDz=f.open_h5file("fDz",h5file::WRITE,0,false);
    h5file * fHx=f.open_h5file("fHx",h5file::WRITE,0,false);    
    h5file * fHy=f.open_h5file("fHy",h5file::WRITE,0,false);
    h5file * fHz=f.open_h5file("fHz",h5file::WRITE,0,false);
    h5file * fBx=f.open_h5file("fBx",h5file::WRITE,0,false);    
    h5file * fBy=f.open_h5file("fBy",h5file::WRITE,0,false);
    h5file * fBz=f.open_h5file("fBz",h5file::WRITE,0,false);
    
    f.output_hdf5(Dielectric,vyz,fEx);
    f.output_hdf5(Permeability,vyz,fEy);
    f.output_hdf5(Ez,vyz,fEz);
    f.output_hdf5(Dx,vyz,fDx);
    f.output_hdf5(Dy,vyz,fDy);
    f.output_hdf5(Dz,vyz,fDz);
    f.output_hdf5(Hx,vyz,fHx);
    f.output_hdf5(Hy,vyz,fHy);
    f.output_hdf5(Hz,vyz,fHz);
    f.output_hdf5(Bx,vyz,fBx);
    f.output_hdf5(By,vyz,fBy);
    f.output_hdf5(Bz,vyz,fBz);*/
    
int stop=0;

    for(int i=1;i<=1000000;i++)
    {
    	if(!stop)
   		   {f.step();}
       if (stop)
       {
         i=1000001;
       }
          //fluxL += f.dt * (left->flux() - right->flux() + bottom->flux() - top->flux());

       
    if ((i%100)==0)
    {
    f.output_hdf5(Hx,vyz);
    f.output_hdf5(Hy,vyz);
    f.output_hdf5(Hz,vyz);
    f.output_hdf5(Bx,vyz);
    f.output_hdf5(By,vyz);
    f.output_hdf5(Bz,vyz);
    f.output_hdf5(Ex,vyz);
    f.output_hdf5(Ey,vyz);
    f.output_hdf5(Ez,vyz);
    f.output_hdf5(Dx,vyz);
    f.output_hdf5(Dy,vyz);
    f.output_hdf5(Dz,vyz);
    f.output_hdf5(Sx,vyz);
    f.output_hdf5(Sy,vyz);
    f.output_hdf5(Sz,vyz);

    }    
       if ((i%100)==0)
      {
        cdouble Vm=compute_Vm(f,ycen);
        cdouble Im=compute_Im(f,ycen);
        cdouble Ve=compute_Ve(f,ycen);
        cdouble Ie=compute_Ie(f,ycen);
        Time<<Im.real()<<" , "<<Im.imag()<<" , "<<Vm.real()<<" , "<<Vm.imag()<<" , "<<Ie.real()<<" , "<<Ie.imag()<<" , "<<Ve.real()<<" , "<<Ve.imag()<<endl;
        //Im , Vm , Ie , Ve 

        monitor_point pin;
        f.get_point(&pin, vec(xcen,ycen,-0.25));
        cdouble E1i = pin.get_component(Ex);
        cdouble E2i = pin.get_component(Ey);
        cdouble E3i = pin.get_component(Ez);
        cdouble D1i = pin.get_component(Dx);
        cdouble D2i = pin.get_component(Dy);
        cdouble D3i = pin.get_component(Dz);
        cdouble H1i = pin.get_component(Hx);
        cdouble H2i = pin.get_component(Hy);
        cdouble H3i = pin.get_component(Hz);
        cdouble B1i = pin.get_component(Bx);
        cdouble B2i = pin.get_component(By);
        cdouble B3i = pin.get_component(Bz);

        monitor_point po;
        f.get_point(&po, vec(xcen,ycen,0.25));
        cdouble E1o = po.get_component(Ex);
        cdouble E2o = po.get_component(Ey);
        cdouble E3o = po.get_component(Ez);
        cdouble D1o = po.get_component(Dx);
        cdouble D2o = po.get_component(Dy);
        cdouble D3o = po.get_component(Dz);
        cdouble H1o = po.get_component(Hx);
        cdouble H2o = po.get_component(Hy);
        cdouble H3o = po.get_component(Hz);
        cdouble B1o = po.get_component(Bx);
        cdouble B2o = po.get_component(By);
        cdouble B3o = po.get_component(Bz);        
        
        FieldsIn<<H1i.real() <<" , "<<H1i.imag()<<" , "<<H2i.real()<<" , "<<H2i.imag()<<" , "<<H3i.real()<<" , "<<H3i.imag()<<" , "<<B1i.real()<<" , "<<B1i.imag()<<" , "<<B2i.real()<<" , "<<B2i.imag()<<" , "<<B3i.real()<<" , "<<B3i.imag()<<" , "<<E1i.real()<<" , "<<E1i.imag()<<" , "<<E2i.real()<<" , "<<E2i.imag()<<" , "<<E3i.real()<<" , "<<E3i.imag()<<" , "<<D1i.real()<<" , "<<D1i.imag()<<" , "<<D2i.real()<<" , "<<D2i.imag()<<" , "<<D3i.real()<<" , "<<D3i.imag()<<endl;
        FieldsOut<<H1o.real()<<" , "<<H1o.imag()<<" , "<<H2o.real()<<" , "<<H2o.imag()<<" , "<<H3o.real()<<" , "<<H3o.imag()<<" , "<<B1o.real()<<" , "<<B1o.imag()<<" , "<<B2o.real()<<" , "<<B2o.imag()<<" , "<<B3o.real()<<" , "<<B3o.imag()<<" , "<<E1o.real()<<" , "<<E1o.imag()<<" , "<<E2o.real()<<" , "<<E2o.imag()<<" , "<<E3o.real()<<" , "<<E3o.imag()<<" , "<<D1o.real()<<" , "<<D1o.imag()<<" , "<<D2o.real()<<" , "<<D2o.imag()<<" , "<<D3o.real()<<" , "<<D3o.imag()<<endl;    
        //Hx , Hy , Hz , Bx , By , Bz , Ex , Ey , Ez , Dx , Dy , Dz 
      }

      if((i%10000)==0) 
      {
      	cout<<"End? (1/0):";
      	cin>>stop;
      }
       
    }

    
    double *fl1 = flux1.flux();
    double *fl2 = flux2.flux();
    cout<<"Flux Harmonics"<<endl;
    for (int i = 0; i < Nfreq; ++i) {
      Fluxes<<(fmin + i * flux1.dfreq)<<" , "<<fl1[i]<<" , "<<fl2[i]<<endl;
      //freq , fluxin , fluxout
    } 
    
    int num_bands=1;
    int * bands=new int [num_bands];
    double * vgrp=new double [num_bands*Nfreq];
    cdouble * coeffs=new cdouble [2*num_bands*Nfreq];
    bands[0]=1;
    for (int i = 1; i <= Nfreq; ++i) {
      vgrp[i]=0.0;
    } 

    //f.get_eigenmode_coefficients(flux1,box1,bands,num_bands,1,divisions,DEFAULT_SUBPIXEL_TOL,coeffs,vgrp);

    cout<<"EigenModes"<<endl;
    for (int i = 0; i < Nfreq; ++i) {
      //cout<<(fmin + i * flux1.dfreq)<<" , "<<bands[i]<<" , "<<coeffs[i]<<endl;
      //freq , fluxin , fluxout
    }

	/*void get_eigenmode_coefficients(dft_flux flux, const volume &eig_vol, int *bands, int num_bands,
                                  int parity, double eig_resolution, double eigensolver_tol,
                                  std::complex<double> *coeffs, double *vgrp,
                                  kpoint_func user_kpoint_func, void *user_kpoint_data,
                                  vec *kpoints, vec *kdom, direction d);
	*/

    
     
    cout<<"Skin Effect"<<endl;
    for (double y=-0.5;y<=0.5;y=y+0.001)  
    {
    	monitor_point pin;
        f.get_point(&pin, vec(xcen,y,zcen));
        cdouble E1i = pin.get_component(Ex);
        cdouble E2i = pin.get_component(Ey);
        cdouble E3i = pin.get_component(Ez);
        cdouble D1i = pin.get_component(Dx);
        cdouble D2i = pin.get_component(Dy);
        cdouble D3i = pin.get_component(Dz);
        cdouble H1i = pin.get_component(Hx);
        cdouble H2i = pin.get_component(Hy);
        cdouble H3i = pin.get_component(Hz);
        cdouble B1i = pin.get_component(Bx);
        cdouble B2i = pin.get_component(By);
        cdouble B3i = pin.get_component(Bz);
        Skin<<H1i.real() <<" , "<<H1i.imag()<<" , "<<H2i.real()<<" , "<<H2i.imag()<<" , "<<H3i.real()<<" , "<<H3i.imag()<<" , "<<B1i.real()<<" , "<<B1i.imag()<<" , "<<B2i.real()<<" , "<<B2i.imag()<<" , "<<B3i.real()<<" , "<<B3i.imag()<<" , "<<E1i.real()<<" , "<<E1i.imag()<<" , "<<E2i.real()<<" , "<<E2i.imag()<<" , "<<E3i.real()<<" , "<<E3i.imag()<<" , "<<D1i.real()<<" , "<<D1i.imag()<<" , "<<D2i.real()<<" , "<<D2i.imag()<<" , "<<D3i.real()<<" , "<<D3i.imag()<<endl;
    }

    cout<<"SpaceEvolution"<<endl;
    //for (double y=(ycen-dymin);y<=(ycen+dymin);y=y+0.001)
    for (double z=-0.5;z<=0.5;z=z+0.01)  
    {
      cdouble Im=compute_Im(f,z);
      cdouble Ie=compute_Ie(f,z);
      cdouble Vm=compute_Vm(f,z);
      cdouble Ve=compute_Ve(f,z);
      Space <<Im.real()<<" , "<<Im.imag()<<" , "<<Vm.real()<<" , "<<Vm.imag()<<" , "<<Ie.real()<<" , "<<Ie.imag()<<" , "<<Ve.real()<<" , "<<Ve.imag()<<endl;
      //Im , Vm , Ie , Ve
    }

/*    f.step();
    f.output_hdf5(Hx,f.gv);
    f.output_hdf5(Hy,f.gv);
    f.output_hdf5(Hz,f.gv);
    f.output_hdf5(Bx,f.gv);
    f.output_hdf5(By,f.gv);
    f.output_hdf5(Bz,f.gv);
    f.output_hdf5(Ex,f.gv);
    f.output_hdf5(Ey,f.gv);
    f.output_hdf5(Ez,f.gv);
    f.output_hdf5(Dx,f.gv);
    f.output_hdf5(Dy,f.gv);
    f.output_hdf5(Dz,f.gv);
*/

    Time.close();
    Space.close();
    FieldsIn.close();
    FieldsOut.close();
    Fluxes.close();
    return 0;
}
