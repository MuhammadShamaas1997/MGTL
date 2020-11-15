/* Copyright (C) 2005-2020 Massachusetts Institute of Technology
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

#include <meep.hpp>
#include <iostream>
#include <fstream>

using namespace meep;
using namespace std;
using std::complex;

typedef std::complex<double> cdouble;

double one(const vec &) { return 1.5; }
double two(const vec &) { return 0.1; }


int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 4;
    std::ofstream FieldsIn;
    FieldsIn.open ("FieldEvolutionIn.txt");
    std::ofstream SourceFFT;
    SourceFFT.open ("SourceFFT.txt");
    std::ofstream Energy;
    Energy.open ("Energy.txt");

  double a = 24.0;
  double ttot = 100.0;

  grid_volume gv = vol3d(3.0,0.0,12.0, a);
  gv.center_origin();
  structure s(gv, one, pml(1.0), identity());
  s.add_susceptibility(two, H_stuff, gyrotropic_susceptibility(vec(0.0,0.0,1.0),1.0, 0.001,0.00001,GYROTROPIC_SATURATED));

  fields f(&s);
  //f.use_real_fields();
  continuous_src_time src(1.0);
  //gaussian_src_time src(1.0,2.0);

for (double fp=0.0;fp<=2.0;fp=fp+0.0001)
    { 
      double yp = 0.0;
      {
        //SourceFFT<<fp<<" "<<src.fourier_transform(fp).real()<<" "<<src.fourier_transform(fp).imag()<<endl;
      }
      
    }

  f.add_point_source(Hx,src,vec(0.0,0.0,-4.5));
    volume vin=volume(vec(-1,-1,-5),vec(1,1,-4));
    volume vout=volume(vec(-1,-1,4),vec(-1,-1,5));
    volume vyz=volume(vec(0,-6,-6),vec(0,6,6));
  int tind=0;
  while (f.time() < ttot) {
	f.step();tind++;
	Energy<<f.magnetic_energy_in_box(vin)<<" , "<<f.magnetic_energy_in_box(vout)<<" , "<<f.magnetic_energy_in_box(gv.surroundings())<<endl;	
	for (double ind=-6.0;ind<=6.0;ind=ind+1)    
		{
    monitor_point pin;
        f.get_point(&pin, vec(0.0,0.0,ind));
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
        
        FieldsIn<<tind<<" , "<<ind<<" , "<<H1i.real() <<" , "<<H1i.imag()<<" , "<<H2i.real()<<" , "<<H2i.imag()<<" , "<<H3i.real()<<" , "<<H3i.imag()<<" , "<<B1i.real()<<" , "<<B1i.imag()<<" , "<<B2i.real()<<" , "<<B2i.imag()<<" , "<<B3i.real()<<" , "<<B3i.imag()<<" , "<<E1i.real()<<" , "<<E1i.imag()<<" , "<<E2i.real()<<" , "<<E2i.imag()<<" , "<<E3i.real()<<" , "<<E3i.imag()<<" , "<<D1i.real()<<" , "<<D1i.imag()<<" , "<<D2i.real()<<" , "<<D2i.imag()<<" , "<<D3i.real()<<" , "<<D3i.imag()<<" , "<<endl;		
}
  //complex<double> E1i = pin.get_component(Ez);  
  //cout<<E1i.real()<<endl;
//f.output_hdf5(Hy,vxz);
//f.output_hdf5(Hx,vxz);

	}
  return 0;
}
