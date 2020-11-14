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

double one(const vec &) { return 1.5; }
double two(const vec &) { return 0.1; }


int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 4;
    std::ofstream FieldsIn;
    FieldsIn.open ("FieldEvolutionIn.txt");

  double a = 24.0;
  double ttot = 100.0;

  grid_volume gv = vol3d(0.0,0.0,12.0, a);
  gv.center_origin();
  structure s(gv, one, no_pml(), identity());
  //s.add_susceptibility(two, H_stuff, gyrotropic_susceptibility(vec(0.0,0.0,1.0),1.0, 0.001,0.00001,GYROTROPIC_SATURATED));

  fields f(&s);
  f.use_real_fields();
  continuous_src_time src(0.8);
  f.add_point_source(Hx,src,vec(0.0,0.0,-4.5));
    volume vxy=volume(vec(-6,-6,0),vec(6,6,0));
    volume vxz=volume(vec(-1,0,-6),vec(1,0,6));
    volume vyz=volume(vec(0,-6,-6),vec(0,6,6));
  int tind=0;
  while (f.time() < ttot) {
	f.step();tind++;	
	for (double ind=-6.0;ind<=6.0;ind++)    
		{
  		monitor_point pin;
  		f.get_point(&pin, vec(0.0,0.0,ind));  
  		complex<double> Hxi = pin.get_component(Hx);
  		complex<double> Hyi = pin.get_component(Hy);
  		FieldsIn<<tind<<" , "<<ind<<" , "<<Hxi.real() <<" , "<<Hxi.imag()<<" , "<<Hyi.real()<<" , "<<Hyi.imag()<<endl;
		}
  //complex<double> E1i = pin.get_component(Ez);  
  //cout<<E1i.real()<<endl;
//f.output_hdf5(Hy,vxz);
//f.output_hdf5(Hx,vxz);

	}
  return 0;
}
