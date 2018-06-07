    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_SFGaugeAction.cc

    Copyright (C) 2015

    Author: David Preti <david.preti@to.infn.it>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace QCD;


int main (int argc, char **argv)
{
 Grid_init(&argc, &argv);

  std::vector<int> latt_size = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  std::vector<int> mpi_layout = GridDefaultMpi();
  GridCartesian Grid(latt_size, simd_layout, mpi_layout);
  GridRedBlackCartesian RBGrid(&Grid);

  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  std::vector<int> seeds({1, 2, 3, 4});
  GridParallelRNG pRNG(&Grid);
  pRNG.SeedFixedIntegers(seeds);

  LatticeGaugeField Umu(&Grid);
  std::vector<LatticeColourMatrix> U(Nd, &Grid);



int STARTCNFG=1000;
int ENDCNFG=1200;
/*
for (int cnfg=STARTCNFG;cnfg<ENDCNFG;cnfg++){
  FieldMetaData header;
  std::string file("./inner_ckpoint_lat." + std::to_string(cnfg));
  std::cout << GridLogMessage << "Reading: " << file << std::endl;
  NerscIO::readConfiguration(Umu,header,file);
*/
  SU3::HotConfiguration(pRNG, Umu);

  for(int mu=0;mu<Nd;mu++){
     U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }


int T   = Umu._grid->GlobalDimensions()[3];
int X   = Umu._grid->GlobalDimensions()[0];
int Y   = Umu._grid->GlobalDimensions()[1];
int Z   = Umu._grid->GlobalDimensions()[2];
int Tmax=T-1;
std::cout << GridLogMessage << "XxYxZxT="<<X<<"x"<<Y<<"x"<<Z<<"x"<<T<<std::endl;



// ##########################################################
// # Here I define the Gauge obs to extract Glueball Masses #
// ##########################################################

std::vector<LatticeColourMatrix> Wl(Nd-1, &Grid), Pl(Nd-1, &Grid);
std::vector<TComplex> Plsum(T), Wlsum(T);
int Nmu=4;
int Mnu=4;
int i=0;
double Plave, Wlave;
RealD N=1./(X*Y*Z*Nc*(Nd-1));
std::vector<RealD> Plres(T),Wlres(T);


for (int r=1;r<=Nmu;r++){
  for (int s=1;s<=Mnu;s++){

i=0;

for (int nu=1;nu<Nd-1;nu++){
  for (int mu=0;mu<nu;mu++){
    WilsonLoops<PeriodicGimplR>::WilsonLoopNxM(Wl[i], U, mu, nu, r, s);
    WilsonLoops<PeriodicGimplR>::dirPlaquette(Pl[i], U, mu, nu);
    i++;
  }
}


for (int j=0;j<Nd-1;j++){
  sliceSum( trace(Pl[j]), Plsum, 3);
  sliceSum( trace(Wl[j]), Wlsum, 3);
  for (int t=0;t<T;t++){
    Plres[t] += N*TensorRemove(Plsum[t]).real();
    Wlres[t] += N*TensorRemove(Wlsum[t]).real();  
  }
}

std::cout << "t Wl(" << r << "," << s << ")" << std::endl;
Plave=0.0;
Wlave=0.0;
for (int t=0; t<T ; t++){
  std::cout << t << " " << Plres[t] << " " << Wlres[t] << std::endl; 
  Plave+=Plres[t];
  Wlave+=Wlres[t];
  Plres[t]=0.0;
  Wlres[t]=0.0;
}
std::cout << "sum:" << " " << double(Wlave/T) << std::endl; 

  }
}



  Grid_finalize();
}



