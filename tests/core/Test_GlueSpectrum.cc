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

int STARTCNFG=50;
int ENDCNFG=10010;

for (int cnfg=STARTCNFG;cnfg<ENDCNFG;cnfg++){
  FieldMetaData header;
  std::string file("./ckpoint_lat." + std::to_string(cnfg));
  std::cout << GridLogMessage << "Reading: " << file << std::endl;
  NerscIO::readConfiguration(Umu,header,file);

  //SU3::ColdConfiguration(pRNG, Umu);

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



// ##########################################################
//   <TrRe[Rect(t0)]TrRe[Rect(t)]>
// ##########################################################
std::vector<LatticeColourMatrix> SpaceRect(Nd-1, &Grid);
int t0=1;
int i=0;

for (int nu=1;nu<Nd-1;nu++){
  for (int mu=0;mu<nu;mu++){
    WilsonLoops<PeriodicGimplR>::dirRectangle(SpaceRect[i], U, mu, nu);
    i++;
  }
}
std::vector<TComplex> SpaceRectSum(T);
std::vector<RealD> Rt(T);
RealD N=0.5/(X*Y*Z*Nc*(Nd-1));  //0.5 is beacuse is Rectangular
for (int j=0;j<Nd-1;j++){
  sliceSum( trace(SpaceRect[j]), SpaceRectSum, 3);
  for (int t=0;t<T;t++){
    Rt[t] += N*TensorRemove(SpaceRectSum[t]).real();
  }
}
std::cout << "Rectangular Wilson Loop Correlator" << std::endl;
for (int t=0;t<T;t++){
  std::cout << t << " " << (Rt[t0]*Rt[t]) << std::endl;
}

// ##########################################################
std::cout << "#############################################" << std::endl;
std::cout << "#############################################" << std::endl;
// ##########################################################

// ##########################################################
//   <TrRe[Plaq(t0)]TrRe[Plaq(t)]>
// ##########################################################
std::vector<LatticeColourMatrix> SpacePlaq0pp(Nd-1, &Grid);
t0=0;
i=0;


//0++
for (int nu=1;nu<Nd-1;nu++){
  for (int mu=0;mu<nu;mu++){
    WilsonLoops<PeriodicGimplR>::dirPlaquette(SpacePlaq0pp[i], U, mu, nu);
    i++;
  }
}
std::vector<TComplex> SpacePlaqSum0pp(T);
std::vector<RealD> Pt0pp(T);
for (int t=0;t<T;t++) Pt0pp[t]=0.0;

N=1.0/(Nc*X*Y*Z); //  /Nd-1/(X*Y*Z) 
for (int j=0;j<Nd-1;j++){
  sliceSum( trace(SpacePlaq0pp[j]), SpacePlaqSum0pp, 3);
  for (int t=0;t<T;t++){
    Pt0pp[t] += N*TensorRemove(SpacePlaqSum0pp[t]).real();
  }
}

//2++
std::vector<LatticeColourMatrix> SpacePlaq2pp(2, &Grid);
i=0;
int mu=0;

for (int nu=1;nu<Nd-1;nu++){
    WilsonLoops<PeriodicGimplR>::dirPlaquette(SpacePlaq2pp[i], U, mu, nu);
    i++;
  }

std::vector<TComplex> SpacePlaqSum2pp(T);
std::vector<RealD> Pt2pp(T);
for (int t=0;t<T;t++) Pt2pp[t]=0.0;

N=1.0/(Nc*X*T*Z);  // /(X*Y*Z)/2.0

for (int j=0;j<2;j++){
  sliceSum( trace(SpacePlaq2pp[j]), SpacePlaqSum2pp, 3);
  for (int t=0;t<T;t++){
    if (j==0){
    Pt2pp[t] += N*TensorRemove(SpacePlaqSum2pp[t]).real();
    } else if (j==1){
    Pt2pp[t] -= N*TensorRemove(SpacePlaqSum2pp[t]).real();
    }
  }
}

std::cout << "Fundamental Wilson Loop Correlator" << std::endl;
for (int t=0;t<T;t++){
  std::cout << t << " " << (Pt0pp[t0]*Pt0pp[t]) << " " << Pt0pp[t0] << " " << Pt0pp[t] << " " << (Pt2pp[t0]*Pt2pp[t]) << " " << Pt2pp[t0] << " " << Pt2pp[t] << std::endl;
  //(Rt[t0]*Rt[t]) << " " << Rt[t0] << " " << Rt[t] <<  std::endl;
}



// ##########################################################
//   <TrRe[Plaq2x2(t0)]TrRe[Plaq2x2(t)]>
// ##########################################################
std::vector<LatticeColourMatrix> SpacePlaq2x20pp(Nd-1, &Grid);
t0=0;
i=0;


//0++
for (int nu=1;nu<Nd-1;nu++){
  for (int mu=0;mu<nu;mu++){
    WilsonLoops<PeriodicGimplR>::dirPlaquette2x2(SpacePlaq2x20pp[i], U, mu, nu);
    i++;
  }
}
std::vector<TComplex> SpacePlaq2x2Sum0pp(T);
std::vector<RealD> P2x2t0pp(T);
for (int t=0;t<T;t++) P2x2t0pp[t]=0.0;

N=1.0/(Nc*X*Y*Z); //  /Nd-1/(X*Y*Z) 
for (int j=0;j<Nd-1;j++){
  sliceSum( trace(SpacePlaq2x20pp[j]), SpacePlaq2x2Sum0pp, 3);
  for (int t=0;t<T;t++){
    P2x2t0pp[t] += N*TensorRemove(SpacePlaqSum0pp[t]).real();
  }
}


//2++
std::vector<LatticeColourMatrix> SpacePlaq2x22pp(2, &Grid);
i=0;
mu=0;

for (int nu=1;nu<Nd-1;nu++){
    WilsonLoops<PeriodicGimplR>::dirPlaquette2x2(SpacePlaq2x22pp[i], U, mu, nu);
    i++;
  }

std::vector<TComplex> SpacePlaq2x2Sum2pp(T);
std::vector<RealD> P2x2t2pp(T);
for (int t=0;t<T;t++) P2x2t2pp[t]=0.0;

N=1.0/(Nc*X*Y*Z);  // /(X*Y*Z)/2.0

for (int j=0;j<2;j++){
  sliceSum( trace(SpacePlaq2x22pp[j]), SpacePlaq2x2Sum2pp, 3);
  for (int t=0;t<T;t++){
    if (j==0){
    P2x2t2pp[t] += N*TensorRemove(SpacePlaq2x2Sum2pp[t]).real();
    } else if (j==1){
    P2x2t2pp[t] -= N*TensorRemove(SpacePlaq2x2Sum2pp[t]).real();
    }
  }
}

std::cout << "2x2 Wilson Loop Correlator" << std::endl;
for (int t=0;t<T;t++){
//  std::cout << t << " " << (P2x2t0pp[t0]*P2x2t0pp[t]) << " " << P2x2t0pp[t0] << " " << P2x2t0pp[t] << " " << (P2x2t2pp[t0]*P2x2t2pp[t]) << " " << P2x2t2pp[t0] << " " << P2x2t2pp[t] << std::endl;
}


} //cnfg loop
  Grid_finalize();
}



