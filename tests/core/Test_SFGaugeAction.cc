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
  SU3::ColdConfiguration(pRNG, Umu);
  std::vector<LatticeColourMatrix> U(Nd, &Grid);

  //FieldMetaData header;
  //std::string file("./ckpoint_lat.4000");
  //NerscIO::readConfiguration(Umu,header,file);

#define NONABELIAN_SF
int T   = Umu._grid->GlobalDimensions()[3];
int X   = Umu._grid->GlobalDimensions()[0];
int Y   = Umu._grid->GlobalDimensions()[1];
int Z   = Umu._grid->GlobalDimensions()[2];
std::cout << GridLogMessage << "XxYxZxT="<<X<<"x"<<Y<<"x"<<Z<<"x"<<T<<std::endl;

assert((X==Y)&&(Y==Z));
int Tmax=T-1;

#ifdef ABELIAN_SF

//SF boundary implementation:
std::vector<double> phi_SF(3),phiprime_SF(3),omega(3);
double eta=0., nuSF=0.;

omega[0] = 1.;
omega[1] = -1./2. + nuSF;
omega[2] = -1./2. - nuSF;

phi_SF[0] = eta*omega[0] - M_PI/3.;
phi_SF[1] = eta*omega[1];
phi_SF[2] = eta*omega[2] + M_PI/3.;

phiprime_SF[0] = -phi_SF[0] - 4./3. * M_PI;
phiprime_SF[1] = -phi_SF[2] + 2./3. * M_PI;
phiprime_SF[2] = -phi_SF[1] + 2./3. * M_PI;

LatticeColourMatrix W_bc(Umu._grid),Wprime_bc(Umu._grid);

W_bc=zero;
Wprime_bc=zero;


for (int i = 0; i < Umu._grid->oSites(); i++){
  W_bc._odata[i]()()(0, 0) = exp(ComplexD(0.0, 1./X)*phi_SF[0]);
  W_bc._odata[i]()()(1, 1) = exp(ComplexD(0.0, 1./X)*phi_SF[1]);
  W_bc._odata[i]()()(2, 2) = exp(ComplexD(0.0, 1./X)*phi_SF[2]);

  Wprime_bc._odata[i]()()(0, 0) = exp(ComplexD(0.0, 1./X)*phiprime_SF[0]);
  Wprime_bc._odata[i]()()(1, 1) = exp(ComplexD(0.0, 1./X)*phiprime_SF[1]);
  Wprime_bc._odata[i]()()(2, 2) = exp(ComplexD(0.0, 1./X)*phiprime_SF[2]);      
}
    
  for (int mu=0;mu<Nd;mu++){
    U[mu] = peekLorentz(Umu,mu);
  }
  Lattice<iScalar<vInteger>> coor(U[3]._grid);
  LatticeCoordinate(coor, 3);  
  U[3] = where(coor==(Tmax), 0.*U[3], U[3]);
  pokeLorentz(Umu, U[3], 3);
  LatticeCoordinate(coor, 3); 
  
  for (int mu=0;mu<Nd-1;mu++){
    U[mu] = where(coor==0, W_bc, U[mu]);
    U[mu] = where(coor==(Tmax), Wprime_bc, U[mu]);
    pokeLorentz(Umu, U[mu], mu);
  }
#endif
 

#ifdef NONABELIAN_SF
//SF non-abelian boundary implementation:
LatticeGaugeField Umu_bc(&Grid);      //input boundary cnfg
SU3::ColdConfiguration(pRNG,Umu_bc);   //this is something coming from an higher lvl hmc
std::vector<LatticeColourMatrix> Ubc(Nd, &Grid);

  for (int mu=0;mu<Nd;mu++){
    U[mu]   = peekLorentz(Umu,mu);
    Ubc[mu] = peekLorentz(Umu_bc,mu);
  }
  Lattice<iScalar<vInteger>> coor(U[3]._grid);
  LatticeCoordinate(coor, 3);  
  U[3] = where(coor==(Tmax), 0.*U[3], U[3]); //<-- HACK
  pokeLorentz(Umu, U[3], 3);
 
  LatticeCoordinate(coor, 3);
  for (int mu=0;mu<Nd-1;mu++){
    U[mu] = where(coor==0, Ubc[mu], U[mu]);
    U[mu] = where(coor==(Tmax), Ubc[mu], U[mu]);
    pokeLorentz(Umu, U[mu], mu);
  }
  
 //std::cout << Umu << std::endl;
#endif

//-----------
//     TEST Action: 

//NOW ANISOTROPIC GAUGE ACTION
std::vector<double> beta={6.0,6.0};
RealD ct=1.0,cs=1.0;

WilsonGaugeSFActionR Action(beta,cs,ct);
std::cout<< Action.LogParameters() << std::endl;
RealD S;
S=Action.S(Umu);
std::cout << GridLogMessage << "S_Wilson_SF= " << S << std::endl;
RealD theta = 1./3. * M_PI * 1./X/X;
RealD norm = 12.*(X*X)*( std::sin(theta) + std::sin(2.*theta) ); 
RealD uSF, dSde;
dSde = ColourWilsonLoops::dSdeta(Umu, beta[1], ct);
uSF= norm/dSde;
std::cout << GridLogMessage << "k   : " << norm << std::endl;
std::cout << GridLogMessage << "dS/deta : " << dSde << std::endl;
std::cout << GridLogMessage << "uSF : " << uSF << std::endl;


RealD plaq;
plaq = ColourWilsonLoops::avgPlaquetteSF(Umu);
std::cout << GridLogMessage << "PlaqSF :" << plaq << std::endl;

//WILSON PLAQUETTE ACTION
//WilsonGaugeActionR Action2(1.0);
//std::cout<< Action2.LogParameters() << std::endl;
//RealD S2;
//S2=Action2.S(Umu);
//std::cout << GridLogMessage << "S_Wilson_periodic= " << S2 << std::endl;



//-----------
//     TEST dSdU: 

// ANISOTROPIC ACTION
//LatticeGaugeField Umu2(&Grid);      //input boundary cnfg
//SU3::HotConfigurationSF(pRNG,Umu2);   //this is something coming from an higher lvl hmc

LatticeGaugeField dSdU(&Grid);   
Action.deriv(Umu, dSdU); 

//std::cout << dSdU << std::endl;

// WILSON PLAQUETTE ACTION
//LatticeGaugeField dSdU2(&Grid);   
//Action2.deriv(Umu, dSdU2); 

LatticeColourMatrix dSdU_mu(&Grid), dSdU2_mu(&Grid);
LatticeColourMatrix diff(&Grid);
for (int i=0;i<Nd;i++){
dSdU_mu   = peekLorentz(dSdU,i);
//dSdU2_mu  = peekLorentz(dSdU2,i);
//std::cout << "dSdU_SF[" << i <<"]:" << dSdU_mu << std::endl;
//std::cout << "dSdU_Wilson[" << i <<"]:" << dSdU2_mu << std::endl;
//diff=dSdU_mu-dSdU2_mu;
//std::cout << GridLogMessage << "DIFF[" << i << "]= " << diff << std::endl;
}

  ////////////////////////////////////
  // Modify the gauge field a little 
  ////////////////////////////////////
  RealD dt = 0.0001;

  LatticeColourMatrix mommu(&Grid); 
  LatticeColourMatrix forcemu(&Grid); 
  LatticeGaugeField mom(&Grid); 
  LatticeGaugeField Uprime(&Grid); 


//fix boundaries in mom
  for(int mu=0;mu<Nd;mu++){
    SU3::GaussianFundamentalLieAlgebraMatrix(pRNG, mommu); // Traceless antihermitian momentum; gaussian in lie alg
    LatticeCoordinate(coor, 3);
    if (mu!=Nd) mommu = where(coor==0 || coor==T-1, 0.*mommu, mommu);
    if (mu==Nd) mommu = where(coor==T-1, 0.*mommu, mommu);
    PokeIndex<LorentzIndex>(mom,mommu,mu);
    // fourth order exponential approx
    parallel_for(auto i=mom.begin();i<mom.end();i++){ // exp(pmu dt) * Umu
      Uprime[i](mu) = Umu[i](mu) + mom[i](mu)*Umu[i](mu)*dt ;
    }
  }
  ComplexD Sprime    = Action.S(Uprime);

  //////////////////////////////////////////////
  // Use derivative to estimate dS
  //////////////////////////////////////////////
  LatticeGaugeField UdSdU(&Grid);
  Action.deriv(Umu,UdSdU);
  LatticeComplex dS(&Grid); dS = zero;

  for(int mu=0;mu<Nd;mu++){
    auto UdSdUmu = PeekIndex<LorentzIndex>(UdSdU,mu);
         mommu   = PeekIndex<LorentzIndex>(mom,mu);
    // Update gauge action density
    // U = exp(p dt) U
    // dU/dt = p U
    // so dSdt = trace( dUdt dSdU) = trace( p UdSdUmu ) 
    dS = dS - trace(mommu*UdSdUmu)*dt*2.0;
   // std::cout<< UdSdUmu << std::endl;
  }
  ComplexD dSpred    = sum(dS);

  std::cout << GridLogMessage << " S      "<<S<<std::endl;
  std::cout << GridLogMessage << " Sprime "<<Sprime<<std::endl;
  std::cout << GridLogMessage << "dS      "<<Sprime-S<<std::endl;
  std::cout << GridLogMessage << "pred dS "<< dSpred <<std::endl;

 // assert( fabs(real(Sprime-S-dSpred)) < 1.0e-2 ) ;

  std::cout<< GridLogMessage << "Done" <<std::endl;

  Grid_finalize();
}



