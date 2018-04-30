/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_WilsonGaugeSF.cc

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>



int main(int argc, char **argv) {
  using namespace Grid;
  using namespace Grid::QCD;

  Grid_init(&argc, &argv);
  GridLogLayout();

   // Typedefs to simplify notation
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
  HMCWrapper TheHMC;

  // Grid from the command line
  std::vector<int> vSIMD({2,1,1,1});
  TheHMC.Resources.AddFourDimGrid("gauge",vSIMD);
  // Possibile to create the module by hand 
  // hardcoding parameters or using a Reader


  // Checkpointer definition
  CheckpointerParameters CPparams;  
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = 10000000;
  CPparams.format = "IEEE64BIG";
  
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);




  std::vector<RealD> beta={9.1544,9.1544};
  RealD ct=1.0 - 6./beta[0] * 0.089; 
  RealD cs=1.0;


  //////////////////////////////////////////////
  //ONLINE MEASUREMENTS:
  
  //Plaquette
  typedef PlaquetteSFMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();         

  //SF COUPLING
  typedef dSdetaMod<HMCWrapper::ImplPolicy> dSdeta;
  dSdetaParameters uSFParms;
  uSFParms.interval = 1;
  uSFParms.betaT = beta[1];  
  uSFParms.ct_SF = ct;   
  TheHMC.Resources.AddObservable<dSdeta>(uSFParms);                
  //////////////////////////////////////////////



  //ACTION: 
  WilsonGaugeSFActionR Waction(beta,cs,ct);
  //WilsonGaugeActionR Waction(beta[0]);

  Waction.isSF = true;
  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&Waction);
  TheHMC.TheAction.push_back(Level1);
  /////////////////////////////////////////////////////////////

  // HMC parameters are serialisable 
  TheHMC.Parameters.MD.MDsteps = 20;
  TheHMC.Parameters.MD.trajL   = 1.0;

  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  TheHMC.Run();  // no smearing

  Grid_finalize();

} // main
