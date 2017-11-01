/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/TBaryon.hpp

Copyright (C) 2015

Author: David Preti <david.preti@csic.es>
Author: Antonin Portelli <antonin.portelli@me.com>

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
directory.
*******************************************************************************/

#ifndef Hadrons_MContraction_Baryon_hpp_
#define Hadrons_MContraction_Baryon_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Baryon contractions
 -----------------------------
 
 * options:
 - q1: input propagator 1 (string)
 - q2: input propagator 2 (string)
 - q3: input propagator 3 (string)
 - gammas: gamma products to insert at sink & source, pairs of gamma matrices 
           (space-separated strings) in angled brackets (i.e. <g_sink g_src>),
           in a sequence (e.g. "<Gamma5 Gamma5><Gamma5 GammaT>").

           Special values: "all" - perform all possible contractions.
 - sink: module to compute the sink to use in contraction (string).
*/

/******************************************************************************
 *                               Baryon                                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)


inline RealD Sign(int a, int b){ 
    return (a==b) ? 0:(b-a)/std::abs(b-a);
}


class BaryonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BaryonPar,
        std::string, q1,
        std::string, q2,
        std::string, q3,
        std::string, output);
};

template <typename FImpl1, typename FImpl2, typename FImpl3>
class TBaryon: public Module<BaryonPar>
{
public:
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    FERM_TYPE_ALIASES(FImpl3, 3);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
            std::vector<std::vector<std::vector<Complex>>>, corr);
    };
public:
    // constructor
    TBaryon(const std::string name);
    // destructor
    virtual ~TBaryon(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(Baryon, ARG(TBaryon<FIMPL, FIMPL, FIMPL>), MContraction);

/******************************************************************************
 *                         TBaryon implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
TBaryon<FImpl1, FImpl2, FImpl3>::TBaryon(const std::string name)
: Module<BaryonPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TBaryon<FImpl1, FImpl2, FImpl3>::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2, par().q3};
    
    return input;
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TBaryon<FImpl1, FImpl2, FImpl3>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TBaryon<FImpl1, FImpl2, FImpl3>::execute(void)
{


  
    LOG(Message) << "Computing baryon contractions '" << getName() << "' using"
    << " quarks '" << par().q1 << "', '" << par().q2 << "', and '"
    << par().q3 << "'" << std::endl;
    
    CorrWriter             writer(par().output);
    PropagatorField1      &q1 = *env().template getObject<PropagatorField1>(par().q1);
    PropagatorField2      &q2 = *env().template getObject<PropagatorField2>(par().q2);
    PropagatorField3      &q3 = *env().template getObject<PropagatorField3>(par().q3);
    LatticeComplex        c(env().getGrid()), NCtrtr(env().getGrid()), NCtr(env().getGrid()), LCtrtr(env().getGrid()), LCtr(env().getGrid()); 
//   SpinMatrix            g[Ns*Ns], g5, Id, Ga[5], Gb[5], P[2];
    std::vector<TComplex> buf;
    Result                Nresult, Lresult;
     
    Gamma::Algebra g[] = {
       Gamma::Algebra::GammaX,
       Gamma::Algebra::GammaY,
       Gamma::Algebra::GammaZ,
       Gamma::Algebra::GammaT,
       Gamma::Algebra::Identity,
       Gamma::Algebra::Gamma5
      };


     // std::cout << "adj g5=" << Gamma::adj[g[5]] << std::endl;                      // Example gamma adjoint
     //  std::cout << "g5=" << Gamma::mul[Gamma::mul[g[5]][g[5]]][g[5]] << std::endl; //Example gammas multiplication
      

// Ga and Gb combinations from Gattringer & Lang pag.130 Eq.(6.19) (and below).
    Gamma::Algebra Ga[] = {
      //spin 1/2
      g[4],                         // Id
      g[5],                         // g5 
      g[4],                         // Id
      //spin 3/2
      g[4],                         // Id
      g[4],                         // Id
      g[4]                          // Id
    };

    Gamma::Algebra Gb[] = {  //TODO
      //spin 1/2
      Gamma::mul[Gamma::mul[g[1]][g[3]]][g[5]],                    // Cg5
      Gamma::mul[g[1]][g[3]],                                      // C
      Gamma::mul[Gamma::mul[Gamma::mul[g[3]][g[1]]][g[3]]][g[5]],  // igTCg5 -> timesI missing!!!
      //spin 3/2
      Gamma::mul[Gamma::mul[g[1]][g[3]]][g[0]],                    // CgX
      Gamma::mul[Gamma::mul[g[1]][g[3]]][g[1]],                    // CgY
      Gamma::mul[Gamma::mul[g[1]][g[3]]][g[2]],                    // CgZ
    };
   
    Gamma::Algebra P[] = {
      g[4],                                                       // Id
      g[4]                                                        // Id
   //   0.5*(g[4] + g[3]),  //Positive Parity Projector           // sum not implemented 
   //   0.5*(g[4] - g[3])   //Negative Parity Projector           // sum not implemented
    };

    RealD four = 4.;
    RealD two = four/2.;
    RealD norm = 1./6.;

  //Prepare antisymmetric tensor epsilon_{ijk} 
  int ep[3][Nc];
  int em[3][Nc];
  int contp=0,contm=0;
  int S1,S2,S3,S4,S5,S6;

  for (int i=0;i<Nc;i++){
    for (int j=0;j<Nc;j++){
      for (int k=0;k<Nc;k++){
        S1=Sign(i,j);
        S2=Sign(j,k);
        S3=Sign(i,k);
        if ((S1*S2*S3)>0){
          ep[contp][0]=i;
          ep[contp][1]=j;
          ep[contp][2]=k;
          contp++;
        } else if ((S1*S2*S3)<0){
          em[contm][0]=i;
          em[contm][1]=j;
          em[contm][2]=k;
          contm++;
        }
      }
    }
  }


// Contractions
LatticeSpinMatrix q1_c(env().getGrid()); q1_c = zero;
LatticeSpinMatrix q2_c(env().getGrid()); q2_c = zero;
LatticeSpinMatrix q3_c(env().getGrid()); q3_c = zero;

Nresult.corr.resize(5);
Lresult.corr.resize(5);

  for (unsigned int iParity=0;iParity<2;iParity++){  // loop over parity + and - 
    for (unsigned int iSrc=0;iSrc<6;iSrc++){         // loop over the 6 possible dirac structures to get spin 1/2 (3 contractions) and spin 3/2 (3 contractions).
      NCtrtr=zero;  // Needed for Nucleon contraction
      NCtr=zero;    // Needed for Nucleon contraction 
      LCtrtr=zero;  // Needed for Lambda contraction
      LCtr=zero;    // Needed for Lambda contraction
 //loop over non-zero permutations     
      for (unsigned int i=0;i<3;i++){
        for (unsigned int j=0;j<3;j++){
    //positive colours permutations                                                                                                                                           
          
    /*   // Uncomment this once the transpose is available for gamma matrices!!!!!!!!!!
    
          q3_c = peekColour(q3,ep[i][1],ep[j][1]);
          q1_c = peekColour(q1,ep[i][0],ep[j][0]);
          q2_c = peekColour(q2,ep[i][2],ep[j][2]);
          
          NCtrtr += trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c ) * trace( transpose(Gamma(Gb[iSrc])*q3_c*transpose(Gamma(Gb[iSrc]))) * q2_c );
          NCtr += trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c * transpose(Gamma(Gb[iSrc])*q3_c*transpose(Gamma(Gb[iSrc]))) * q2_c );

          LCtrtr += four * trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q3_c ) * trace( transpose(Gamma(Gb[iSrc])*q2_c*transpose(Gamma(Gb[iSrc]))) * q1_c );
          LCtrtr += trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q2_c ) * trace( transpose(Gamma(Gb[iSrc])*q3_c*transpose(Gamma(Gb[iSrc]))) * q1_c );
          LCtrtr += trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c ) * trace( transpose(Gamma(Gb[iSrc])*q3_c*transpose(Gamma(Gb[iSrc]))) * q2_c );
          LCtr -= two*trace( q3_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c * Gamma(Gb[iSrc])*transpose(q2_c)*Gamma(Gb[iSrc]) );
          LCtr += two*trace( q2_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q3_c * transpose(Gamma(Gb[iSrc]))*transpose(q1_c)*Gamma(Gb[iSrc]) );
          LCtr -= two*trace( q1_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q3_c * transpose(Gamma(Gb[iSrc]))*transpose(q2_c)*transpose(Gamma(Gb[iSrc])) );
          LCtr += two*trace( q3_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q2_c * transpose(Gamma(Gb[iSrc]))*transpose(q1_c)*Gamma(Gb[iSrc]) );
          LCtr -= trace( q1_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q2_c * Gamma(Gb[iSrc])*transpose(q3_c)*transpose(Gamma(Gb[iSrc])) );
          LCtr -= trace( q2_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c * Gamma(Gb[iSrc])*transpose(q3_c)*transpose(Gamma(Gb[iSrc])) );
          q3_c = zero; q1_c = zero; q2_c = zero;
          
          q3_c = peekColour(q3,em[i][1],em[j][1]);
          q1_c = peekColour(q1,em[i][0],em[j][0]);
          q2_c = peekColour(q2,em[i][2],em[j][2]);

          NCtrtr += trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c ) * trace( transpose(Gamma(Gb[iSrc])*q3_c*transpose(Gamma(Gb[iSrc]))) * q2_c );
          NCtr += trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c * transpose(Gamma(Gb[iSrc])*q3_c*transpose(Gamma(Gb[iSrc]))) * q2_c );
          
          LCtrtr += four * trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q3_c ) * trace( transpose(Gamma(Gb[iSrc])*q2_c*transpose(Gamma(Gb[iSrc]))) * q1_c );
          LCtrtr += trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q2_c ) * trace( transpose(Gamma(Gb[iSrc])*q3_c*transpose(Gamma(Gb[iSrc]))) * q1_c );
          LCtrtr += trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c ) * trace( transpose(Gamma(Gb[iSrc])*q3_c*transpose(Gamma(Gb[iSrc]))) * q2_c );
          LCtr -= two*trace( q3_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c * Gamma(Gb[iSrc])*transpose(q2_c)*Gamma(Gb[iSrc]) );
          LCtr += two*trace( q2_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q3_c * transpose(Gamma(Gb[iSrc]))*transpose(q1_c)*Gamma(Gb[iSrc]) );
          LCtr -= two*trace( q1_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q3_c * transpose(Gamma(Gb[iSrc]))*transpose(q2_c)*transpose(Gamma(Gb[iSrc])) );
          LCtr += two*trace( q3_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q2_c * transpose(Gamma(Gb[iSrc]))*transpose(q1_c)*Gamma(Gb[iSrc]) );
          LCtr -= trace( q1_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q2_c * Gamma(Gb[iSrc])*transpose(q3_c)*transpose(Gamma(Gb[iSrc])) );
          LCtr -= trace( q2_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c * Gamma(Gb[iSrc])*transpose(q3_c)*transpose(Gamma(Gb[iSrc])) );
          q3_c = zero; q1_c = zero; q2_c = zero;

    //negative colours permutations                                                                                                                                           
          q3_c = peekColour(q3,ep[i][1],em[j][1]);
          q1_c = peekColour(q1,ep[i][0],em[j][0]);
          q2_c = peekColour(q2,ep[i][2],em[j][2]);
          
          NCtrtr -= trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c ) * trace( transpose(Gamma(Gb[iSrc])*q3_c*transpose(Gamma(Gb[iSrc]))) * q2_c );
          NCtr -= trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c * transpose(Gamma(Gb[iSrc])*q3_c*transpose(Gamma(Gb[iSrc]))) * q2_c );
          
          LCtrtr -= four * trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q3_c ) * trace( transpose(Gamma(Gb[iSrc])*q2_c*transpose(Gamma(Gb[iSrc]))) * q1_c );
          LCtrtr -= trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q2_c ) * trace( transpose(Gamma(Gb[iSrc])*q3_c*transpose(Gamma(Gb[iSrc]))) * q1_c );
          LCtrtr -= trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c ) * trace( transpose(Gamma(Gb[iSrc])*q3_c*transpose(Gamma(Gb[iSrc]))) * q2_c );
          LCtr += two*trace( q3_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c * Gamma(Gb[iSrc])*transpose(q2_c)*Gamma(Gb[iSrc]) );
          LCtr -= two*trace( q2_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q3_c * transpose(Gamma(Gb[iSrc]))*transpose(q1_c)*Gamma(Gb[iSrc]) );
          LCtr += two*trace( q1_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q3_c * transpose(Gamma(Gb[iSrc]))*transpose(q2_c)*transpose(Gamma(Gb[iSrc])) );
          LCtr -= two*trace( q3_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q2_c * transpose(Gamma(Gb[iSrc]))*transpose(q1_c)*Gamma(Gb[iSrc]) );
          LCtr += trace( q1_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q2_c * Gamma(Gb[iSrc])*transpose(q3_c)*transpose(Gamma(Gb[iSrc])) );
          LCtr += trace( q2_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c * Gamma(Gb[iSrc])*transpose(q3_c)*transpose(Gamma(Gb[iSrc])) );
          q3_c = zero; q1_c = zero; q2_c = zero;

          q3_c = peekColour(q3,em[i][1],ep[j][1]);
          q1_c = peekColour(q1,em[i][0],ep[j][0]);
          q2_c = peekColour(q2,em[i][2],ep[j][2]);

          NCtrtr -= trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c ) * trace( transpose(Gamma(Gb[iSrc])*q3_c*transpose(Gamma(Gb[iSrc]))) * q2_c );
          NCtr -= trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c * transpose(Gamma(Gb[iSrc])*q3_c*transpose(Gamma(Gb[iSrc]))) * q2_c );
        
          LCtrtr -= four * trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q3_c ) * trace( transpose(Gamma(Gb[iSrc])*q2_c*transpose(Gamma(Gb[iSrc]))) * q1_c );
          LCtrtr -= trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q2_c ) * trace( transpose(Gamma(Gb[iSrc])*q3_c*transpose(Gamma(Gb[iSrc]))) * q1_c );
          LCtrtr -= trace( (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c ) * trace( transpose(Gamma(Gb[iSrc])*q3_c*transpose(Gamma(Gb[iSrc]))) * q2_c );       
          LCtr += two*trace( q3_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c * Gamma(Gb[iSrc])*transpose(q2_c)*Gamma(Gb[iSrc]) );
          LCtr -= two*trace( q2_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q3_c * transpose(Gamma(Gb[iSrc]))*transpose(q1_c)*Gamma(Gb[iSrc]) );
          LCtr += two*trace( q1_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q3_c * transpose(Gamma(Gb[iSrc]))*transpose(q2_c)*transpose(Gamma(Gb[iSrc])) );
          LCtr -= two*trace( q3_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q2_c * transpose(Gamma(Gb[iSrc]))*transpose(q1_c)*Gamma(Gb[iSrc]) );        
          LCtr += trace( q1_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q2_c * Gamma(Gb[iSrc])*transpose(q3_c)*transpose(Gamma(Gb[iSrc])) );
          LCtr += trace( q2_c * (Gamma(Ga[iSrc])*Gamma(P[iParity])*Gamma(Ga[iSrc])) * q1_c * Gamma(Gb[iSrc])*transpose(q3_c)*transpose(Gamma(Gb[iSrc])) );
          q3_c = zero; q1_c = zero; q2_c = zero;
          */


        }
      }
      c = (NCtrtr + NCtr);
      sliceSum(c, buf, Tp);
      Nresult.corr[iSrc].resize(2);
      Nresult.corr[iSrc][iParity].resize(buf.size());

      c = norm * (LCtrtr + LCtr); 
      sliceSum(c, buf, Tp);
      Lresult.corr[iSrc].resize(2);
      Lresult.corr[iSrc][iParity].resize(buf.size());

      for (unsigned int t=0; t<buf.size(); t++){
        Nresult.corr[iSrc][iParity][t] =  TensorRemove(buf[t]);
        Lresult.corr[iSrc][iParity][t] =  TensorRemove(buf[t]);
      }
    }
  }

    //----------------------------
 
  write(writer, "baryon", Nresult);
  write(writer, "baryon", Lresult);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_Baryon_hpp_
