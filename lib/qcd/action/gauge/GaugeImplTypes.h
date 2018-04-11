/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/gauge/GaugeImpl.h

Copyright (C) 2015

Author: paboyle <paboyle@ph.ed.ac.uk>

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
#ifndef GRID_GAUGE_IMPL_TYPES_H
#define GRID_GAUGE_IMPL_TYPES_H

namespace Grid {
namespace QCD {

////////////////////////////////////////////////////////////////////////
// Implementation dependent gauge types
////////////////////////////////////////////////////////////////////////

#define INHERIT_GIMPL_TYPES(GImpl)                  \
  typedef typename GImpl::Simd Simd;                \
  typedef typename GImpl::LinkField GaugeLinkField; \
  typedef typename GImpl::Field GaugeField;         \
  typedef typename GImpl::ComplexField ComplexField;\
  typedef typename GImpl::SiteField SiteGaugeField; \
  typedef typename GImpl::SiteComplex SiteComplex;  \
  typedef typename GImpl::SiteLink SiteGaugeLink;

#define INHERIT_FIELD_TYPES(Impl)		    \
  typedef typename Impl::Simd Simd;		    \
  typedef typename Impl::ComplexField ComplexField; \
  typedef typename Impl::SiteField SiteField;	    \
  typedef typename Impl::Field Field;

// hardcodes the exponential approximation in the template
template <class S, int Nrepresentation = Nc, int Nexp = 12 > class GaugeImplTypes {
public:
  typedef S Simd;

  template <typename vtype> using iImplScalar     = iScalar<iScalar<iScalar<vtype> > >;
  template <typename vtype> using iImplGaugeLink  = iScalar<iScalar<iMatrix<vtype, Nrepresentation> > >;
  template <typename vtype> using iImplGaugeField = iVector<iScalar<iMatrix<vtype, Nrepresentation> >, Nd>;

  typedef iImplScalar<Simd>     SiteComplex;
  typedef iImplGaugeLink<Simd>  SiteLink;
  typedef iImplGaugeField<Simd> SiteField;

  typedef Lattice<SiteComplex> ComplexField;
  typedef Lattice<SiteLink>    LinkField; 
  typedef Lattice<SiteField>   Field;

  // Guido: we can probably separate the types from the HMC functions
  // this will create 2 kind of implementations
  // probably confusing the users
  // Now keeping only one class


  // Move this elsewhere? FIXME
  static inline void AddLink(Field &U, LinkField &W,
                                  int mu) { // U[mu] += W
    PARALLEL_FOR_LOOP
    for (auto ss = 0; ss < U._grid->oSites(); ss++) {
      U._odata[ss]._internal[mu] =
          U._odata[ss]._internal[mu] + W._odata[ss]._internal;
    }
  }

  ///////////////////////////////////////////////////////////
  // Move these to another class
  // HMC auxiliary functions
  static inline void generate_momenta(Field &P, GridParallelRNG &pRNG) {
    // specific for SU gauge fields
    LinkField Pmu(P._grid);
    Pmu = zero;
    for (int mu = 0; mu < Nd; mu++) {
      SU<Nrepresentation>::GaussianFundamentalLieAlgebraMatrix(pRNG, Pmu);
      PokeIndex<LorentzIndex>(P, Pmu, mu);
    }
  }

  static inline void generate_momentaSF(Field &P, GridParallelRNG &pRNG) {
    // specific P0_{mu}(x,t) for SF where momenta @ t=0,T-1 are identically 0. 
    LinkField Pmu(P._grid);
    Pmu = zero;
    Lattice<iScalar<vInteger>> coor(Pmu._grid);
    LatticeCoordinate(coor, Nd-1);
    int T = Pmu._grid->GlobalDimensions()[Nd-1];
    for (int mu = 0; mu < Nd; mu++) {
      SU<Nrepresentation>::GaussianFundamentalLieAlgebraMatrix(pRNG, Pmu);
        if (mu!=Nd-1) Pmu = where((coor==0 || coor==T-1), 0.*Pmu, Pmu); 
        if (mu==Nd-1) Pmu = where((coor==T-1), 0.*Pmu, Pmu); 
      PokeIndex<LorentzIndex>(P, Pmu, mu);
    }
  }

  static inline Field projectForce(Field &P) { return Ta(P); }

  static inline void update_field(Field& P, Field& U, double ep){
    parallel_for(int ss=0;ss<P._grid->oSites();ss++){
      for (int mu = 0; mu < Nd; mu++) 
        U[ss]._internal[mu] = ProjectOnGroup(Exponentiate(P[ss]._internal[mu], ep, Nexp) * U[ss]._internal[mu]); 
    }
  }

  static inline void update_fieldSF(Field& P, Field& U, double ep){
    //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl; 
    //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl; 
    //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl; 
    //std::cout << "Uold=" << U << std::endl; 
  LatticeColourMatrix Umu(U._grid), Pmu(P._grid);
  int lvol = U._grid->lSites();
  std::vector<int> lcoor;
  ColourMatrix Umux, Pmux;
  int T   = U._grid->GlobalDimensions()[3];

  for (int mu = 0; mu < Nd; mu++){ 
      Umu = peekLorentz(U,mu);
      Pmu = peekLorentz(P,mu);
    for(int site=0;site<lvol;site++){
      U._grid->LocalIndexToLocalCoor(site, lcoor);
//      std::cout << lcoor << " mu= " << mu << std::endl;
      peekLocalSite(Umux, Umu, lcoor);
      peekLocalSite(Pmux, Pmu, lcoor);  
      if (mu!=3 && (lcoor[3]!=0 && lcoor[3]!=T-1)){            
        Umux = ProjectOnGroup(Exponentiate(Pmux, ep, Nexp) * Umux);    
      }
      if (mu==3 && lcoor[3]!=T-1){        
        Umux = ProjectOnGroup(Exponentiate(Pmux, ep, Nexp) * Umux);
      }
      pokeLocalSite(Umux, Umu, lcoor); 
      }
      pokeLorentz(U, Umu, mu);
    }
    //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl; 
    //std::cout << "MOMENTI=" << P << std::endl; 
    //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl; 
    //std::cout << "Unew=" << U << std::endl; 
  }


  static inline RealD FieldSquareNorm(Field& U){
    LatticeComplex Hloc(U._grid);
    Hloc = zero;
    for (int mu = 0; mu < Nd; mu++) {
      auto Umu = PeekIndex<LorentzIndex>(U, mu);
      Hloc += trace(Umu * Umu);
    }
    Complex Hsum = sum(Hloc);
    return Hsum.real();
  }

  static inline void HotConfiguration(GridParallelRNG &pRNG, Field &U) {
    SU<Nc>::HotConfiguration(pRNG, U);
  }

  static inline void TepidConfiguration(GridParallelRNG &pRNG, Field &U) {
    SU<Nc>::TepidConfiguration(pRNG, U);
  }

  static inline void ColdConfiguration(GridParallelRNG &pRNG, Field &U) {
    SU<Nc>::ColdConfiguration(pRNG, U);
  }

  static inline void HotConfigurationSF(GridParallelRNG &pRNG, Field &U) {
    SU<Nc>::HotConfigurationSF(pRNG, U);
  }

  static inline void TepidConfigurationSF(GridParallelRNG &pRNG, Field &U) {
    SU<Nc>::TepidConfigurationSF(pRNG, U);
  }

  static inline void ColdConfigurationSF(GridParallelRNG &pRNG, Field &U) {
    SU<Nc>::ColdConfigurationSF(pRNG, U);
  }

    static inline void HotConfigurationNonAbelianSF(GridParallelRNG &pRNG, Field &U, Field &U_bc) {
    SU<Nc>::HotConfigurationNonAbelianSF(pRNG, U, U_bc);
  }

  static inline void TepidConfigurationNonAbelianSF(GridParallelRNG &pRNG, Field &U, Field &U_bc) {
    SU<Nc>::TepidConfigurationNonAbelianSF(pRNG, U, U_bc);
  }

  static inline void ColdConfigurationNonAbelianSF(GridParallelRNG &pRNG, Field &U, Field &U_bc) {
    SU<Nc>::ColdConfigurationNonAbelianSF(pRNG, U, U_bc);
  }

};


typedef GaugeImplTypes<vComplex, Nc> GimplTypesR;
typedef GaugeImplTypes<vComplexF, Nc> GimplTypesF;
typedef GaugeImplTypes<vComplexD, Nc> GimplTypesD;

typedef GaugeImplTypes<vComplex, SU<Nc>::AdjointDimension> GimplAdjointTypesR;
typedef GaugeImplTypes<vComplexF, SU<Nc>::AdjointDimension> GimplAdjointTypesF;
typedef GaugeImplTypes<vComplexD, SU<Nc>::AdjointDimension> GimplAdjointTypesD;


} // QCD
} // Grid

#endif // GRID_GAUGE_IMPL_TYPES_H
