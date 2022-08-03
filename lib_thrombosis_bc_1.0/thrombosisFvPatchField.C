/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"
#include "localEulerDdtScheme.H"

#include "thrombosisFvPatchField.H"

#define TMPL template<class Type>
#define CLASS Foam::thrombosisFvPatchField<Type>

#define DUMPHERE Info << "DUMPHERE " << __FILE__ <<": "<< __LINE__ << nl;
#define DUMP(x)  Info << "DUMP " << #x <<": "<< x << nl;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TMPL
//------------------------------------------------------------
CLASS::thrombosisFvPatchField
//------------------------------------------------------------
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    _reactive_species("undefined"),
    _reactive_zone("all"),
    _nonreactive_type("none"),
    _adp_injection(0),
    _is_init(false)
{
    _reactive_aabb        = List<point>();
    this->refValue()      = Zero;
    this->refGrad()       = Zero;
    this->valueFraction() = 0.0;
}


TMPL
//------------------------------------------------------------
CLASS::thrombosisFvPatchField
//------------------------------------------------------------
(
    const thrombosisFvPatchField& pF,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(pF, p, iF, mapper),
    _reactive_species(pF._reactive_species),
    _reactive_zone(pF._reactive_zone),
    _reactive_aabb(pF._reactive_aabb),
    _nonreactive_type(pF._nonreactive_type),
    _adp_injection(pF._adp_injection),
    _is_init(false)
{}


TMPL
//------------------------------------------------------------
CLASS::thrombosisFvPatchField
//------------------------------------------------------------
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p,iF),
    _reactive_species(dict.lookupOrDefault<word>("reactive_species", "undefined")),
    _reactive_zone(dict.lookupOrDefault<word>("reactive_zone", "all")),
    _nonreactive_type(dict.lookupOrDefault<word>("nonreactive_type", "none")),
    _adp_injection(dict.lookupOrDefault<scalar>("adp_injection", 0.0)),
    _is_init(false)
{
    _reactive_aabb = List<point>();

    if( dict.found("reactive_aabb") )
    {
        dict.lookup("reactive_aabb") >> _reactive_aabb;
    }

    if( dict.found("value") )
    {
      fvPatchField<Type>::operator=( Field<Type>("value", dict, p.size() ) );
    }
    else
    {
      //fvPatchField<Type>::operator=( this->patchInternalField() );
      FatalIOErrorInFunction
      (
        dict
      ) << "Issue:" << nl
        << " - On patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file: " << this->internalField().objectPath() << nl
        << " - BC value field is not defined."
        << exit(FatalIOError);
    }

    this->refValue()      = *this;
    this->refGrad()       = Zero;
    this->valueFraction() = 0.0;

    if( _nonreactive_type == "refValue" )
    {
      fvPatchField<Type>::operator= ( Field<Type>("nonreactive_refValue", dict, p.size()) );
      this->refGrad()       = Zero;
      this->valueFraction() = 1.0;
    }
    if( _nonreactive_type == "refGradient" )
    {
      this->refGrad() = Field<Type>("nonreactive_refGradient", dict, p.size());
      this->valueFraction() = 0.0;
    }    
    if( _nonreactive_type == "mixed" )
    {
      fvPatchField<Type>::operator= ( Field<Type>("nonreactive_refValue", dict, p.size()) );
      this->refGrad() = Field<Type>("nonreactive_refGradient", dict, p.size());
      this->valueFraction() = scalarField("nonreactive_valueFraction", dict, p.size());
    }
    if( _nonreactive_type == "none" ) // set default values
    {
      bool is_zero_gradient = false;
      bool is_zero_value    = false;
      bool is_zero_mixed    = false;
  
      bool handled = 0;
      if( _reactive_species == "AP"  ) { is_zero_mixed = true; handled=1; }
      if( _reactive_species == "RP"  ) { is_zero_mixed = true; handled=1; }
      if( _reactive_species == "II"  ) { is_zero_mixed = true; handled=1; }
      if( _reactive_species == "AT"  ) { is_zero_gradient = true; handled=1; }
      if( _reactive_species == "IIa"  ) { is_zero_gradient = true; handled=1; }
      if( _reactive_species == "apr" ) { is_zero_gradient = true; handled=1; }
      if( _reactive_species == "aps" ) { is_zero_gradient = true; handled=1; }
      if( _reactive_species == "APd" ) { is_zero_value = true; handled=1; }
      if( _reactive_species == "RPd" ) { is_zero_value = true; handled=1; }
      if( _reactive_species == "APs" ) { is_zero_value = true; handled=1; }
      if( _reactive_species == "XII" ) { is_zero_gradient = true; handled=1; }
      if( _reactive_species == "XIIa" ) { is_zero_gradient = true; handled=1; }
      if( _reactive_species == "Fg" ) { is_zero_gradient = true; handled=1; }
      if( _reactive_species == "Fgd" ) { is_zero_value = true; handled=1; }
      if( _reactive_species == "Fn" ) { is_zero_gradient = true; handled=1; }
      if( _reactive_species == "Fnd" ) { is_zero_value = true; handled=1; }

      if( is_zero_mixed ) this->valueFraction() = 0.0;
      if( is_zero_mixed ) this->refValue()      = Zero;
      if( is_zero_mixed ) this->refGrad()       = Zero;
  
      if( is_zero_gradient ) this->valueFraction() = 0.0;
      if( is_zero_gradient ) this->refGrad()       = Zero;

      if( is_zero_value    ) this->valueFraction() = 1.0;
      if( is_zero_value    ) this->refValue()      = Zero;

      _nonreactive_type = "mixed";

      if( ! handled )
      FatalIOErrorInFunction
      (
        dict
      ) << "Issue:" << nl
        << " - On patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << ", nonreactive_type is defined as: " << _nonreactive_type << nl
        << " - Since reactive_zone is not 'all', nonreactive_type must be set to either: refValue, refGradient, or mixed. "
        << exit(FatalIOError);
    }

    /* debug
    //label patchi = this->patch().index();
    _init();
    auto field = this->data();
    forAll( _reactive_flag, i )
    {
      //Info << i <<" "<< this->operator[](patchi) << endl;
      Info << i <<" "<< this->operator[](i) << endl;
      Info << i <<". "<< field[i] << endl;
      //if( _reactive_flag[i] ) field[i] = this->value[i];
    }
    */
}


TMPL
//------------------------------------------------------------
CLASS::thrombosisFvPatchField
//------------------------------------------------------------
(
    const thrombosisFvPatchField& pF
)
:
    mixedFvPatchField<Type>(pF),
    _reactive_species(pF._reactive_species),
    _reactive_zone(pF._reactive_zone),
    _reactive_aabb(pF._reactive_aabb),
    _nonreactive_type(pF._nonreactive_type),
    _adp_injection(pF._adp_injection),
    _is_init(false)
{}


TMPL
//------------------------------------------------------------
CLASS::thrombosisFvPatchField
//------------------------------------------------------------
(
    const thrombosisFvPatchField& pF,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(pF,iF),
    _reactive_species(pF._reactive_species),
    _reactive_zone(pF._reactive_zone),
    _reactive_aabb(pF._reactive_aabb),
    _nonreactive_type(pF._nonreactive_type),
    _adp_injection(pF._adp_injection),
    _is_init(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

TMPL void
//---------------------------------------------------
CLASS::_init() 
//---------------------------------------------------
{  
  // 1. flag which patch faces are reactive

  bool init_flag = true; if( _reactive_zone == "inside_aabb" ) init_flag = false;

  // burgreen fails for v70 _reactive_flag  = boolField( this->patch().size(), init_flag );
  _reactive_flag  = Field<bool>( this->patch().size(), init_flag );

  int num_aabb = _reactive_aabb.size()/2;

  const vectorField& fc = this->patch().Cf();

  forAll( fc, i )
  {
    bool inside = true;
    for( label bb=0; bb < num_aabb; bb++)
    {
      int bb0 = (2*bb)+0;
      int bb1 = (2*bb)+1;
      if( fc[i][0] < _reactive_aabb[bb0][0] ) inside = false;
      if( fc[i][0] > _reactive_aabb[bb1][0] ) inside = false;
      if( fc[i][1] < _reactive_aabb[bb0][1] ) inside = false;
      if( fc[i][1] > _reactive_aabb[bb1][1] ) inside = false;
      if( fc[i][2] < _reactive_aabb[bb0][2] ) inside = false;
      if( fc[i][2] > _reactive_aabb[bb1][2] ) inside = false;
    }
    if( _reactive_zone == "inside_aabb"  && inside ) _reactive_flag[i] = true;
    if( _reactive_zone == "outside_aabb" && inside ) _reactive_flag[i] = false;
  }

  // 2. finales

  _is_init = true;
}

TMPL void
//---------------------------------------------------
CLASS::updateCoeffs() 
//---------------------------------------------------
{
//DUMP( _reactive_species )
  if( ! _is_init ) _init();

//DUMPHERE
  label patchi = this->patch().index();
//DUMP( patchi )
//Info << " - On patch " << this->patch().name() << nl;

  //const vectorField& delta = this->patch().delta().cref();
  const scalarField& deltaCoeff = this->patch().deltaCoeffs();

  Field<Type> field; this->patchInternalField( field );

  Field<Type>& refValue      = this->refValue();
  Field<Type>& refGrad       = this->refGrad();
  scalarField& valueFraction = this->valueFraction();

  #define RETRIEVE(v) \
  const volScalarField& vol_##v = this->db().objectRegistry::template lookupObject<volScalarField>(#v); \
  const auto& v = vol_##v.boundaryField()[patchi];

  RETRIEVE(DpltIn)
  RETRIEVE(DIIIn)
  RETRIEVE(DaprIn)
  RETRIEVE(DapsIn)
  RETRIEVE(DIIaIn)
  RETRIEVE(DXIIIn)
  RETRIEVE(SwIn)
  RETRIEVE(krpdBIn)
  RETRIEVE(kapdBIn)
  RETRIEVE(kapa)
  RETRIEVE(kspa)
  RETRIEVE(ksXIIIn)
  RETRIEVE(phiAIn)
  RETRIEVE(phiRIn)
  //RETRIEVE(betaThrIn)
  RETRIEVE(lambdajIn)
  RETRIEVE(spjIn)
  RETRIEVE(fstbIn)
  RETRIEVE(thetaIn)
  RETRIEVE(fembB)
  RETRIEVE(dTime)
  RETRIEVE(APd)
  RETRIEVE(RPd)
  RETRIEVE(AP)
  RETRIEVE(APs)
  RETRIEVE(RP)
  RETRIEVE(II)
  RETRIEVE(IIa)
  RETRIEVE(XII)
  RETRIEVE(kFgdBIn)
  RETRIEVE(kFndBIn)
  RETRIEVE(DFnIn)
  RETRIEVE(DFgIn)
  RETRIEVE(Fg)
  RETRIEVE(Fgd)
  RETRIEVE(Fn)
  RETRIEVE(Fnd)
  RETRIEVE(kcat6In)
  RETRIEVE(km6In)
  //RETRIEVE(XIIa)
  //RETRIEVE(apr)
  //RETRIEVE(aps)
  //RETRIEVE(AT)

  //scalar smalln = 1e-3;
  scalar smalln = 1e-20;
  scalar eta_efect = 0.05;

  forAll( _reactive_flag, i )
  {
    if( _reactive_flag[i] )
    {
      scalar delta = 1./deltaCoeff[i];
      if( _reactive_species == "AP" )
      {
        //valueExpression "0.0";
        //gradientExpression "(xp>InjP1 && xp<InjP2) ? fembB*APd/DpltIn : 0.0";
        //fractionExpression "(xp>InjP1 && xp<InjP2) ? 1.0/(1.0 + DpltIn/(mag(delta())*(SwIn+smalln)*kapdBIn)) : 0.0";

        refGrad[i] = fembB[i]*APd[i]/DpltIn[i];
        valueFraction[i] =  1.0/(1.0 + DpltIn[i]/(delta*(SwIn[i]+smalln)*kapdBIn[i]));
      }
      if( _reactive_species == "RP" )
      {
        //valueExpression "0.0";
        //gradientExpression "0.0";
        //fractionExpression "(xp>InjP1 && xp<InjP2) ? 1.0/(1.0 + DpltIn/(mag(delta())*(SwIn+smalln)*krpdBIn)) : 0.0";

        valueFraction[i] = 1.0/(1.0 + DpltIn[i]/(delta*(SwIn[i]+smalln)*krpdBIn[i]));
      }
      if( _reactive_species == "II" )
      {
       // valueFraction[i] = 1.0/(1.0 + DIIIn[i]/(delta*betaThrIn[i]*(phiAIn[i]*APd[i]+phiRIn[i]*RPd[i]+smalln)));
        valueFraction[i] = 1.0/(1.0 + DIIIn[i]/(delta*(phiAIn[i]*APd[i]+phiRIn[i]*RPd[i]+smalln)));
      }
      if( _reactive_species == "APd" ) 
      {
        //valueExpression	"(xp>InjP1 && xp<InjP2) ? (APd + (SwIn*kapdBIn*AP+thetaIn*SwIn*krpdBIn*RP+kapa*RPd+kspa*RPd-(fstbIn+fembB)*APd)*dTime) : 0.0";
        //gradientExpression "0";
        //fractionExpression "1.0";

        scalar value = SwIn[i]*kapdBIn[i]*AP[i] 
                     + thetaIn[i]*SwIn[i]*krpdBIn[i]*RP[i] 
                     + kapa[i]*RPd[i]
                     + kspa[i]*RPd[i]
                     - (fstbIn[i]+fembB[i])*APd[i];
        refValue[i] = APd[i] + value*dTime[i];
      }
      if( _reactive_species == "RPd" ) 
      {
        //valueExpression	"(xp>InjP1 && xp<InjP2) ? (RPd + ((1-thetaIn)*SwIn*krpdBIn*RP-kapa*RPd-kspa*RPd-fembB*RPd)*dTime) : 0.0";
        //gradientExpression "0";
        //fractionExpression "1.0";

        scalar value = (1.0-thetaIn[i])*SwIn[i]*krpdBIn[i]*RP[i]
                     - kapa[i]*RPd[i]
                     - kspa[i]*RPd[i]
                     - fembB[i]*RPd[i];
        refValue[i] = RPd[i] + value*dTime[i];
      }      
      if( _reactive_species == "AT" ) 
      {
        scalar grad = 0.0;
        //refValue[i] = field[i] + grad/deltaCoeff[i];
        refGrad[i] = grad;
      }
      if( _reactive_species == "IIa" ) 
      {
	//valueExpression	"0.0";
	//gradientExpression "(xp>InjP1 && xp<InjP2) ? (phiAIn*APd+phiRIn*RPd)*II/DtbIn : 0.0";
	//fractionExpression "0";

        scalar grad = ( phiAIn[i]*APd[i] + phiRIn[i]*RPd[i] ) * II[i] / DIIaIn[i];
        //refValue[i] = field[i] + grad/deltaCoeff[i];
        refGrad[i] = grad;
      }
      if( _reactive_species == "apr" ) 
      {
        //valueExpression	"0.0";
        //gradientExpression "(xp>InjP1 && xp<InjP2) ? (apdInj+lambdajIn*(kapa*RPd + kspa*RPd + thetaIn*SwIn*krpdBIn*RP))/DaprIn : 0.0";
        //fractionExpression "0";

        scalar grad = _adp_injection
                    + lambdajIn[i]*( kapa[i]*RPd[i] )
                    + lambdajIn[i]*( kspa[i]*RPd[i] )
                    + lambdajIn[i]*( thetaIn[i]*SwIn[i]*krpdBIn[i]*RP[i] );
        grad /= DaprIn[i];
        refGrad[i] = grad;
      }
      if( _reactive_species == "aps" ) 
      {
        //valueExpression "0.0";
        //gradientExpression "(xp>InjP1 && xp<InjP2) ? APd*spjIn/DapsIn : 0.0";
        //fractionExpression "0";

        scalar grad = APd[i]*spjIn[i]/DapsIn[i];
        refGrad[i] = grad;
      }
      if( _reactive_species == "APs" ) 
      {
        //valueExpression  "(xp>InjP1 && xp<InjP2) ? (APs+(fstbIn*APd)*dTime) : 0.0";
        //gradientExpression "0";
        //fractionExpression "1.0";
        scalar value = fstbIn[i]*APd[i];
        refValue[i] = APs[i] + value*dTime[i];
      }
      if( _reactive_species == "XII" ) 
      {
        valueFraction[i] = 1.0-(1.0/(1.0 + delta*ksXIIIn[i]/DXIIIn[i]));
      }
      if( _reactive_species == "XIIa" ) 
      {
        scalar grad = ksXIIIn[i]*XII[i]/DXIIIn[i];
        refGrad[i] = grad;
      }
      if( _reactive_species == "Fg" )
      {
        refGrad[i] = fembB[i]*Fgd[i]/DFgIn[i];
        valueFraction[i] =  1.0/(1.0 + DFgIn[i]/(delta*(SwIn[i]+smalln)*kFgdBIn[i]));
      }

      if( _reactive_species == "Fgd" )
      {

        scalar value = -fembB[i]*Fgd[i]+kFgdBIn[i]*Fg[i]*SwIn[i]-kcat6In[i]*Fgd[i]*IIa[i]/(km6In[i]+Fgd[i])*eta_efect;
        refValue[i] = Fgd[i] + value*dTime[i];
      }
      if( _reactive_species == "Fn" )
      {
        refGrad[i] = fembB[i]*Fnd[i]/DFnIn[i];
        valueFraction[i] =  1.0/(1.0 + DFnIn[i]/(delta*(SwIn[i]+smalln)*kFndBIn[i]));
      }

      if( _reactive_species == "Fnd" )
      {
        scalar value = -fembB[i]*Fnd[i]+kFndBIn[i]*Fn[i]*SwIn[i]+kcat6In[i]*Fgd[i]*IIa[i]/(km6In[i]+Fgd[i])*eta_efect;
        refValue[i] = Fnd[i] + value*dTime[i];
      }
    }
  }

  mixedFvPatchField<Type>::updateCoeffs();
}

TMPL void
//---------------------------------------------------
CLASS::write(Ostream& os) const
//---------------------------------------------------
{
  #if defined (OPENFOAM_ORG_60) 
    #include "write_OF_org_60.ipp"
  #endif

  #if defined (OPENFOAM_ORG) 
    #include "write_OF_org.ipp"
  #endif

  #if defined (OPENFOAM_COM) 
    #include "write_OF_com.ipp"
  #endif
}

#undef TMPL 
#undef CLASS 

// ************************************************************************* //
