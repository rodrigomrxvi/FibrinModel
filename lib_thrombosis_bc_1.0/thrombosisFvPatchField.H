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

Class
    Foam::thrombosisFvPatchField

Description
    Sets reactive boundary conditions for the Cornell thrombosis model.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type  thrombosis;      // required
        reactive_species   RP; // required
        value uniform 0;       // required
        reactive_zone   inside_aabb; // optional, default=all
        reactive_aabb // optional
        (
          (150e-6 -1e6 -1e6) (170e-6 1e6 1e6) // min/max of axis-aligned bounding box (aabb)
        );
        adp_injection 1e3; // optional & only for apr equation, default = 0
        nonreactive_type  mixed; // optional, default = none
        nonreactive_refValue  uniform 0; // optional
        nonreactive_refGradient  uniform 0; // optional
        nonreactive_valueFraction  uniform 0; // optional
    }
    \endverbatim

SourceFiles
    thrombosisFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef thrombosisFvPatchField_H
#define thrombosisFvPatchField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class thrombosisFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class thrombosisFvPatchField
:
    public mixedFvPatchField<Type>
{
protected:

    // Private data

    word         _reactive_species; // RP|TB|etc.

    word         _reactive_zone;  // optional: all|inside_aabb|outside_aabb (default:all)

    List<point>  _reactive_aabb; // list of min/max aabb points

    word         _nonreactive_type;  // optional: none|refValue|refGradient|mixed (default:none)

    boolList     _reactive_flag; // sets whether face is reactive or not

    scalar       _adp_injection; // optionally used only for apr equation

    bool         _is_init; 

public:

    //- Runtime type information
    TypeName("thrombosis");


    // Constructors

        //- Construct from patch and internal field
        thrombosisFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        thrombosisFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given thrombosisFvPatchField
        //  onto a new patch
        thrombosisFvPatchField
        (
            const thrombosisFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        thrombosisFvPatchField
        (
            const thrombosisFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new thrombosisFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        thrombosisFvPatchField
        (
            const thrombosisFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new thrombosisFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

        // Evaluation functions
       
            virtual void updateCoeffs();

        //- Write
            virtual void write(Ostream&) const;

private: 

        void  _init();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "thrombosisFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
