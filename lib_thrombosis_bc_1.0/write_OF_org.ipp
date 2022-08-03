// _write_OF_org.ipp
// supports openfoam-org v70 and higher 

    fvPatchField<Type>::write(os);

    writeEntryIfDifferent<word>(os, "reactive_species", "undefined", _reactive_species);
    writeEntryIfDifferent<word>(os, "reactive_zone", "all", _reactive_zone);

    if( _reactive_aabb.size() > 0 )
    {
        int nbbox = _reactive_aabb.size()/2;
        os.writeKeyword("reactive_aabb");
        os << nl;
        os << tab << "(" << nl;
        for( label i=0; i < nbbox; i++)
        {
          os << tab << "  " << _reactive_aabb[(2*i)+0] << " " << _reactive_aabb[(2*i)+1] << nl;
        }
        os << tab << ");" << nl;
    }

    writeEntryIfDifferent<word>(os, "nonreactive_type", "none", _nonreactive_type);

    if( _adp_injection > 0 )
    {
      writeEntry( os, "adp_injection", _adp_injection );
    }
    if( _nonreactive_type == "refValue" )
    {
      writeEntry(os, "nonreactive_refValue", this->refValue() );
    }
    if( _nonreactive_type == "refGradient" )
    {
      writeEntry(os, "nonreactive_refGradient", this->refGrad() );
    }
    if( _nonreactive_type == "mixed" )
    {
      writeEntry(os, "nonreactive_refValue", this->refValue() );
      writeEntry(os, "nonreactive_refGradient", this->refGrad() );
      writeEntry(os, "nonreactive_valueFraction", this->valueFraction() );
    }

    writeEntry(os, "value", *this);

