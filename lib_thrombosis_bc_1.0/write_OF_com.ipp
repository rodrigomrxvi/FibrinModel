// _write_OF_com.ipp
// supports openfoam-com v1906 and higher

    fvPatchField<Type>::write(os);

    os.writeEntryIfDifferent<word>("reactive_species", "undefined", _reactive_species);
    os.writeEntryIfDifferent<word>("reactive_zone", "all", _reactive_zone);

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

    os.writeEntryIfDifferent<word>("nonreactive_type", "none", _nonreactive_type);

    if( _adp_injection > 0 )
    {
      os.writeEntry( "adp_injection", _adp_injection );
    }
    if( _nonreactive_type == "refValue" )
    {
      this->refValue().writeEntry("nonreactive_refValue", os);
    }
    if( _nonreactive_type == "refGradient" )
    {
      this->refGrad().writeEntry("nonreactive_refGradient", os);
    }
    if( _nonreactive_type == "mixed" )
    {
      this->refValue().writeEntry("nonreactive_refValue", os);
      this->refGrad().writeEntry("nonreactive_refGradient", os);
      this->valueFraction().writeEntry("nonreactive_valueFraction", os);
    }

    this->writeEntry("value", os);

