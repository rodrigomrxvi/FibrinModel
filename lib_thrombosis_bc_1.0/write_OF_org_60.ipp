// _write_OF_org_60.ipp
// supports openfoam-org v60 

    fvPatchField<Type>::write(os);

    this->writeEntryIfDifferent( os, "reactive_species", word("undefined"), _reactive_species );
    this->writeEntryIfDifferent(os, "reactive_zone", word("all"), _reactive_zone);

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

    this->writeEntryIfDifferent(os, "nonreactive_type", word("none"), _nonreactive_type);

    if( _adp_injection > 0 )
    {
      os.writeKeyword( word("adp_injection") );
      os <<" "<< _adp_injection << token::END_STATEMENT << endl;
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
