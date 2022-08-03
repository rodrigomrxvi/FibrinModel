
---------------------------------------------------------------
2020.07.03

A reactive wall thrombosis bc that is independent of swak4Foam 

Author: greg.burgreen@gmail.com

Supports the WeiTaoWu-2016 thrombosis model, namely, the 10 species:

  reactive_species == "AP"
  reactive_species == "PT"
  reactive_species == "RP"
  reactive_species == "TB"
  reactive_species == "apr"
  reactive_species == "aps"
  reactive_species == "AT"
  reactive_species == "RPd"
  reactive_species == "APd"
  reactive_species == "APs"

---------------------------------------------------------------
Example of the boundary condition specification:

  <patchName>
  {
    type  thrombosis;      // required
    reactive_species   RP; // required
    value uniform 0;       // required
    reactive_zone   inside_aabb; // optional, default=all
    reactive_aabb // optional
    (
      // list of min/max of axis-aligned bounding boxs (aabb)
      (150e-6 -1e6 -1e6) (170e-6 1e6 1e6) 
    );
    adp_injection 1e3; // optional & only for apr equation, default = 0
    nonreactive_type  mixed; // optional, default = none
    nonreactive_refValue  uniform 0; // optional
    nonreactive_refGradient  uniform 0; // optional
    nonreactive_valueFraction  uniform 0; // optional
  }

---------------------------------------------------------------
Example BC spec:

  WALL_BOTTOM
  {
    type thrombosis;
    reactive_species AP;
    reactive_zone inside_aabb;
    reactive_aabb
    (
      (150e-06 -1e+06 -1e+06) (170e-06 1e+06 1e+06)
      (190e-06 -1e+06 -1e+06) (210e-06 1e+06 1e+06)
    );
    value uniform 0;
  }

---------------------------------------------------------------
Other sample usages:

  See ../weitaoFoamOF6/sample_cases.
  See ./sample_bcs/hx09-snappyhex

---------------------------------------------------------------
To compile and use:

  1. Edit Make/options to set the version of OF with one of:
     -DOPENFOAM_COM (supports OF v1906 and higher)
     -DOPENFOAM_ORG (supports OF v7.0 and higher)
     -DOPENFOAM_ORG_60 (supports OF v6.0)
  2. wmake
     - this will create lib_thrombosis_bc_1.0 inside FOAM_USER_LIBBIN
  3. In system/controldict, replace the swak4Foam libs with lib_thrombosis_bc_1.0
  4. Replace species bcs following the examples above
     