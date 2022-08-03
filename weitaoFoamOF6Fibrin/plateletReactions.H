//Griffith's model:

    Gamma=k1T*HeCon*ATH/(alphaHep*Kat*Kt+alphaHep*Kat*IIa+ATH*IIa);

    forAll(mesh.cells(), celli)
    {
		if (Gamma[celli] < scalar(0.0))
		{
			Gamma[celli] = scalar(1.0e-30);
		} 
    }

    forAll(mesh.boundary(), celli)
    {
		forAll(mesh.boundary()[celli],facei)
		{
			if (Gamma.boundaryField()[celli][facei] < scalar(0.0))
			{
			Gamma.boundaryFieldRef()[celli][facei] =  scalar(1.0e-30);//changed to Ref
			} 
		}
    }

//PLT activation due to agonists:

    Omega = wapr*apr/aprC + waps*aps/apsC + wTB*IIa/betaThr/tbC;

    forAll(mesh.cells(), celli)
    {
		if (Omega[celli] < scalar(1.0))
		{
			Omega[celli] = scalar(0.0);
		} 
    }

    forAll(mesh.boundary(), celli)
    {
		forAll(mesh.boundary()[celli],facei)
		{
			if (Omega.boundaryField()[celli][facei] < scalar(1.0))
			{
				Omega.boundaryFieldRef()[celli][facei] =  scalar(0.0);//changed to Ref in OF6
			} 
		}
    }
    
    kapa = Omega/tct;
    kapa = min(1/tact,kapa);

//PLT diffusion coefficient:

    gU = fvc::grad(U);
    shearRateU = sqrt(scalar(2.0)*tr(symm(gU)&symm(gU)));
    Dplt = DpltB + ShearID*shearRateU;//shear-induced diffusion
    DpltIn = Dplt;//DpltIn is used in boundary conditions

//PLT activation due to shear stress:

    tau = nu*rho*shearRateU;//in Pa, converted to dyne/cm2 when used in tctspa (multiplied by 10.0)
	  tctspa = 4.0e6*pow(10.0*tau/unitk*unitm*units*units + smalln,-2.3)*units;
	  tctspa = max(tact,tctspa);//ERROR? and should tact be different from agonist tact

	  kspa = scalar(1.0)/tctspa;

//---Prevent high activation at inlet of atrium Cannula simulation
 
    const label cellZoneID = mesh.cellZones().findZoneID("ENTRY");
    const cellZone& zone = mesh.cellZones()[cellZoneID];
    const cellZoneMesh& zoneMesh = zone.zoneMesh();
    const labelList& cellsZone = zoneMesh[cellZoneID]; // list of ZoneID cells
 
    forAll(cellsZone, cI)
    {
        kspa[cellsZone[cI]] = 0;
        krpd[cellsZone[cI]] = 0;
        kapd[cellsZone[cI]] = 0;
        kFnd[cellsZone[cI]] = 0;
    }

//Thrombus volume calculation

    forAll(mesh.cells(), celli)
    {
		if (THVF[celli] > THVFC.value())	//only takes into account cells with THVF > 0.1
		{
			cellThrombVol[celli] = THVF[celli]*meshVolumeF[celli];
		} 
    }

    thrombVol = gSum(cellThrombVol)*unitm*unitm*unitm;//in [m3]


