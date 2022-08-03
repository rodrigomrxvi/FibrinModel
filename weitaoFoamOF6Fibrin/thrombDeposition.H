
    dTime=runTime.deltaT();

    pltVF = (RPd+APd+APs)/PLTmax;// [PLTmax2D at wall] & [PLTmax3D at internalField] 
    pltVF = min(pltVF,1.0);
	pltVF.correctBoundaryConditions();

    FnVF = (Fgd +Fnd)/Fnmax;// [PLTmax2D at wall] & [PLTmax3D at internalField] 
    FnVF = min(FnVF,1.0);
	FnVF.correctBoundaryConditions();

    THVF = (RPd+APd+APs)/PLTmax + (Fgd +Fnd)/Fnmax;// [PLTmax2D at wall] & [PLTmax3D at internalField] 
    THVF = min(THVF,1.0);
	THVF.correctBoundaryConditions();

    hinFunc = THVF*(scalar(1.0)+scalar(6.5)*THVF);

    SwIn = 1.0-((RPd+APd+APs)/PLTmax2D + (Fgd+Fnd)/Fnmax2D);

//-----------------------------------------------------------------------------------------------------------------------

    shearEmb = exp((-0.0095)*10.0*tau/tauEmb);
    shearEmb = min(1.0,shearEmb);
    shearEmb = max(0.0,shearEmb);

    fembB = (1.0-exp((-0.0095)*10.0*tau/tauEmbB))/units;//in boundary
    fembB = min(1.0/units,fembB);
    fembB = max(0.0/units,fembB);
    fembB = fembB*switchEmb;

    kembf = (0.0*fvc::interpolate(THVF));
    kpdf = (0.0*fvc::interpolate(THVF));
    kpdf1 = (1.0*fvc::interpolate(THVF));

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
		if(THVF[owner[facei]]>THVFC.value())
		{
			kpdf[facei] = kpdf1[facei];
			kembf[facei] = 1.0;
		}
    }
    forAll(neighbour, facei)
    {
		if(THVF[neighbour[facei]]>THVFC.value())
		{
			kpdf[facei] = kpdf1[facei];
			kembf[facei] = 1.0;
		}
    }

    forAll(mesh.boundaryMesh(), patchI)
    {
		// Coupled boundaries (processor, cyclic, etc)
		if (mesh.boundaryMesh()[patchI].coupled())
		{
			forAll(mesh.boundary()[patchI], faceI)
			{
				if(THVF.boundaryField()[patchI][faceI]>THVFC.value())
				{
					kpdf.boundaryFieldRef()[patchI][faceI] = kpdf1.boundaryField()[patchI][faceI];//changed to Ref
					kembf.boundaryFieldRef()[patchI][faceI] = 1.0;//changed to Ref
				}
			} 
		}
		//- Standard boundaries - if you want to include them
		if (mesh.boundaryMesh()[patchI].type()=="wall")
		{
			forAll(mesh.boundary()[patchI], faceI)
			{
				kpdf.boundaryFieldRef()[patchI][faceI] = kpdf1.boundaryField()[patchI][faceI];//changed to Ref
				kembf.boundaryFieldRef()[patchI][faceI] = 1.0;//changed to Ref
			}
		}
    }

    kemb=(fvc::surfaceSum(kembf*mesh.magSf()))/unitm/unitm;//this is div(kembf.n)
    kapd=(fvc::surfaceSum(kpdf*mesh.magSf())*kaa)/unitm/unitm/unitm;
    krpd=(fvc::surfaceSum(kpdf*mesh.magSf())*kra)/unitm/unitm/unitm;
    kFnd=(fvc::surfaceSum(kpdf*mesh.magSf())*kfn)/unitm/unitm/unitm;
    forAll(mesh.cells(),celli)
    {
		kapd[celli] /=meshVolumeF[celli];
		krpd[celli] /=meshVolumeF[celli];
		kemb[celli] /=meshVolumeF[celli];
		kFnd[celli] /=meshVolumeF[celli];
    };

    femb = Lplt*kemb*(1.0-shearEmb)/unitm/units*switchEmb;

    forAll(mesh.cells(), celli)
    {
		if(THVF[celli] >= 0.999)//turn off deposition after a cell fills up 
		{
			kapd[celli] = 0;
			krpd[celli] = 0;
			kFnd[celli] = 0;
		}
    }

//---------------------------------------- To avoid THVF drop in thrombus core ----------------------------------------
/* flow through the porous body of the thrombus can create artificiially high shear leading to errosion and instability */ 

    forAll(mesh.cells(), celli)
    {
		if(THVF[celli]>0.9) 
		{
			femb[celli] = 0;
			fembB[celli] = 0;
		}
    }
    forAll(mesh.boundary(), celli)
    {
		forAll(mesh.boundary()[celli],facei)
		{
			if(THVF.boundaryField()[celli][facei]>0.9) 
			{
				femb.boundaryFieldRef()[celli][facei] = 0;//changed to Ref
				fembB.boundaryFieldRef()[celli][facei] = 0;//changed to Ref
			}
		}
    }

