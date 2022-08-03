    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );	

    scalar switchEmb = readScalar(transportProperties.lookup("switchEmb")); 		//switch shear cleaning on-off
    scalar depositionNo = readScalar(transportProperties.lookup("depositionNo"));	//update velocity field per depositionNo loop
    scalar velLoopNo = readScalar(transportProperties.lookup("velLoopNo"));			//deltaT for pimple loop = deltaT/velLoopNo
    scalar thrombosisSW = readScalar(transportProperties.lookup("thrombosisSW"));	//switch speciesEquations.H and thrombusDeposition.H on-off


//Diffusion coef:

    dimensionedScalar nu 
    (
        transportProperties.lookup("nu")
    );

    dimensionedScalar rho 
    (
        transportProperties.lookup("rho")
    );

    dimensionedScalar ShearID 	//coef in platelet diffusivity equation 6.0e-13
    (
        transportProperties.lookup("ShearID")
    );

    dimensionedScalar DpltB 	//Brownian diffusion coef for RP & AP
    (
        transportProperties.lookup("DpltB")
    );

    dimensionedScalar Dpt 	//diffusion coefficient for prothrombin
    (
        transportProperties.lookup("Dpt")
    );

    dimensionedScalar Dtb 	//diffusion coefficient for thrombin
    (
        transportProperties.lookup("Dtb")
    );

    dimensionedScalar Dat 	//diffusion coefficient for ATIII
    (
        transportProperties.lookup("Dat")
    );

    dimensionedScalar Dapr 	//PLT-released agonist (ADP) diffusivity
    (
        transportProperties.lookup("Dapr")
    );

    dimensionedScalar Daps 	//PLT-synthesized agonist (TxA2) diffusivity
    (
        transportProperties.lookup("Daps")
    );

    dimensionedScalar DXII 	//diffusion coefficient for XII
    (
        transportProperties.lookup("DXII")
    );

    dimensionedScalar DXIIa 	//diffusion coefficient for XIIa
    (
        transportProperties.lookup("DXIIa")
    );

    dimensionedScalar DII 	//diffusion coefficient for II
    (
        transportProperties.lookup("DII")
    );

    dimensionedScalar DIIa 	//diffusion coefficient for IIa
    (
        transportProperties.lookup("DIIa")
    );

    dimensionedScalar DV 	//diffusion coefficient for V
    (
        transportProperties.lookup("DV")
    );

    dimensionedScalar DVa 	//diffusion coefficient for Va
    (
        transportProperties.lookup("DVa")
    );

    dimensionedScalar DATH 	//diffusion coefficient for ATH
    (
        transportProperties.lookup("DATH")
    );

    dimensionedScalar DATIIa 	//diffusion coefficient for ATIIa
    (
        transportProperties.lookup("DATIIa")
    );

    dimensionedScalar DFg    	//diffusion coefficient for Fg
    (
        transportProperties.lookup("DFg")
    );

    dimensionedScalar DFn    	//diffusion coefficient for Fn
    (
        transportProperties.lookup("DFn")
    );

//Reaction coef:

    dimensionedScalar wapr 	//weight of apr in Omega equation
    (
        transportProperties.lookup("wapr")
    );

    dimensionedScalar waps 	//weight of aps in Omega equation
    (
        transportProperties.lookup("waps")
    );

    dimensionedScalar wTB 	//weight of TB in Omega equation
    (
        transportProperties.lookup("wTB")
    );

    dimensionedScalar tct 	//characteristic time constant for PLT activation
    (
        transportProperties.lookup("tct")
    );

    dimensionedScalar tact 	//physical PLT activation characteristic time
    (
        transportProperties.lookup("tact")
    );

    dimensionedScalar aprC 	//threshold concentration of ADP for PLT activation
    (
        transportProperties.lookup("aprC")
    );

    dimensionedScalar apsC 	//threshold concentration of TxA2 for PLT activation
    (
        transportProperties.lookup("apsC")
    );

    dimensionedScalar tbC 	//threshold concentration of TB for PLT activation
    (
        transportProperties.lookup("tbC")
    );

    dimensionedScalar lambdaj	//amount of agonist (ADP) released per platelet
    (
        transportProperties.lookup("lambdaj")
    );

    dimensionedScalar spj 	//rate of synthesis of agonist TxA2
    (
        transportProperties.lookup("spj")
    );

    dimensionedScalar k1j 	//1st-order reaction rate constant for inhibition of TxA2
    (
        transportProperties.lookup("k1j")
    );

    dimensionedScalar k1T 	//first-order rate constant in Griffith's model
    (
        transportProperties.lookup("k1T")
    );

    dimensionedScalar k2 	//second-order rate constant (used in absence of heparin)
    (
        transportProperties.lookup("k2")
    );

    dimensionedScalar Kt 	//dissociation constant for heparin/thrombin in Griffith's model
    (
        transportProperties.lookup("Kt")
    );

    dimensionedScalar Kat 	//dissociation constant for heparin/ATIII in Griffith's model
    (
        transportProperties.lookup("Kat")
    );

    dimensionedScalar alphaHep 	//factor for change in affinity of heparin for ATIII or thrombin
    (
        transportProperties.lookup("alphaHep")
    );

    dimensionedScalar betaThr 	//conversion factor to convert thrombin concentration from U.ml-1 to uM
    (
        transportProperties.lookup("betaThr")
    );

    dimensionedScalar phiR 	//rate constant for TB generation from PT at the surface of RP
    (
        transportProperties.lookup("phiR")
    );

    dimensionedScalar phiA 	//rate constant for TB generation from PT at the surface of AP
    (
        transportProperties.lookup("phiA")
    );

    dimensionedScalar HeCon 	//Heparin concentration
    (
        transportProperties.lookup("HeCon")
    );

    dimensionedScalar k1c  // reaction rate XII --> XIIa
    (
        "k1c",
        transportProperties.lookup("k1c")
    );

    dimensionedScalar km2  // reaction rate XIIa --> IIa
    (
        "km2",
        transportProperties.lookup("km2")
    );

    dimensionedScalar kcat2  // reaction rate XIIa --> IIa
    (
        "kcat2",
        transportProperties.lookup("kcat2")
    );

    dimensionedScalar km3  // reaction rate II --> Va
    (
        "km3",
        transportProperties.lookup("km3")
    );

    dimensionedScalar kcat3  // reaction rate II --> Va
    (
        "kcat3",
        transportProperties.lookup("kcat3")
    );

    dimensionedScalar km4  // reaction rate Va --> IIa
    (
        "km4",
        transportProperties.lookup("km4")
    );
   
    dimensionedScalar kcat4  // reaction rate Va --> IIa
    (
        "kcat4",
        transportProperties.lookup("kcat4")
    );
   
    dimensionedScalar k5c  // reaction rate AT + IIa --> ATIIa
    (
        "k5c",
        transportProperties.lookup("k5c")
    );

    dimensionedScalar km6  // reaction rate Fg + IIa --> Fn + IIa
    (
        "km6",
        transportProperties.lookup("km6")
    );

    dimensionedScalar kcat6  // reaction rate Fg + IIa --> Fn + IIa
    (
        "kcat6",
        transportProperties.lookup("kcat6")
    );

    dimensionedScalar k6  // reaction rate Fg + IIa --> Fn + IIa From Govindarajan 
    (
        "k6",
        transportProperties.lookup("k6")
    );

    dimensionedScalar ksXII  // reaction rate XII --> XIIa
    (
        "ksXII",
        transportProperties.lookup("ksXII")
    );

//Deposition & Embolization coef:

    dimensionedScalar theta 	//fraction of PLTs activating upon contact
    (
        transportProperties.lookup("theta")
    );

    dimensionedScalar C2 	//thrombus resistance coef
    (
        transportProperties.lookup("C2")
    );

    dimensionedScalar THVFC 	//critical PLT volume fraction (maximum packing fraction)
    (
        transportProperties.lookup("THVFC")
    );

    dimensionedScalar PLTmax3D 	//max conc of platelets in space; PLTmax2D/DiaPLT
    (
        transportProperties.lookup("PLTmax3D")
    );

    dimensionedScalar PLTmax2D 	//max conc of platelets on surface;
    (
        transportProperties.lookup("PLTmax2D")
    );

    dimensionedScalar Fnmax3D 	//max conc of fibrin in space; Fnmax2D/DiaFn
    (
        transportProperties.lookup("Fnmax3D")
    );

    dimensionedScalar Fnmax2D 	//max conc of fibrin on surface;
    (
        transportProperties.lookup("Fnmax2D")
    );

    dimensionedScalar kra 	//related to krpd
    (
        transportProperties.lookup("kra")
    );

    dimensionedScalar kaa	//related to kapd
    (
        transportProperties.lookup("kaa")
    );

    dimensionedScalar kfn	//related to kfnd
    (
        transportProperties.lookup("kfn")
    );

    dimensionedScalar krpdB 	//RP-Boundary(Wall) deposition rate 
    (
        transportProperties.lookup("krpdB")
    );

    dimensionedScalar kapdB	//AP-Boundary(Wall) deposition rate 
    (
        transportProperties.lookup("kapdB")
    );

    dimensionedScalar kFgdB	//Fg-Boundary(Wall) deposition rate
    (
        transportProperties.lookup("kFgdB")
    );

    dimensionedScalar kFndB	//Fn-Boundary(Wall) deposition rate
    (
        transportProperties.lookup("kFndB")
    );

    dimensionedScalar tauEmb 	//related to shear embolization
    (
        transportProperties.lookup("tauEmb")
    );

    dimensionedScalar tauEmbB 	//related to shear embolization
    (
        transportProperties.lookup("tauEmbB")
    );

    dimensionedScalar fstb 	
    (
        transportProperties.lookup("fstb")
    );

    dimensionedScalar Lplt 	//platelet diameter
    (
        transportProperties.lookup("Lplt")
    );


//---------- For manipulating units --------------------------------

    dimensionedScalar unitk
    (
        transportProperties.lookup("unitk")
    );

    dimensionedScalar unitm
    (
        transportProperties.lookup("unitm")
    );

    dimensionedScalar units
    (
        transportProperties.lookup("units")
    );

    dimensionedScalar smalln
    (
        transportProperties.lookup("smalln")
    );

//---------- For bounding the species concentration --------------------------------

    dimensionedScalar minRP("minRP",dimensionSet(0, -3, 0, 0, 0, 1, 0),1.0e5);
    dimensionedScalar minAP("minAP",dimensionSet(0, -3, 0, 0, 0, 1, 0),1.0e5);
    dimensionedScalar minAT("minAT",dimensionSet(0, -3, 0, 0, 1, 0, 0),0.0);
    dimensionedScalar minapr("minapr",dimensionSet(0, -3, 0, 0, 1, 0, 0),0.0);
    dimensionedScalar minaps("minaps",dimensionSet(0, -3, 0, 0, 1, 0, 0),0.0);
    dimensionedScalar minRPd("minRPd",dimensionSet(0, -3, 0, 0, 0, 1, 0),0.0);
    dimensionedScalar minAPd("minAPd",dimensionSet(0, -3, 0, 0, 0, 1, 0),0.0);
    dimensionedScalar minAPs("minAPs",dimensionSet(0, -3, 0, 0, 0, 1, 0),0.0);
     dimensionedScalar minXII("minXII",dimensionSet(0, -3, 0, 0, 1, 0, 0),0.0);
    dimensionedScalar minXIIa("minXIIa",dimensionSet(0, -3, 0, 0, 1, 0, 0),0.0);
    dimensionedScalar minII("minII",dimensionSet(0, -3, 0, 0, 1, 0, 0),0.0);
    dimensionedScalar minIIa("minIIa",dimensionSet(0, -3, 0, 0, 1, 0, 0),0.0);
    dimensionedScalar minV("minV",dimensionSet(0, -3, 0, 0, 1, 0, 0),0.0);
    dimensionedScalar minVa("minVa",dimensionSet(0, -3, 0, 0, 1, 0, 0),0.0);
    dimensionedScalar minATH("minATH",dimensionSet(0, -3, 0, 0, 1, 0, 0),0.0);
    dimensionedScalar minATIIa("minATIIa",dimensionSet(0, -3, 0, 0, 1, 0, 0),0.0);
    dimensionedScalar minFg("minFg",dimensionSet(0, -3, 0, 0, 1, 0, 0),0.0);
    dimensionedScalar minFn("minFn",dimensionSet(0, -3, 0, 0, 1, 0, 0),0.0);
    dimensionedScalar minFgd("minFgd",dimensionSet(0, -3, 0, 0, 1, 0, 0),0.0);
    dimensionedScalar minFnd("minFnd",dimensionSet(0, -3, 0, 0, 1, 0, 0),0.0);

//*******************************************************************************//

//Species:

 Info<< "Reading field RP\n" << endl;		//Resting Platelets
    volScalarField RP
    (
        IOobject
        (
            "RP",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field AP\n" << endl;		//Activated Platelets
    volScalarField AP
    (
        IOobject
        (
            "AP",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field apr\n" << endl;		//PLT-released agonist (ADP)
    volScalarField apr
    (
        IOobject
        (
            "apr",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field aps\n" << endl;		//PLT-synthesized agonist (TxA2)
    volScalarField aps
    (
        IOobject
        (
            "aps",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field RPd\n" << endl;		//Deposited Resting Platelets
    volScalarField RPd
    (
        IOobject
        (
            "RPd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field APd\n" << endl;		//Deposited Activated Platelets
    volScalarField APd
    (
        IOobject
        (
            "APd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field APs\n" << endl;		//Deposited and Stabilized Platelets
    volScalarField APs
    (
        IOobject
        (
            "APs",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field XII\n" << endl;       //Coagulation factor XII
    volScalarField XII
    (
        IOobject
        (
            "XII",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field XIIa\n" << endl;       //Coagulation factor XIIa
    volScalarField XIIa
    (
        IOobject
        (
            "XIIa",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field II\n" << endl;       //Coagulation factor II
    volScalarField II
    (
        IOobject
        (
            "II",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field IIa\n" << endl;       //Coagulation factor IIa
    volScalarField IIa
    (
        IOobject
        (
            "IIa",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field V\n" << endl;          //Coagulation factor V
    volScalarField V
    (
        IOobject
        (
            "V",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field Va\n" << endl;          //Coagulation factor Va
    volScalarField Va
    (
        IOobject
        (
            "Va",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field ATH\n" << endl;         //Anti thrombin for coagulation
    volScalarField ATH
    (
        IOobject
        (
            "ATH",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field ATIIa\n" << endl;         //Complex AT=thrombin
    volScalarField ATIIa
    (
        IOobject
        (
            "ATIIa",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field Fg\n" << endl;
    volScalarField Fg
    (
        IOobject
        (
            "Fg",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field Fn\n" << endl;
    volScalarField Fn
    (
        IOobject
        (
            "Fn",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field Fgd\n" << endl;
    volScalarField Fgd
    (
        IOobject
        (
            "Fgd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Reading field Fnd\n" << endl;
    volScalarField Fnd
    (
        IOobject
        (
            "Fnd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField PLTmax
    (
        IOobject
        (
            "PLTmax",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField Fnmax
    (
        IOobject
        (
            "Fnmax",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


//*******************************************************************************//

//PLT diffusion coefficient:

   volTensorField gU(fvc::grad(U));

   volScalarField shearRateU
    (
        IOobject
        (
            "shearRateU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		sqrt(scalar(2.0)*tr(symm(gU)&symm(gU)))
    );

    volScalarField Dplt(DpltB + ShearID*shearRateU);

//THVF volume fraction and hindrance function:

    volScalarField THVF
    (
        IOobject
        (
            "THVF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		(RPd+APd+APs)/PLTmax + (Fgd+Fnd)/Fnmax
    );

    volScalarField pltVF
    (
        IOobject
        (
            "pltVF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		(RPd+APd+APs)/PLTmax 
    );

    volScalarField FnVF
    (
        IOobject
        (
            "FnVF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		(Fgd+Fnd)/Fnmax
    );

    volScalarField hinFunc	
    (
        IOobject
        (
            "hinFunc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		THVF*(scalar(1.0)+scalar(6.5)*THVF)
    ); 

    volScalarField cellThrombVol
    (
        IOobject
        (
            "cellThrombVol",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("cellThrombVol",dimensionSet(0, 3, 0, 0, 0, 0, 0),0.0)
    );

    volScalarField thrombVol
    (
        IOobject
        (
            "thrombVol",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("thrombVol",dimensionSet(0, 3, 0, 0, 0, 0, 0),0.0)
    );

    volScalarField thrombVolPlot //in [10^-6 cm3]
    (
        IOobject
        (
            "thrombVolPlot",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		thrombVol*1.0e12 
    );


//Platelet activation:

    volScalarField Omega	//PLT activation function
    (
        IOobject
        (
            "Omega",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	scalar(1.0)*apr/aprC + scalar(3.3)*aps/apsC + scalar(30.0)*IIa/tbC/betaThr
    ); 

    volScalarField kapa 	//PLT activation due to agonists
    (
        IOobject
        (
            "kapa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
	dimensionedScalar("kapa",dimensionSet(0, 0, -1, 0, 0, 0, 0),0.0)
    );

//PLT activation due to shear stress:

   volScalarField tau
    (
        IOobject
        (
            "tau",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	nu*rho*shearRateU //in [dyne/cm2]
    );

    volScalarField tctspa 	
    (
        IOobject
        (
            "tctspa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	4.0e6*pow(tau/unitk*unitm*units*units + smalln,-2.3)*units
    );

    volScalarField kspa
    (
        IOobject
        (
            "kspa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
	dimensionedScalar("kspa",dimensionSet(0, 0, -1, 0, 0, 0, 0),0.0)
    );

//Deposition rates:

    scalarField meshVolumeF(mesh.V());

    surfaceScalarField kembf(0.0*fvc::interpolate(THVF));
    surfaceScalarField kpdf(0.0*fvc::interpolate(THVF));
    surfaceScalarField kpdf1(fvc::interpolate(THVF));

   volScalarField krpd 		//rate of deposition of RP to APd
    (
        IOobject
        (
            "krpd",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	fvc::surfaceSum(fvc::interpolate(THVF/THVFC)*mesh.magSf())*kra/unitm/unitm/unitm
    );

    volScalarField kapd 	//rate of deposition of AP to APd
    (
        IOobject
        (
            "kapd",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	fvc::surfaceSum(fvc::interpolate(THVF/THVFC)*mesh.magSf())*kaa/unitm/unitm/unitm
    );

    volScalarField kFnd 	//rate of deposition of Fn to FnB
    (
        IOobject
        (
            "kFnd",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	fvc::surfaceSum(fvc::interpolate(THVF/THVFC)*mesh.magSf())*kfn/unitm/unitm/unitm
    );

//Shear embolization:

    volScalarField kemb 	
    (
        IOobject
        (
            "kemb",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	fvc::surfaceSum(fvc::interpolate(THVF/THVFC)*mesh.magSf())/unitm/unitm
    );

    volScalarField femb 	
    (
        IOobject
        (
            "femb",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	fvc::surfaceSum(fvc::interpolate(THVF/THVFC)*mesh.magSf())*Lplt/unitm/unitm/unitm/units
    );

    volScalarField shearEmb 	
    (
        IOobject
        (
            "shearEmb",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	exp((-0.0095)*tau/tauEmb)
    );


//Griffith's template model:

    volScalarField Gamma
    (
        IOobject
        (
            "Gamma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	k1T*HeCon*ATH/(alphaHep*Kat*Kt+alphaHep*Kat*IIa+ATH*IIa)
    );

