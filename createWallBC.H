    volScalarField krpdBIn
    (
        IOobject
        (
            "krpdBIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("krpdBIn",krpdB)
    );

    volScalarField kapdBIn
    (
        IOobject
        (
            "kapdBIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("kapdBIn",kapdB)
    );

    volScalarField fembB
    (
        IOobject
        (
            "fembB",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		(1.0-exp((-0.0095)*tau/tauEmbB))/units
    );

     volScalarField kFgdBIn
    (
        IOobject
        (
            "kFgdBIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("kFgdBIn",kFgdB)
    );

    volScalarField kFndBIn
    (
        IOobject
        (
            "kFndBIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("kFndBIn",kFndB)
    );

    volScalarField SwIn
    (
        IOobject
        (
            "SwIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		1.0-((RPd+APd+APs)/PLTmax2D + (Fgd+Fnd)/Fnmax2D)
    );

     volScalarField DpltBIn
    (
        IOobject
        (
            "DpltBIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		mesh,
		dimensionedScalar("DpltBIn",DpltB)
    );

     volScalarField ShearIDIn
    (
        IOobject
        (
            "ShearIDIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("ShearIDIn",ShearID)
    );

    volScalarField DpltIn	
    (
        IOobject
        (
            "DpltIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("DpltIn",dimensionSet(0, 2, -1, 0, 0, 0, 0),1.58e-13)//recalculated in plateletReactions.H
    ); 

     volScalarField DaprIn
    (
        IOobject
        (
            "DaprIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("DaprIn",Dapr)
    );

	volScalarField DapsIn
    (
        IOobject
        (
            "DapsIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("DapsIn",Daps)
    );

     volScalarField DIIIn
    (
        IOobject
        (
            "DIIIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("DIIIn",DII)
    );

     volScalarField DIIaIn
    (
        IOobject
        (
            "DIIaIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("DIIaIn",DIIa)
    );

     volScalarField DatIn
    (
        IOobject
        (
            "DatIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("DatIn",Dat)
    );

    volScalarField thetaIn
    (
        IOobject
        (
            "thetaIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("thetaIn",theta)
    );

    volScalarField fstbIn
    (
        IOobject
        (
            "fstbIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("fstbIn",fstb)
    );


    volScalarField lambdajIn
    (
        IOobject
        (
            "lambdajIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("lambdajIn",lambdaj)
    );

    volScalarField spjIn
    (
        IOobject
        (
            "spjIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("spjIn",spj)
    );

    volScalarField phiAIn
    (
        IOobject
        (
            "phiAIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("phiAIn",phiA)
    );

    volScalarField phiRIn
    (
        IOobject
        (
            "phiRIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("phiRIn",phiR)
    );

     volScalarField betaThrIn
    (
        IOobject
        (
            "betaThrIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("betaThrIn",betaThr)
    );

    volScalarField DXIIIn
    (
        IOobject
        (
            "DXIIIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("DXIIIn",DXII)
    );

     volScalarField ksXIIIn
    (
        IOobject
        (
            "ksXIIIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("ksXIIIn",ksXII)
    );

    volScalarField dTime	
    (
        IOobject
        (
            "dTime",
	    runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "dTime",
            dimensionSet(0, 0, 1, 0, 0, 0, 0),
            runTime.deltaT().value()
        )
    );

    volScalarField DFgIn	
    (
        IOobject
        (
            "DFgIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("DFgIn",dimensionSet(0, 2, -1, 0, 0, 0, 0),3.1e-11)//recalculated in plateletReactions.H
    ); 

    volScalarField DFnIn	
    (
        IOobject
        (
            "DFnIn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("DFnIn",dimensionSet(0, 2, -1, 0, 0, 0, 0),2.47e-11)//recalculated in plateletReactions.H
    ); 

    volScalarField kcat6In	
    (
        IOobject
        (
            "kcat6In",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("kcat6In",dimensionSet(0, 0, -1, 0, 0, 0, 0),80.0)//
    ); 

    volScalarField km6In	
    (
        IOobject
        (
            "km6In",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("km6In",dimensionSet(0, -3, 0, 0, 1, 0, 0),6500.0e3)//
    ); 

