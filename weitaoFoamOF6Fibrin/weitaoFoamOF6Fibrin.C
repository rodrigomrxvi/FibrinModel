/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    pimpleFoam

Description
    Large time-step transient solver for incompressible, flow using the PIMPLE
    (merged PISO-SIMPLE) algorithm.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable finite volume options, e.g. MRF, explicit porosity

    weitaoFoamVWFcsOF6: SciRep model + vWFcs model adapted for OF6

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"	//was #include "turbulenceModel.H" 
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"					// was #include "fvIOoptionList.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createSpecies.H" //*
    #include "createWallBC.H"  //*
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
    scalar iDeltaT(runTime.deltaTValue());	//*//save initial deltaT value
    scalar runTimeNoSet(0);

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

	runTimeNoSet += 1;
        Info<< "runTimeNoSet = " << runTimeNoSet << nl << endl;	

	if (runTimeNoSet>=depositionNo)
	{

	    runTimeNoSet = 0;
	    runTime.setDeltaT(iDeltaT/velLoopNo);	//*// set new deltaT for pimple loop
	    Info<< "deltaT = " <<  runTime.deltaT().value() << endl;

	    // --- Pressure-velocity PIMPLE corrector loop
	    while (pimple.loop())
	    {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

			#include "UEqn.H"

			// --- Pressure corrector loop
			while (pimple.correct())
			{
				#include "pEqn.H"
			}

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
	    }

	}

	runTime.setDeltaT(iDeltaT);	//*// restore initial deltaT for thrombosis equations 
	
    #include "plateletReactions.H"
	if(thrombosisSW==1.0)
	{
		#include "thrombDeposition.H"
		#include "speciesEquations.H"	
	}

	Info<< "MaxP" << max(p) << endl;
	Info<< "MinP" << min(p) << endl;
	Info<< "MaxU" << max(mag(U)) << endl;
	Info<< "MaxShR" << max(shearRateU) << endl;
  Info<< "MaxpltVF" << max(pltVF) << endl;
  Info<< "MaxFnVF" << max(FnVF) << endl;
  Info<< "MaxTHVF" << max(THVF) << endl;
	Info<< "ThrombVol" << max(thrombVol) << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
