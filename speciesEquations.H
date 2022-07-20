
//Convection-Diffusion-Reaction equations:

        fvScalarMatrix RPEqn
        (
            fvm::ddt(RP)
           +fvm::div(phi, RP)
           -fvm::laplacian(Dplt, RP)	
          == 
           -fvm::Sp(kapa, RP)
           -fvm::Sp(kspa, RP)	
           -fvm::Sp(krpd, RP)	
		   +femb*RPd
        );
        RPEqn.relax();            
        RPEqn.solve();


        fvScalarMatrix APEqn
        (
            fvm::ddt(AP)
           +fvm::div(phi, AP)
           -fvm::laplacian(Dplt, AP)
          == 
	    	kapa*RP
           +kspa*RP	
		   +femb*APd
           -fvm::Sp(kapd, AP)
        );
        APEqn.relax();            
        APEqn.solve();


        fvScalarMatrix aprEqn		
        (
            fvm::ddt(apr)
           +fvm::div(phi, apr)
           -fvm::laplacian(Dapr, apr)	
          == 
            lambdaj*(kapa*RP
           +kspa*RP+kapa*RPd
           +kspa*RPd+theta*krpd*RP)
        );
        aprEqn.relax();            
        aprEqn.solve();


        fvScalarMatrix apsEqn
        (
            fvm::ddt(aps)
           +fvm::div(phi, aps)
           -fvm::laplacian(Daps, aps)
          ==
            spj*(AP+APd)			
           -fvm::Sp(k1j, aps)		
        );
        apsEqn.relax();            
        apsEqn.solve();


        fvScalarMatrix RPdEqn
        (
            fvm::ddt(RPd)
          == 
	    (1-theta)*krpd*RP
           -fvm::Sp(kapa, RPd)
           -fvm::Sp(kspa, RPd)
           -fvm::Sp(femb, RPd)
        );
        RPdEqn.relax();            
        RPdEqn.solve();


        fvScalarMatrix APdEqn
        (
            fvm::ddt(APd)
          == 
			theta*krpd*RP
           +kapd*AP
           +kapa*RPd
           +kspa*RPd
           -fvm::Sp(femb+fstb, APd)
        );
        APdEqn.relax();            
        APdEqn.solve();


        fvScalarMatrix APsEqn
        (
            fvm::ddt(APs)
          == 
            fvm::Sp(0.0/units, APs)
            +fstb*APd
        );
        APsEqn.relax();            
        APsEqn.solve();


        // Scalar transport equations for coagulation Factors
        fvScalarMatrix XIIEqn
            (
               fvm::ddt(XII)
            +  fvm::div(phi, XII)
            ==
              fvm::laplacian(DXII, XII)
            );

        XIIEqn.relax();
        XIIEqn.solve();

        fvScalarMatrix XIIaEqn
            (
               fvm::ddt(XIIa)
            +  fvm::div(phi, XIIa)
            ==
              fvm::laplacian(DXIIa, XIIa)
            );

        XIIaEqn.relax();
        XIIaEqn.solve();

        fvScalarMatrix IIEqn
            (
               fvm::ddt(II)
            +  fvm::div(phi, II)
            -  fvm::laplacian(DII, II)
            ==
               -(kcat4*Va*II/(km4+II))-(kcat2*XIIa*II/(km2+II))
               -betaThr*(phiA*(AP+APd)+phiR*(RP+RPd))*II
            );

        IIEqn.relax();
        IIEqn.solve();
        fvScalarMatrix IIaEqn
            (
               fvm::ddt(IIa)
            +  fvm::div(phi, IIa)
            -  fvm::laplacian(DIIa, IIa)
            ==
               kcat2*XIIa*II/(km2+II)+(kcat4*Va*II/(km4+II))
               -(k5c*IIa*ATH)
              +II*(phiA*(AP+APd)+phiR*(RP+RPd))*betaThr
              //-fvm::Sp(Gamma,IIa)  // this term is when Heparing is taken into account, TB should be II
            );

        IIaEqn.relax();
        IIaEqn.solve();

        fvScalarMatrix VEqn
            (
               fvm::ddt(V)
            +  fvm::div(phi, V)
            -  fvm::laplacian(DV, V)
            ==
               -kcat3*V*IIa/(km3+V)
            );

        VEqn.relax();
        VEqn.solve();

        fvScalarMatrix VaEqn
            (
               fvm::ddt(Va)
            +  fvm::div(phi, Va)
            -  fvm::laplacian(DVa, Va)
            ==
               kcat3*V*IIa/(km3+V)
            );

        VaEqn.relax();
        VaEqn.solve();

        fvScalarMatrix ATHEqn
        (
            fvm::ddt(ATH)
           +fvm::div(phi, ATH)
           -fvm::laplacian(DATH, ATH)
          ==
           -k5c*IIa*ATH
           //-Gamma*IIa  // this term is when Heparing is taken into account, TB should be II
        );

        ATHEqn.relax();
        ATHEqn.solve();

        fvScalarMatrix ATIIaEqn
            (
               fvm::ddt(ATIIa)
            +  fvm::div(phi, ATIIa)
            -  fvm::laplacian(DATIIa, ATIIa)
            ==
               k5c*IIa*ATH
            );

        ATIIaEqn.relax();
        ATIIaEqn.solve();

        fvScalarMatrix FgEqn
            (
               fvm::ddt(Fg)
            +  fvm::div(phi, Fg)
            -  fvm::laplacian(DFg, Fg)
            ==
            - kcat6*Fg*IIa/(km6+Fg)*0.05 // 0.05 from Chen et al Plos One 2019
            );

        FgEqn.relax();
        FgEqn.solve();

        fvScalarMatrix FnEqn
            (
               fvm::ddt(Fn)
            +  fvm::div(phi, Fn)
            -  fvm::laplacian(DFn, Fn)
            ==
               kcat6*Fg*IIa/(km6+Fg)*0.05 // from Chen et al Plos One 2019
               -kFnd*Fn
               +femb*Fnd  // embolized Fibrin
            );

        FnEqn.relax();
        FnEqn.solve();

        fvScalarMatrix FgdEqn
            (
               fvm::ddt(Fgd)
            //==
            //   -femb*Fgd
            );
        FgdEqn.relax();
        FgdEqn.solve();

        fvScalarMatrix FndEqn
            (
                fvm::ddt(Fnd)
              ==
                kFnd*Fn    // Depositon
                -femb*Fnd  // embolized Fibrin
            );
        FndEqn.relax();
        FndEqn.solve();


//*******************************************************************************//

//Bound species concentration fields for stability:

RP = max(minRP,RP);
AP = max(minAP,AP);
apr = max(minapr,apr);
aps = max(minaps,aps);

RPd = max(minRPd,RPd);
APd = max(minAPd,APd);
APs = max(minAPs,APs);

RPd = min(PLTmax,RPd);
APd = min(PLTmax,APd);
APs = min(PLTmax,APs);

XII = max(minXII,XII);
XIIa = max(minXIIa,XIIa);
II = max(minII,II);
IIa = max(minIIa,IIa);
V = max(minV,V);
Va = max(minVa,Va);
ATH = max(minATH,ATH);
ATIIa = max(minATIIa,ATIIa);
Fg = max(minFg,Fg);
Fn = max(minFn,Fn);
Fgd = max(minFgd,Fgd);
Fnd = max(minFnd,Fnd);
Fnd = min(Fnmax,Fnd);


forAll(mesh.boundary(), patchi)
{
   forAll(mesh.boundary()[patchi], facei)
   {
		if (RP.boundaryField()[patchi][facei]<0.0){RP.boundaryFieldRef()[patchi][facei]=1e-10;}//changed to boundaryFieldRef() for OF6
		if (AP.boundaryField()[patchi][facei]<0.0){AP.boundaryFieldRef()[patchi][facei]=1e-10;} 
		if (apr.boundaryField()[patchi][facei]<0.0){apr.boundaryFieldRef()[patchi][facei]=1e-10;} 
		if (aps.boundaryField()[patchi][facei]<0.0){aps.boundaryFieldRef()[patchi][facei]=1e-10;} 

    if (XII.boundaryField()[patchi][facei]<0.0){XII.boundaryFieldRef()[patchi][facei]=1e-10;}//changed to boundaryFieldRef() for OF6
		if (XIIa.boundaryField()[patchi][facei]<0.0){XIIa.boundaryFieldRef()[patchi][facei]=1e-10;}
		if (II.boundaryField()[patchi][facei]<0.0){II.boundaryFieldRef()[patchi][facei]=1e-10;}
		if (IIa.boundaryField()[patchi][facei]<0.0){IIa.boundaryFieldRef()[patchi][facei]=1e-10;}
		if (V.boundaryField()[patchi][facei]<0.0){V.boundaryFieldRef()[patchi][facei]=1e-10;}
		if (Va.boundaryField()[patchi][facei]<0.0){Va.boundaryFieldRef()[patchi][facei]=1e-10;}
		if (ATH.boundaryField()[patchi][facei]<0.0){ATH.boundaryFieldRef()[patchi][facei]=1e-10;}
		if (ATIIa.boundaryField()[patchi][facei]<0.0){ATIIa.boundaryFieldRef()[patchi][facei]=1e-10;}
		if (Fg.boundaryField()[patchi][facei]<0.0){Fg.boundaryFieldRef()[patchi][facei]=1e-10;}
		if (Fn.boundaryField()[patchi][facei]<0.0){Fn.boundaryFieldRef()[patchi][facei]=1e-10;}

		if (RPd.boundaryField()[patchi][facei]<0.0){RPd.boundaryFieldRef()[patchi][facei]=1e-10;} 
		if (APd.boundaryField()[patchi][facei]<0.0){APd.boundaryFieldRef()[patchi][facei]=1e-10;} 
		if (APs.boundaryField()[patchi][facei]<0.0){APs.boundaryFieldRef()[patchi][facei]=1e-10;} 
    if (Fnd.boundaryField()[patchi][facei]<0.0){Fnd.boundaryFieldRef()[patchi][facei]=1e-10;}
		if (Fgd.boundaryField()[patchi][facei]<0.0){Fgd.boundaryFieldRef()[patchi][facei]=1e-10;}

		if (RPd.boundaryField()[patchi][facei]>PLTmax2D.value()){RPd.boundaryFieldRef()[patchi][facei]=PLTmax2D.value();} 
		if (APd.boundaryField()[patchi][facei]>PLTmax2D.value()){APd.boundaryFieldRef()[patchi][facei]=PLTmax2D.value();} 
		if (APs.boundaryField()[patchi][facei]>PLTmax2D.value()){APs.boundaryFieldRef()[patchi][facei]=PLTmax2D.value();} 
    if (Fgd.boundaryField()[patchi][facei]>Fnmax2D.value()){Fgd.boundaryFieldRef()[patchi][facei]=Fnmax2D.value();}
		if (Fnd.boundaryField()[patchi][facei]>Fnmax2D.value()){Fnd.boundaryFieldRef()[patchi][facei]=Fnmax2D.value();}
   }
}


