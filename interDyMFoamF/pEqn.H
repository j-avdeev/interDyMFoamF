{
    rAU = 1.0/UEqn.A();

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, Uf)
    );

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    surfaceScalarField phig
    (
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryField(),
        (
            phiHbyA.boundaryField()
          - (mesh.Sf().boundaryField() & U.boundaryField())
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );
	
	scalarField& D = p_rghEqn.diag(); //*

	// Calculate the sum-mag off-diagonal from the interior faces
	scalarField sumOff(D.size(), 0.0);

	// Handle the boundary contributions to the diagonal
	forAll(p_rgh.boundaryField(), patchI)
	{
	    const fvPatchField<scalar>& ptf = p_rgh.boundaryField()[patchI];

	    if (ptf.size())
	    {
		const labelUList& pa = p_rghEqn.lduAddr().patchAddr(patchI);
		scalarField& iCoeffs = p_rghEqn.internalCoeffs()[patchI];

		if (ptf.coupled())
		{
		    const scalarField& pCoeffs = p_rghEqn.boundaryCoeffs()[patchI];

		    // For coupled boundaries add the diagonal and
		    // off-diagonal contributions
		    forAll(pa, face)
		    {
			D[pa[face]] += component(iCoeffs[face], 0);
			sumOff[pa[face]] += mag(component(pCoeffs[face], 0));
		    }
		}
		else
		{
		    // For non-coupled boundaries add the maximum magnitude diagonal
		    // contribution to ensure stability
		    forAll(pa, face)
		    {
			D[pa[face]] += cmptMax(cmptMag(iCoeffs[face]));
		    }
		}
	    }
	}
	
	volScalarField F0
	(
//             "F0",
            mag(fvc::grad(p_rgh))		// forDel!
        );
// 	F0.internalField() = mag(D)+mag(sumOff);
// // 	F0.internalField() = mag(D)+mag(sumOff);
// 	F0.internalField() = mag(D)+mag(sumOff);
// // //         volScalarField normalisedF
// // //         (
// // //             "normalisedF",
// // //             F0/max(F0)				// normalisation of F0
// // //         );
// // //         normalisedF.writeOpt() = IOobject::AUTO_WRITE;
// // //         F0.clear();
	
	normalisedF.internalField() = mag(D)+mag(sumOff);
// 	Info << "F = " << F << endl;
// 	Info << "F-mag(D) = " << F-mag(D) << endl;
// 	Info << "mag(D)-mag(sumOff) = " << mag(D)-mag(sumOff) << endl;
	
        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
        
    }

    #include "continuityErrs.H"

    {
        Uf = fvc::interpolate(U);
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        Uf += n*(phi/mesh.magSf() - (n & Uf));
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
    
    
}
