// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevSigma(U)
  - (1.0 - scalar(bodyForceStabilization)) * turbines.force()
    
 ==
    fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

//volVectorField debugForce
//(
//    "debugForce",
 //   turbines.force()
//);

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - fvc::snGrad(p)
                  + scalar(bodyForceStabilization) * fvc::flux(turbines.force())/mesh.magSf()
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
