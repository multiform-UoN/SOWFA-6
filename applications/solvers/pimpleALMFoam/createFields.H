#include "createRDeltaT.H"


Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, laminarTransport)
);

// Create source term object (not needed - ABL?)
//drivingForce<vector> momentumSourceTerm(ABLProperties, "momentum",U);

#include "createMRF.H"
#include "createFvOptions.H"
bool bodyForceStabilization = mesh.schemesDict().lookupOrDefault<bool>("bodyForceStabilization",true);

// Create an object of the horizontalWindTurbineArray class if there is to be a turbine array
turbineModels::propALMAdvanced turbines(U);
