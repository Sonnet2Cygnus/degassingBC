/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  10
     \\/     M anipulation  |

  degassingBC class

    This class implements the degassing boundary condition model。
    It is used to handle degassing processes in multiphase flow simulations by adding source terms to phase equations.
\*---------------------------------------------------------------------------*/



#ifndef degassingBC_H
#define degassingBC_H


#include "fvModel.H"
#include "fvCellSet.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "phaseSystem.H"


namespace Foam
{
namespace fv
{
class degassingBC:public fvModel
{
    // Private member data

        // List of names corresponding to degassing phases
        wordList degassingBC_;

        // Name of the continuous phase 
        word continuousPhase_;

        // Name of the phase fraction field 
        word alphaName_;

        // Constant reference to the phaseSystem, providing access to phase-related properties
        const phaseSystem& fluid_;

        //  List of names for alpha (phase fraction) fields
        wordList fieldNames_;

        // List of names for outlet boundary patches
        wordList patchNames_;

        // List of cell labels where source terms will be added
        labelList cells_;

        // Total volume of cells considered
        scalar V_;


    // Private member functions

        // Checks if a given field name corresponds to a degassing phase
        bool isDegassingPhases(const word& fieldName) const;

        // Calculates the total volume of cells attached to specified patches
        scalar patchCellVolume(const wordList& patches) const;

        // Template function to add source terms of a specific type to the transport equation
        template<class Type>

        void addSupType                         // Add source term type
        (
            const volScalarField& alpha,        // Read the phase fraction field
            fvMatrix<Type>& eqn,                // Read the FVM matrix representing the transport equation
            const word& fieldName               // Name of the field to which the source is added
        ) const;
        
        // Air limination
        void clampAirFraction() const;
 


public:
    // Runtime type information macro, enabling runtime type identification

        TypeName("degassingBC");

    // Constructors

        // Main constructor initializing the degassingBC model with necessary parameters

        degassingBC     // In OpenFOAM, the order of parameter has to be name; modelType; dict; mesh
        (
            const word& name,               // Name of boundary condition
            const word& modelType,          // Type of model
            const dictionary& dict,         // Dictionary containing model parameters
            const fvMesh& mesh              // Reference to the finite volume mesh
        );

        degassingBC(const degassingBC&) = delete;

        virtual ~degassingBC();


    // Member fucntions

        // Check functions

            // Retrun true if the option adds a source term to the given field's transport equation
            virtual bool addsSupToField(const word& fieldName) const;

            // Return the list of fields for which the option adds source term to the transport equation
            virtual wordList addSupFields() const;

        // Source term addition functions

            // Add a source term to a phase equation with constant density
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_RHO_SUP);

        
        // mesh motion

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();
    
    // member operators 

        // Disallow default bitwise assignment
        void operator=(const degassingBC&) = delete;
};
}
}
#endif