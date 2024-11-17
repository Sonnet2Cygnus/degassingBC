#include "degassingBC.H"
#include "fvMatrix.H"
#include "fvmSup.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"
#include "sizeGroup.H"


namespace Foam
{
namespace fv
{
    // Static Data Members

        // define type name and debug level
        defineTypeNameAndDebug(degassingBC, 0);

        // add degassingBC cloass to the runtime
        addToRunTimeSelectionTable
        (
            fvModel,
            degassingBC,
            dictionary
        );
}
}

// check if the phase is a degassing phase
bool Foam::fv::degassingBC::isDegassingPhases(const word& fieldName) const
{
    return pos0(findIndex(degassingBC_, IOobject::group(fieldName)));
}

// calculate the volume of cells attached to patches
Foam::scalar
Foam::fv::degassingBC::patchCellVolume(const wordList& patches) const
{
    scalar V(0.0);

    forAll(patches, patchi)
    {
        const fvPatch& patchField = mesh().boundary()[patches[patchi]];

        forAll(patchField.faceCells(), i)
        {
            V += mesh().V()[patchField.faceCells()[i]];
        }
    }

    reduce(V, sumOp<scalar>());

    return V;
}

// add a source term to a phase equation with constant density
template<class Type>
void Foam::fv::degassingBC::addSupType
(
    const volScalarField& alpha,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    const VolField<Type>& psi = eqn.psi();
    const scalar rDeltaT = 1.0/mesh().time().deltaT().value();

    // 创建内部源项场 create internal source term field
    volScalarField::Internal Sp
    (
        IOobject
        (
            name() + ":" + fieldName + ":Sp",
            mesh().time().IOobject::name(),  // 这个部分似乎有多重歧义,保留
            mesh()
        ),
        mesh(),
        dimensionedScalar(dimless/dimTime, Zero),
        false
    );

    if (isDegassingPhases(fieldName))
    {
        // 对每个单元格添加源项  add source term for each cell
        forAll(cells_, i)
        {
            Sp[cells_[i]] -= alpha[cells_[i]]*rDeltaT;
        }
    }

    // 调试信息输出 debug information output
    if (debug)
    {
        Info << name() << ": " << fieldName << ": Sp min, max: "
             << min(Sp).value() << ", " << max(Sp).value() << endl;

        if (mesh().time().writeTime()) Sp.write();
    }

    // add source term to the equation 
    eqn += fvm::SuSp(Sp, psi);

    clampAirFraction();

}

void Foam::fv::degassingBC::clampAirFraction() const
{
    // 使用 const 来保证从 lookupObject 得到的是常量引用
      volScalarField& alphaAirConst = const_cast<volScalarField&>(mesh().lookupObject<volScalarField>(alphaName_));

    // creat  a copt
    volScalarField alphaAir = alphaAirConst;

    // limit the alphaair in 0-1
    alphaAir = max(alphaAir, 0.0);
    alphaAir = min(alphaAir, 1.0);

    // 输出限制后的最小值和最大值以进行验证
    alphaAirConst = alphaAir;

}



// Constructors 构造函数
namespace Foam
{
namespace fv
{
    degassingBC::degassingBC
    (
        const word& name,
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh
    )
   :
        fvModel(name, modelType, dict, mesh),
        degassingBC_(dict.lookup<wordList>("degassingBC")),
        continuousPhase_(dict.lookupOrDefault<word>("continuousPhase", "default")),
        alphaName_(dict.lookupOrDefault<word>("alphaName", "alpha")),
        fluid_(mesh.lookupObject<phaseSystem>(phaseSystem::propertiesName)),
        patchNames_(dict.lookup<wordList>("patchNames")),
        cells_(dict.lookupOrDefault<labelList>("cells", labelList())),
        V_(patchCellVolume(patchNames_))
    {
        read(dict);

        if (continuousPhase_ == "default")
        {
            // 使用正确的作用域引用 phaseModelList
            const phaseSystem::phaseModelList& allPhases = fluid_.phases();

            // 确保系统是两相系统
            if (allPhases.size() != 2)
            {
                FatalErrorInFunction
                    << "Exactly two phases are required when continuousPhase is set to default."
                    << endl << exit(FatalError);
            }

            // 获取退化相
            const phaseModel& degassingPhase = allPhases[0];
            const phaseModel& continuousPhase = fluid_.otherPhase(degassingPhase);

            // 确保 continuousPhase 存在
            fieldNames_.append(continuousPhase.rho()().name());
            continuousPhase_ = continuousPhase.name();
        }
        else
        {
            // 查找连续相
            const phaseSystem::phaseModelList& allPhases = fluid_.phases();
            label continuousPhaseIndex = -1;

            forAll(allPhases, i)
            {
                if (allPhases[i].name() == continuousPhase_)
                {
                    continuousPhaseIndex = i;
                    break;
                }
            }

            if (continuousPhaseIndex == -1)
            {
                FatalErrorInFunction
                    << "Specified continuous phase '" << continuousPhase_
                    << "' not found in phase system."
                    << endl << exit(FatalError);
            }

            const phaseModel& continuousPhaseModel = allPhases[continuousPhaseIndex];
            fieldNames_.append(continuousPhaseModel.rho()().name());
        }

        // list of solution variable to apply correction source term 
        const wordList solutionVariables
        (
            wordList{"U", "k", "R", "epsilon", "omega"}
        );

        forAll(degassingBC_, phasei)
        {
            const phaseModel& phase = fluid_.phases()[degassingBC_[phasei]];

            // Add phase density field
            fieldNames_.append(phase.rho()().name());
            forAll(solutionVariables, variablei)
            {
                const word variable = solutionVariables[variablei];
                const word fieldName =
                    IOobject::groupName(variable, phase.name());

                if (pos0(findIndex(mesh.names(), fieldName)))
                {
                    fieldNames_.append(fieldName);
                }
            }

            //  Add field of solution variable for energy
            if (!phase.isothermal())
            {
                fieldNames_.append(phase.thermo().he().name());
            }

            // Add size group fields
            if (isA<diameterModels::velocityGroup>(phase.dPtr()()))
            {
                const diameterModels::velocityGroup& velGroup=
                    refCast<const diameterModels::velocityGroup>
                    (
                        phase.dPtr()()
                    );

                forAll(velGroup.sizeGroups(), i)
                {
                    const diameterModels::sizeGroup* fi =
                        velGroup.sizeGroups()(i);

                    fieldNames_.append(fi->name());
                }
            }
        }

        if (debug)
        {
            word phases = degassingBC_[0];
            for (Foam::label i=1; i<degassingBC_.size(); i++)
            {
                phases += ", " + degassingBC_[i];
            }
            Info<< "    Degassing Phases: " << phases << nl
                << "    Continuous phase: " << continuousPhase_ << endl;

        }

        forAll(patchNames_, patchi)
        {
            const fvPatch& patchField = mesh.boundary()[patchNames_[patchi]];

            forAll(patchField.faceCells(), i)
            {
                cells_.append(patchField.faceCells()[i]);
            }
        }

        V_ = patchCellVolume(patchNames_);
    }
}
}

// destuctor
Foam::fv::degassingBC::~degassingBC()
{}

// member functions
bool Foam::fv::degassingBC::addsSupToField(const word& fieldName) const
{
    return findIndex(fieldNames_, fieldName) != -1;
}

Foam::wordList Foam::fv::degassingBC::addSupFields() const
{
    return fieldNames_;
}

FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_RHO_SUP,
    fv::degassingBC
);

void Foam::fv::degassingBC::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::degassingBC::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::degassingBC::distribute(const polyDistributionMap&)
{}


bool Foam::fv::degassingBC::movePoints()
{
    V_ = patchCellVolume(patchNames_);

    return true;
}
