#ifndef material_h
#define material_h

#include "femcmpnn.h"

#include "matconst.h"
#include "matstatus.h"
#include "materialmode.h"
#include "timestep.h"
#include "internalstatetype.h"
#include "internalstatevaluetype.h"
#include "matresponsemode.h"
#include "dictionary.h"
#include "chartype.h"

///@name Input fields for Material
//@{
#define _IFT_Material_density "d"
#define _IFT_Material_castingtime "castingtime"
#define _IFT_Material_preCastingTimeMat "precastingtimemat"
//@}

namespace fem 
{
#define STRAIN_STEPS 10.0

class GaussPoint;
class Dictionary;
class FloatArray;
class FloatMatrix;
class Element;
class ProcessCommunicator;

/**
 * Abstract base class for all material models. Declares the basic common interface
 * to all material models. Derived classes should expand this interface, because they are
 * assumed to be  base classes for analysis specific tasks (for example mechanical or
 * thermal analysis).
 *
 * Instance of integration point class is assumed to be implicit argument to
 * all method, depending on internal state in point of consideration.
 * To provide opportunity for storing arbitrary material model related history variables
 * in integration points, associated material status class is introduced.
 * Each new material model class should be declared together with its associated status class
 * (derived from MaterialStatus class). This status can be seen as simple container,
 * storing necessary history variables and providing some access and modification methods.
 * Each integration point can contain material status. Material model should create
 * unique copy of its associated status in each integration point.
 * Because integration point is parameter of all messages to material model
 * class, material model therefore can easily access  all history variables it needs.
 *
 * The attribute 'propertyDictionary' contains all the properties of a material
 * like its Young modulus, its mass density or Poisson ratio.
 *
 * Its task is to indicate whether there required material mode is valid for receiver
 * (method hasMaterialModeCapability). Note: for some material models and linear materials
 * there need not exist support for assembling material char matrix at material level,
 * all is handled properly at crossSection level (_2dBeam mode, 3dShellMode, ...).
 * But this function must indicate whether mode is valid or not for real stress computation.
 *
 * @see MaterialStatus class
 * @see GaussPoint class
 */

/**
* Task:
* give the material constitutive matrix which related the the stress and strain status
* 实现材料本构模型，材料刚度，应力更新
* 历史变量的存储
* 同一材料，在不同维数和单元下，格式不一样，这个如何处理？  材料行为是由单元决定的
*/
class FEM_EXPORT Material : public FEMComponent
{
protected:
   
   
public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    Material(int n, Domain *d);
    /// Destructor.
    virtual ~Material() = default;

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    void printYourself() override;

};
} // end namespace fem
#endif // material_h
