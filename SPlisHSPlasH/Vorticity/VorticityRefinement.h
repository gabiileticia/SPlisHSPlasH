#ifndef __VorticityRefinement_h__
#define __VorticityRefinement_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "VorticityBase.h"

namespace SPH
{
	/** \brief This class implements the vorticity refinement model introduced
	* by Liu et al. [XXXXX]
	*
	* References:
	* - 
	*/
	class VorticityRefinement : public VorticityBase
	{    
	protected:
        std::vector<Vector3r> m_vorticity;
		// Real m_vorticityCoeff;

		// virtual void initParameters();

	public:
        static int IDEAL_VORTICITY_REFINEMENT_ALPHA;

		VorticityRefinement(FluidModel *model);
		virtual ~VorticityRefinement(void);
	};
}

#endif
