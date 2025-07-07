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
        std::vector<Vector3r> m_vorticity_linear_field;
		std::vector<Vector3r> m_vorticity_current;
		std::vector<Vector3r> m_vorticity_derivative;
		std::vector<Vector3r> m_vorticity_dissipation;
		std::vector<Vector3r> m_stream;
		std::vector<Vector3r> m_delta_velocity;
		Real m_vorticity_refinement_alpha;

		virtual void initParameters();

	public:
        static int IDEAL_VORTICITY_REFINEMENT_ALPHA;
		static int VISCOSITY_VORT_REF;

		VorticityRefinement(FluidModel *model);
		virtual ~VorticityRefinement(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new VorticityRefinement(model); }

		virtual void step();
		virtual void reset();
	};
}

#endif
