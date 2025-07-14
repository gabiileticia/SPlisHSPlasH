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
		std::vector<Vector3r> m_vorticity_equation;
		std::vector<Vector3r> m_v_pre;
		std::vector<Vector3r> m_v_post;
		std::vector<Vector3r> m_x_pre;
		std::vector<Vector3r> m_vorticity_rate_laplacian;
		std::vector<Vector3r> m_vorticity_rate_gradient;
		std::vector<Vector3r> m_gradV_x;
		std::vector<Vector3r> m_gradV_y;
		std::vector<Vector3r> m_gradV_z;

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

		virtual void performNeighborhoodSearchSort();
	};
}

#endif
