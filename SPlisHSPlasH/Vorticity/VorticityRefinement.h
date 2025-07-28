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
		std::vector<Vector3r> m_vorticity_advected;
		std::vector<Vector3r> m_vorticity_corrected_end;
		std::vector<Vector3r> m_vorticity_derivative;
		std::vector<Vector3r> m_vorticity_dissipation;
		std::vector<Vector3r> m_stream;
		std::vector<Vector3r> m_delta_velocity;
		std::vector<Vector3r> m_vorticity_equation;
		std::vector<Vector3r> m_a_adv;
		std::vector<Vector3r> m_v_adv;
		std::vector<Vector3r> m_velocity_from_dfsph;
		std::vector<Vector3r> m_velocity_corrected_end;
		std::vector<Vector3r> m_position_from_dfsph;
		std::vector<Vector3r> m_vorticity_rate_laplacian;
		std::vector<Vector3r> m_vorticity_rate_gradient;
		std::vector<Vector3r> m_gradV_x;
		std::vector<Vector3r> m_gradV_y;
		std::vector<Vector3r> m_gradV_z;
		std::vector<Vector3r> m_last_acceleration;
		std::vector<Real> m_direction_velocity;

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
