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
		std::vector<Vector3r> m_vorticity_init;
		std::vector<Vector3r> m_vorticity_corrected_end;
		std::vector<Vector3r> m_vorticity_derivative;
		std::vector<Vector3r> m_vorticity_dissipation;
		std::vector<Vector3r> m_vorticity_final_dissipation;
		std::vector<Real> m_lambdatwo;
		std::vector<Real> m_is_vortex;
		std::vector<Vector3r> m_stream;
		std::vector<Vector3r> m_stream2d;
		std::vector<Vector3r> m_stream3d;
		std::vector<Vector3r> m_stream2d_vol;
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
		std::vector<Real> m_direction_vort_dev;
		std::vector<Vector3r> m_velocity_advected;
		std::vector<Vector3r> m_acceleration_vr;
		std::vector<Vector3r> m_position_advected;
		Real m_vorticity_refinement_alpha;
		Real m_v_v;
		Real m_flag;

		virtual void initParameters();

	public:
		static int IDEAL_VORTICITY_REFINEMENT_ALPHA;
		static int VISCOSITY_VORT_REF;
		static int VORTICITY_FLAG;

		VorticityRefinement(FluidModel *model);
		virtual ~VorticityRefinement(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new VorticityRefinement(model); }

		virtual void step();
		virtual void reset();

		virtual void performNeighborhoodSearchSort();
	};
}

#endif
