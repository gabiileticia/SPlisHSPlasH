#include "VorticityRefinement.h"
#include "../TimeManager.h"
#include "../Simulation.h"

using namespace SPH;

int MicropolarModel_Bender2017::IDEAL_VORTICITY_REFINEMENT_ALPHA = -1;

VorticityRefinement::VorticityRefinement(FluidModel *model) :
	VorticityBase(model)
{
    //////////////////////////////////////////////////
    // essa aqui eh um constructor
    //////////////////////////////////////////////////
	m_vorticity.resize(model->numParticles(), Vector3r::Zero());
	m_stream.resize(model->numParticles(), Vector3r::Zero());
    last_vorticity.resize(model->numParticles(), Vector3r::Zero());
    last_gradVelocity.resize(model->numParticles(), Vector3r::Zero());
	m_vorticityRefinementAlpha = static_cast<Real>(1.0);
    m_viscosityKinematic = static_cast<Real>(0.05); // try to use the one from the simulation 

	model->addField({ "vorticity", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity[i][0]; }, true });
    model->addField({ "stream", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_stream[i][0]; }, true });
}

VorticityRefinement::~VorticityRefinement(void)
{
    //////////////////////////////////////////////////
    // essa aqui eh um destructor
    //////////////////////////////////////////////////
	m_model->removeFieldByName("vorticity");
    m_model->removeFieldByName("stream");

	m_vorticity.clear();
    m_stream.clear();
    last_vorticity.clear();
    last_gradVelocity.clear();
}

void VorticityRefinement::initParameters()
{
    //////////////////////////////////////////////////
    // essa aqui eh para setup de parametros inicias (acho q da interface)
    //////////////////////////////////////////////////
	VorticityBase::initParameters();

 	IDEAL_VORTICITY_REFINEMENT_ALPHA = createNumericParameter("vorticityRefinementAlpha", "Ideal Vorticity Refinement", &m_vorticityRefinementAlpha);
 	setGroup(IDEAL_VORTICITY_REFINEMENT_ALPHA, "Fluid Model|Vorticity");
 	setDescription(IDEAL_VORTICITY_REFINEMENT_ALPHA, "Ideal voricity refinment (alpha). Controls the amount of turbulence added to every simulation time step.");
 	RealParameter* rparam = static_cast<RealParameter*>(getParameter(IDEAL_VORTICITY_REFINEMENT_ALPHA));
 	// rparam->setMinValue(0.0); they say alpha is in R
}

void VorticityRefinement::step()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	if (numParticles == 0)
		return;

	FluidModel *model = m_model;
	const Real density0 = model->getDensity0();

	const Real dt = TimeManager::getCurrent()->getTimeStepSize();
	// const Real invDt = static_cast<Real>(1.0) / dt;

	const Real h = sim->getSupportRadius();
	const Real h2 = h*h;

    Real d = 3.0;
	if (sim->is2DSimulation())
		d = 2.0;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
            
			const Vector3r &xi = m_model->getPosition(i);
			const Vector3r &vi = m_model->getVelocity(i);

            Vector3r &last_vorticity = last_vorticity[i];
            Vector3r &vorticity = m_vorticity[i];
			vorticity.setZero();

            Vector3r &stream = m_stream[i];
			stream.setZero();

            Vector3r gradV;
            gradV.setZero();

            Vector3r vortEquationGradient;
            vortEquationGradient.setZero();

            Vector3r vortEquationLaplacian;
            vortEquationLaplacian.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
            forall_fluid_neighbors_in_same_phase(
                const Vector3r &vj = m_model->getVelocity(neighborIndex);
				const Real density_j = m_model->getDensity(neighborIndex);
                const Real mass_j = m_model->getMass(neighborIndex);
                //vorticity through linear field
                vorticity += (mass_j /density_j) * (vi  - vj).cross(gradW);
            );

			forall_fluid_neighbors_in_same_phase(
                
                const Vector3r &xj = m_model->getPosition(neighborIndex);
                const Vector3r &vj = m_model->getVelocity(neighborIndex);
				const Real density_j = m_model->getDensity(neighborIndex);
                const Real mass_j = m_model->getMass(neighborIndex);

                const Vector3r xij = xi - xj;
                const Vector3r gradW = sim->gradW(xij);
                
                // compute vorticity equation - computar vorticity ideal
                //  vorticity * gradV + v_v * laplacian(vorticity)
                gradV += (mass_j/density_j) * (vj - vi) * (gradW);

                vortEquationGradient[0] = last_vorticity.dot(gradV[0]);
                vortEquationGradient[1] = last_vorticity.dot(gradV[1]);
                vortEquationGradient[2] = last_vorticity.dot(gradV[2]);
                
                const Vector3r vort_ij = last_vorticity[i] - last_vorticity[neighborIndex];
                vortEquationLaplacian += 2 * (d + 2) * m_viscosityKinematic * (mass_j/density_j) *  ((vort_ij * xij) / (xij * xij + 0.01 * h2)) * gradW;

                const Vector3r dVorticity = vortEquationGradient + vortEquationLaplacian;
                
                const Vector3r vorticity_ideal = last_vorticity + dt * dVorticity;

                // dissipation of vorticity - diff entre ideal e linear field
                const Vector3r dissipation = vorticity_ideal - vorticity;

			);
            // update do q eh last vorticity
            last_vorticity[i] = vorticity;

            forall_fluid_neighbors_in_same_phase(

                const Vector3r &xj = m_model->getVelocity(neighborIndex);
                const Vector3r xij = xi - xj;
                const Real vol_j = m_model->getVolume(neighborIndex)ç

                // compute stream function
                stream += (1 / (4 * M_PI)) *  (dissipation * vol_j) / std::abs(xij);             
            );     

		}

        #pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{

            Vector3r &vi = m_model->getVelocity(i);
            const Vector3r &xi = m_model->getPosition(i);
            Vector3r &stream_i = stream[i];

            Vector3r deltaVelocity;
            deltaVelocity.setZero();

            forall_fluid_neighbors_in_same_phase(

                const Real density_j = m_model->getDensity(neighborIndex);
                const Real mass_j = m_model->getMass(neighborIndex);

                const Vector3r &xj = m_model->getPosition(neighborIndex);

                const Vector3r xij = xi - xj;
                const Vector3r gradW = sim->gradW(xij);

                // refinement of linear velocity - delta v
                deltaVelocity += (mass_j/density_j) * (stream_i - stream[neighborIndex]).cross(gradW);
            );

            // refine linear velocity - nova velocidade
            vi += m_vorticityRefinementAlpha * deltaVelocity;
        }
	}
}


void VorticityRefinement::reset()
{
	for (unsigned int i = 0; i < m_model->numParticles(); i++)
        m_vorticity[i].setZero();
        m_stream[i].setZero();
        last_vorticity[i].setZero();
        last_gradVelocity[i].setZero();
}


// do i need to think of neighboorhod search


