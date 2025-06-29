#include "VorticityRefinement.h"
#include <iostream>
#include "../TimeManager.h"
#include "../Simulation.h"
#include "../Viscosity/ViscosityBase.h"

using namespace SPH;
using namespace GenParam;

int VorticityRefinement::IDEAL_VORTICITY_REFINEMENT_ALPHA = -1;

VorticityRefinement::VorticityRefinement(FluidModel *model) :
	VorticityBase(model)
{
	m_vorticity_linear_field.resize(model->numParticles(), Vector3r::Zero());
    m_vorticity_current.resize(model->numParticles(), Vector3r::Zero());
    m_vorticity_derivative.resize(model->numParticles(), Vector3r::Zero());
    m_vorticity_dissipation.resize(model->numParticles(), Vector3r::Zero());
	m_stream.resize(model->numParticles(), Vector3r::Zero());
	m_vorticity_refinement_alpha = static_cast<Real>(1.0);

	model->addField({ "vorticity_current", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_current[i][0]; }, true });
    model->addField({ "stream", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_stream[i][0]; }, true });
    model->addField({ "vorticity_linear_field", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_linear_field[i][0]; }, true });
    model->addField({ "vorticity_derivative", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_derivative[i][0]; }, true });
    model->addField({ "vorticity_dissipation", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_dissipation[i][0]; }, true });
}

VorticityRefinement::~VorticityRefinement(void)
{
	m_model->removeFieldByName("vorticity_current");
    m_model->removeFieldByName("stream");
    m_model->removeFieldByName("vorticity_linear_field");
    m_model->removeFieldByName("vorticity_derivative");
    m_model->removeFieldByName("vorticity_dissipation");

	m_vorticity_linear_field.clear();
    m_vorticity_current.clear();
    m_vorticity_derivative.clear();
    m_vorticity_dissipation.clear();
    m_stream.clear();
}

void VorticityRefinement::initParameters()
{
	VorticityBase::initParameters();

 	IDEAL_VORTICITY_REFINEMENT_ALPHA = createNumericParameter("vorticityRefinementAlpha", "Ideal Vorticity Refinement", &m_vorticity_refinement_alpha);
 	setGroup(IDEAL_VORTICITY_REFINEMENT_ALPHA, "Fluid Model|Vorticity");
 	setDescription(IDEAL_VORTICITY_REFINEMENT_ALPHA, "Ideal voricity refinment (alpha). Controls the amount of turbulence added to every simulation time step.");
 	RealParameter* rparam = static_cast<RealParameter*>(getParameter(IDEAL_VORTICITY_REFINEMENT_ALPHA));

}

void VorticityRefinement::step()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	FluidModel *model = m_model;

    Real m_v_v;
    ViscosityBase *m_visc = m_model->getViscosityBase();
    if (!m_visc)
        m_v_v = 0.0;
    else
        m_v_v = m_visc->model_viscosity();

	const Real dt = TimeManager::getCurrent()->getTimeStepSize();

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
            // first loop: compute vorticity though linear field and compute dissipation
            
			const Vector3r &xi = m_model->getPosition(i);
			const Vector3r &vi = m_model->getVelocity(i);

            Vector3r &vorticity_linear_field = m_vorticity_linear_field[i];
			vorticity_linear_field.setZero();
            Vector3r &vorticity_dissipation = m_vorticity_dissipation[i];
            vorticity_dissipation.setZero();
            
            Vector3r &vorticity_current = m_vorticity_current[i];
            Vector3r &vorticity_derivative = m_vorticity_derivative[i];
            Vector3r vorticity_equation; //compute vorticity equation
            vorticity_equation = vorticity_current + dt * vorticity_derivative;

            forall_fluid_neighbors_in_same_phase(
                Vector3r &vj = m_model->getVelocity(neighborIndex);
				Real density_j = m_model->getDensity(neighborIndex);
                Real mass_j = m_model->getMass(neighborIndex);

                Vector3r xij = xi - xj;
                Vector3r gradW = sim->gradW(xij);
                //vorticity through linear field
                vorticity_linear_field += (mass_j /density_j) * (vi  - vj).cross(gradW);
			);
            // dissipation of vorticity - diff entre ideal e linear field
            vorticity_dissipation = vorticity_equation - vorticity_linear_field;
		}
        
        #pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
            // 2 loop: compute stream
            const Vector3r &xi = m_model->getPosition(i);

            Vector3r &stream = m_stream[i];
			stream.setZero();

            forall_fluid_neighbors_in_same_phase(

                const Vector3r xij = xi - xj;
                const Real vol_j = m_model->getVolume(neighborIndex);
                // compute stream function
                stream += (1 / (4 * M_PI)) *  (m_vorticity_dissipation[neighborIndex] * vol_j) / xij.norm();             
            ); 
        }    

        #pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
            // 3 loop: delta v u update v
            Vector3r &vi = m_model->getVelocity(i);
            const Vector3r &xi = m_model->getPosition(i);
            Vector3r &stream_i = m_stream[i];

            Vector3r delta_velocity;
            delta_velocity.setZero();

            forall_fluid_neighbors_in_same_phase(
                const Real density_j = m_model->getDensity(neighborIndex);
                const Real mass_j = m_model->getMass(neighborIndex);

                const Vector3r xij = xi - xj;
                const Vector3r gradW = sim->gradW(xij);

                // refinement of linear velocity - delta v
                delta_velocity += (mass_j/density_j) * (stream_i - m_stream[neighborIndex]).cross(gradW);
            );

            // refine linear velocity - nova velocidade
            vi += m_vorticity_refinement_alpha * delta_velocity;
        }

        #pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
            //4 loop: compute final vorticity
            const Vector3r &xi = m_model->getPosition(i);
            Vector3r &vi = m_model->getVelocity(i);
            Vector3r &vorticity_current = m_vorticity_current[i];
            vorticity_current.setZero();


            // update do q eh last vorticity   
            forall_fluid_neighbors_in_same_phase(
                Vector3r &vj = m_model->getVelocity(neighborIndex);
				Real density_j = m_model->getDensity(neighborIndex);
                Real mass_j = m_model->getMass(neighborIndex);

                Vector3r xij = xi - xj;
                Vector3r gradW = sim->gradW(xij);
                //vorticity through linear field
                vorticity_current += (mass_j /density_j) * (vi  - vj).cross(gradW);
            );
        }

        #pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
            //5 loop: compute final delta vorticity

            const Vector3r &xi = m_model->getPosition(i);
            Vector3r &vi = m_model->getVelocity(i);
            Vector3r &vorticity_current = m_vorticity_current[i];
            Vector3r &vorticity_derivative = m_vorticity_derivative[i];
            vorticity_derivative.setZero();

            Vector3r gradV_x;
            gradV_x.setZero();
            Vector3r gradV_y;
            gradV_y.setZero();
            Vector3r gradV_z;
            gradV_z.setZero();

            Vector3r vorticity_rate_gradient;
            vorticity_rate_gradient.setZero();
            Vector3r vorticity_rate_laplacian;
            vorticity_rate_laplacian.setZero();

            // update do q eh last vorticity   
            forall_fluid_neighbors_in_same_phase(
                Vector3r &vj = m_model->getVelocity(neighborIndex);
                Real density_j = m_model->getDensity(neighborIndex);
                Real mass_j = m_model->getMass(neighborIndex);

                Vector3r xij = xi - xj;
                Vector3r gradW = sim->gradW(xij);
                const Vector3r vort_ij = vorticity_current - m_vorticity_current[neighborIndex];
                
                // compute vorticity equation - computar vorticity ideal
                //  vorticity * gradV + v_v * laplacian(vorticity)
                gradV_x += (mass_j / density_j) * (vj.x() - vi.x()) * (gradW);
                gradV_y += (mass_j / density_j) * (vj.y() - vi.y()) * (gradW);
                gradV_z += (mass_j / density_j) * (vj.z() - vi.z()) * (gradW);

                vorticity_rate_laplacian += 2 * (d + 2) * m_v_v * (mass_j / density_j) * ((vort_ij.dot(xij)) / (xij.dot(xij) + 0.01 * h2)) * gradW;
            );

            vorticity_rate_gradient.x() = vorticity_current.dot(gradV_x);
            vorticity_rate_gradient.y() = vorticity_current.dot(gradV_y);
            vorticity_rate_gradient.z() = vorticity_current.dot(gradV_z);

            vorticity_derivative = vorticity_rate_gradient + vorticity_rate_laplacian;  

        }
	}
}


void VorticityRefinement::reset()
{
	for (unsigned int i = 0; i < m_model->numParticles(); i++){
        m_vorticity_linear_field[i].setZero();
        m_vorticity_current[i].setZero();
        m_vorticity_derivative[i].setZero();
        m_vorticity_dissipation[i].setZero();
        m_stream[i].setZero();
    }
}