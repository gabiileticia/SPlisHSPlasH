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
    m_vorticity_init.resize(model->numParticles(), Vector3r::Zero());
    m_vorticity_corrected_end.resize(model->numParticles(), Vector3r::Zero());
    m_vorticity_derivative.resize(model->numParticles(), Vector3r::Zero());
    m_vorticity_dissipation.resize(model->numParticles(), Vector3r::Zero());
    m_stream.resize(model->numParticles(), Vector3r::Zero());
    m_stream2d.resize(model->numParticles(), Vector3r::Zero());
    m_stream3d.resize(model->numParticles(), Vector3r::Zero());
    m_stream2d_vol.resize(model->numParticles(), Vector3r::Zero());
    m_delta_velocity.resize(model->numParticles(), Vector3r::Zero());
    m_vorticity_equation.resize(model->numParticles(), Vector3r::Zero());
    m_v_adv.resize(model->numParticles(), Vector3r::Zero());
    m_velocity_from_dfsph.resize(model->numParticles(), Vector3r::Zero());
    m_velocity_corrected_end.resize(model->numParticles(), Vector3r::Zero());
    m_position_from_dfsph.resize(model->numParticles(), Vector3r::Zero());
    m_vorticity_rate_laplacian.resize(model->numParticles(), Vector3r::Zero());
    m_vorticity_rate_gradient.resize(model->numParticles(), Vector3r::Zero());
    m_gradV_x.resize(model->numParticles(), Vector3r::Zero());
    m_gradV_y.resize(model->numParticles(), Vector3r::Zero());
    m_gradV_z.resize(model->numParticles(), Vector3r::Zero());
    m_direction_velocity.resize(model->numParticles(), 0.0);

    m_vorticity_refinement_alpha = static_cast<Real>(1.0);
    
    model->addField({ "m_vorticity_init", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_init[i][0]; }, true });
    model->addField({ "vorticity_corrected_end", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_corrected_end[i][0]; }, true });
    model->addField({ "stream", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_stream[i][0]; }, true });
    model->addField({ "m_stream2d", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_stream2d[i][0]; }, true });
    model->addField({ "m_stream3d", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_stream3d[i][0]; }, true });
    model->addField({ "m_stream2d_vol", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_stream2d_vol[i][0]; }, true });
    model->addField({ "vorticity_linear_field", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_linear_field[i][0]; }, true });
    model->addField({ "vorticity_derivative", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_derivative[i][0]; }, true });
    model->addField({ "vorticity_dissipation", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_dissipation[i][0]; }, true });
    model->addField({ "delta_velocity", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_delta_velocity[i][0]; }, true });
    model->addField({ "m_vorticity_equation", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_equation[i][0]; }, true });
    model->addField({ "m_velocity_from_dfsph", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_velocity_from_dfsph[i][0]; }, true });
    model->addField({ "m_velocity_corrected_end", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_velocity_corrected_end[i][0]; }, true });
    model->addField({ "m_position_from_dfsph", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_position_from_dfsph[i][0]; }, true });
    model->addField({ "m_vorticity_rate_laplacian", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_rate_laplacian[i][0]; }, true });
    model->addField({ "m_vorticity_rate_gradient", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_rate_gradient[i][0]; }, true });
    model->addField({ "m_gradV_x", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_gradV_x[i][0]; }, true });
    model->addField({ "m_gradV_y", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_gradV_y[i][0]; }, true });
    model->addField({ "m_gradV_z", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_gradV_z[i][0]; }, true });
	model->addField({ "m_direction_velocity", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_direction_velocity[i]; }, true });
}

VorticityRefinement::~VorticityRefinement(void)
{
    m_model->removeFieldByName("m_vorticity_init");
    m_model->removeFieldByName("vorticity_corrected_end");
    m_model->removeFieldByName("stream");
    m_model->removeFieldByName("vorticity_linear_field");
    m_model->removeFieldByName("vorticity_derivative");
    m_model->removeFieldByName("vorticity_dissipation");
    m_model->removeFieldByName("delta_velocity");
    m_model->removeFieldByName("m_vorticity_equation");
    m_model->removeFieldByName("m_v_adv");
    m_model->removeFieldByName("m_velocity_from_dfsph");
    m_model->removeFieldByName("m_velocity_corrected_end");
    m_model->removeFieldByName("m_position_from_dfsph");
    m_model->removeFieldByName("m_vorticity_rate_laplacian");
    m_model->removeFieldByName("m_vorticity_rate_gradient");
    m_model->removeFieldByName("m_gradV_x");
    m_model->removeFieldByName("m_gradV_y");
    m_model->removeFieldByName("m_gradV_z");
    m_model->removeFieldByName("m_direction_velocity");


    m_vorticity_linear_field.clear();
    m_vorticity_init.clear();
    m_vorticity_corrected_end.clear();
    m_vorticity_derivative.clear();
    m_vorticity_dissipation.clear();
    m_stream.clear();
    m_delta_velocity.clear();
    m_vorticity_equation.clear();
    m_v_adv.clear();
    m_velocity_from_dfsph.clear();
    m_velocity_corrected_end.clear();
    m_position_from_dfsph.clear();
    m_vorticity_rate_laplacian.clear();
    m_vorticity_rate_gradient.clear();
    m_gradV_x.clear();
    m_gradV_y.clear();
    m_gradV_z.clear();
    m_direction_velocity.clear();
}

void VorticityRefinement::initParameters()
{
    VorticityBase::initParameters();
    IDEAL_VORTICITY_REFINEMENT_ALPHA = createNumericParameter("vorticityRefinementAlpha", "Ideal Vorticity Refinement", &m_vorticity_refinement_alpha);
    setGroup(IDEAL_VORTICITY_REFINEMENT_ALPHA, "Fluid Model|Vorticity");
    setDescription(IDEAL_VORTICITY_REFINEMENT_ALPHA, "Ideal voricity refinment (alpha). Controls the amount of turbulence added to every simulation time step.");
    RealParameter* rparam = static_cast<RealParameter*>(getParameter(IDEAL_VORTICITY_REFINEMENT_ALPHA));

}
/*
#ifdef USE_AVX

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
    Real vorticity_refinement_alpha = m_vorticity_refinement_alpha;


    Real d = 3.0;
    if (sim->is2DSimulation()) {
        d = 2.0;
    }

    #pragma omp parallel default(shared)
    {
        
        #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            // 1st loop: compute vorticity though linear field and compute dissipation
            //compute linear field vorticity based on velocity from dfsph on current timestep

            Vector3r &xi = m_model->getPosition(i);
            Vector3r &vi = m_model->getVelocity(i);

            const Vector3f8 xi_avx(xi);
			const Vector3f8 vi_avx(vi);
            
            //saving initial velocity for analysis////////
            m_velocity_from_dfsph[i] = vi;
            /////////////////////////////////////////////

            Vector3r& vorticity_linear_field = m_vorticity_linear_field[i];
            vorticity_linear_field.setZero();
            Vector3f8 vorticity_linear_field_avx;
			vorticity_linear_field_avx.setZero();
            
            forall_fluid_neighbors_in_same_phase_avx(
                const Scalarf8 Vj_avx = convert_zero(model->getVolume(0), count);
				compute_Vj_gradW_samephase();

                const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &model->getVelocity(0), count);
                
                //vorticity through linear field
                vorticity_linear_field_avx -= (vi_avx - vj_avx) % V_gradW;
            );

            vorticity_linear_field[0] += vorticity_linear_field_avx.x().reduce();
			vorticity_linear_field[1] += vorticity_linear_field_avx.y().reduce();
			vorticity_linear_field[2] += vorticity_linear_field_avx.z().reduce();

            m_vorticity_equation[i] = m_vorticity_corrected_end[i] + dt * m_vorticity_derivative[i];
            m_vorticity_dissipation[i] = m_vorticity_equation[i] - vorticity_linear_field;
        }

        #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            // 2nd loop: compute stream
            Vector3r& xi = m_model->getPosition(i);
            const Vector3f8 xi_avx(xi);
    
            Vector3r &stream = m_stream[i];
            stream.setZero();

            Vector3f8 stream_avx;
			stream_avx.setZero();
    
            forall_fluid_neighbors_in_same_phase_avx(
                const Scalarf8 Vj_avx = convert_zero(model->getVolume(0), count);
				
                Vector3f8 xij = xi_avx - xj_avx;

                const Vector3f8 vort_diss_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &m_vorticity_dissipation[0], count);
                // compute stream function
                stream_avx += (vort_diss_avx) * (Vj_avx / (xij.norm()));
            ); 

            stream[0] += stream_avx.x().reduce();
			stream[1] += stream_avx.y().reduce();
			stream[2] += stream_avx.z().reduce();

            stream = stream * 0.25 / M_PI;
        }  

        #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            // compute delta v and update v
            Vector3r &vi = m_model->getVelocity(i);
            Vector3r &xi = m_model->getPosition(i);
            Vector3r &stream_i = m_stream[i];

            const Vector3f8 xi_avx(xi);
			const Vector3f8 vi_avx(vi);
            const Vector3f8 stream_i_avx(m_stream[i]);
    
            Vector3r &delta_velocity = m_delta_velocity[i];
            delta_velocity.setZero();

            Vector3f8 delta_velocity_avx;
			delta_velocity_avx.setZero();
    
            forall_fluid_neighbors_in_same_phase_avx(
                const Scalarf8 Vj_avx = convert_zero(model->getVolume(0), count);
				compute_Vj_gradW_samephase();

                const Vector3f8 stream_j_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &m_stream[0], count);
    
                // compute delta velocity
                delta_velocity_avx -= (stream_i_avx - stream_j_avx) % V_gradW;
            );
            delta_velocity[0] += delta_velocity_avx.x().reduce();
			delta_velocity[1] += delta_velocity_avx.y().reduce();
			delta_velocity[2] += delta_velocity_avx.z().reduce();
    
            m_direction_velocity[i] = vi.dot(delta_velocity);
            // refine linear velocity
            vi += vorticity_refinement_alpha * delta_velocity;

            // saving final velocity
            m_velocity_corrected_end[i] = vi;
        }

        #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            //4th loop: compute final vorticity
            Vector3r &xi = m_model->getPosition(i);
            Vector3r &vi = m_model->getVelocity(i);

            const Vector3f8 xi_avx(xi);
			const Vector3f8 vi_avx(vi);

            Vector3r& vorticity_corrected_end = m_vorticity_corrected_end[i];
            vorticity_corrected_end.setZero();
            Vector3f8 vorticity_corrected_end_avx;
			vorticity_corrected_end_avx.setZero();
            
            forall_fluid_neighbors_in_same_phase_avx(
                const Scalarf8 Vj_avx = convert_zero(model->getVolume(0), count);
				compute_Vj_gradW_samephase();

                const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &model->getVelocity(0), count);
                
                //vorticity through linear field
                vorticity_corrected_end_avx -= (vi_avx - vj_avx) % V_gradW;
            );

            vorticity_corrected_end[0] += vorticity_corrected_end_avx.x().reduce();
			vorticity_corrected_end[1] += vorticity_corrected_end_avx.y().reduce();
			vorticity_corrected_end[2] += vorticity_corrected_end_avx.z().reduce();
        }

        #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            //5th loop: compute vorticity derivative based on new vorticity

            Vector3r &xi = m_model->getPosition(i);
            Vector3r &vi = m_model->getVelocity(i);
            
            const Vector3f8 xi_avx(xi);
			const Vector3f8 vi_avx(vi);
            
            Vector3r &vorticity_corrected_end = m_vorticity_corrected_end[i];
            const Vector3f8 vorticity_corrected_end_i_avx(m_vorticity_corrected_end[i]);

            Vector3r &vorticity_derivative = m_vorticity_derivative[i];
            vorticity_derivative.setZero();

            Vector3f8 vorticity_derivative_avx;
			vorticity_derivative_avx.setZero();
    
            Vector3r &gradV_x = m_gradV_x[i];
            gradV_x.setZero();
            Vector3r &gradV_y = m_gradV_y[i];
            gradV_y.setZero();
            Vector3r &gradV_z = m_gradV_z[i];
            gradV_z.setZero();

            Vector3f8 gradV_x_avx;
			gradV_x_avx.setZero();
            Vector3f8 gradV_y_avx;
			gradV_y_avx.setZero();
            Vector3f8 gradV_z_avx;
			gradV_z_avx.setZero();
    
            Vector3r &vorticity_rate_laplacian = m_vorticity_rate_laplacian[i];
            vorticity_rate_laplacian.setZero();

            Vector3f8 vorticity_rate_laplacian_avx;
			vorticity_rate_laplacian_avx.setZero();
    
            forall_fluid_neighbors_in_same_phase_avx(
                const Scalarf8 Vj_avx = convert_zero(model->getVolume(0), count);
				compute_Vj_gradW_samephase();

				const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &model->getVelocity(0), count);
                Vector3f8 xij = xi_avx - xj_avx;
                
                const Vector3f8 vorticity_corrected_end_j_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &m_vorticity_corrected_end[0], count);
                Vector3f8 vort_ij = vorticity_corrected_end_i_avx - vorticity_corrected_end_j_avx;
                    
                //  vorticity * gradV + v_v * laplacian(vorticity)
                gradV_x_avx += V_gradW * (vj_avx.x() - vi_avx.x());
                gradV_y_avx += V_gradW * (vj_avx.y() - vi_avx.y());
                gradV_z_avx += V_gradW * (vj_avx.z() - vi_avx.z());
    
                vorticity_rate_laplacian_avx += V_gradW * ((vort_ij.dot(xij)) / (xij.squaredNorm() + 0.01 * h2));
            );

            gradV_x[0] += gradV_x_avx.x().reduce();
			gradV_x[1] += gradV_x_avx.y().reduce();
			gradV_x[2] += gradV_x_avx.z().reduce();

            gradV_y[0] += gradV_y_avx.x().reduce();
			gradV_y[1] += gradV_y_avx.y().reduce();
			gradV_y[2] += gradV_y_avx.z().reduce();

            gradV_z[0] += gradV_z_avx.x().reduce();
			gradV_z[1] += gradV_z_avx.y().reduce();
			gradV_z[2] += gradV_z_avx.z().reduce();

            vorticity_rate_laplacian[0] += vorticity_rate_laplacian_avx.x().reduce();
			vorticity_rate_laplacian[1] += vorticity_rate_laplacian_avx.y().reduce();
			vorticity_rate_laplacian[2] += vorticity_rate_laplacian_avx.z().reduce();
            vorticity_rate_laplacian = vorticity_rate_laplacian * 2 * (d + 2) * m_v_v;

            Vector3r &vorticity_rate_gradient = m_vorticity_rate_gradient[i];
            vorticity_rate_gradient.setZero();

            vorticity_rate_gradient.x() = vorticity_corrected_end.dot(gradV_x);
            vorticity_rate_gradient.y() = vorticity_corrected_end.dot(gradV_y);
            vorticity_rate_gradient.z() = vorticity_corrected_end.dot(gradV_z);
    
            vorticity_derivative = vorticity_rate_gradient + vorticity_rate_laplacian;
        }
    }
}


#else
*/

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

    TimeManager* tm = TimeManager::getCurrent();
    const Real dt = tm->getTimeStepSize();
    const Real t = tm->getTime();

    const Real h = sim->getSupportRadius();
    const Real h2 = h*h;
    const Real r = sim->getParticleRadius();
    const Real r2 = r * r;
    const Real vorticity_refinement_alpha = m_vorticity_refinement_alpha;
    


    Real d = 3.0;
    if (sim->is2DSimulation()) {
        d = 2.0;
    }

  /*  Vector3r a;
    Vector3r b;
    Vector3r dd;

    dd.x() = -1;
    dd.y() = -1;
    dd.z() = -2;

    a.x() = 2;
    a.y() = 2;
    a.z() = 2;

    b.x() = 3;
    b.y() = 4;
    b.z() = 3;
    std::cout << a.cross(b) << "\n";

    Vector3r c;

    c.x() = 0.02;
    c.y() = 0.02;
    c.z() = 0.02;

    std::cout << sim->gradW(c) << "\n";
    std::cout << dd.cross(sim->gradW(c)) << "\n";*/
    
    /*if (sim->is2DSimulation())
    {
         #pragma omp parallel default(shared)
         {
             #pragma omp for schedule(static)  
            for (int i = 0; i < (int)numParticles; i++)
            {
                Vector3r &vi = m_model->getVelocity(i);
                Vector3r& xi = m_model->getPosition(i);
                vi.z() = 0.0;
                xi.z() = 0.0;
            }
         }
    }*/

     #pragma omp parallel default(shared)
     {
        #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            // 1st loop: compute vorticity though linear field and compute dissipation
            
            //compute linear field vorticity based on velocity from dfsph on current timestep
            const Vector3r &xi = m_model->getPosition(i);
            const Vector3r &vi = m_model->getVelocity(i);
            
            //saving initial velocity for analysis////////
            m_velocity_from_dfsph[i] = vi;
            m_vorticity_init[i] = m_vorticity_corrected_end[i];
            /////////////////////////////////////////////

            Vector3r &vorticity_linear_field = m_vorticity_linear_field[i];
            vorticity_linear_field.setZero();
            
            forall_fluid_neighbors_in_same_phase(
                const Vector3r &vj = m_model->getVelocity(neighborIndex);
                const Real density_j = m_model->getDensity(neighborIndex);
                const Real mass_j = m_model->getMass(neighborIndex);
                
                const Vector3r xij = xi - xj;
                const Vector3r gradW = sim->gradW(xij);
                //vorticity through linear field
                vorticity_linear_field += (mass_j / density_j) * (vi - vj).cross(gradW);
            );

            /*if (sim->is2DSimulation()) {
                vorticity_linear_field.x() = 0.0;
                vorticity_linear_field.y() = 0.0;
            }*/

            m_vorticity_equation[i] = m_vorticity_corrected_end[i] + dt * m_vorticity_derivative[i];
            m_vorticity_dissipation[i] = m_vorticity_equation[i] - vorticity_linear_field;

            if (t == 0) {
                m_vorticity_dissipation[i].setZero();
            }
        }
     }

     #pragma omp parallel default(shared)
     {
        #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            // 2nd loop: compute stream
            const Vector3r &xi = m_model->getPosition(i);
    
            Vector3r &stream = m_stream[i];
            stream.setZero();
    
            forall_fluid_neighbors_in_same_phase(
    
               const Vector3r xij = xi - xj;
               const Real density_j = m_model->getDensity(neighborIndex);
               const Real mass_j = m_model->getMass(neighborIndex);
    
               // compute stream function
               if (sim->is2DSimulation()) {
                   stream.z() -= (0.5 / M_PI) * (m_vorticity_dissipation[neighborIndex].z() * M_PI * r2) * log(xij.norm() + 0.01 * h2);
               }
               else {
                   stream += (0.25 / M_PI) * (m_vorticity_dissipation[neighborIndex] * mass_j / density_j) / (xij.norm() + 0.01 * h2);
               }
            ); 
        }  
     }

     #pragma omp parallel default(shared)
     {
         #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            // compute delta v and update v
            Vector3r &vi = m_model->getVelocity(i);
            const Vector3r &xi = m_model->getPosition(i);
            const Vector3r &stream_i = m_stream[i];
    
            Vector3r &delta_velocity = m_delta_velocity[i];
            delta_velocity.setZero();
    
            forall_fluid_neighbors_in_same_phase(
                const Real density_j = m_model->getDensity(neighborIndex);
                const Real mass_j = m_model->getMass(neighborIndex);
    
                const Vector3r xij = xi - xj;
                const Vector3r gradW = sim->gradW(xij);
                const Vector3r &stream_j = m_stream[neighborIndex];
    
                // compute delta velocity
                delta_velocity += (mass_j/ density_j) * (stream_i - stream_j).cross(gradW);
            );

        /*    if (sim->is2DSimulation()) {
                delta_velocity.z() = 0.0;
            }*/
    
            // refine linear velocity
            if ((vorticity_refinement_alpha * delta_velocity).norm() < vi.norm()){
                vi += vorticity_refinement_alpha * delta_velocity;
            }

            // saving final velocity
            m_velocity_corrected_end[i] = vi;
        }
     }

     #pragma omp parallel default(shared)
     {
        # pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {           
            //4th loop: compute final vorticity
            const Vector3r &vi = m_model->getVelocity(i);
            const Vector3r &xi = m_model->getPosition(i);

            // compute final corrected vorticity
            Vector3r& vorticity_corrected_end = m_vorticity_corrected_end[i];
            vorticity_corrected_end.setZero();

            forall_fluid_neighbors_in_same_phase(
                const Vector3r &vj = m_model->getVelocity(neighborIndex);
                const Real density_j = m_model->getDensity(neighborIndex);
                const Real mass_j = m_model->getMass(neighborIndex);

                const Vector3r xij = xi - xj;
                const Vector3r gradW = sim->gradW(xij);
                //vorticity through linear field
                vorticity_corrected_end += (mass_j / density_j) * (vi - vj).cross(gradW);
            );

            if (sim->is2DSimulation()) {
                vorticity_corrected_end.x() = 0.0;
                vorticity_corrected_end.y() = 0.0;
            }
        }

    }

    #pragma omp parallel default(shared)
     {
        #pragma omp for schedule(static)  
         for (int i = 0; i < (int)numParticles; i++)
         {
             //5th loop: compute vorticity derivative based on new vorticity

             const Vector3r& xi = m_model->getPosition(i);
             const Vector3r& vi = m_model->getVelocity(i);

             const Vector3r& vort_i = m_vorticity_corrected_end[i];
             Vector3r& vorticity_derivative = m_vorticity_derivative[i];
             vorticity_derivative.setZero();

             Vector3r& gradV_x = m_gradV_x[i];
             gradV_x.setZero();
             Vector3r& gradV_y = m_gradV_y[i];
             gradV_y.setZero();
             Vector3r& gradV_z = m_gradV_z[i];
             gradV_z.setZero();

             Vector3r& vorticity_rate_gradient = m_vorticity_rate_gradient[i];
             vorticity_rate_gradient.setZero();
             Vector3r& vorticity_rate_laplacian = m_vorticity_rate_laplacian[i];
             vorticity_rate_laplacian.setZero();

             forall_fluid_neighbors_in_same_phase(
                 const Vector3r & vj = m_model->getVelocity(neighborIndex);
                 const Real density_j = m_model->getDensity(neighborIndex);
                 const Real mass_j = m_model->getMass(neighborIndex);
                 const Vector3r & vort_j = m_vorticity_corrected_end[neighborIndex];

                 const Vector3r xij = xi - xj;
                 const Vector3r gradW = sim->gradW(xij);
                 const Vector3r vort_ij = vort_i - vort_j;

                 gradV_x += (mass_j / density_j) * (vj.x() - vi.x()) * (gradW);
                 gradV_y += (mass_j / density_j) * (vj.y() - vi.y()) * (gradW);
                 gradV_z += (mass_j / density_j) * (vj.z() - vi.z()) * (gradW);

                 if (sim->is2DSimulation()) {
                     vorticity_rate_laplacian.z() += 2 * (d + 2) * m_v_v * (mass_j / density_j) * vort_ij.z() * xij.dot(gradW) / (xij.squaredNorm() + 0.01 * h2);
                 }
                 else {
                     vorticity_rate_laplacian += 2 * (d + 2) * m_v_v * (mass_j / density_j) * ((vort_ij.dot(xij)) / (xij.squaredNorm() + 0.01 * h2)) * gradW;
                 }
            );

            vorticity_rate_gradient.x() = vort_i.dot(gradV_x);
            vorticity_rate_gradient.y() = vort_i.dot(gradV_y);
            vorticity_rate_gradient.z() = vort_i.dot(gradV_z);

            vorticity_derivative = vorticity_rate_gradient + vorticity_rate_laplacian;
        }
     }

 }


/*
#endif
*/
void VorticityRefinement::reset()
{
	for (unsigned int i = 0; i < m_model->numParticles(); i++){
        m_vorticity_linear_field[i].setZero();
        m_vorticity_init[i].setZero();
        m_vorticity_corrected_end[i].setZero();
        m_vorticity_derivative[i].setZero();
        m_vorticity_dissipation[i].setZero();
        m_stream[i].setZero();
        m_delta_velocity[i].setZero();
        m_vorticity_equation[i].setZero();
        m_v_adv[i].setZero();
        m_velocity_from_dfsph[i].setZero();
        m_velocity_corrected_end[i].setZero();
        m_position_from_dfsph[i].setZero();
        m_vorticity_rate_laplacian[i].setZero();
        m_vorticity_rate_gradient[i].setZero();
        m_gradV_x[i].setZero();
        m_gradV_y[i].setZero();
        m_gradV_z[i].setZero();
        m_direction_velocity[i] = 0.0;
    }
}


void SPH::VorticityRefinement::performNeighborhoodSearchSort()
{
    const unsigned int numPart = m_model->numActiveParticles();
    if (numPart == 0)
        return;

    Simulation* sim = Simulation::getCurrent();
    auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
    // d.sort_field(&m_vorticity_linear_field[0]);
    // d.sort_field(&m_vorticity_advected[0]);
    d.sort_field(&m_vorticity_corrected_end[0]);
    d.sort_field(&m_vorticity_derivative[0]);
    // d.sort_field(&m_vorticity_dissipation[0]);
    // d.sort_field(&m_stream[0]);
    // d.sort_field(&m_delta_velocity[0]);
    //d.sort_field(&m_vorticity_equation[0]);
    // d.sort_field(&m_v_adv[0]);
    // d.sort_field(&m_velocity_from_dfsph[0]);
     //d.sort_field(&m_velocity_corrected_end[0]);
    // d.sort_field(&m_position_from_dfsph[0]);
    // d.sort_field(&m_vorticity_rate_laplacian[0]);
    // d.sort_field(&m_vorticity_rate_gradient[0]);
    // d.sort_field(&m_gradV_x[0]);
    // d.sort_field(&m_gradV_y[0]);
    // d.sort_field(&m_gradV_z[0]);
    // d.sort_field(&m_direction_velocity[0]);
}