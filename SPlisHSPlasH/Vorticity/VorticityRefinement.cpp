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
    m_vorticity_next.resize(model->numParticles(), Vector3r::Zero());
    m_vorticity_derivative.resize(model->numParticles(), Vector3r::Zero());
    m_vorticity_dissipation.resize(model->numParticles(), Vector3r::Zero());
    m_stream.resize(model->numParticles(), Vector3r::Zero());
    m_delta_velocity.resize(model->numParticles(), Vector3r::Zero());
    m_vorticity_equation.resize(model->numParticles(), Vector3r::Zero());
    m_a_adv.resize(model->numParticles(), Vector3r(0.0, -9.81, 0.0));
    m_v_adv.resize(model->numParticles(), Vector3r::Zero());
    m_v_pre.resize(model->numParticles(), Vector3r::Zero());
    m_v_post.resize(model->numParticles(), Vector3r::Zero());
    m_x_pre.resize(model->numParticles(), Vector3r::Zero());
    m_vorticity_rate_laplacian.resize(model->numParticles(), Vector3r::Zero());
    m_vorticity_rate_gradient.resize(model->numParticles(), Vector3r::Zero());
    m_gradV_x.resize(model->numParticles(), Vector3r::Zero());
    m_gradV_y.resize(model->numParticles(), Vector3r::Zero());
    m_gradV_z.resize(model->numParticles(), Vector3r::Zero());

    m_vorticity_refinement_alpha = static_cast<Real>(1.0);
    
    model->addField({ "vorticity_current", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_current[i][0]; }, true });
    model->addField({ "vorticity_next", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_next[i][0]; }, true });
    model->addField({ "stream", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_stream[i][0]; }, true });
    model->addField({ "vorticity_linear_field", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_linear_field[i][0]; }, true });
    model->addField({ "vorticity_derivative", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_derivative[i][0]; }, true });
    model->addField({ "vorticity_dissipation", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_dissipation[i][0]; }, true });
    model->addField({ "delta_velocity", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_delta_velocity[i][0]; }, true });
    model->addField({ "m_vorticity_equation", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_equation[i][0]; }, true });
    model->addField({ "m_a_adv", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_a_adv[i][0]; }, true });
    model->addField({ "m_v_pre", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_v_pre[i][0]; }, true });
    model->addField({ "m_v_post", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_v_post[i][0]; }, true });
    model->addField({ "m_x_pre", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_x_pre[i][0]; }, true });
    model->addField({ "m_vorticity_rate_laplacian", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_rate_laplacian[i][0]; }, true });
    model->addField({ "m_vorticity_rate_gradient", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vorticity_rate_gradient[i][0]; }, true });
    model->addField({ "m_gradV_x", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_gradV_x[i][0]; }, true });
    model->addField({ "m_gradV_y", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_gradV_y[i][0]; }, true });
    model->addField({ "m_gradV_z", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_gradV_z[i][0]; }, true });
}

VorticityRefinement::~VorticityRefinement(void)
{
    m_model->removeFieldByName("vorticity_current");
    m_model->removeFieldByName("vorticity_next");
    m_model->removeFieldByName("stream");
    m_model->removeFieldByName("vorticity_linear_field");
    m_model->removeFieldByName("vorticity_derivative");
    m_model->removeFieldByName("vorticity_dissipation");
    m_model->removeFieldByName("delta_velocity");
    m_model->removeFieldByName("m_vorticity_equation");
    m_model->removeFieldByName("m_a_adv");
    m_model->removeFieldByName("m_v_adv");
    m_model->removeFieldByName("m_v_pre");
    m_model->removeFieldByName("m_v_post");
    m_model->removeFieldByName("m_x_pre");
    m_model->removeFieldByName("m_vorticity_rate_laplacian");
    m_model->removeFieldByName("m_vorticity_rate_gradient");
    m_model->removeFieldByName("m_gradV_x");
    m_model->removeFieldByName("m_gradV_y");
    m_model->removeFieldByName("m_gradV_z");


    m_vorticity_linear_field.clear();
    m_vorticity_current.clear();
    m_vorticity_next.clear();
    m_vorticity_derivative.clear();
    m_vorticity_dissipation.clear();
    m_stream.clear();
    m_delta_velocity.clear();
    m_vorticity_equation.clear();
    m_a_adv.clear();
    m_v_adv.clear();
    m_v_pre.clear();
    m_v_post.clear();
    m_x_pre.clear();
    m_vorticity_rate_laplacian.clear();
    m_vorticity_rate_gradient.clear();
    m_gradV_x.clear();
    m_gradV_y.clear();
    m_gradV_z.clear();
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

            m_v_pre[i] = m_model->getVelocity(i);
            //advected velocity
            m_v_adv[i] = m_v_post[i] + dt*m_a_adv[i];
            // m_vorticity_current[i] = m_vorticity_next[i];
        }

        #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {

            //4th loop: compute linear field vorticity
            Vector3r &xi = m_model->getPosition(i);
            Vector3r &vi = m_model->getVelocity(i);

            Vector3r& xi_pre = m_x_pre[i];
            
            Vector3r& vorticity_linear_field = m_vorticity_linear_field[i];
            vorticity_linear_field.setZero();

            Vector3r& vorticity_current = m_vorticity_current[i];
            vorticity_current.setZero();
            
            forall_fluid_neighbors_in_same_phase(
                Vector3r &vj = m_model->getVelocity(neighborIndex);
                Real density_j = m_model->getDensity(neighborIndex);
                Real mass_j = m_model->getMass(neighborIndex);
                
                Vector3r xij = xi - xj;
                Vector3r gradW = sim->gradW(xij);
                //vorticity through linear field
                vorticity_linear_field += (mass_j / density_j) * (vi - vj).cross(gradW);

                Vector3r xij_pre = xi_pre - m_x_pre[neighborIndex];
                Vector3r gradW_pre = sim->gradW(xij_pre);
                vorticity_current += (mass_j / density_j) * (m_v_adv[i] - m_v_adv[j]).cross(gradW_pre);
            );
        }


        #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            //5th loop: compute vorticity derivative based on new vorticity and next velocity
    
            Vector3r& xi = m_x_pre[i];
            Vector3r& vi = m_v_adv[i];
            Vector3r &vorticity_current = m_vorticity_current[i];
            Vector3r& vorticity_linear_field = m_vorticity_linear_field[i];
            Vector3r &vorticity_derivative = m_vorticity_derivative[i];
            vorticity_derivative.setZero();
    
            Vector3r &gradV_x = m_gradV_x[i];
            gradV_x.setZero();
            Vector3r &gradV_y = m_gradV_y[i];
            gradV_y.setZero();
            Vector3r &gradV_z = m_gradV_z[i];
            gradV_z.setZero();
    
            Vector3r &vorticity_rate_gradient = m_vorticity_rate_gradient[i];
            vorticity_rate_gradient.setZero();
            Vector3r &vorticity_rate_laplacian = m_vorticity_rate_laplacian[i];
            vorticity_rate_laplacian.setZero();
    
            forall_fluid_neighbors_in_same_phase(
                Vector3r & vj = m_v_adv[neighborIndex];
                Real density_j = m_model->getDensity(neighborIndex);
                Real mass_j = m_model->getMass(neighborIndex);
    
                Vector3r xij = xi - m_x_pre[neighborIndex];
                Vector3r gradW = sim->gradW(xij);
                Vector3r vort_ij = vorticity_current - m_vorticity_current[neighborIndex];
                    
                //  vorticity * gradV + v_v * laplacian(vorticity)
                gradV_x += (mass_j / density_j) * (vj.x() - vi.x()) * (gradW);
                gradV_y += (mass_j / density_j) * (vj.y() - vi.y()) * (gradW);
                gradV_z += (mass_j / density_j) * (vj.z() - vi.z()) * (gradW);
    
                vorticity_rate_laplacian += 2 * (d + 2) * m_v_v * (mass_j / density_j) * ((vort_ij.dot(xij)) / (xij.squaredNorm() + 0.01 * h2)) * gradW;
            );
    
             //vorticity_rate_gradient = vorticity_current.x() * gradV_x + vorticity_current.y() * gradV_y + vorticity_current.z() * gradV_z;

            vorticity_rate_gradient.x() = vorticity_current.dot(gradV_x);
            vorticity_rate_gradient.y() = vorticity_current.dot(gradV_y);
            vorticity_rate_gradient.z() = vorticity_current.dot(gradV_z);
    
            vorticity_derivative = vorticity_rate_gradient + vorticity_rate_laplacian;
    
            Vector3r& vorticity_equation = m_vorticity_equation[i]; //compute vorticity equation
            vorticity_equation = m_vorticity_next[i] + dt * vorticity_derivative;

            Vector3r& vorticity_dissipation = m_vorticity_dissipation[i]; //compute vorticity equation
            vorticity_dissipation = vorticity_equation - vorticity_linear_field;
        }

        #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            // 2nd loop: compute stream
            Vector3r &xi = m_model->getPosition(i);
    
            Vector3r &stream = m_stream[i];
            stream.setZero();
    
            forall_fluid_neighbors_in_same_phase(
    
                Vector3r xij = xi - xj;
                Real density_j = m_model->getDensity(neighborIndex);
                Real mass_j = m_model->getMass(neighborIndex);
    
                // compute stream function
                stream += (0.25 / M_PI) * (m_vorticity_dissipation[neighborIndex] * mass_j / density_j) / xij.norm();
            ); 
        }  

        #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            // 3rd loop: compute delta v and update v
            Vector3r &vi = m_model->getVelocity(i);
            Vector3r &xi = m_model->getPosition(i);
            Vector3r &stream_i = m_stream[i];
    
            Vector3r &delta_velocity = m_delta_velocity[i];
            delta_velocity.setZero();
    
            forall_fluid_neighbors_in_same_phase(
                Real density_j = m_model->getDensity(neighborIndex);
                Real mass_j = m_model->getMass(neighborIndex);
    
                Vector3r xij = xi - xj;
                Vector3r gradW = sim->gradW(xij);
    
                // compute delta velocity
                delta_velocity += (mass_j/ density_j) * (stream_i - m_stream[neighborIndex]).cross(gradW);
            );
    
            // refine linear velocity
            vi += vorticity_refinement_alpha * delta_velocity;
        }


        #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            
            //updating stuff for next loop

            m_a_adv[i] = m_model->getAcceleration(i);
            m_v_post[i] = m_model->getVelocity(i);
            m_x_pre[i] = m_model->getPosition(i);


            //4th loop: compute linear field vorticity
            Vector3r& xi = m_model->getPosition(i);
            Vector3r& vi = m_model->getVelocity(i);


            Vector3r& vorticity_next = m_vorticity_next[i];
            vorticity_next.setZero();

            forall_fluid_neighbors_in_same_phase(
                Vector3r& vj = m_model->getVelocity(neighborIndex);
                Real density_j = m_model->getDensity(neighborIndex);
                Real mass_j = m_model->getMass(neighborIndex);

                Vector3r xij = xi - xj;
                Vector3r gradW = sim->gradW(xij);
                //vorticity through linear field
                vorticity_next += (mass_j / density_j) * (vi - vj).cross(gradW);
            );
        }



    }
}


void VorticityRefinement::reset()
{
	for (unsigned int i = 0; i < m_model->numParticles(); i++){
        m_vorticity_linear_field[i].setZero();
        m_vorticity_current[i].setZero();
        m_vorticity_next[i].setZero();
        m_vorticity_derivative[i].setZero();
        m_vorticity_dissipation[i].setZero();
        m_stream[i].setZero();
        m_delta_velocity[i].setZero();
        m_vorticity_equation[i].setZero();
        m_a_adv[i] = Vector3r(0.0, -9.81, 0.0);
        m_v_adv[i].setZero();
        m_v_pre[i].setZero();
        m_v_post[i].setZero();
        m_x_pre[i].setZero();
        m_vorticity_rate_laplacian[i].setZero();
        m_vorticity_rate_gradient[i].setZero();
        m_gradV_x[i].setZero();
        m_gradV_y[i].setZero();
        m_gradV_z[i].setZero();
    }
}


void SPH::VorticityRefinement::performNeighborhoodSearchSort()
{
    const unsigned int numPart = m_model->numActiveParticles();
    if (numPart == 0)
        return;

    Simulation* sim = Simulation::getCurrent();
    auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
    d.sort_field(&m_vorticity_linear_field[0]);
    d.sort_field(&m_vorticity_current[0]);
    d.sort_field(&m_vorticity_next[0]);
    d.sort_field(&m_vorticity_derivative[0]);
    d.sort_field(&m_vorticity_dissipation[0]);
    d.sort_field(&m_stream[0]);
    d.sort_field(&m_delta_velocity[0]);
    d.sort_field(&m_vorticity_equation[0]);
    d.sort_field(&m_a_adv[0]);
    d.sort_field(&m_v_adv[0]);
    d.sort_field(&m_v_pre[0]);
    d.sort_field(&m_v_post[0]);
    d.sort_field(&m_x_pre[0]);
    d.sort_field(&m_vorticity_rate_laplacian[0]);
    d.sort_field(&m_vorticity_rate_gradient[0]);
    d.sort_field(&m_gradV_x[0]);
    d.sort_field(&m_gradV_y[0]);
    d.sort_field(&m_gradV_z[0]);
}