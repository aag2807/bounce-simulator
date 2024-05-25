#include "ParticleSimulator.h"
#include <iostream>
#include <cmath>

void ParticleSimulator::addParticle(const Particle& particle) {
    particles.push_back(particle);
}

void ParticleSimulator::update(double dt) {
    computeDensityAndPressure();
    computeForces();

    for (auto& particle : particles) {
        particle.update(dt);
    }
}

void ParticleSimulator::render(sf::RenderWindow& window) const {
    for (const auto& particle : particles) {
        particle.draw(window);
    }
}

void ParticleSimulator::createParticle(double x, double y) {
    Particle newParticle(x, y);
    addParticle(newParticle);
}

void ParticleSimulator::clearParticles() {
    particles.clear();
}

void ParticleSimulator::computeDensityAndPressure() {
    for (auto& pi : particles) {
        double density = 0.0;
        for (auto& pj : particles) {
            double dx = pj.getX() - pi.getX();
            double dy = pj.getY() - pi.getY();
            double rSquared = dx * dx + dy * dy;
            if (rSquared < kernelRadius * kernelRadius) {
                density += mass * poly6Kernel(rSquared);
            }
        }
        pi.setDensity(density);
        pi.setPressure(gasConstant * (density - restDensity));
    }
}

void ParticleSimulator::computeForces() {
    for (auto& pi : particles) {
        double fx = 0.0;
        double fy = 0.0;

        for (auto& pj : particles) {
            if (&pi == &pj) continue;

            double dx = pj.getX() - pi.getX();
            double dy = pj.getY() - pi.getY();
            double r = std::sqrt(dx * dx + dy * dy);

            if (r < kernelRadius) {
                // Pressure force
                double pressureTerm = -0.5 * (pi.getPressure() + pj.getPressure()) * spikyKernelGradient(r);
                fx += pressureTerm * (dx / r);
                fy += pressureTerm * (dy / r);

                // Viscosity force
                double viscosityTerm = viscosity * (pj.getVX() - pi.getVX()) * viscosityKernelLaplacian(r);
                fx += viscosityTerm;
                fy += viscosityTerm * (dy / dx);
            }
        }

        // Apply computed forces
        pi.applyForce(fx, fy);
    }

    // Reset forces for next iteration
    for (auto& particle : particles) {
        particle.resetForce();
    }
}

double ParticleSimulator::poly6Kernel(double rSquared) const {
    static const double factor = 315.0 / (64.0 * M_PI * std::pow(kernelRadius, 9));
    return factor * std::pow(kernelRadius * kernelRadius - rSquared, 3);
}

double ParticleSimulator::spikyKernelGradient(double r) const {
    static const double factor = 15.0 / (M_PI * std::pow(kernelRadius, 6));
    return factor * std::pow(kernelRadius - r, 2);
}

double ParticleSimulator::viscosityKernelLaplacian(double r) const {
    static const double factor = 45.0 / (M_PI * std::pow(kernelRadius, 6));
    return factor * (kernelRadius - r);
}
