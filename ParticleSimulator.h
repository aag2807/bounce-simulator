#ifndef PARTICLESIMULATOR_H
#define PARTICLESIMULATOR_H

#include "Particle.h"
#include <vector>
#include <SFML/Graphics.hpp>

class ParticleSimulator {
public:
    void addParticle(const Particle& particle);
    void update(double dt);
    void render(sf::RenderWindow& window) const;
    void createParticle(double x, double y);
    void clearParticles();

private:
    std::vector<Particle> particles;

    // SPH constants
    static constexpr double restDensity = 1000.0;
    static constexpr double gasConstant = 2000.0;
    static constexpr double viscosity = 250.0;
    static constexpr double kernelRadius = 10.0;
    static constexpr double mass = 1.0;

    void computeDensityAndPressure();
    void computeForces();
    double poly6Kernel(double rSquared) const;
    double spikyKernelGradient(double r) const;
    double viscosityKernelLaplacian(double r) const;
};

#endif // PARTICLESIMULATOR_H
