#ifndef PARTICLE_H
#define PARTICLE_H

#include <SFML/Graphics.hpp>

class Particle {
public:
    Particle(double x, double y);
    void update(double dt);
    void draw(sf::RenderWindow& window) const;
    void applyForce(double fx, double fy);
    void resetForce();
    double getX() const;
    double getY() const;
    double getVX() const;
    double getVY() const;
    double getDensity() const;
    double getPressure() const;
    void setPosition(double newX, double newY);
    void setVelocity(double newVX, double newVY);
    void setDensity(double density);
    void setPressure(double pressure);

private:
    double x, y;  // Position
    double vx, vy;  // Velocity
    double fx, fy;  // Force
    double density, pressure;
    sf::CircleShape shape;
    static constexpr double radius = 5.0;
    static constexpr double mass = 1.0;
    static constexpr double gravity = 200.0;
};

#endif // PARTICLE_H
