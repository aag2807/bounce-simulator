#include "Particle.h"
#include <cmath>

Particle::Particle(double x, double y)
        : x(x), y(y), vx(0), vy(0), fx(0), fy(0), density(0), pressure(0), shape(radius) {
    shape.setFillColor(sf::Color::Blue);
    shape.setPosition(x, y);
}

void Particle::update(double dt) {
    vx += fx * dt;
    vy += (fy + gravity * mass) * dt;
    x += vx * dt;
    y += vy * dt;

    // Collision with window boundaries
    if (x < radius) { x = radius; vx = -vx; }
    if (x > 800 - radius) { x = 800 - radius; vx = -vx; }
    if (y < radius) { y = radius; vy = -vy; }
    if (y > 600 - radius) { y = 600 - radius; vy = -vy; }

    shape.setPosition(x, y);
}

void Particle::draw(sf::RenderWindow& window) const {
    window.draw(shape);
}

void Particle::applyForce(double fx, double fy) {
    this->fx += fx;
    this->fy += fy;
}

void Particle::resetForce() {
    fx = 0;
    fy = 0;
}

double Particle::getX() const {
    return x;
}

double Particle::getY() const {
    return y;
}

double Particle::getVX() const {
    return vx;
}

double Particle::getVY() const {
    return vy;
}

double Particle::getDensity() const {
    return density;
}

double Particle::getPressure() const {
    return pressure;
}

void Particle::setPosition(double newX, double newY) {
    x = newX;
    y = newY;
    shape.setPosition(x, y);
}

void Particle::setVelocity(double newVX, double newVY) {
    vx = newVX;
    vy = newVY;
}

void Particle::setDensity(double density) {
    this->density = density;
}

void Particle::setPressure(double pressure) {
    this->pressure = pressure;
}
