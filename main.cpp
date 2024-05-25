#include <SFML/Graphics.hpp>
#include "Particle.h"
#include "ParticleSimulator.h"

int main() {
    sf::RenderWindow window(sf::VideoMode(800, 600), "Bounce Simulator");

    ParticleSimulator simulator;

    sf::Clock clock;

    while (window.isOpen()) {
        sf::Event event{};
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            } else if (event.type == sf::Event::MouseButtonPressed) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    simulator.createParticle(event.mouseButton.x, event.mouseButton.y);
                }
            } else if (event.type == sf::Event::MouseMoved) {
                if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
                    simulator.createParticle(event.mouseMove.x, event.mouseMove.y);
                }
            } else if (event.type == sf::Event::KeyPressed) {
                auto keycode = event.key.code;
                switch (keycode) {

                    // pressing space bar
                    case 39: // space mac
                    case 32: // space windows
                    case 57: //space linux
                        simulator.clearParticles();
                        break;
                    default:
                        break;
                }
            }
        }

        double dt = clock.restart().asSeconds();

        simulator.update(dt);

        window.clear();
        simulator.render(window);
        window.display();
    }

    return 0;
}
