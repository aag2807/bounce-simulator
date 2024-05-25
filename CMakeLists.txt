cmake_minimum_required(VERSION 3.28)
project(LiquidSimulator)

set(CMAKE_CXX_STANDARD 17)

# Find SFML package
find_package(SFML 2.5 COMPONENTS graphics window system REQUIRED)

# Add the source files
add_executable(LiquidSimulator main.cpp
        Particle.cpp
        Particle.h
        ParticleSimulator.cpp
        ParticleSimulator.h)

# Link SFML libraries
target_link_libraries(LiquidSimulator sfml-graphics sfml-window sfml-system)
