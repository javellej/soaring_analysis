# Soaring analyser

This project intends to gather open source tracking data from gliders around the world to gain knowledge about areas of lift in space an time.
The intention is to help glider pilots discover patterns and establish a strategy to select their trajectory.


## How to build

The build process is based on the `meson` build system:
   - `CC=gcc CXX=g++ meson build` -> this creates a directory called `build/`
   - `cd build/ && ninja`


## Dependencies

Tools:
   - `meson` and `ninja`
   - `locales`
   - `libeigen3-dev`
   - `libpng-dev`
   - `libpng++-dev`
   - `libboost-dev`

Code:
   - `cpp-lib` (https://gitlab.com/gewesp/cpp-lib) integrated as an external dependency
   - based on an example integration from http://wiki.glidernet.org/dev-cpp
