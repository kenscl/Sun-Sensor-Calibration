# Sun-Sensor Simulator
This repository contains the code for my Bachelors thesis

# General-Setup
```
git clone --recurse-submodules  https://github.com/kenscl/Sun-Sensor-Calibration.git sun-sensor
mkdir sun-sensor/sun-sensor-simulator/build
cd sun-sensor/sun-sensor-simulator/build
cmake ..
make -j 8
```

# Executing the code
There are two compilation targets:
- sun-sensor, execute with `./sun-sensor`
- alignment, execute with `./alignment`

# Error handling 
You might encounter issues with matplotlib-cpp.
If this occurs make sure gtk and python gtk are correctly installed.
On arch:
```
sudo pacman -S gtk4 python-gobject
```
