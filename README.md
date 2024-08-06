# Bachelorarbeit-Code
This is the Python-Code for my Bachelors Thesis.

# General-Setup
```
git clone --recurse-submodules  https://github.com/kenscl/Bachelorarbeit-Code.git sun-sensor
cd sun-sensor
python3 -m venv python/.venv 
source python/.venv/bin/activate
pip install -r python/requirements.txt
cmake -S sun-sensor-simulator/ -B sun-sensor-simulator/build
cmake --build sun-sensor-simulator/build
```

# Error handling 
You might encounter issues with matplotlib-cpp.
If this occurs make sure gtk and python gtk are correctly installed.
On arch:
```
sudo pacman -S gtk4 python-gobject
```
