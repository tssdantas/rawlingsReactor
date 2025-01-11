# Detect the operating system
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
  # macOS specific commands
  OS := macos
  COMMAND = g++-13 -std=c++14 -I/opt/homebrew/Cellar/boost/1.86.0_2/include -L/opt/homebrew/Cellar/boost/1.86.0_2/lib -lboost_system -Wall -g etano_d.cpp -o etano_d

else ifeq ($(findstring MINGW,$(UNAME_S)),MINGW)
  # Windows specific commands (using MinGW environment)
  OS := windows
  COMMAND = g++-13 -std=c++14 -I/opt/homebrew/Cellar/boost/1.86.0_2/include -L/opt/homebrew/Cellar/boost/1.86.0_2/lib -lboost_system -Wall -g etano.cpp -o etano
else ifeq ($(UNAME_S),Linux)
  # Linux specific commands
  OS := linux
  #COMMAND = g++-10 -std=c++14 -I/.include -Wall -g etano.cpp -o etano
  #COMMAND = g++-10 -std=c++14 -I/mnt/e/apps/MATLAB/extern/include \
       -I/mnt/e/apps/MATLAB/simulink/include \
       -fPIC -shared -o benzeno.mexa64 benzeno_m.cpp \
       -L/mnt/e/apps/MATLAB/bin/glnxa64 -lmx -lmex -lmat

  COMMAND = g++-10 -std=c++17 benzeno_m.cpp -I/mnt/e/apps/MATLAB/extern/include -L/mnt/e/apps/MATLAB/bin/glnxa64 -lMatlabEngine -lMatlabDataArray -o benzeno_m


else
  $(error Unsupported OS)
endif

# Default target
all:
	$(COMMAND)
