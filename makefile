# Detect the operating system
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
  # macOS specific commands
  OS := macos
  COMMAND = g++-13 -std=c++14 -I/opt/homebrew/Cellar/boost/1.86.0_2/include -L/opt/homebrew/Cellar/boost/1.86.0_2/lib -lboost_system -Wall -g etano_d.cpp -o etano_d

else ifeq ($(findstring MINGW,$(UNAME_S)),MINGW)
  # Windows specific commands (using MinGW environment)
  OS := windows
  COMMAND = g++ -std=c++14 -g etano.cpp -o etano
else ifeq ($(UNAME_S),Linux)
  # Linux specific commands
  OS := linux
  COMMAND = g++-10 -std=c++14 -I/usr/include -L/usr/lib -lboost_system -Wall -g etano_d.cpp -o etano


else
  $(error Unsupported OS)
endif

# Default target
all:
	$(COMMAND)
