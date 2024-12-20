# Detect the operating system
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
  # macOS specific commands
  OS := macos
  COMMAND = echo "This is macOS"
else ifeq ($(findstring MINGW,$(UNAME_S)),MINGW)
  # Windows specific commands (using MinGW environment)
  OS := windows
  COMMAND = g++-13 -std=c++14 -I/opt/homebrew/Cellar/boost/1.86.0_2/include -L/opt/homebrew/Cellar/boost/1.86.0_2/lib -lboost_system -Wall -g etano.cpp -o etano
else ifeq ($(UNAME_S),Linux)
  # Linux specific commands
  OS := linux
  #COMMAND = g++-10 -std=c++14 -I/.include -Wall -g etano.cpp -o etano
  COMMAND = g++-10 -std=c++14 -I/.include -Wall -g benzeno.cpp -o benzeno
else
  $(error Unsupported OS)
endif

# Default target
all:
	$(COMMAND)
