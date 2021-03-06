CPPFLAGS += -std=c++11 -W -Wall -g -O3 -I include
CPPFLAGS += -Wno-unused-parameter

SHELL	= /bin/bash

LDLIBS += -ljpeg
LDLIBS += -ltbb

all: reference_tools user_simulator

bin/% : src/%.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

bin/user/user_simulator: include/user_simulator.hpp include/graphs/user_heat.hpp

reference_tools : bin/ref/simulator bin/tools/generate_heat_rect

user_simulator : bin/user/user_simulator
