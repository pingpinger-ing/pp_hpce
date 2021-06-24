CPPFLAGS += -std=c++11 -W -Wall -g -O3 -I include
CPPFLAGS += -Wno-unused-parameter

SHELL	= /bin/bash

LDLIBS += -ljpeg
LDLIBS += -ltbb

all: reference_tools user_simulator

bin/% : src/%.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

bin/user/simulator: include/user_simulator.hpp include/graphs/user_heat.hpp

reference_tools : bin/ref/simulator bin/tools/generate_heat_rect

user_simulator : bin/user/simulator

bin w records results:
	mkdir -p $@

.SECONDARY:
.DELETE_ON_ERROR:

UPD	?= 1024
SUB	?= 4
SPA	?= 8
LOG	?= 0

w/heat%_rect.graph: | bin/tools/generate_heat_rect w
	bin/tools/generate_heat_rect $* $(UPD) $(SUB) $(SPA) > $@

w/heat_%.graph: | w
	octave src/tools/generate_heat_$*.m > $@

w/%.ref.mjpeg: w/%.graph | bin/ref/simulator
	(time $| --log-level $(LOG) $< $(@:.mjpeg=.stats) $@) 2>&1 | tee $(@:.mjpeg=.log)

w/%.user.mjpeg: w/%.graph bin/user/simulator
	(time bin/user/simulator --log-level $(LOG) $< $(@:.mjpeg=.stats) $@) 2>&1 | tee $(@:.mjpeg=_$(LOG).log)

results/%.pass: w/%.ref.mjpeg w/%.user.mjpeg | results
	-diff -q $^
	-diff -q w/$*.ref.stats w/$*.user.stats

test:	results/heat128_rect.pass \
