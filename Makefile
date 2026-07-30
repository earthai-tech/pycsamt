# pyCSAMT -- developer convenience Makefile
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
#
# Thin forwarders to pycsamt/models/_solver_build/*.sh for people working
# from a git checkout who prefer `make <target>` over spelling out that
# path. This is a convenience layer only -- it has no logic of its own and
# ships nothing extra; the scripts it calls are the single source of truth
# (see pycsamt/models/_solver_build/README.md and
# docs/source/user_guide/models/compilation.rst).
#
# Anyone using an installed (non-checkout) pycsamt should prefer the
# packaged CLI instead, which works the same way from any directory and
# needs no `make`:
#
#   pycsamt build modem3d --auto-install -y
#   pycsamt build occam2d --clean
#   pycsamt build mare2dem
#
# Requires a bash-compatible shell (see BASH below to override).
#
# Usage:
#   make modem2d ARGS="--auto-install -y"
#   make modem3d ARGS="--clean"
#   make occam2d
#   make mare2dem
#   make build-help                 # list solver targets
#   make BASH="/path/to/bash" modem3d

BASH        ?= bash
BUILD_DIR   := pycsamt/models/_solver_build
ARGS        ?=

.PHONY: modem2d modem3d occam2d mare2dem build-help

modem2d:
	$(BASH) $(BUILD_DIR)/modem2d.sh $(ARGS)

modem3d:
	$(BASH) $(BUILD_DIR)/modem3d.sh $(ARGS)

occam2d:
	$(BASH) $(BUILD_DIR)/occam2d.sh $(ARGS)

mare2dem:
	$(BASH) $(BUILD_DIR)/mare2dem.sh $(ARGS)

build-help:
	$(BASH) $(BUILD_DIR)/build.sh --help
