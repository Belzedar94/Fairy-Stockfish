# Forward all top-level make targets to the engine build system in src/
.PHONY: default

default:
	$(MAKE) -C src

%:
	$(MAKE) -C src $@
