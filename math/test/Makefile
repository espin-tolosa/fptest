# Check if the compiler supports FMA
HAS_FMA := $(shell echo | $(CC) -mfma -E - > /dev/null 2>&1 && echo 1 || echo 0)

ifeq ($(HAS_FMA),1)
    CFLAGS += -mfma
    FMA_ENABLED := 1
else
    FMA_ENABLED := 0
endif

all:
	@echo "FMA enabled: $(FMA_ENABLED)"
