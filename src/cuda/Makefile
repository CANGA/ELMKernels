OBJECT        = ../


default: all

include $(OBJECT)config/Makefile.config


# Order is important due to .mod file creation
OBJS = CanopyHydrology_Interception.cu.o \
	CanopyHydrology_FracWet.cu.o \
	CanopyHydrology_FracH2OSfc.cu.o

all: $(OBJS)
	$(AR) cr libelm.a $(OBJS)

clean:
	@$(ELM_CLEAN); 
	@cd $(ELM_UTILS_DIR); $(ELM_CLEAN)
