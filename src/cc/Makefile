OBJECT        = ../


default: all

include $(OBJECT)config/Makefile.config


# Order is important due to .mod file creation
OBJS =CanopyHydrology_Interception.cc.o \
	CanopyHydrology_FracWet.cc.o \
	CanopyHydrology_FracH2OSfc.cc.o

all: $(OBJS)
	$(AR) cr libelm.a $(OBJS)

clean:
	@$(ELM_CLEAN); 
