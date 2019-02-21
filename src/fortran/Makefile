OBJECT        = ../


default: all

include $(OBJECT)/config/Makefile.config
FCFLAGS += -I$(ELM_UTILS_DIR)

# Order is important due to .mod file creation
OBJS = \
	$(ELM_UTILS_DIR)/shr_kind_mod.F90.o \
	$(ELM_UTILS_DIR)/shr_sys_mod.F90.o \
	$(ELM_UTILS_DIR)/shr_infnan_mod.F90.o \
	$(ELM_UTILS_DIR)/shr_strconvert_mod.F90.o \
	$(ELM_UTILS_DIR)/shr_assert_mod.F90.o \
	landunit_varcon.F90.o \
	column_varcon.F90.o \
	clm_varpar.F90.o \
	clm_varctl.F90.o \
	CanopyHydrology_Interception.F90.o \
	CanopyHydrology_SnowWater.F90.o \
	CanopyHydrology_FracWet.F90.o \
	CanopyHydrology_FracH2OSfc.F90.o

all: $(OBJS)
	$(AR) cr libelm.a $(OBJS)

clean:
	@$(ELM_CLEAN); 
	@cd $(ELM_UTILS_DIR); $(ELM_CLEAN)
