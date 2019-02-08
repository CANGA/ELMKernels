SRCDIR	      =	.
OBJECT        = .


default: all

include $(OBJECT)/config/Makefile.config
FCFLAGS += -I $(ELM_UTILS_DIR) -fbounds-check -g3

# Order is important due to .mod file creation
OBJS = \
	$(ELM_UTILS_DIR)/shr_kind_mod.o \
	$(ELM_UTILS_DIR)/shr_sys_mod.o \
	$(ELM_UTILS_DIR)/shr_infnan_mod.o \
	$(ELM_UTILS_DIR)/shr_strconvert_mod.o \
	$(ELM_UTILS_DIR)/shr_assert_mod.o \
	landunit_varcon.o \
	column_varcon.o \
	CanopyHydrologyKern1.o

all: $(OBJS)
	$(AR) cr libelm.a $(OBJS)

clean:
	@$(ELM_CLEAN)

allclean: clean
	@cd $(ELM_UTILS_DIR); $(ELM_CLEAN)


