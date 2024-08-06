###########################
# Misc stuff for pretty-fing things
###########################
# Output using the 216-color rules mode
rule_file = $(notdir $(1))
rule_path = $(patsubst %/,%,$(dir $(1)))
last_path = $(notdir $(patsubst %/,%,$(dir $(1))))
ansicolor = $(shell echo $(call last_path,$(1)) | cksum | cut -b1-2 | xargs -IS expr 2 \* S + 17)
emacs_out = @printf "  %10s %s/%s\n" $(1) $(call rule_path,$(2)) $(call rule_file,$(2))
color_out = @if [ -t 1 ]; then \
				printf "  %10s \033[38;5;%d;1m%s\033[m/%s\n" \
					$(1) $(call ansicolor,$(2)) \
					$(call rule_path,$(2)) $(call rule_file,$(2)); else \
				printf "  %10s %s\n" $(1) $(2); fi
# if TERM=dumb, use it, otherwise switch to the term one
output = $(if $(TERM:dumb=),$(call color_out,$1,$2),$(call emacs_out,$1,$2))

# if V is set to non-nil, turn the verbose mode
quiet = $(if $(V),$($(1)),$(call output,$1,$@);$($(1)))

# make-4.3 allows string literals like "#include" in variables, but older versions need "\#include". Specifically, the following code:
#
# X := $(shell echo "#foo")
#
# works with make-4.3, but fails with previous versions:
#
# Makefile:1: *** unterminated call to function 'shell': missing ')'.  Stop.
#
# Older versions work if you spell it "\#foo", but 4.3 will include the backslash. We define $(HASH), which works consistently across versions.
HASH := \#
###########################
# END Misc stuff for pretty-fing things
###########################

PETSc.pc := $(PETSC_DIR)/$(PETSC_ARCH)/lib/pkgconfig/PETSc.pc

# ASAN must be left empty if you don't want to use it
ASAN ?=

CC = $(call pkgconf, --variable=ccompiler $(PETSc.pc))
CFLAGS = -std=c99 \
  $(call pkgconf, --variable=cflags_extra $(PETSc.pc)) \
  $(call pkgconf, --cflags-only-other $(PETSc.pc)) \
  $(OPT) $(OPT_EXAMPLES)
CPPFLAGS = $(call pkgconf, --cflags-only-I $(PETSc.pc)) \
  $(call pkgconf, --variable=cflags_dep $(PETSc.pc))
LDFLAGS = $(call pkgconf, --libs-only-L --libs-only-other $(PETSc.pc))
LDFLAGS += $(patsubst -L%, $(call pkgconf, --variable=ldflag_rpath $(PETSc.pc))%, $(call pkgconf, --libs-only-L $(PETSc.pc)))
LDLIBS = $(call pkgconf, --libs-only-l $(PETSc.pc)) -lm

AFLAGS ?= -fsanitize=address #-fsanitize=undefined -fno-omit-frame-pointer
CFLAGS += $(if $(ASAN),$(AFLAGS))
FFLAGS += $(if $(ASAN),$(AFLAGS))
LDFLAGS += $(if $(ASAN),$(AFLAGS))
CPPFLAGS += -I./include

OBJDIR := build

src.c := $(sort $(wildcard ./*.c))
src.o = $(src.c:%.c=$(OBJDIR)/%.o)

all: exec

exec: $(src.o) | $(PETSc.pc)
	$(call quiet,LINK.o) $^ $(LOADLIBES) $(LDLIBS) -o $@

.SECONDEXPANSION: # to expand $$(@D)/.DIR
%/.DIR :
	@mkdir -p $(@D)
	@touch $@

# Quiet, color output
quiet ?= $($(1))

$(OBJDIR)/%.o : %.c | $$(@D)/.DIR
	$(call quiet,CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $(abspath $<)

print: $(PETSc.pc)
	$(info CC      : $(CC))
	$(info CFLAGS  : $(CFLAGS))
	$(info CPPFLAGS: $(CPPFLAGS))
	$(info LDFLAGS : $(LDFLAGS))
	$(info LDLIBS  : $(LDLIBS))
	$(info OPT     : $(OPT))
	@true

print-% :
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [        flavor]: $(flavor $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info )
	@true

$(PETSc.pc):
	$(if $(wildcard $@),,$(error \
	  PETSc config not found. Please set PETSC_DIR and PETSC_ARCH))

.PHONY: all print clean

pkgconf = $(shell pkg-config $(if $(STATIC),--static) $1 | sed -e 's/^"//g' -e 's/"$$//g')

-include $(src.o:%.o=%.d)
