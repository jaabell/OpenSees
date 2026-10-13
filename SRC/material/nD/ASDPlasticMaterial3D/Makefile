include ../../../../Makefile.def

OBJS       = OPS_AllASDPlasticMaterial3Ds.o \
	ASD_material_definitions_0.o \
	ASD_material_definitions_1.o \
	ASD_material_definitions_2.o \
	ASD_material_definitions_3.o \
	ASD_material_definitions_4.o \
	ASD_material_definitions_5.o \
	ASD_material_definitions_6.o \
	ASD_material_definitions_7.o

all:    $(OBJS)

# Miscellaneous

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.



