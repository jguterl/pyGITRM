# ******************************************************************
#  This document is confidential information.
#  Copyright 1997-2023 Simmetrix Inc. All rights reserved. This 
#  document is an unpublished work fully protected by the United 
#  States copyright laws and is considered a trade secret belonging 
#  to the copyright holder. Disclosure, use, or reproduction without 
#  the written authorization of Simmetrix Inc. is prohibited. 
# ******************************************************************/ 
FNAME = simTranslate
PROJECTROOT = /home/nathd/simmetrix
PLATFORM = x64_rhel8_gcc83
PVERS := 360
BUILD := 2023.1-231028
EXECNAME := $(FNAME)

  INCLUDES := -I$(PROJECTROOT)/$(BUILD)/include
  LIBPATH := \
    -L$(PROJECTROOT)/$(BUILD)/lib/$(PLATFORM) \
    -L$(PROJECTROOT)/$(BUILD)/lib/$(PLATFORM)/psKrnl
  LIBS := \
    $(LIBPATH) \
    -lSimMeshing \
    -lSimMeshTools \
    -lSimParasolid$(PVERS) \
    -lSimModel \
    -lSimLicense \
    -lpskernel 

  COMPFLAG := -c
  OBJNAME := $(FNAME).o
  OBJTARGET := -o $(OBJNAME)
  LINKTARGET := -o $(EXECNAME)

  ifneq ($(findstring x64_rhel,$(PLATFORM)),)
    CC = g++
    LINK = g++
    CCFLAGS := -O2 -std=c++11
    LINKFLAGS :=
    SOLIBPATH := -Xlinker -rpath -Xlinker \
      $(PROJECTROOT)/$(BUILD)/lib/$(PLATFORM)/psKrnl
    LIBS := $(LIBS) $(SOLIBPATH) \
      -lpthread \
      -lm 
  endif

  ifneq ($(findstring x64_rhel8,$(PLATFORM)),)
    LIBS := $(LIBS) \
      -ltirpc 
  endif



$(EXECNAME): $(OBJNAME)
	-mkdir -p bin/$(PLATFORM)
	$(LINK) $(LINKFLAGS) $(OBJNAME) $(LINKTARGET) $(LIBS)

$(OBJNAME): $(FNAME).cc
	-mkdir -p obj/$(PLATFORM)
	$(CC) $(CCFLAGS) $(INCLUDES) $(COMPFLAG) $(FNAME).cc $(OBJTARGET)

clean:
	rm -f $(EXECNAME) $(OBJNAME)

