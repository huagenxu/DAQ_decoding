
CC=g++
CXX=gcc
RM=rm -f
CFLAGS=-c -Wall

SRCS= madc_decoding.c
OBJS=$(subst .c, .o, $(SRCS))

ROOTCFLAGS	= $(shell root-config --cflags)
ROOTLIBS	= $(shell root-config --libs) -lNew
ROOTGLIBS       = $(shell root-config --glibs)

INCLUDESROOT    := 'root-config --incdir'
INCLUDES        := $PWD
LIBROOT         := 'root-config --libs'
LIBS            = $(ROOTLIBS) $(ROOTGLIBS)
CPPFLAGS       += $(ROOTCFLAGS)


all: decoding

#decoding: $(OBJS)
#	$(CXX) $(CPPFLAGS) $(LDFLAGS) -o decoding $(OBJS) $(LDLIBS)



decoding: madc_decoding.o
	$(CC) madc_decoding.o -o decoding -I$(INCLUDESROOT) $(LIBS)



madc_decoding.o: madc_decoding.c
	$(CC) $(CFLAGS) $(CPPFLAGS) madc_decoding.c -I$(INCLUDESROOT) $(LIBS)


clean:
	rm -rf *o decoding
