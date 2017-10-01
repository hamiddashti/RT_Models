CC = cc
INCLUDE = 
CFLAGS = ${INCLUDE} -g  -DVERBOSEjunk
LDFLAGS = -g
LIBS = -lm
PGM=gort1

SRCS = main.c ${PGM}.c crown_proj.c \
	S.c PD_s.c K_open.c Lk.c  misc.c  \
	EPgap.c  gort_brdf.c pgap.c gort_apply.c

OBJS = main.o  ${PGM}.o  crown_proj.o \
	S.o PD_s.o K_open.o Lk.o  misc.o \
	EPgap.o gort_brdf.o pgap.o gort_apply.o

${PGM}: ${OBJS}
	${CC} ${LDFLAGS} -o ${PGM} ${OBJS} ${LIBS}

main.o: main.c
	${CC} ${CFLAGS} -c main.c

crown_proj.o: crown_proj.c
	${CC} ${CFLAGS} -c crown_proj.c

S.o: S.c
	${CC} ${CFLAGS} -c S.c

PD_s.o: PD_s.c
	${CC} ${CFLAGS} -c PD_s.c

K_open.o: K_open.c
	${CC} ${CFLAGS} -c K_open.c

Lk.o: Lk.c
	${CC} ${CFLAGS} -c Lk.c
 
misc.o: misc.c
	${CC} ${CFLAGS} -c misc.c
   
EPgap.o: EPgap.c
	${CC} ${CFLAGS} -c EPgap.c
 
gort_brdf.o: gort_brdf.c
	${CC} ${CFLAGS} -c gort_brdf.c
 
pgap.o: pgap.c
	${CC} ${CFLAGS} -c pgap.c
 
gort_apply.o: gort_apply.c
	${CC} ${CFLAGS} -c gort_apply.c
 
${PGM}.o: ${PGM}.c
	${CC} ${CFLAGS} -c ${PGM}.c

clean:
	rm -f core ${OBJS} *.o

archive: ${SRCS} Makefile
	shar ${SRCS} Makefile > ${PGM}.shar


