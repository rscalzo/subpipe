#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <cir_memory.h>

static void send_err(size_t);

extern void *cir_malloc(size_t numbytes) {
    void *new;

    /* Malloc the space and see if it succeeded. */

    new = malloc(numbytes);
    if (new == NULL)
	send_err(numbytes);
    return(new);
}

extern void *cir_calloc(size_t nelm, size_t nb) {
    void *new;

    /* Calloc the space and see if it succeeded */

    new = cir_malloc(nelm*nb);
    bzero(new,nelm*nb);
    return(new);
}

extern void *cir_realloc(void *old, size_t numbytes) {
    void *new;

    /* Check that the old pointer had a value. If not, then just do a malloc */

    if (old == NULL)
	return(cir_malloc(numbytes));

    /* Do the realloc and check that it worked */

    new = realloc(old,numbytes);
    if (new == NULL)
	send_err(numbytes);
    return(new);
}

static void send_err(size_t numbytes) {
    fprintf(stderr,"Failed to allocate %d bytes of memory\n",numbytes);
    exit(2);
}
