#include <stdio.h>
#include <stdlib.h>
#include "vertices.h"

/* local functions */
vertices *new_vert();
void free_vert();

/*------------------------------------------------------------------------------
  Allocate memory for vertices structure of nvmax vertices.

   Input: nvmax = desired number of vertices.
  Return value: pointer to a new vertices structure,
		or null if failed to allocate memory.
*/ 
vertices *new_vert(nvmax)
int nvmax;
{
    vertices *vert = 0x0;

    /* allocate memory for new vertices structure */
    vert = (vertices *) malloc(sizeof(vertices));
    if (!vert) return(0x0);

    /* allocate memory for array of nvmax az-el structures */
    vert->v = (azel *) malloc(sizeof(azel) * nvmax);
    if (!vert->v) return(0x0);

    vert->nvmax = nvmax;

    return(vert);
}

/*------------------------------------------------------------------------------
  Free vertices memory.
*/ 
void free_vert(vert)
vertices *vert;
{
    if (vert) {
	if (vert->v) {
	    free(vert->v);
	    vert->v = 0x0;
	}
	free(vert);
    }
}
