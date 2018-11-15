#ifndef __MATALLOC_H__
extern void *matalloc(size_t size, void * data, int ndim, ...);
extern void matfree(void *p);
#define __MATALLOC_H__
#endif

