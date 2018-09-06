#ifndef __SORT_H__
#define __SORT_H__

#include "Constant.h"

#include <stddef.h>

void QuickSort(void *base, size_t nelem, size_t width,
               int (*cmp)(const void *, const void *));

void QuickSort2(void *base, void * base2, size_t nelem, size_t width,
               int (*cmp)(const void *, const void *));

void * BinarySearch(const void *key, const void *base,
               size_t nelem, size_t width,
               int (*cmp)(const void *, const void *));

#endif
