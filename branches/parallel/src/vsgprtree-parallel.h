/* LIBVSG - Visaurin Geometric Library
 * Copyright (C) 2006-2007 Pierre Gay
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __VSGPRTREE_PARALLEL_H__
#define __VSGPRTREE_PARALLEL_H__

#include <vsg/vsgmpi.h>
#ifdef VSG_HAVE_MPI
#include <vsg/vsgpackedmsg.h>
#endif

#include <vsg/vsgprtree-common.h>

G_BEGIN_DECLS;


/* typedefs */

typedef struct _VsgPRTreeParallelConfig VsgPRTreeParallelConfig;
typedef struct _VsgParallelMigrateVTable VsgParallelMigrateVTable;
typedef struct _VsgParallelVTable VsgParallelVTable;

typedef gpointer (*VsgMigrableAllocDataFunc) (gboolean resident,
                                              gpointer user_data);

typedef void (*VsgMigrableDestroyDataFunc) (gpointer data, gboolean resident,
                                            gpointer user_data);

typedef void (*VsgMigrableReductionDataFunc) (gpointer a, gpointer b,
                                              gpointer user_data);

#ifdef VSG_HAVE_MPI
typedef void (*VsgMigrablePackDataFunc) (gpointer var, VsgPackedMsg *pm,
                                         gpointer user_data);

typedef void (*VsgMigrableUnpackDataFunc) (gpointer var, VsgPackedMsg *pm,
                                           gpointer user_data);

/* structs */

/**
 * VsgParallelMigrateVTable:
 * @pack: a packing function
 * @pack_data: user provided arbitrary data to be used as a @pack argument
 * @unpack: a unpacking function
 * @unpack_data: user provided arbitrary data to be used as an @unpack argument
 * @reduce: a reduction function
 * @reduce_data: user provided arbitrary data to be used as a @reduce argument
 *
 * a set of functions to manipulate a datatype through
 * pack, unpack and reduce operations.
 */
struct _VsgParallelMigrateVTable {
  VsgMigrablePackDataFunc pack;
  gpointer pack_data;

  VsgMigrableUnpackDataFunc unpack;
  gpointer unpack_data;

  VsgMigrableReductionDataFunc reduce;
  gpointer reduce_data;
};

#endif

/**
 * VsgParallelVTable:
 * @alloc: an allocation function
 * @alloc_data: user provided arbitrary data to be used as a @alloc argument
 * @destroy: a deallocation function
 * @destroy_data: user provided arbitrary data to be used as a @destroy argument
 * @migrate: a #VsgParallelMigrateVTable. WARNING: Only available when
 * VSG_HAVE_MPI is defined.
 * @visit_forward: a #VsgParallelMigrateVTable. WARNING: Only available when
 * VSG_HAVE_MPI is defined.
 * @visit_backward: a #VsgParallelMigrateVTable. WARNING: Only available when
 * VSG_HAVE_MPI is defined.
 *
 * a set of functions for manipulate some data type with allocate,
 * deallocate and migration operations.
 */
struct _VsgParallelVTable {

  /* memory management functions */
  VsgMigrableAllocDataFunc alloc;
  gpointer alloc_data;

  VsgMigrableDestroyDataFunc destroy;
  gpointer destroy_data;

#ifdef VSG_HAVE_MPI
  /* migration functions */
  VsgParallelMigrateVTable migrate;

  /* functions for volatile migrations (visits) */
  VsgParallelMigrateVTable visit_forward;
  VsgParallelMigrateVTable visit_backward;
#endif
};

/**
 * VsgPRTreeParallelConfig:
 * @point: a #VsgParallelVTable for #Vsgpoint
 * @region:  a #VsgParallelVTable for #VsgRegion
 * @node_data:  a #VsgParallelVTable for node's user data.
 * @communicator:  a #MPI_Comm. WARNING: Only available when
 * VSG_HAVE_MPI is defined.
 *
 * a configuration for parallel operations on #VsgPRTree**. 
 */
struct _VsgPRTreeParallelConfig {
  VsgParallelVTable point;
  VsgParallelVTable region;
  VsgParallelVTable node_data;

#ifdef VSG_HAVE_MPI
  MPI_Comm communicator;
#endif
};

G_END_DECLS;

#endif
