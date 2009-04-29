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

#ifdef VSG_HAVE_MPI
typedef void (*VsgMigrablePackDataFunc) (gpointer var, VsgPackedMsg *pm,
                                         gpointer user_data);

typedef void (*VsgMigrableUnpackDataFunc) (gpointer var, VsgPackedMsg *pm,
                                           gpointer user_data);
/* structs */
struct _VsgParallelMigrateVTable {
  VsgMigrablePackDataFunc pack;
  gpointer pack_data;

  VsgMigrableUnpackDataFunc unpack;
  gpointer unpack_data;
};

#endif

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

struct _VsgPRTreeParallelConfig {
#ifdef VSG_HAVE_MPI
  MPI_Comm communicator;
#else
  gpointer padding;
#endif

  VsgParallelVTable point;
  VsgParallelVTable region;
  VsgParallelVTable node_data;
};

G_END_DECLS;

#endif
