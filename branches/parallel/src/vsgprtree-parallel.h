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
#include <vsg/vsgpackedmsg.h>
#include <vsg/vsgprtree-common.h>

G_BEGIN_DECLS;


/* typedefs */

typedef struct _VsgPRTreeParallelConfig VsgPRTreeParallelConfig;
typedef struct _VsgParallelVTable VsgParallelVTable;

typedef gpointer (*VsgMigrableAllocDataFunc) (gboolean resident,
                                              gpointer user_data);

typedef void (*VsgMigrableDestroyDataFunc) (gpointer data, gboolean resident,
                                            gpointer user_data);

typedef void (*VsgMigrablePackDataFunc) (gpointer var, VsgPackedMsg *pm,
                                         gpointer user_data);

typedef void (*VsgMigrableUnpackDataFunc) (gpointer var, VsgPackedMsg *pm,
                                           gpointer user_data);


/* structs */
struct _VsgParallelVTable {

  /* memory management functions */
  VsgMigrableAllocDataFunc alloc;
  gpointer alloc_data;

  VsgMigrableDestroyDataFunc destroy;
  gpointer destroy_data;

  /* migration functions */
  VsgMigrablePackDataFunc migrate_pack;
  gpointer migrate_pack_data;

  VsgMigrableUnpackDataFunc migrate_unpack;
  gpointer migrate_unpack_data;

  /* functions for volatile migrations (visits) */
  VsgMigrablePackDataFunc visit_forward_pack;
  gpointer visit_forward_pack_data;

  VsgMigrableUnpackDataFunc visit_forward_unpack;
  gpointer visit_forward_unpack_data;

  VsgMigrablePackDataFunc visit_backward_pack;
  gpointer visit_backward_pack_data;

  VsgMigrableUnpackDataFunc visit_backward_unpack;
  gpointer visit_backward_unpack_data;

};

struct _VsgPRTreeParallelConfig {
  MPI_Comm communicator;

  VsgParallelVTable point;
  VsgParallelVTable region;
  VsgParallelVTable node_data;
};

G_END_DECLS;

#endif