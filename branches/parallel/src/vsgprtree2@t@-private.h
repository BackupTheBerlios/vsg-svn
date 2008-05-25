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

#ifndef __VSGPRTREE2@T@_PRIVATE_H__
#define __VSGPRTREE2@T@_PRIVATE_H__

#include "vsgprtree-parallel.h"
#include "vsgprtree2@t@.h"

G_BEGIN_DECLS;

/* binary localization macros: */
#define VSG_LOC2_COMP(loc,axis) ( \
 (vsgloc2) VSG_LOC2_MASK & ((loc) & (axis)) \
)

#define VSG_LOC2_ALONG_X(loc) ( VSG_LOC2_COMP (loc, VSG_LOC2_X) )
#define VSG_LOC2_ALONG_Y(loc) ( VSG_LOC2_COMP (loc, VSG_LOC2_Y) )

/* 
 * Convenience macros.
 */
#define CALL_POINT2@T@_LOC(config, one, other) ( \
(config)->point_loc_func ((one), (other), (config)->point_loc_data) \
)

#define CALL_REGION2@T@_LOC(config, one, other) ( \
(config)->region_loc_func ((one), (other), (config)->region_loc_data) \
)

#define CALL_POINT2@T@_DIST(config, one, other) ( \
(config)->point_dist_func ((one), (other), (config)->point_dist_data) \
)

/* forward typedefs */
typedef struct _VsgPRTree2@t@Leaf VsgPRTree2@t@Leaf;
typedef struct _VsgPRTree2@t@Int VsgPRTree2@t@Int;
typedef struct _VsgPRTree2@t@Node VsgPRTree2@t@Node;

typedef struct _VsgPRTree2@t@Config VsgPRTree2@t@Config;

typedef void (*VsgPRTree2@t@InternalFunc) (VsgPRTree2@t@Node *node,
                                           VsgPRTree2@t@NodeInfo *info,
                                           gpointer user_data);


/* private structs */
struct _VsgPRTree2@t@Leaf {

  gpointer isint;
  GSList *point;
};

struct _VsgPRTree2@t@Int {

  VsgPRTree2@t@Node *children[4];
};

struct _VsgPRTree2@t@Node {

  union {
    gpointer isint; /* first pointer of struct is 0 if leaf, not 0 otherwise */
    VsgPRTree2@t@Int interior;
    VsgPRTree2@t@Leaf leaf;
  } variable;

  /* bounding box limits */
  VsgVector2@t@ lbound;
  VsgVector2@t@ ubound;

  /* center of the box */
  VsgVector2@t@ center;

  /* counts */
  guint point_count;
  guint region_count;

  /* VsgRegion2 list */
  GSList *region_list;

  /* user data */
  gpointer user_data;

  /* parallel status */
  VsgParallelStatus parallel_status;
};

#define PRTREE2@T@NODE_ISINT(node) ( \
((node) != NULL) && ((node)->variable.isint != 0) \
)
#define PRTREE2@T@NODE_ISLEAF(node) ( \
((node) != NULL) && ((node)->variable.isint == 0) \
)

#define PRTREE2@T@NODE_INT(node) ( \
(node)->variable.interior \
)

#define PRTREE2@T@NODE_CHILD(node, i) ( \
PRTREE2@T@NODE_INT(node).children[i] \
)

#define PRTREE2@T@NODE_LEAF(node) ( \
(node)->variable.leaf \
)

#define PRTREE2@T@NODE_IS_REMOTE(node) ( \
VSG_PARALLEL_STATUS_IS_REMOTE (node->parallel_status) \
)

#define PRTREE2@T@NODE_IS_LOCAL(node) ( \
VSG_PARALLEL_STATUS_IS_LOCAL (node->parallel_status) \
)

#define PRTREE2@T@NODE_IS_SHARED(node) ( \
VSG_PARALLEL_STATUS_IS_SHARED (node->parallel_status) \
)

#define PRTREE2@T@NODE_PROC(node) ( \
VSG_PARALLEL_STATUS_PROC (node->parallel_status) \
)

struct _VsgPRTree2@t@Config {

  /* localization methods */
  VsgPoint2@t@LocDataFunc point_loc_func;
  gpointer point_loc_data;

  VsgRegion2@t@LocDataFunc region_loc_func;
  gpointer region_loc_data;

  /* point distance func */
  VsgPoint2@t@DistDataFunc point_dist_func;
  gpointer point_dist_data;

  /* config */
  guint max_point;

  /* spatial tolerance */
  @type@ tolerance;

  /* user node data */
  gpointer user_data_model;
  GType user_data_type;

  /* children order in traversals */
  VsgChildrenOrderDataFunc children_order;
  gpointer children_order_data;
  gpointer root_key;

  /* parallel tree configuration */
  VsgPRTreeParallelConfig parallel_config;
};

struct _VsgPRTree2@t@ {

  /* node part */
  VsgPRTree2@t@Node *node;

  /* tree configuration */
  VsgPRTree2@t@Config config;

  /* place to store pending message of inter processor VsgRegion */
  GSList *pending_shared_regions;
};

/* private functions */

void vsg_prtree2@t@node_free (VsgPRTree2@t@Node *node,
                              const VsgPRTree2@t@Config *config);

VsgPRTree2@t@Node *
vsg_prtree2@t@node_alloc (const VsgVector2@t@ *lbound,
                          const VsgVector2@t@ *ubound,
                          const VsgPRTree2@t@Config *config);

void _vsg_prtree2@t@node_get_info (VsgPRTree2@t@Node *node,
                                   VsgPRTree2@t@NodeInfo *node_info,
                                   VsgPRTree2@t@NodeInfo *father_info);

VsgPRTree2@t@Node *_vsg_prtree2@t@node_get_child_at (VsgPRTree2@t@Node *node,
                                                     const VsgVector2@t@ *pos,
                                                     gint depth);

void vsg_prtree2@t@_traverse_custom_internal (VsgPRTree2@t@ *prtree2@t@,
                                              GTraverseType order,
                                              VsgRegion2@t@LocDataFunc sel_func,
                                              VsgRegion2 selector,
                                              gpointer sel_data,
                                              VsgPRTree2@t@InternalFunc func,
                                              gpointer user_data);

void vsg_prtree2@t@node_make_int (VsgPRTree2@t@Node *node,
                                  const VsgPRTree2@t@Config *config);

G_END_DECLS;

#endif /* __VSGPRTREE2@T@_PRIVATE_H__ */
