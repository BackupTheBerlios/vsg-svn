#include "vsgprtree2@t@-parallel.h"
#include "vsgprtree2@t@-private.h"

#include "vsgcommbuffer.h"
#include "vsgcommbuffer-private.h"

#include "string.h"

typedef struct _VTableAndMsg VTableAndMsg;

struct _VTableAndMsg
{
  VsgParallelVTable *vtable;
  VsgPackedMsg *msg;
};

static void _send_point (VsgPoint2 pt, VTableAndMsg *vm)
{
  vm->vtable->migrate.pack (pt, vm->msg, vm->vtable->migrate.pack_data);
  vm->vtable->destroy (pt, TRUE, vm->vtable->destroy_data);
}

static void _send_region (VsgRegion2 rg, VTableAndMsg *vm)
{
  vm->vtable->migrate.pack (rg, vm->msg, vm->vtable->migrate.pack_data);
  vm->vtable->destroy (rg, TRUE, vm->vtable->destroy_data);
}

static void _send_shared_region (VsgRegion2 rg, VTableAndMsg *vm)
{
  vm->vtable->migrate.pack (rg, vm->msg, vm->vtable->migrate.pack_data);
}

void vsg_prtree2@t@_set_parallel (VsgPRTree2@t@ *tree,
                                  VsgPRTreeParallelConfig *pconfig)
{
  MPI_Comm comm;
  VsgParallelVTable *pt_vtable;
  VsgParallelVTable *rg_vtable;
  
  VsgVector2@t@ bounds[2];
  gint rk, sz;

  g_return_if_fail (tree != NULL);
  g_return_if_fail (pconfig != NULL);

  memcpy (&tree->config.parallel_config, pconfig,
          sizeof (VsgPRTreeParallelConfig));

  comm = pconfig->communicator;
  pt_vtable = &pconfig->point;
  rg_vtable = &pconfig->region;
  
  MPI_Comm_rank (comm, &rk);
  MPI_Comm_size (comm, &sz);


  if (rk != 0)
    {
      gint cnt;
      VsgPackedMsg pt_send = VSG_PACKED_MSG_STATIC_INIT (comm);
      VsgPackedMsg rg_send = VSG_PACKED_MSG_STATIC_INIT (comm);
      VTableAndMsg pt_vm = {&tree->config.parallel_config.point, &pt_send};
      VTableAndMsg rg_vm = {&tree->config.parallel_config.region, &rg_send};

      /* send points to 0 */
      cnt = vsg_prtree2@t@_point_count (tree);
      vsg_packed_msg_send_append (&pt_send, &cnt, 1, MPI_INT);
      vsg_prtree2@t@_foreach_point (tree, (GFunc) _send_point, &pt_vm);

      vsg_packed_msg_send (&pt_send, 0, 1001);
      vsg_packed_msg_drop_buffer (&pt_send);

      /* get 0's new bounds */
      MPI_Bcast (bounds, 2, VSG_MPI_TYPE_VECTOR2@T@, 0, comm);

/*       g_printerr ("%d: bounds", rk); */
/*       vsg_vector2d_write (&bounds[0], stderr); */
/*       vsg_vector2d_write (&bounds[1], stderr); */
/*       g_printerr ("\n"); */

      /* send regions to 0 */
      cnt = vsg_prtree2@t@_region_count (tree);
      vsg_packed_msg_send_append (&rg_send, &cnt, 1, MPI_INT);
      vsg_prtree2@t@_foreach_region (tree, (GFunc) _send_region,
                                     &rg_vm);

      vsg_packed_msg_send (&rg_send, 0, 1002);
      vsg_packed_msg_drop_buffer (&rg_send);

      /* Transform root node to an empty remote */
      vsg_prtree2@t@node_free (tree->node, &tree->config);

      tree->node = vsg_prtree2@t@node_alloc_no_data (&bounds[0], &bounds[1],
                                                     &tree->config);

      tree->node->parallel_status.storage = VSG_PARALLEL_REMOTE;
      tree->node->parallel_status.proc = 0;
    }
  else
    {
      gint src;
      VsgPackedMsg pt_recv = VSG_PACKED_MSG_STATIC_INIT (comm);
      VsgPackedMsg rg_recv = VSG_PACKED_MSG_STATIC_INIT (comm);

      /* receive other procs points */
      for (src=1; src<sz; src++)
        {
          gint i, cnt;

          vsg_packed_msg_recv (&pt_recv, src, 1001);
          vsg_packed_msg_recv_read (&pt_recv, &cnt, 1, MPI_INT);

          for (i=0; i<cnt; i++)
            {
              VsgPoint2 *pt = pt_vtable->alloc (TRUE, pt_vtable->alloc_data);

              pt_vtable->migrate.unpack (pt, &pt_recv,
                                         pt_vtable->migrate.unpack_data);
              vsg_prtree2@t@_insert_point (tree, pt);
            }

          vsg_packed_msg_drop_buffer (&pt_recv);
        }

      /* communicate new tree bounds to other processors */
      vsg_prtree2@t@_get_bounds (tree, &bounds[0], &bounds[1]);
      MPI_Bcast (bounds, 2, VSG_MPI_TYPE_VECTOR2@T@, 0, comm);

/*       g_printerr ("%d: bounds", rk); */
/*       vsg_vector2d_write (&bounds[0], stderr); */
/*       vsg_vector2d_write (&bounds[1], stderr); */
/*       g_printerr ("\n"); */

      /* receive other procs regions */
      for (src=1; src<sz; src++)
        {
          gint i, cnt;

          vsg_packed_msg_recv (&rg_recv, src, 1002);
          vsg_packed_msg_recv_read (&rg_recv, &cnt, 1, MPI_INT);

          for (i=0; i<cnt; i++)
            {
              VsgRegion2 rg = rg_vtable->alloc (TRUE, rg_vtable->alloc_data);

              rg_vtable->migrate.unpack (rg, &rg_recv,
                                         rg_vtable->migrate.unpack_data);
              vsg_prtree2@t@_insert_region (tree, rg);
            }

          vsg_packed_msg_drop_buffer (&rg_recv);
        }

    }
}

/* selector function used in traverse_custom_internal to avoid
   traversing all local nodes */
static vsgrloc2 _selector_skip_local_nodes (VsgRegion2 *selector,
                                            VsgPRTree2@t@NodeInfo *node_info,
                                            gpointer data)
{
  if (VSG_PRTREE2@T@_NODE_INFO_IS_LOCAL (node_info)) return 0x0;

  return VSG_RLOC2_MASK;
}

typedef struct _MigrateData MigrateData;
struct _MigrateData
{
  VsgParallelVTable *vtable;
  gint rk;
  VsgCommBuffer *cb;
};

static void _migrate_traverse_point_send (VsgPRTree2@t@Node *node,
                                          VsgPRTree2@t@NodeInfo *node_info,
                                          MigrateData *md)
{
  if (PRTREE2@T@NODE_IS_REMOTE (node) &&
      node_info->parallel_status.proc != md->rk)
    {
      VsgPackedMsg *pm = 
        vsg_comm_buffer_get_send_buffer (md->cb, node->parallel_status.proc);
      VTableAndMsg vm = {md->vtable, pm};

      g_slist_foreach (PRTREE2@T@NODE_LEAF (node).point,
                       (GFunc) _send_point, &vm);
      g_slist_free (PRTREE2@T@NODE_LEAF (node).point);
      PRTREE2@T@NODE_LEAF (node).point = NULL;
    }
}

static void _migrate_traverse_region_send (VsgPRTree2@t@Node *node,
                                           VsgPRTree2@t@NodeInfo *node_info,
                                           MigrateData *md)
{
  if (PRTREE2@T@NODE_IS_REMOTE (node) &&
      node_info->parallel_status.proc != md->rk)
    {
      VsgPackedMsg *pm = 
        vsg_comm_buffer_get_send_buffer (md->cb, node->parallel_status.proc);
      VTableAndMsg vm = {md->vtable, pm};

      g_slist_foreach (node->region_list, (GFunc) _send_region, &vm);
      g_slist_free (node->region_list);
      node->region_list = NULL;
    }
}

/* static void _print_region (gpointer rg, FILE *file) */
/* { */
/*   fprintf (file, "%p ", rg); */
/* } */

void vsg_prtree2@t@_migrate_flush (VsgPRTree2@t@ *tree)
{
  MPI_Comm comm;
  VsgCommBuffer *cb;
  VsgPackedMsg pm = VSG_PACKED_MSG_NULL;
  MigrateData md;
  VsgParallelVTable *pt_vtable;
  VsgParallelVTable *rg_vtable;
  VTableAndMsg vm;
  gint src, rk, sz;

  g_return_if_fail (tree != NULL);
  comm = tree->config.parallel_config.communicator;
  g_return_if_fail (comm != MPI_COMM_NULL);
  
  MPI_Comm_rank (comm, &rk);
  MPI_Comm_size (comm, &sz);

  cb = vsg_comm_buffer_new (comm);
  vsg_packed_msg_init (&pm, comm);

  pt_vtable = &tree->config.parallel_config.point;
  rg_vtable = &tree->config.parallel_config.region;

  md.vtable = pt_vtable;
  md.rk = rk;
  md.cb = cb;

  /* send the pending VsgPoint to their remote location */

  vsg_prtree2@t@_traverse_custom_internal (tree, G_PRE_ORDER,
                                           (VsgRegion2@t@InternalLocDataFunc)
                                           _selector_skip_local_nodes,
                                           NULL, NULL,
                                           (VsgPRTree2@t@InternalFunc)
                                           _migrate_traverse_point_send,
                                           &md);

  vsg_comm_buffer_share (cb);

  /* get local VsgPoint from remote processors */

  for (src=0; src<sz; src++)
    {
      VsgPackedMsg *recv;

      if (src == rk) continue;

      recv = vsg_comm_buffer_get_recv_buffer (cb, src);

      /* Reception buffers are allocated at the exact size of the message.
       * We can therefore use the advance of successive unpack to control the
       * number of received VsgPoint.
       */
      while (recv->position < recv->size)
        {
          VsgPoint2 *pt = pt_vtable->alloc (TRUE, pt_vtable->alloc_data);

          pt_vtable->migrate.unpack (pt, recv, pt_vtable->migrate.unpack_data);
          vsg_prtree2@t@_insert_point (tree, pt);
          
        }

      vsg_comm_buffer_drop_recv_buffer (cb, src);
    }

  if (rg_vtable->migrate.pack != NULL)
    {
      /* send the pending shared VsgRegion to every other processor */
      vsg_packed_msg_init (&pm, comm);

      vm.vtable = rg_vtable;
      vm.msg = &pm;

/*       g_printerr ("%d: pending shared regions: ", rk); */
/*       g_slist_foreach (tree->pending_shared_regions, */
/*                        (GFunc) _print_region, stderr); */
/*       g_printerr ("\n"); */

      g_slist_foreach (tree->pending_shared_regions,
                       (GFunc) _send_shared_region, &vm);

      g_slist_free (tree->pending_shared_regions);
      tree->pending_shared_regions = NULL;

      vsg_comm_buffer_set_bcast (cb, &pm);

      vsg_comm_buffer_share (cb);

      vsg_packed_msg_drop_buffer (&pm);

      /* send the pending (in "remote" nodes) VsgRegion to their destination */
      md.vtable = rg_vtable;
      md.rk = rk;
      md.cb = cb;

      vsg_prtree2@t@_traverse_custom_internal (tree, G_PRE_ORDER,
                                               (VsgRegion2@t@InternalLocDataFunc)
                                               _selector_skip_local_nodes,
                                               NULL, NULL,
                                               (VsgPRTree2@t@InternalFunc)
                                               _migrate_traverse_region_send,
                                               &md);

      vsg_comm_buffer_share (cb);

      /* get local VsgRegion from remote processors */

      for (src=0; src<sz; src++)
        {
          VsgPackedMsg *recv;

          if (src == rk) continue;

          recv = vsg_comm_buffer_get_recv_buffer (cb, src);

          /* Reception buffers are allocated at the exact size of the message.
           * We can therefore use the advance of successive unpack to control
           * the number of received VsgRegion.
           */
          while (recv->position < recv->size)
            {
              VsgRegion2 *rg = rg_vtable->alloc (TRUE, rg_vtable->alloc_data);

              rg_vtable->migrate.unpack (rg, recv,
                                         rg_vtable->migrate.unpack_data);
              vsg_prtree2@t@_insert_region (tree, rg);
            }

          vsg_comm_buffer_drop_recv_buffer (cb, src);
        }
    }

  vsg_comm_buffer_free (cb);

  /* remote trees depth may have changed */
  tree->config.remote_depth_dirty = TRUE;
}

static void _prtree2@t@node_fix_counts_local (VsgPRTree2@t@Node *node)
{
  /* *warning*: not a recursive function. Children counts are assumed valid */

  gint pcnt = 0;
  gint rcnt = 0;

  if (PRTREE2@T@NODE_ISINT (node))
    {
      vsgloc2 i;
      for (i=0; i<4; i++)
        {
          pcnt += PRTREE2@T@NODE_INT (node).children[i]->point_count;
          rcnt += PRTREE2@T@NODE_INT (node).children[i]->region_count;
        }
    }
  else
    {
      pcnt = g_slist_length (PRTREE2@T@NODE_LEAF (node).point);;
    }

  node->point_count = pcnt;
  node->region_count = rcnt + g_slist_length (node->region_list);
}

static void _migrate_pack_node_header (VsgPRTree2@t@NodeInfo *info,
                                       VsgPackedMsg *msg, gint dst,
                                       VsgPRTree2@t@Config *config)
{
  vsg_packed_msg_send_append (msg, &info->id, 1, VSG_MPI_TYPE_PRTREE_KEY2@T@);
  vsg_packed_msg_send_append (msg, &dst, 1, MPI_INT);
}

static void _migrate_pack_node (VsgPRTree2@t@Node *node, VsgPackedMsg *msg,
                                VsgPRTree2@t@Config *config,
                                gboolean shared)
{
  VsgPRTreeParallelConfig *pc = &config->parallel_config;
  gint datapresent, point_count, region_count;

  datapresent = pc->node_data.migrate.pack != NULL && node->user_data != NULL;
  point_count = 0;
  if (PRTREE2@T@NODE_ISLEAF (node))
    point_count = node->point_count;
  region_count = g_slist_length (node->region_list); /* use region number
                                                      * at this level only.
                                                      */

  /* pack message definition: */
  vsg_packed_msg_send_append (msg, &datapresent, 1, MPI_INT);
  vsg_packed_msg_send_append (msg, &point_count, 1, MPI_INT);
  vsg_packed_msg_send_append (msg, &region_count, 1, MPI_INT);

  /* pack the node data */
  if (datapresent)
    {
      pc->node_data.migrate.pack (node->user_data, msg,
                                  pc->node_data.migrate.pack_data);

      if (! shared)
        {
          if (pc->node_data.destroy)
            pc->node_data.destroy (node->user_data, TRUE,
                                   pc->node_data.destroy_data);
          else
            g_boxed_free (config->user_data_type, node->user_data);

          node->user_data = NULL;
        }
    }

  /* pack the points */
  if (point_count > 0)
    {
      VTableAndMsg vm = {&pc->point, msg};
      g_slist_foreach (PRTREE2@T@NODE_LEAF (node).point, (GFunc) _send_point,
                       &vm);

      g_slist_free (PRTREE2@T@NODE_LEAF (node).point);
      PRTREE2@T@NODE_LEAF (node).point = NULL;
    }
  node->point_count = 0; /* remove non local points count */

  /* pack the regions */
  if (region_count > 0)
    {
      VTableAndMsg vm = {&pc->region, msg};

      if (shared)
        g_slist_foreach (node->region_list, (GFunc) _send_shared_region, &vm);
      else
        {
          g_slist_foreach (node->region_list, (GFunc) _send_region, &vm);

          g_slist_free (node->region_list);
          node->region_list = NULL;
          node->region_count = 0;
        }
    }

  if (shared) _prtree2@t@node_fix_counts_local (node);
    
}

static gint _prtree2@t@node_get_children_proc (VsgPRTree2@t@Node *node)
{
  vsgloc2 i;
  gint children_proc;
  VsgPRTree2@t@Node *children[4];

  g_assert (PRTREE2@T@NODE_ISINT (node));

  memcpy (children, PRTREE2@T@NODE_INT (node).children,
          4 * sizeof (VsgPRTree2@t@Node *));

  if (PRTREE2@T@NODE_IS_SHARED (children[0]))
    return -1;
  else
    children_proc =
      PRTREE2@T@NODE_PROC (children[0]);

  for (i=1; i<4; i++)
    {
      gint proc =
        PRTREE2@T@NODE_PROC (children[i]);

      if (PRTREE2@T@NODE_IS_SHARED (children[i]) || proc != children_proc)
        return -1;
    }

  return children_proc;
}

static void _flatten_remote (VsgPRTree2@t@Node *node,
                             const VsgPRTree2@t@Config *config)
{
  /* destroy remaining children */
  if (PRTREE2@T@NODE_ISINT (node))
    {
      gint i;

      for (i=0; i<4; i++)
        {
          vsg_prtree2@t@node_free (PRTREE2@T@NODE_CHILD (node, i), config);
          PRTREE2@T@NODE_CHILD (node, i) = NULL;
        }
    }

  /* remote nodes aren't aware of points and regions stored on another
   * processor */
  g_slist_foreach (node->region_list, 
                   (GFunc) config->parallel_config.region.destroy,
                   config->parallel_config.region.destroy_data);
  g_slist_free (node->region_list);
  node->region_list = NULL;

  /* remove precedent user_data: not needed for a remote node */
  if (node->user_data != NULL)
    {
      const VsgPRTreeParallelConfig *pc = &config->parallel_config;

      if (pc->node_data.destroy)
        pc->node_data.destroy (node->user_data, TRUE,
                               pc->node_data.destroy_data);
      else
        g_boxed_free (config->user_data_type, node->user_data);

      node->user_data = NULL;
    }

  node->point_count = 0;
  node->region_count = 0;
}

static void _traverse_flatten_remote (VsgPRTree2@t@Node *node,
				      VsgPRTree2@t@NodeInfo *node_info,
				      const VsgPRTree2@t@Config *config)
{
  if (PRTREE2@T@NODE_IS_REMOTE (node))
    {
      _flatten_remote (node, config);
    }
}

typedef struct _DistributeData DistributeData;
struct _DistributeData {
  VsgPRTree2@t@DistributionFunc func;
  gpointer user_data;
  gint rk;
  gint sz;
  VsgCommBuffer *cb;
  VsgPackedMsg *bcast;
  VsgPRTree2@t@Config *config;
};

static void _traverse_distribute_nodes (VsgPRTree2@t@Node *node,
                                        VsgPRTree2@t@NodeInfo *node_info,
                                        DistributeData *dd)
{
  VsgParallelStorage old_storage = node->parallel_status.storage;
  VsgParallelStorage new_storage;
  gint dst;

  if (PRTREE2@T@NODE_ISLEAF (node))
    {
      g_assert (old_storage != VSG_PARALLEL_SHARED);

      dst = dd->func (node_info, dd->user_data);

      if (dst<0 || dst>=dd->sz) return; /* don't know: do nothing */
    }
  else
    {
      dst = _prtree2@t@node_get_children_proc (node);
    }

  if (dst < 0)
    new_storage = VSG_PARALLEL_SHARED;
  else
    new_storage = (dst == dd->rk) ? VSG_PARALLEL_LOCAL : VSG_PARALLEL_REMOTE;

  if (dst != dd->rk)
    {
      if (old_storage == VSG_PARALLEL_LOCAL)
        {
          VsgPackedMsg *msg;

          /* tell all processors about the migration. */
          _migrate_pack_node_header (node_info, dd->bcast, dst, dd->config);

          /* choose between ptp communication or bcast */
          if (new_storage == VSG_PARALLEL_REMOTE)
            msg = vsg_comm_buffer_get_send_buffer (dd->cb, dst);
          else
            msg = dd->bcast;

          /* write node to the corresponding message */
          _migrate_pack_node (node, msg, dd->config, dst < 0);
        }

    }

  if (old_storage == VSG_PARALLEL_REMOTE) return; /* check migrations later */

  if (new_storage == VSG_PARALLEL_REMOTE)
    {
     _flatten_remote (node, dd->config);
    }

  /* update node's parallel status */
  node->parallel_status.storage = new_storage;
  node->parallel_status.proc = dst;
}

void vsg_prtree2@t@node_insert_child (VsgPRTree2@t@Node *node,
                                      const VsgPRTree2@t@Config *config,
                                      VsgParallelStorage storage,
                                      gint dst,
                                      VsgPRTreeKey2@t@ key,
                                      gpointer user_data,
                                      GSList *points,
                                      GSList *regions)
{
  if (key.depth > 0)
    {
      vsgloc2 child = vsg_prtree_key2@t@_child (&key);
      gint children_proc;

      if (PRTREE2@T@NODE_ISLEAF (node))
        {
          /* formerly remote nodes need a new user_data to be allocated */
          if (PRTREE2@T@NODE_IS_REMOTE (node))
            {
              const VsgPRTreeParallelConfig *pc = &config->parallel_config;

              g_assert (node->user_data == NULL);

              if (pc->node_data.alloc != NULL)
                {
                  node->user_data =
                    pc->node_data.alloc (TRUE, pc->node_data.alloc_data);
                }
              else if (config->user_data_type != G_TYPE_NONE)
                {
                  node->user_data = g_boxed_copy (config->user_data_type,
                                                  config->user_data_model);
                }
            }

          /* transform into interior node */
          vsg_prtree2@t@node_make_int (node, config);
        }

      key.depth --;

      vsg_prtree2@t@node_insert_child (PRTREE2@T@NODE_CHILD (node, child),
                                       config, storage, dst, key, user_data,
                                       points, regions);

      _prtree2@t@node_fix_counts_local (node);

      children_proc = _prtree2@t@node_get_children_proc (node);

      if (children_proc < 0)
        {
          node->parallel_status.storage = VSG_PARALLEL_SHARED;
          node->parallel_status.proc = -1; /* set an invalid value */
        }
      else
        {
          g_assert (children_proc == dst);

          if (storage == VSG_PARALLEL_REMOTE)
            {
              _flatten_remote (node, config);
            }

          node->parallel_status.storage = storage;
          node->parallel_status.proc = dst;
        }

      return;
    }

  if (storage == VSG_PARALLEL_REMOTE)
    {
      _flatten_remote (node, config);
    }

  if (user_data != NULL)
    {
      if (node->user_data != NULL)
        {
          const VsgPRTreeParallelConfig *pc = &config->parallel_config;

          /* destroy precedent user_data */
          if (pc->node_data.destroy)
            pc->node_data.destroy (node->user_data, TRUE,
                                   pc->node_data.destroy_data);
          else
            g_boxed_free (config->user_data_type, node->user_data);
        }
      else
        {
          g_assert (node->parallel_status.storage == VSG_PARALLEL_REMOTE);
        }

      node->user_data = user_data;
    }

  if (points != NULL)
    {
      g_assert (PRTREE2@T@NODE_ISLEAF (node));

      PRTREE2@T@NODE_LEAF (node).point =
        g_slist_concat (points, PRTREE2@T@NODE_LEAF (node).point);
    }

  if (regions != NULL)
    node->region_list = g_slist_concat (regions, node->region_list);

  _prtree2@t@node_fix_counts_local (node);

  node->parallel_status.storage = storage;
  node->parallel_status.proc = dst;
}

void vsg_prtree2@t@_distribute_nodes (VsgPRTree2@t@ *tree,
                                      VsgPRTree2@t@DistributionFunc func,
                                      gpointer user_data)
{
  MPI_Comm comm;
  VsgCommBuffer *cb;
  VsgCommBuffer *bcastcb;
  VsgPackedMsg bcast = VSG_PACKED_MSG_NULL;
  VsgPRTreeParallelConfig *pc;
  DistributeData dd;
  gint rk, sz, src;

  g_return_if_fail (tree != NULL);
  pc = &tree->config.parallel_config;
  comm = pc->communicator;
  g_return_if_fail (comm != MPI_COMM_NULL);
  
  MPI_Comm_rank (comm, &rk);
  MPI_Comm_size (comm, &sz);

  cb = vsg_comm_buffer_new (comm);
  bcastcb = vsg_comm_buffer_new (comm);
  vsg_packed_msg_init (&bcast, comm);

  dd.func = func;
  dd.user_data = user_data;
  dd.rk = rk;
  dd.sz = sz;
  dd.cb = cb;
  dd.bcast = &bcast;
  dd.config = &tree->config;

/*   g_printerr ("%d: before gather\n", rk); */

  /* gather all migration messages */
  vsg_prtree2@t@_traverse_custom_internal (tree, G_POST_ORDER, NULL, NULL, NULL,
					   (VsgPRTree2@t@InternalFunc)
					   _traverse_distribute_nodes, &dd);

/*   g_printerr ("%d: after gather\n", rk); */

  /* send/receive pending messages */
  vsg_comm_buffer_set_bcast (bcastcb, &bcast);

/*   g_printerr ("%d: before share (bcast = %dBytes)\n", rk, bcast.position); */

  vsg_comm_buffer_share (bcastcb);

/*   g_printerr ("%d: bcast ok\n%d: before private msg\n", rk, rk); */

  vsg_comm_buffer_share (cb);

  vsg_packed_msg_drop_buffer (&bcast);

/*   g_printerr ("%d: after share\n", rk); */

  /* decode received messages */
  for (src=0; src<sz; src++)
    {
      VsgPackedMsg *hdrmsg, *msg;

      if (src == rk) continue;

      hdrmsg = vsg_comm_buffer_get_recv_buffer (bcastcb, src);
      msg = vsg_comm_buffer_get_recv_buffer (cb, src);

      while (hdrmsg->position < hdrmsg->size)
        {
          VsgPRTreeKey2@t@ id;
          gint dst;

          /* unpack the header */
          vsg_packed_msg_recv_read (hdrmsg, &id, 1,
                                    VSG_MPI_TYPE_PRTREE_KEY2@T@);
          vsg_packed_msg_recv_read (hdrmsg, &dst, 1, MPI_INT);

/*           g_printerr ("%d: unpacking src=%d id=[", rk, src); */
/*           vsg_prtree_key2@t@_write (&id, stderr); */
/*           g_printerr ("] dst=%d\n", dst); */

          if (dst == rk || dst < 0)
            { /* we are destination of thi message (bcast or not) */
              VsgParallelStorage storage =
                (dst >= 0) ? VSG_PARALLEL_LOCAL : VSG_PARALLEL_SHARED;
              VsgPackedMsg *unpack = (dst >= 0) ? msg : hdrmsg;
              gpointer data = NULL;
              GSList *points = NULL;
              GSList *regions = NULL;
              gint datapresent, point_count, region_count, i;

              /* load node contents */

              vsg_packed_msg_recv_read (unpack, &datapresent, 1, MPI_INT);
              vsg_packed_msg_recv_read (unpack, &point_count, 1, MPI_INT);
              vsg_packed_msg_recv_read (unpack, &region_count, 1, MPI_INT);
/*               g_printerr ("%d: unpacking src=%d : d=%d p=%d r=%d\n", */
/*                           rk, src, datapresent, point_count, region_count); */

              if (datapresent)
                {
                  /* create a temporary node_data */
                  data = pc->node_data.alloc (FALSE, pc->node_data.alloc_data);

                  pc->node_data.migrate.unpack (data, unpack,
                                                pc->node_data.migrate.unpack_data);
                }

              g_assert (point_count == 0 || storage == VSG_PARALLEL_LOCAL);

              for (i=0; i<point_count; i++)
                {
                  VsgPoint2 pt =
                    pc->point.alloc (TRUE, pc->point.alloc_data);

                  pc->point.migrate.unpack (pt, unpack,
                                            pc->point.migrate.unpack_data);

                  points = g_slist_prepend (points, pt);
                }

              for (i=0; i<region_count; i++)
                {
                  VsgRegion2 pt =
                    pc->region.alloc (TRUE, pc->region.alloc_data);

                  pc->region.migrate.unpack (pt, unpack,
                                             pc->region.migrate.unpack_data);

                  regions = g_slist_prepend (regions, pt);
                }

              /* insert the node in the tree */
              vsg_prtree2@t@node_insert_child (tree->node, &tree->config,
                                               storage, dst, id, data,
                                               points, regions);
            }
          else
            {
              /* we just witness the migration but we have to keep track of
               * remote nodes location.
               */

              /* insert the node in remote mode */
              vsg_prtree2@t@node_insert_child (tree->node, &tree->config,
                                               VSG_PARALLEL_REMOTE, dst,
                                               id, NULL, NULL, NULL);
              
            }

/*           g_printerr ("%d: unpacking done\n", rk); */
        }
    }

  vsg_comm_buffer_free (bcastcb);
  vsg_comm_buffer_free (cb);

  /* fix all remote nodes: remove remaining subtree and locally stored
   * regions */
  vsg_prtree2@t@_traverse_custom_internal (tree, G_POST_ORDER,
                                           (VsgRegion2@t@InternalLocDataFunc)
                                           _selector_skip_local_nodes,
                                           NULL, NULL,
					   (VsgPRTree2@t@InternalFunc)
					   _traverse_flatten_remote,
					   &tree->config);

  /* remote trees depth may have changed */
  tree->config.remote_depth_dirty = TRUE;
}

/* #include <sys/types.h> */
/* #include <unistd.h> */

#define DIRECTION_FORWARD (0)
#define DIRECTION_BACKWARD (1)

/* messages tags */
#define VISIT_FORWARD_TAG (100)
#define VISIT_BACKWARD_TAG (101)
#define END_FW_TAG (102)
#define VISIT_SHARED_TAG (103)


/*
 * Holds a visiting node while it waits to be sent back to its original
 * processor.
 */
typedef struct _WaitingVisitor WaitingVisitor;
struct _WaitingVisitor {
  VsgPRTree2@t@Node *node;
  VsgPRTreeKey2@t@ id;
  gint8 flag;
  gint src;
};

static WaitingVisitor *_waiting_visitor_new (VsgPRTree2@t@Node *node,
                                             VsgPRTreeKey2@t@ *id,
                                             gint src)
{
  WaitingVisitor *wv = g_slice_new (WaitingVisitor);

  wv->node = node;
  wv->id = *id;
  wv->src = src;
  wv->flag = 0;

  return wv;
}

static void _waiting_visitor_free (WaitingVisitor *wv)
{
  g_slice_free (WaitingVisitor, wv);
}

static void _wv_free_gfunc (WaitingVisitor *wv, gpointer data)
{
  _waiting_visitor_free (wv);
}

/*
 * Holds a message and the request associated with a particular transfer (send
 * or receive).
 */
typedef struct _VsgNFProcMsg VsgNFProcMsg;
struct _VsgNFProcMsg {
  VsgPackedMsg send_pm;

  GSList *forward_pending;
  GSList *backward_pending;

  MPI_Request request;
};

static void vsg_nf_proc_msg_init (VsgNFProcMsg *nfpm, MPI_Comm comm)
{
  vsg_packed_msg_init (&nfpm->send_pm, comm);

  nfpm->forward_pending = NULL;
  nfpm->backward_pending = NULL;

  nfpm->request = MPI_REQUEST_NULL;
}

static VsgNFProcMsg *vsg_nf_proc_msg_new (MPI_Comm comm)
{
  VsgNFProcMsg *nfpm = g_slice_new0 (VsgNFProcMsg);

  vsg_nf_proc_msg_init(nfpm, comm);

  return nfpm;
}

static void vsg_nf_proc_msg_free (VsgNFProcMsg *nfpm)
{
  vsg_packed_msg_drop_buffer (&nfpm->send_pm);

  g_slist_foreach (nfpm->forward_pending, (GFunc) _wv_free_gfunc, NULL);
  g_slist_free (nfpm->forward_pending);
  g_slist_foreach (nfpm->backward_pending, (GFunc) _wv_free_gfunc, NULL);
  g_slist_free (nfpm->backward_pending);

  g_slice_free (VsgNFProcMsg, nfpm);
}

/*
 * Initializes a VsgNFConfig2@t@.
 */
void vsg_nf_config2@t@_init (VsgNFConfig2@t@ *nfc,
                             MPI_Comm comm,
                             VsgPRTree2@t@FarInteractionFunc far_func,
                             VsgPRTree2@t@InteractionFunc near_func,
                             gpointer user_data)
{
  nfc->far_func = far_func;
  nfc->near_func = near_func;
  nfc->user_data = user_data;

  vsg_packed_msg_init (&nfc->recv, comm);
  nfc->procs_msgs = NULL;

   if (comm != MPI_COMM_NULL)
    {
      MPI_Comm_rank (comm, &nfc->rk);
      MPI_Comm_size (comm, &nfc->sz);

      nfc->procs_msgs =
        g_hash_table_new_full (g_direct_hash, g_direct_equal, NULL,
                               (GDestroyNotify) vsg_nf_proc_msg_free);
    }
  else
    {
      nfc->rk = 0;
      nfc->sz = 1;
    }

   nfc->tmp_node_data = NULL;

   nfc->forward_pending_nb = 0;
   nfc->backward_pending_nb = 0;

   nfc->pending_end_forward = nfc->sz-1;
   nfc->pending_backward_msgs = 0;

   nfc->all_fw_sends = 0; nfc->all_fw_recvs = 0;
   nfc->all_bw_sends = 0; nfc->all_bw_recvs = 0;

   nfc->shared_far_interaction_counter = 0;
}

/*
 * Called before every shared/shared far interaction, it computes wether the
 * current processor has to do it or no.
 */
gboolean vsg_nf_config2@t@_shared_far_interaction_skip (VsgNFConfig2@t@ *nfc)
{
  /* to dispatch far interaction between shared nodes on all processors */
  gboolean doit = nfc->shared_far_interaction_counter % nfc->sz != nfc->rk;

  nfc->shared_far_interaction_counter ++;

  return doit;
}

/*
 * Frees all memory allocated inside a VsgNFConfig2@t@. (doesn't free @nfc
 * itself).
 */
void vsg_nf_config2@t@_clean (VsgNFConfig2@t@ *nfc)
{
  vsg_packed_msg_drop_buffer (&nfc->recv);

  if (nfc->procs_msgs != NULL)
    g_hash_table_destroy (nfc->procs_msgs);

}

/*
 * Searchs for the VsgNFProcMsg associated with a particular processor.
 * Allocates a new one if necessary.
 */
static VsgNFProcMsg * vsg_nf_config2@t@_proc_msgs_lookup (VsgNFConfig2@t@ *nfc,
                                                          gint proc)
{
  VsgNFProcMsg *nfpm = g_hash_table_lookup (nfc->procs_msgs,
                                            GINT_TO_POINTER (proc));

  if (nfpm == NULL)
    {
      nfpm = vsg_nf_proc_msg_new (nfc->recv.communicator);

      g_hash_table_insert (nfc->procs_msgs, GINT_TO_POINTER (proc), nfpm);
    }

  return nfpm;
}

static void _visit_forward_pack (gpointer data, VTableAndMsg *vm)
{
  vm->vtable->visit_forward.pack (data, vm->msg,
                                  vm->vtable->visit_forward.pack_data);
}

static void _visit_backward_pack (gpointer data, VTableAndMsg *vm)
{
  vm->vtable->visit_backward.pack (data, vm->msg,
                                   vm->vtable->visit_backward.pack_data);
}

static void _visit_destroy_point (VsgPoint2 pt, VsgParallelVTable *vtable)
{
  vtable->destroy (pt, FALSE, vtable->destroy_data);
}

static void _visit_destroy_region (VsgRegion2 rg, VsgParallelVTable *vtable)
{
  vtable->destroy (rg, FALSE, vtable->destroy_data);
}

/*
 * Allocates a floating VsgPRTree2@t@Node to hold a visiting node from another
 * processor.
 */
static VsgPRTree2@t@Node *_new_visiting_node (VsgPRTree2@t@ *tree,
                                              VsgPRTreeKey2@t@ *id,
                                              gint src)
{
  VsgPRTree2@t@Config *config = &tree->config;
  VsgPRTreeParallelConfig *pc = &config->parallel_config;
  VsgPRTree2@t@Node *node;
  VsgVector2@t@ *lb = &tree->node->lbound;
  VsgVector2@t@ *ub = &tree->node->ubound;
  VsgVector2@t@ newlb, newub;

  /* compute the bounds of the visiting node */
  if (id->depth == 0)
    {
      vsg_vector2@t@_copy (lb, &newlb);
      vsg_vector2@t@_copy (ub, &newub);
    }
  else
    {
      VsgVector2@t@ delta;

      vsg_vector2@t@_sub (ub, lb, &delta);

      vsg_vector2@t@_scalp (&delta, 1. / (1 << id->depth), &delta);

      newlb.x = lb->x + id->x * delta.x;
      newlb.y = lb->y + id->y * delta.y;

      vsg_vector2@t@_add (&newlb, &delta, &newub);
    }

  /* allocate new node without node's user_data */
  node = vsg_prtree2@t@node_alloc (&newlb, &newub, NULL);

  /* separately allocate node's user_data */
  if (pc->node_data.alloc != NULL)
    {
      /* allocate non resident data */
      node->user_data = pc->node_data.alloc (FALSE, pc->node_data.alloc_data);
    }
  else
    {
      /* fallback with the GBoxed method */
      if (config->user_data_type != G_TYPE_NONE)
        node->user_data = g_boxed_copy (config->user_data_type,
                                     config->user_data_model);
    }

  node->parallel_status.storage = VSG_PARALLEL_REMOTE;
  node->parallel_status.proc = src;

  return node;
}

/*
 * Deallocates the VsgPRTree2@t@Node structure used by a visiting node.
 */
static void _visiting_node_free (VsgPRTree2@t@Node *node,
                                VsgPRTree2@t@Config *config)
{
  VsgPRTreeParallelConfig *pc = &config->parallel_config;

  if (node->user_data != NULL)
    {
      if (pc->node_data.destroy)
        pc->node_data.destroy (node->user_data, FALSE,
                               pc->node_data.destroy_data);
      else
        g_boxed_free (config->user_data_type, node->user_data);

      node->user_data = NULL;
    }

  if (PRTREE2@T@NODE_ISLEAF (node))
    {
      g_slist_foreach (PRTREE2@T@NODE_LEAF (node).point,
                       (GFunc) _visit_destroy_point,
                       &pc->point);

      g_slist_free (PRTREE2@T@NODE_LEAF (node).point);
      PRTREE2@T@NODE_LEAF (node).point = NULL;
    }

  g_slist_foreach (node->region_list, (GFunc) _visit_destroy_region,
                   &pc->region);

  g_slist_free (node->region_list);
  node->region_list = NULL;

  vsg_prtree2@t@node_free (node, config);
}

/*
 * Packs a node into a VsgPackedMsg according to the diirection specified and
 * wether the node is to be deallocated when done or not.
 */
static void _visit_pack_node (VsgPRTree2@t@Node *node, VsgPackedMsg *msg,
                              VsgPRTree2@t@Config *config,
                              gint8 direction, gboolean shared)
{
  VsgPRTreeParallelConfig *pc = &config->parallel_config;
  VsgParallelMigrateVTable *node_data_migrate;
  GFunc send;

  gint8 datapresent;
  gint point_count;
  gint region_count;
  gboolean do_data, do_point, do_region;

  if (direction == DIRECTION_FORWARD)
    {
      /* configure for forward sends and check availability of functions */
      node_data_migrate = &pc->node_data.visit_forward;
      send = (GFunc) _visit_forward_pack;

      do_data = node_data_migrate->pack != NULL &&
        node_data_migrate->unpack != NULL;
      do_point = pc->point.visit_forward.pack != NULL &&
        pc->point.visit_forward.unpack != NULL;
      do_region = pc->region.visit_forward.pack != NULL &&
        pc->region.visit_forward.unpack != NULL;
    }
  else
    {
      /* configure for backward sends and check availability of functions */
      node_data_migrate = &pc->node_data.visit_backward;
      send = (GFunc) _visit_backward_pack;

      do_data = node_data_migrate->pack != NULL &&
        node_data_migrate->unpack != NULL;
      do_point = pc->point.visit_backward.pack != NULL &&
        pc->point.visit_backward.unpack != NULL;
      do_region = pc->region.visit_backward.pack != NULL &&
        pc->region.visit_backward.unpack != NULL;
    }

  datapresent = do_data && node->user_data != NULL;

  point_count = 0;
  if (do_point && PRTREE2@T@NODE_ISLEAF (node))
    point_count = node->point_count;

  region_count = 0;
  if (do_region)
    region_count = g_slist_length (node->region_list); /* use region number
                                                        * at this level only.
                                                        */

  /* pack message definition: */
  vsg_packed_msg_send_append (msg, &datapresent, 1, MPI_BYTE);
  vsg_packed_msg_send_append (msg, &point_count, 1, MPI_INT);
  vsg_packed_msg_send_append (msg, &region_count, 1, MPI_INT);

  /* pack the node data */
  if (datapresent)
    {
      node_data_migrate->pack (node->user_data, msg,
                               node_data_migrate->pack_data);
    }

  /* pack the points */
  if (point_count > 0)
    {
      VTableAndMsg vm = {&pc->point, msg};
      g_slist_foreach (PRTREE2@T@NODE_LEAF (node).point, send,
                       &vm);

    }

  /* pack the regions */
  if (region_count > 0)
    {
      VTableAndMsg vm = {&pc->region, msg};

      g_slist_foreach (node->region_list, (GFunc) send, &vm);
    }

  if (direction == DIRECTION_BACKWARD && ! shared)
    _visiting_node_free (node, config);

}

/*
 * Reads a node from a VsgPackedMsg and stores it into a VsgPRTree2@t@Node
 * according to the value of the @direction argument.
 */
static void _visit_unpack_node (VsgPRTree2@t@Config *config,
                                VsgPRTree2@t@Node *node,
                                VsgNFConfig2@t@ *nfc,
                                gint8 direction)
{
  VsgPRTreeParallelConfig *pc = &config->parallel_config;
  gint8 datapresent;
  gint i, point_count, region_count;
  VsgParallelMigrateVTable *point;
  VsgParallelMigrateVTable *region;
  VsgParallelMigrateVTable *node_data;

  if (direction == DIRECTION_FORWARD)
    {
      point = &pc->point.visit_forward;
      region = &pc->region.visit_forward;
      node_data = &pc->node_data.visit_forward;
    }
  else
    {
      point = &pc->point.visit_backward;
      region = &pc->region.visit_backward;
      node_data = &pc->node_data.visit_backward;
    }

  /* unpack message definition */
  vsg_packed_msg_recv_read (&nfc->recv, &datapresent, 1, MPI_BYTE);
  vsg_packed_msg_recv_read (&nfc->recv, &point_count, 1, MPI_INT);
  vsg_packed_msg_recv_read (&nfc->recv, &region_count, 1, MPI_INT);

  /* unpack user_data */
  if (datapresent)
    {
      if (direction == DIRECTION_FORWARD)
        {
          node_data->unpack (node->user_data, &nfc->recv,
                             node_data->unpack_data);
        }
      else
        {
           node_data->unpack (nfc->tmp_node_data, &nfc->recv,
                              node_data->unpack_data);
           node_data->reduce (nfc->tmp_node_data, node->user_data,
                              node_data->reduce_data);
        }
    }

  /* unpack the points */
  if (PRTREE2@T@NODE_ISLEAF(node) && point_count > 0)
    {
      GSList *points = NULL;
      VsgPoint2 pt;

      if (direction == DIRECTION_FORWARD)
        {
          /* create a new list of points */
          for (i=0; i<point_count; i++)
            {
              pt = pc->point.alloc (FALSE, pc->point.alloc_data);

              point->unpack (pt, &nfc->recv, point->unpack_data);

              points = g_slist_append (points, pt);
            }

          PRTREE2@T@NODE_LEAF(node).point = points;
          node->point_count = point_count;
        }
      else
        {
          g_assert (point_count == node->point_count);

          /* unpack on existing points */
          points = PRTREE2@T@NODE_LEAF(node).point;
          while (points != NULL)
            {
              pt = (VsgPoint2) points->data;

              point->unpack (pt, &nfc->recv, point->unpack_data);

              points = g_slist_next (points);
            }
        }
    }

  /* unpack the regions */
  if (region_count > 0)
    {
      GSList *regions = NULL;
      VsgRegion2 pt;

      if (direction == DIRECTION_FORWARD)
        {
          /* create a new list of regions */
          for (i=0; i<region_count; i++)
            {
              pt = pc->region.alloc (FALSE, pc->region.alloc_data);

              region->unpack (pt, &nfc->recv, region->unpack_data);

              regions = g_slist_append (regions, pt);
            }

          node->region_list = regions;
          node->region_count = region_count;
        }
      else
        {
          g_assert (region_count == g_slist_length (node->region_list));

          /* unpack on existing regions */
          regions = node->region_list;
          while (regions != NULL)
            {
              pt = (VsgRegion2) regions->data;

              region->unpack (pt, &nfc->recv, region->unpack_data);

              regions = g_slist_next (regions);
            }
        }
    }
}

static void _do_send_forward_node (VsgPRTree2@t@ *tree,
                                   VsgNFConfig2@t@ *nfc,
                                   VsgNFProcMsg *nfpm,
                                   VsgPRTree2@t@Node *node,
                                   VsgPRTreeKey2@t@ *id,
                                   gint proc)
{
  VsgPackedMsg *msg = &nfpm->send_pm;

  msg->position = 0;

  vsg_packed_msg_send_append (msg, id, 1, VSG_MPI_TYPE_PRTREE_KEY2@T@);
  _visit_pack_node (node, msg, &tree->config, DIRECTION_FORWARD, FALSE);

  vsg_packed_msg_isend (msg, proc, VISIT_FORWARD_TAG, &nfpm->request);

  nfc->all_fw_sends ++;
}

static void _send_pending_forward_node (VsgPRTree2@t@ *tree,
                                        VsgNFConfig2@t@ *nfc,
                                        VsgNFProcMsg *nfpm,
                                        gint proc)
{
  GSList *first = nfpm->forward_pending;
  WaitingVisitor *wv = (WaitingVisitor *) first->data;

  nfpm->forward_pending = g_slist_next (nfpm->forward_pending);
  g_slist_free1 (first);

  _do_send_forward_node (tree, nfc, nfpm, wv->node, &wv->id, proc);

  _waiting_visitor_free (wv);

  nfc->forward_pending_nb --;
  nfc->pending_backward_msgs ++;
}

/* static gint _bw_pending_max = 0; */
/* static gint _bw_pending_sum = 0; */
/* static gint _bw_pending_calls = 0; */

static void _do_send_backward_node (VsgPRTree2@t@ *tree,
                                    VsgNFConfig2@t@ *nfc,
                                    VsgNFProcMsg *nfpm,
                                    VsgPRTree2@t@Node *node,
                                    VsgPRTreeKey2@t@ *id,
                                    gint proc)
{
  VsgPackedMsg *msg = &nfpm->send_pm;

/*   gint _bw_pending_length = g_slist_length (nfpm->backward_pending); */

/*   _bw_pending_calls ++; */
/*   _bw_pending_sum += _bw_pending_length; */
/*   _bw_pending_max = MAX (_bw_pending_max, _bw_pending_length); */

  msg->position = 0;

  vsg_packed_msg_send_append (msg, id, 1, VSG_MPI_TYPE_PRTREE_KEY2@T@);
  _visit_pack_node (node, msg, &tree->config, DIRECTION_BACKWARD, FALSE);

  vsg_packed_msg_isend (msg, proc, VISIT_BACKWARD_TAG, &nfpm->request);

  nfc->all_bw_sends ++;
}

static void _send_pending_backward_node (VsgPRTree2@t@ *tree,
                                         VsgNFConfig2@t@ *nfc,
                                         VsgNFProcMsg *nfpm,
                                         gint proc)
{
  /* check for backward messages */
  GSList *first = nfpm->backward_pending;
  WaitingVisitor *wv = (WaitingVisitor *) first->data;

  nfpm->backward_pending = g_slist_next (nfpm->backward_pending);
  g_slist_free1 (first);

  _do_send_backward_node (tree, nfc, nfpm, wv->node,
                          &wv->id, proc);

  _waiting_visitor_free (wv);

  nfc->backward_pending_nb --;
}

static void _propose_node_forward (VsgPRTree2@t@ *tree,
                                   VsgNFConfig2@t@ *nfc,
                                   gint proc,
                                   VsgPRTree2@t@Node *node,
                                   VsgPRTreeKey2@t@ *id)
{
  gint flag;
  VsgNFProcMsg *nfpm;

  nfpm = vsg_nf_config2@t@_proc_msgs_lookup (nfc, proc);

  MPI_Test (&nfpm->request, &flag, MPI_STATUS_IGNORE);

  if (flag && nfpm->backward_pending == NULL)
    {
      _do_send_forward_node (tree, nfc, nfpm, node, id, proc);

/*       g_printerr ("%d : propose fw direct to %d - ", nfc->rk, proc); */
/*       vsg_prtree_key2@t@_write (id, stderr); */
/*       g_printerr ("\n"); */

      nfc->pending_backward_msgs ++;
    }
  else
    {
      WaitingVisitor *wv = _waiting_visitor_new (node, id, proc);

      nfpm->forward_pending = g_slist_append (nfpm->forward_pending, wv);

/*       g_printerr ("%d : propose fw delay to %d - ", nfc->rk, proc); */
/*       vsg_prtree_key2@t@_write (id, stderr); */
/*       g_printerr ("\n"); */

      if (flag)
        {
          _send_pending_backward_node (tree, nfc, nfpm, proc);
        }

      nfc->forward_pending_nb ++;
    }

}

static void _propose_node_backward (VsgPRTree2@t@ *tree,
                                    VsgNFConfig2@t@ *nfc,
                                    gint proc,
                                    WaitingVisitor *wv)
{
  gint flag;
  VsgNFProcMsg *nfpm;

  nfpm = vsg_nf_config2@t@_proc_msgs_lookup (nfc, proc);

  MPI_Test (&nfpm->request, &flag, MPI_STATUS_IGNORE);

  if (flag != 0)
    {
      _do_send_backward_node (tree, nfc, nfpm, wv->node, &wv->id, proc);

/*       g_printerr ("%d : propose bw direct to %d - ", nfc->rk, proc); */
/*       vsg_prtree_key2@t@_write (&wv->id, stderr); */
/*       g_printerr ("\n"); */

      _waiting_visitor_free (wv);
    }
  else
    {
      nfpm->backward_pending = g_slist_append (nfpm->backward_pending, wv);

/*       g_printerr ("%d : propose bw delay to %d - ", nfc->rk, proc); */
/*       vsg_prtree_key2@t@_write (&wv->id, stderr); */
/*       g_printerr ("\n"); */

      nfc->backward_pending_nb ++;
    }

}

/*
 * Stores visitor's node info and near/far functions for a
 * vsg_prtree2@t@_traverse operation.
 */
typedef struct _NIAndFuncs NIAndFuncs;
struct _NIAndFuncs {
  VsgPRTree2@t@Node *ref;
  VsgPRTree2@t@NodeInfo ref_info;
  VsgPRTree2@t@FarInteractionFunc far_func;
  VsgPRTree2@t@InteractionFunc near_func;
  gpointer user_data;
  VsgPRTreeKey2@t@ *ref_ancestry_ids;
  gint8 done_flag;
};

/*
 * computes near/far interactions for a visiting node during a tree traversal.
 */
static void _traverse_visiting_nf (VsgPRTree2@t@Node *node,
                                   VsgPRTree2@t@NodeInfo *node_info,
                                   NIAndFuncs *niaf)
{
  VsgPRTree2@t@NodeInfo *ref_info = &niaf->ref_info;
  guint8 node_depth;
  gboolean fardone;
/*       gint rk; */

/*       MPI_Comm_rank (MPI_COMM_WORLD, &rk); */

  node_depth = node_info->id.depth;

  if (node_depth <= ref_info->id.depth)
    {
      VsgPRTreeKey2@t@ *ref_id = &niaf->ref_ancestry_ids[node_depth];
      if (vsg_prtree_key2@t@_is_neighbour (ref_id, &node_info->id))
        {
          if (ref_info->point_count == 0 || node->point_count == 0) return;

          if (PRTREE2@T@NODE_ISLEAF (node) ||
              node_depth == ref_info->id.depth)
            {
/*         g_printerr ("%d : near interaction [", rk); */
/*         vsg_prtree_key2@t@_write (&ref_info->id, stderr); */
/*         g_printerr ("] ["); */
/*         vsg_prtree_key2@t@_write (&node_info->id, stderr); */
/*         g_printerr ("]\n"); */

              vsg_prtree2@t@node_recursive_near_func (niaf->ref, ref_info,
                                                      node, node_info,
                                                      niaf->near_func,
                                                      niaf->user_data);
              niaf->done_flag |= 1;
            }
        }
      else if (node_depth == ref_info->id.depth &&
               PRTREE2@T@NODE_IS_LOCAL (node) &&
               node_info->point_count != 0)
        {
/*         g_printerr ("%d : far interaction [", rk); */
/*         vsg_prtree_key2@t@_write (&ref_info->id, stderr); */
/*         g_printerr ("] ["); */
/*         vsg_prtree_key2@t@_write (&node_info->id, stderr); */
/*         g_printerr ("]\n"); */

          fardone = niaf->far_func (ref_info, node_info, niaf->user_data);
          if (! fardone)
            g_critical ("far_func() -> FALSE not handled in \"%s\"",
                        __PRETTY_FUNCTION__);
          niaf->done_flag |= 1<<1;
        }
    }
}

static vsgrloc2 _selector_nf_visitor (VsgPRTree2@t@NodeInfo *ref_info,
                                      VsgPRTree2@t@NodeInfo *node_info,
                                      VsgPRTreeKey2@t@ *ref_ancestry_ids)
{
  static const vsgrloc2 ancestor_order_rlocs[] = {
    VSG_RLOC2_MASK - (VSG_RLOC2_COMP (0)-1),
    VSG_RLOC2_MASK - (VSG_RLOC2_COMP (1)-1),
    VSG_RLOC2_MASK - (VSG_RLOC2_COMP (2)-1),
    VSG_RLOC2_MASK - (VSG_RLOC2_COMP (3)-1),
  };
  VsgPRTreeKey2@t@ *node_key = &node_info->id;
  gint d = node_key->depth;

  if (d >= ref_info->id.depth)
    {
/*       if (ref_info->point_count == 0/\*  || d > ref_info->id.depth *\/) */
        return 0x0;

      return VSG_RLOC2_MASK;
    }

  if (vsg_prtree_key2@t@_is_neighbour (&ref_ancestry_ids[d], node_key))
    {
      if (vsg_prtree_key2@t@_equals (&ref_ancestry_ids[d], node_key))
        {
          vsgloc2 loc = (ref_ancestry_ids[d+1].x & 1) |
            ((ref_ancestry_ids[d+1].y & 1) << 1);

          return ancestor_order_rlocs[loc];
        }
      return VSG_RLOC2_MASK;
    }

  return 0x0;
}

/* static gint _visitors = 0; */
/* static gint _unused_visitors = 0; */

/*
 * operates a traversal of near/far interactions for a visiting node.
 */
static gboolean _compute_visiting_node (VsgPRTree2@t@ *tree,
                                        VsgNFConfig2@t@ *nfc,
                                        WaitingVisitor *wv)
{
  NIAndFuncs niaf;
  VsgPRTreeKey2@t@ ref_ancestry_ids[sizeof (@key_type@) * 8];
  gint i;

  niaf.ref = wv->node;
  niaf.far_func = nfc->far_func;
  niaf.near_func = nfc->near_func;
  niaf.user_data = nfc->user_data;

  _vsg_prtree2@t@node_get_info (wv->node, &niaf.ref_info, NULL, 0);
  memcpy (&niaf.ref_info.id, &wv->id, sizeof (VsgPRTreeKey2@t@));

  vsg_prtree_key2@t@_copy (&ref_ancestry_ids[niaf.ref_info.id.depth],
                           &niaf.ref_info.id);
  for (i = niaf.ref_info.id.depth-1; i >= 0; i --)
    {
      vsg_prtree_key2@t@_get_father (&ref_ancestry_ids[i+1],
                                     &ref_ancestry_ids[i]);
    }

  niaf.ref_ancestry_ids = ref_ancestry_ids;

  niaf.done_flag = 0;

  vsg_prtree2@t@_traverse_custom_internal (tree, G_POST_ORDER,
                                           (VsgRegion2@t@InternalLocDataFunc)
                                           _selector_nf_visitor,
                                           &niaf.ref_info, ref_ancestry_ids,
                                           (VsgPRTree2@t@InternalFunc)
                                           _traverse_visiting_nf,
                                           &niaf);

/*   if (niaf.done_flag == 0) _unused_visitors ++; */
/*   _visitors ++; */

/*   if (niaf.done_flag == 0) */
/*     { */
/*       g_printerr ("%d : reject node [", nfc->rk); */
/*       vsg_prtree_key2@t@_write (&wv->id, stderr); */
/*       g_printerr ("]\n"); */
/* /\*       _visiting_node_free (wv->node, &tree->config); *\/ */
/* /\*       _waiting_visitor_free (wv); *\/ */
/* /\*       g_slist_free_1 (first); *\/ */

/* /\*       return FALSE; *\/ */
/*     } */
/*   else */
/*     { */
/*       g_printerr ("%d : accepted node [", nfc->rk); */
/*       vsg_prtree_key2@t@_write (&wv->id, stderr); */
/*       g_printerr ("]\n"); */
     
/*     } */

  return TRUE;
}



typedef struct _NFSendData NFSendData;
struct _NFSendData {
  VsgPRTree2@t@ *tree;
  VsgNFConfig2@t@ *nfc;
};

/*
 * checks for a specific VsgNFPocMsg request completion and fills the
 * request with pending msgs.
 */
static void _fill_procs_msgs (gpointer key, VsgNFProcMsg *nfpm,
                              NFSendData *nfsd)
{
  gint flag;

  MPI_Test (&nfpm->request, &flag, MPI_STATUS_IGNORE);

  if (flag)
    {
      gint proc = GPOINTER_TO_INT (key);

      if (nfpm->backward_pending != NULL)
        {
          /* First, check for backward messages */
          _send_pending_backward_node (nfsd->tree, nfsd->nfc, nfpm, proc);
        }
      else if (nfpm->forward_pending != NULL)
        {
          /* Fallback into forward message */
          _send_pending_forward_node (nfsd->tree, nfsd->nfc, nfpm, proc);
        }
    }
}


static void vsg_prtree2@t@_nf_check_send (VsgPRTree2@t@ *tree,
                                          VsgNFConfig2@t@ *nfc)
{
  NFSendData nfsd = {tree, nfc};

  g_hash_table_foreach (nfc->procs_msgs, (GHFunc) _fill_procs_msgs, &nfsd);
}

/*
 * Probes for incoing messages. When called with @blocking == TRUE, will wait
 * until any receive happens. The time spent waiting will be used in
 * pending computations (of visitors) or in sending pending messages.
 */
gboolean vsg_prtree2@t@_nf_check_receive (VsgPRTree2@t@ *tree,
                                          VsgNFConfig2@t@ *nfc, gint tag,
                                          gboolean blocking)
{
  gint flag = FALSE;
  VsgPRTree2@t@Config *config = &tree->config;
  VsgPRTreeParallelConfig * pc = &config->parallel_config;
  MPI_Status status;
  gint rk;
  MPI_Comm_rank (MPI_COMM_WORLD, &rk);

  MPI_Iprobe (MPI_ANY_SOURCE, tag, pc->communicator, &flag, &status);

  if (blocking)
    {
      while (!flag)
        {
          if ((nfc->forward_pending_nb + nfc->backward_pending_nb) > 0)
            {
              vsg_prtree2@t@_nf_check_send (tree, nfc);
            }
          else
            {
              /* fallback asleep for just a moment before rechecking */
              g_usleep (1);
            }

          MPI_Iprobe (MPI_ANY_SOURCE, tag, pc->communicator, &flag, &status);
        }
    }

  if (flag)
    {
      VsgPRTreeKey2@t@ id;
      VsgPRTree2@t@Node *node;
      WaitingVisitor *wv;

/*       g_printerr ("%d : check_recv begin from %d tag=%d\n", nfc->rk, */
/*                   status.MPI_SOURCE, status.MPI_TAG); */

      vsg_packed_msg_recv (&nfc->recv, status.MPI_SOURCE, status.MPI_TAG);

/*       g_printerr ("%d : check_recv recv completed\n", nfc->rk); */

      switch (status.MPI_TAG) {
      case VISIT_FORWARD_TAG:
        vsg_packed_msg_recv_read (&nfc->recv, &id, 1,
                                  VSG_MPI_TYPE_PRTREE_KEY2@T@);

/*         g_printerr ("%d(%d) : fw recv from %d - ", nfc->rk, getpid (), */
/*                     status.MPI_SOURCE); */
/*         vsg_prtree_key2@t@_write (&id, stderr); */
/*         g_printerr ("\n"); */
/*         fflush (stderr); */

        node = _new_visiting_node (tree, &id, status.MPI_SOURCE);

        wv = _waiting_visitor_new (node, &id, status.MPI_SOURCE);

        _visit_unpack_node (config, node, nfc, DIRECTION_FORWARD);

        /* compute nf interactions with local tree */
        if (_compute_visiting_node (tree, nfc, wv))
          {
            _propose_node_backward (tree, nfc, status.MPI_SOURCE, wv);
          }
        else /* TODO: implement dropped visitors */
          {
            _propose_node_backward (tree, nfc, status.MPI_SOURCE, wv);
          }

        nfc->all_fw_recvs ++;

        break;
      case VISIT_BACKWARD_TAG:
        vsg_packed_msg_recv_read (&nfc->recv, &id, 1,
                                  VSG_MPI_TYPE_PRTREE_KEY2@T@);

        node = vsg_prtree2@t@node_key_lookup (tree->node, id);

/*         g_printerr ("%d(%d) : bw recv from %d - ", nfc->rk, getpid (), */
/*                     status.MPI_SOURCE); */
/*         vsg_prtree_key2@t@_write (&id, stderr); */
/*         g_printerr ("\n"); */
/*         fflush (stderr); */

        g_assert (PRTREE2@T@NODE_IS_LOCAL (node));

        _visit_unpack_node (config, node, nfc, DIRECTION_BACKWARD);

        nfc->pending_backward_msgs --;
/*             g_printerr ("bw recv(%d)\n", nfc->rk); */
/*         fflush (stderr); */

        nfc->all_bw_recvs ++;
        break;
      case END_FW_TAG:
        nfc->pending_end_forward --;
        break;
      default:
        g_error ("Unrecognized Message TAG %d.", status.MPI_TAG);
        break;

      }

/*       g_printerr ("%d : check_recv end %d %d %d\n", nfc->rk, */
/*                   nfc->pending_backward_msgs, */
/*                   nfc->forward_pending_nb, */
/*                   nfc->backward_pending_nb); */

    }
  return flag != FALSE;
}

typedef struct _NodeRemoteData NodeRemoteData;
struct _NodeRemoteData
{
  VsgPRTree2@t@ *tree;
  VsgNFConfig2@t@ *nfc;
  VsgPRTree2@t@Node *ref_node;
  VsgPRTree2@t@NodeInfo *ref_info;
  VsgPRTreeKey2@t@ *ref_ancestry_ids;
  gboolean *procs;
  gboolean sent;
};

static vsgrloc2 _selector_nf_remote (VsgPRTree2@t@NodeInfo *ref_info,
                                     VsgPRTree2@t@NodeInfo *node_info,
                                     VsgPRTreeKey2@t@ *ref_ancestry_ids)
{
  static const vsgrloc2 ancestor_order_rlocs[] = {
    VSG_RLOC2_MASK - (VSG_RLOC2_COMP (0)-1),
    VSG_RLOC2_MASK - (VSG_RLOC2_COMP (1)-1),
    VSG_RLOC2_MASK - (VSG_RLOC2_COMP (2)-1),
    VSG_RLOC2_MASK - (VSG_RLOC2_COMP (3)-1),
  };
  VsgPRTreeKey2@t@ *node_key = &node_info->id;
  gint d = node_key->depth;

  if (d >= ref_info->id.depth)
    {
      /* we can't skip non leaf nodes like in the following because
         all ref_info's subtree could be skipped (see parallel_check in
         vsgprrtee2@t@-extras.c:vsg_prtree2@t@node_near_far_traversal()) */
      /* if (!ref_info->isleaf) return 0x0; */

      return VSG_RLOC2_MASK;
    }

  if (vsg_prtree_key2@t@_is_neighbour (&ref_ancestry_ids[d], node_key))
    {
      if (VSG_PRTREE2@T@_NODE_INFO_IS_LOCAL (node_info)) return 0x0;

      if (vsg_prtree_key2@t@_equals (&ref_ancestry_ids[d], node_key))
        {
          vsgloc2 loc = (ref_ancestry_ids[d+1].x & 1) |
            ((ref_ancestry_ids[d+1].y & 1) << 1);

          return ancestor_order_rlocs[loc];
        }

      return VSG_RLOC2_MASK;
    }

  return 0x0;
}

static inline @key_type@ _key_coord_dist (@key_type@ ref, @key_type@ node,
                                          @key_type@ clamped, guint8 height)
{
  if (node == ref)
    return 0;

  if (node > ref)
    return (1<<(height)) - clamped;

  return clamped;
}

/*
 * used in vsg_prtree2@t@_traverse to register which processors a node is to
 * be sent to for a near/far interaction.
 */
static void _traverse_check_remote_neighbours (VsgPRTree2@t@Node *node,
                                               VsgPRTree2@t@NodeInfo *node_info,
                                               NodeRemoteData *data)
{
  if (PRTREE2@T@NODE_IS_REMOTE (node))
    {
      gint proc = PRTREE2@T@NODE_PROC (node);
      VsgPRTree2@t@NodeInfo *ref_info = data->ref_info;
      guint8 node_depth;

     if (data->procs[proc]) return;

/*       gint rk; */

/*       MPI_Comm_rank (MPI_COMM_WORLD, &rk); */

       node_depth = node_info->id.depth;

      if (node_depth <= ref_info->id.depth)
        {
          VsgPRTreeKey2@t@ *ref_id = &data->ref_ancestry_ids[node_depth];

          if (! vsg_prtree_key2@t@_is_neighbour (ref_id, &node_info->id))
            {
              if (node_depth < ref_info->id.depth)
                return;
            }
          else
            {
              if (node_depth < ref_info->id.depth &&
                  PRTREE2@T@NODE_LEAF (node).remote_depth > 0)
                {
                  /* check for mindepth */
                  VsgPRTreeKey2@t@ clamped_ref;
                  guint8 mindepth;
                  guint8 shift;
                  guint8 height;

                  /* use remote depth to avoid unwanted fw sends */
                  mindepth = node_info->depth +
                    PRTREE2@T@NODE_LEAF (node).remote_depth;

                  /* don't check deeper than node's depth */
                  mindepth = MIN (mindepth, data->ref_info->depth);

                  shift = ref_info->id.depth - mindepth;
                  height = mindepth - node_info->depth;

                  vsg_prtree_key2@t@_truncate (&ref_info->id, shift, &clamped_ref);

                  vsg_prtree_key2@t@_sever (&clamped_ref, height, &clamped_ref);

                  if (_key_coord_dist (ref_id->x, node_info->id.x,
                                       clamped_ref.x, height) > 3)
                    return;

                  if (_key_coord_dist (ref_id->y, node_info->id.y,
                                       clamped_ref.y, height) > 3)
                    return;
                }
            }
        }
      else
        if (ref_info->point_count == 0) return;

      _propose_node_forward (data->tree, data->nfc, proc,
                             data->ref_node, &data->ref_info->id);

      data->procs[proc] = TRUE;
      data->sent = TRUE;

      return;

    }
}

/*
 * checks wether some specified node is to be sent to distant processors in
 * order to complete a near/far traversal. Returns TRUE if node id
 * LOCAL and was sent to another processor.
 */
gboolean
vsg_prtree2@t@_node_check_parallel_near_far (VsgPRTree2@t@ *tree,
                                             VsgNFConfig2@t@ *nfc,
                                             VsgPRTree2@t@Node *node,
                                             VsgPRTree2@t@NodeInfo *info)
{
  gboolean ret = TRUE;
  gint rk, sz;

  if (tree->config.parallel_config.communicator == MPI_COMM_NULL)
    return ret;

  MPI_Comm_size (tree->config.parallel_config.communicator, &sz);
  if (sz < 2) return ret;

  MPI_Comm_rank (tree->config.parallel_config.communicator, &rk);

  vsg_prtree2@t@_nf_check_receive (tree, nfc, MPI_ANY_TAG, FALSE);

  if (VSG_PRTREE2@T@_NODE_INFO_IS_LOCAL (info))
    {
      gint i;
      NodeRemoteData nrd;
      VsgPRTreeKey2@t@ ref_ancestry_ids[sizeof (@key_type@) * 8];

      nrd.tree = tree;
      nrd.nfc = nfc;
      nrd.ref_node = node;
      nrd.ref_info = info;

      vsg_prtree_key2@t@_copy (&ref_ancestry_ids[info->id.depth], &info->id);
      for (i = info->id.depth-1; i >= 0; i --)
        {
          vsg_prtree_key2@t@_get_father (&ref_ancestry_ids[i+1],
                                         &ref_ancestry_ids[i]);
        }

      nrd.ref_ancestry_ids = ref_ancestry_ids;

      nrd.procs = g_alloca (nfc->sz * sizeof (gboolean));
      memset (nrd.procs, 0, nfc->sz * sizeof (gboolean));
      nrd.sent = FALSE;

      vsg_prtree2@t@_traverse_custom_internal (tree, G_PRE_ORDER,
                                               (VsgRegion2@t@InternalLocDataFunc)
                                               _selector_nf_remote,
                                               info, ref_ancestry_ids,
                                               (VsgPRTree2@t@InternalFunc)
                                               _traverse_check_remote_neighbours,
                                               &nrd);

      ret = nrd.sent;
    }

  vsg_prtree2@t@_nf_check_send (tree, nfc);

  return ret;
}

typedef struct _NodeDataReduce NodeDataReduce;

struct _NodeDataReduce {
  VsgPackedMsg *msg;
  VsgParallelMigrateVTable *data_vtable;
  gpointer tmp_node_data;
};

/*
 * used in a traversal to unpack all shared nodes of a tree from a VsgPackedMsg
 * in a backward direction (that is: supposedly accumulating/reducing with the
 * output part of the near/far interaction).
 */
static void _unpack_shared (VsgPRTree2@t@Node *node,
                            VsgPRTree2@t@NodeInfo *node_info,
                            NodeDataReduce *ndr)
{
  if (PRTREE2@T@NODE_IS_SHARED (node))
    {
/*       gint rk; */
/*       MPI_Comm_rank (MPI_COMM_WORLD, &rk); */

/*       g_printerr ("%d : unpacking shared ", rk); */
/*       vsg_prtree_key2@t@_write (&node_info->id, stderr); */
/*       g_printerr ("\n"); */
      ndr->data_vtable->unpack (ndr->tmp_node_data, ndr->msg,
                                ndr->data_vtable->unpack_data);

      ndr->data_vtable->reduce (ndr->tmp_node_data, node->user_data,
                                ndr->data_vtable->reduce_data);
    }
}

/*
 * used in a traversal to pack all shared nodes of a tree into a VsgPackedMsg
 * in a backward direction (that is: supposedly accumulating/reducing with the
 * output part of the near/far interaction).
 */
static void _pack_shared (VsgPRTree2@t@Node *node,
                          VsgPRTree2@t@NodeInfo *node_info,
                          NodeDataReduce *ndr)
{
  if (PRTREE2@T@NODE_IS_SHARED (node))
    {
/*       gint rk; */
/*       MPI_Comm_rank (MPI_COMM_WORLD, &rk); */

/*       g_printerr ("%d : packing shared ", rk); */
/*       vsg_prtree_key2@t@_write (&node_info->id, stderr); */
/*       g_printerr ("\n"); */
      ndr->data_vtable->pack (node->user_data, ndr->msg,
                              ndr->data_vtable->pack_data);

    }
}

static void
_shared_nodes_allreduce_internal (VsgPRTree2@t@ *tree,
                                  VsgParallelMigrateVTable *data_vtable,
                                  gpointer tmp_node_data)
{
  gint rk, sz;
  VsgPackedMsg send_msg =
    VSG_PACKED_MSG_STATIC_INIT (tree->config.parallel_config.communicator);
  VsgPackedMsg recv_msg =
    VSG_PACKED_MSG_STATIC_INIT (tree->config.parallel_config.communicator);

  NodeDataReduce ndr_send = {&send_msg, data_vtable, tmp_node_data};
  NodeDataReduce ndr_recv = {&recv_msg, data_vtable, tmp_node_data};
  MPI_Request requests[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
  gint step = 0;
  gint maxoffset, offset;
  gint quo;
  gint mod, div, src, dst;

  MPI_Comm_rank (tree->config.parallel_config.communicator, &rk);
  MPI_Comm_size (tree->config.parallel_config.communicator, &sz);

  while ((1<<step) < sz)
    step ++;

  maxoffset = 1<<(step-1);

/*   g_printerr ("%d : allreduce steps=%d maxoffset=%d\n", rk, step, */
/*               maxoffset); */

  while (step > 0)
    {
      step --;
      offset = 1 << step;
      quo = offset << 1;
      mod = (rk + offset) % quo;
      div = rk / quo;

      dst = mod + div * quo;
      src = dst;
      if (src >= sz) src = src % maxoffset;

      if (dst < sz && dst != rk)
        {
          /* reinit msg for packing a new message */
          send_msg.position = 0;

          /* fill ndr_send with shared nodes */
          vsg_prtree2@t@_traverse_custom_internal (tree, G_PRE_ORDER,
                                                   (VsgRegion2@t@InternalLocDataFunc)
                                                   _selector_skip_local_nodes,
                                                   NULL, NULL,
                                                   (VsgPRTree2@t@InternalFunc)
                                                   _pack_shared,
                                                   &ndr_send);

          /* send to first destination */
/*           g_printerr ("%d : allreduce sending (step %d) to %d\n", rk, */
/*                       step, dst); */
          vsg_packed_msg_isend (&send_msg, dst, VISIT_SHARED_TAG, &requests[0]);

          requests[1] = MPI_REQUEST_NULL;

          /* try alternate destination */
          dst += maxoffset;
          if (dst < sz && dst !=rk)
            {
              int dstmod = (dst + offset) % quo;
              int dstdiv = dst / quo;
              int dstsrc = dstmod + dstdiv * quo;

              /* send to alternate destination only if it has no
                 communiocation this step */
              if (dstsrc >= sz)
                {
/*                   g_printerr ("%d : allreduce sending2 (step %d) to %d\n", */
/*                               rk, step, dst); */
                  vsg_packed_msg_isend (&send_msg, dst, VISIT_SHARED_TAG,
                                        &requests[1]);
                }
            }
        }

      if (src != rk)
        {
          /* receive from source */
/*           g_printerr ("%d : allreduce recv (step %d) from %d\n", rk, */
/*                       step, src); */
          vsg_packed_msg_recv (&recv_msg, src, VISIT_SHARED_TAG);

          /* add results to shared nodes */
          vsg_prtree2@t@_traverse_custom_internal (tree, G_PRE_ORDER,
                                                   (VsgRegion2@t@InternalLocDataFunc)
                                                   _selector_skip_local_nodes,
                                                   NULL, NULL,
                                                   (VsgPRTree2@t@InternalFunc)
                                                   _unpack_shared,
                                                   &ndr_recv);
        }
/*       else */
/*           g_printerr ("%d : allreduce no recv\n", rk); */

      MPI_Waitall (2, requests, MPI_STATUSES_IGNORE);
/*       g_printerr ("%d : allreduce (step %d) ok\n", rk, step); */
    }

  vsg_packed_msg_drop_buffer (&send_msg);
  vsg_packed_msg_drop_buffer (&recv_msg);
}

/**
 * vsg_prtree2@t@_shared_nodes_allreduce:
 * @tree: A #VsgPRTree2@t@.
 * @data_vtable: A #VsgParallelMigrateVTable specifying som reduction operator.
 *
 * copies the semantic of the MPI_Allreduce operation on user_data for all
 * shared nodes of a tree. Reduction operator is given by the user
 * through @data_vtable.
 */
void
vsg_prtree2@t@_shared_nodes_allreduce (VsgPRTree2@t@ *tree,
                                       VsgParallelMigrateVTable *data_vtable)
{
  VsgPRTreeParallelConfig *pconfig = &tree->config.parallel_config;
  gpointer tmp_node_data = NULL;

  if (pconfig->node_data.alloc != NULL)
    {
      tmp_node_data =
        pconfig->node_data.alloc (FALSE, pconfig->node_data.alloc_data);
    }
  else if (tree->config.user_data_type != G_TYPE_NONE)
    {
      tmp_node_data = g_boxed_copy (tree->config.user_data_type,
                                    tree->config.user_data_model);
    }

  _shared_nodes_allreduce_internal (tree, data_vtable, tmp_node_data);

  if (pconfig->node_data.destroy != NULL)
    {
      pconfig->node_data.destroy (tmp_node_data, FALSE,
                                  pconfig->node_data.destroy_data);
    }
  else if (tree->config.user_data_type != G_TYPE_NONE)
    {
      g_boxed_free (tree->config.user_data_type, tmp_node_data);
    }
}

/*
 * waits for a specific VsgNFPocMsg request completion.
 */
static void _wait_procs_msgs (gpointer key, VsgNFProcMsg *nfpm, gpointer data)
{
  MPI_Wait (&nfpm->request, MPI_STATUS_IGNORE);
}

/*
 * finishes all communication and computations involved in a parallel near/far
 * interaction.
 */
void
vsg_prtree2@t@_nf_check_parallel_end (VsgPRTree2@t@ *tree,
                                      VsgNFConfig2@t@ *nfc)
{
  MPI_Comm comm = tree->config.parallel_config.communicator;
  gint i, dst;
  gint msg = 0;
  GTimer *timer = g_timer_new ();


/*   g_printerr ("%d(%d) : parallel_end begin (fw pending wv=%d) (bw pending wv=%d)\n", */
/*               nfc->rk, getpid (), */
/*               nfc->forward_pending_nb, */
/*               nfc->backward_pending_nb); */

  while (nfc->forward_pending_nb > 0)
    {
      vsg_prtree2@t@_nf_check_receive (tree, nfc, MPI_ANY_TAG, FALSE);
      vsg_prtree2@t@_nf_check_send (tree, nfc);
    }

  for (i=1; i<nfc->sz; i++)
    {
      dst = (nfc->rk+i) % nfc->sz;
      vsg_prtree2@t@_nf_check_send (tree, nfc);                                                                              
      MPI_Send (&msg, 0, MPI_INT, dst, END_FW_TAG, comm);
    }

/*   g_printerr ("%d(%d) : end fw sent\n", nfc->rk, getpid ()); */

  /* check all remaining messages */
  while (nfc->pending_end_forward > 0)
    {
/*     g_printerr ("%d : check %d\n", nfc->rk,  nfc->end_forward_received); */
      vsg_prtree2@t@_nf_check_send (tree, nfc);
      vsg_prtree2@t@_nf_check_receive (tree, nfc, MPI_ANY_TAG, TRUE);
    }

/*   g_printerr ("%d(%d) : end fw received\n", nfc->rk, getpid ()); */

  /* now, no forward visitor should be left incoming */

  /* do all remaining stuff */
  while ((nfc->forward_pending_nb + nfc->backward_pending_nb) > 0)
    {
/*       g_printerr ("%d(%d) : pending bw %d\n", */
/*                   nfc->rk, getpid (), nfc->pending_backward_msgs); */
      vsg_prtree2@t@_nf_check_send (tree, nfc);
      vsg_prtree2@t@_nf_check_receive (tree, nfc, MPI_ANY_TAG, FALSE);
    }

  while (nfc->pending_backward_msgs > 0)
    {
      vsg_prtree2@t@_nf_check_receive (tree, nfc, MPI_ANY_TAG, TRUE);
    }

/*   g_printerr ("%d(%d) : pending bw recv ok\n", nfc->rk, getpid ()); */

  g_hash_table_foreach (nfc->procs_msgs, (GHFunc) _wait_procs_msgs, NULL);

  MPI_Barrier (comm);

/*   g_printerr ("%d : allreduce begin\n", nfc->rk); */


  _shared_nodes_allreduce_internal (tree,
                                    &tree->config.parallel_config.node_data.visit_backward, nfc->tmp_node_data);

/*   g_printerr ("%d : parallel_end done (%f seconds)\n", nfc->rk, */
/*               g_timer_elapsed (timer, NULL)); */

/*   g_printerr ("%d : unused visitors %d (%d total)\n", nfc->rk, */
/*               _unused_visitors, _visitors); */

/*   g_printerr ("%d : nodes msgs (fw send=%d recv=%d) (bw send=%d recv=%d)\n", */
/*               nfc->rk, nfc->all_fw_sends, nfc->all_fw_recvs, nfc->all_bw_sends, */
/*               nfc->all_bw_recvs); */

/*   g_printerr ("%d : pending bw stats: max=%d avg=%f nb=%d\n", nfc->rk, */
/*               _bw_pending_max, */
/*               _bw_pending_sum / (gdouble) _bw_pending_calls, */
/*               _bw_pending_calls); */

  g_timer_destroy (timer);
}

static guint8 _prtree2@t@node_mindepth (const VsgPRTree2@t@Node *node)
{
  guint8 res = 8 * sizeof (@type@); /* max depth is the number of bits of
                                       @type@ */
  vsgloc2 i;

  /* empty nodes can be omitted here: return /infinite/ depth */
  if (node->point_count == 0)
    return res;

  if (PRTREE2@T@NODE_ISLEAF (node)) return 0;

  for (i=0; i<4; i++)
    {
      guint8 tmp = _prtree2@t@node_mindepth (PRTREE2@T@NODE_CHILD (node, i));
      if (tmp < res) res = tmp;
    }

  return (res == G_MAXUINT8) ? res : res + 1;
}

static void _remote_depths_array_build (VsgPRTree2@t@Node *node,
                                        VsgPRTree2@t@NodeInfo *node_info,
                                        GArray *array)
{
  if (! PRTREE2@T@NODE_IS_SHARED (node))
    {
      /* detect root of local (remote) subtrees */
      if (node_info->father_info == NULL ||
          VSG_PRTREE2@T@_NODE_INFO_IS_SHARED (node_info->father_info))
        {
          guint8 i = 0;

          if (PRTREE2@T@NODE_IS_LOCAL (node))
            i = _prtree2@t@node_mindepth (node);

          g_array_append_val (array, i);
        }
    }
}



static void _remote_depths_store (VsgPRTree2@t@Node *node,
                                  VsgPRTree2@t@NodeInfo *node_info,
                                  guint8 ** depths)
{
  if (! PRTREE2@T@NODE_IS_SHARED (node))
    {
      /* detect root of local (remote) subtrees */
      if (node_info->father_info == NULL ||
          VSG_PRTREE2@T@_NODE_INFO_IS_SHARED (node_info->father_info))
        {
          if (PRTREE2@T@NODE_IS_REMOTE (node))
            PRTREE2@T@NODE_LEAF (node).remote_depth = **depths;

          /* go to the next depth entry */
          (*depths) ++;
        }
    }
}

void vsg_prtree2@t@_update_remote_depths (VsgPRTree2@t@ *tree)
{
  GArray *array;
  GArray *reduced;
  guint8 *depths;

  g_return_if_fail (tree != NULL);

  array = g_array_sized_new (FALSE, FALSE, sizeof (guint8), 1024);
  vsg_prtree2@t@_traverse_custom_internal (tree, G_PRE_ORDER, NULL, NULL, NULL,
                                           (VsgPRTree2@t@InternalFunc)
                                           _remote_depths_array_build,
                                           array);

  /* prepare reduced for storing the results of the Allreduce */
  reduced = g_array_sized_new (FALSE, TRUE, sizeof (guint8), array->len);
  reduced = g_array_set_size (reduced, array->len);

  MPI_Allreduce (array->data, reduced->data, array->len, MPI_UNSIGNED_CHAR,
                 MPI_MAX, tree->config.parallel_config.communicator);

  g_array_free (array, TRUE);

  depths = (guint8 *) reduced->data;

  vsg_prtree2@t@_traverse_custom_internal (tree, G_PRE_ORDER, NULL, NULL, NULL,
                                           (VsgPRTree2@t@InternalFunc)
                                           _remote_depths_store,
                                           &depths);

  g_array_free (reduced, TRUE);

  tree->config.remote_depth_dirty = FALSE;
}

static void _local_leaves_count (VsgPRTree2@t@NodeInfo *node_info,
                                 GArray *array)
{
  if (! VSG_PRTREE2@T@_NODE_INFO_IS_SHARED (node_info))
    {
      /* detect root of local (remote) subtrees */
      if (node_info->father_info == NULL ||
          VSG_PRTREE2@T@_NODE_INFO_IS_SHARED (node_info->father_info))
        {
          gint i = 0;
          g_array_append_val (array, i);
        }

      /* increment local leaves */
      if (node_info->isleaf && VSG_PRTREE2@T@_NODE_INFO_IS_LOCAL (node_info) &&
          node_info->point_count != 0)
        g_array_index (array, gint, array->len-1) ++;
    }
}

typedef struct _ContiguousDistData ContiguousDistData;
struct _ContiguousDistData {
  GArray *array;       /* array of number of leaves on each local subtree */
  gint total_leaves;   /* total number of leaves in the tree (sum of array) */
  gint q, r, m;        /* distribution variables */
  gint current_index;  /* number of already checked array entries */
  gint current_lcount; /* number of already checked leaves */
  VsgPRTreeKey2@t@ last_subtree_id; /* key of the last local subtree */
};

static gint contiguous_dist (VsgPRTree2@t@NodeInfo *node_info,
                             ContiguousDistData *cda)
{
  if (VSG_PRTREE2@T@_NODE_INFO_IS_LOCAL (node_info))
    {
      gint ret = 0;

      if (cda->current_lcount >= cda->m)
        ret = (cda->current_lcount - cda->r) / cda->q;
      else
        ret = cda->current_lcount / (cda->q + 1);

      if (node_info->point_count == 0) return ret;

      cda->current_lcount ++;

      if (! vsg_prtree_key2@t@_is_ancestor (&cda->last_subtree_id,
                                          &node_info->id))
        {
          VsgPRTree2@t@NodeInfo *ancestor = node_info;

          while (ancestor->father_info &&
                 VSG_PRTREE2@T@_NODE_INFO_IS_LOCAL (ancestor->father_info))
            ancestor = ancestor->father_info;

          cda->current_index ++;
          vsg_prtree_key2@t@_copy (&cda->last_subtree_id, &ancestor->id);
        }

/*       g_printerr ("%d: sending node ", rk); */
/*       vsg_vector2@t@_write (&node_info->center, stderr); */
/*       g_printerr (" to %d\n", ret); */

      return ret;
    }

  g_assert (VSG_PRTREE2@T@_NODE_INFO_IS_REMOTE (node_info));

  cda->current_lcount += g_array_index (cda->array, gint, cda->current_index);
  cda->current_index ++;

  return -1;
}

static const VsgPRTreeKey2@t@ _dummy_key = {0, 0, 255};

/**
 * vsg_prtree2@t@_distribute_contiguous_leaves:
 * @tree: a #VsgPRTree2@t@.
 *
 * Performs parallel distribution of the nodes in @tree so as leaves
 * are sent to processors in order to have contiguous segments of (non
 * empty) leaves of roughly the same size in each
 * processor. "Contiguous segments" is here meant as if leaves were
 * numbered successively in an ordinary traversal of @tree. This means
 * that this distribution is subject to the current children ordering
 * of the tree.
 */
void vsg_prtree2@t@_distribute_contiguous_leaves (VsgPRTree2@t@ *tree)
{
  MPI_Comm comm;
  ContiguousDistData cda;
  GArray *array;
  GArray *reduced;
  gint i, sz;

  g_return_if_fail (tree != NULL);

  array = g_array_sized_new (FALSE, FALSE, sizeof (gint), 1024);
  comm = tree->config.parallel_config.communicator;
  
/*   MPI_Comm_rank (comm, &rk); */
  MPI_Comm_size (comm, &sz);

  /* accumulate number of local leaves in the local subtrees' roots */
  vsg_prtree2@t@_traverse (tree, G_PRE_ORDER,
                         (VsgPRTree2@t@Func) _local_leaves_count,
                         array);

  /* prepare reduced for storing the results of the Allreduce */
  reduced = g_array_sized_new (FALSE, TRUE, sizeof (gint), array->len);
  reduced = g_array_set_size (reduced, array->len);

  MPI_Allreduce (array->data, reduced->data, array->len, MPI_INT, MPI_MAX,
                 comm);

/*   g_printerr ("%d: contiguous allreduce size=%d\n", rk, array->len); */

  g_array_free (array, TRUE);

  cda.array = reduced;
  cda.total_leaves = 0;
  for (i=0; i<reduced->len; i++)
    cda.total_leaves += g_array_index (reduced, gint, i);

  if (cda.total_leaves > 0)
    {
      cda.q = cda.total_leaves / sz;
      cda.r = cda.total_leaves % sz;
      cda.m = (cda.q+1) * cda.r;
      cda.current_lcount = 0;
      cda.current_index = 0;
      cda.last_subtree_id = _dummy_key;

      vsg_prtree2@t@_distribute_nodes (tree,
                                     (VsgPRTree2@t@DistributionFunc) contiguous_dist,
                                     &cda);

/*   g_printerr ("%d : reduced-len=%d current-index=%d\n", rk, reduced->len, cda.current_index); */
    }

  g_array_free (reduced, TRUE);
}

typedef struct _ScatterData ScatterData;
struct _ScatterData {
  gint cptr;
  gint sz;
};

static gint scatter_dist (VsgPRTree2@t@NodeInfo *node_info, ScatterData *sd)
{
  if (VSG_PRTREE2@T@_NODE_INFO_IS_LOCAL (node_info))
    {
      gint dst = sd->cptr;

      sd->cptr ++;
      sd->cptr = sd->cptr % sd->sz;

      return dst;
    }

  return -1;
}

/**
 * vsg_prtree2@t@_distribute_scatter_leaves:
 * @tree: a #VsgPRTree2@t@.
 *
 * Performs parallel distribution of the nodes in @tree so as all
 * leaves will be scattered across all processors, in the order they
 * are encountered in a regular traversal (ie. first leaf for
 * processor 0, second leaf for processor 1, etc... and cycling when
 * the number of processors is reached).
 */
void vsg_prtree2@t@_distribute_scatter_leaves (VsgPRTree2@t@ *tree)
{
  MPI_Comm comm;
  ScatterData sd;

  g_return_if_fail (tree != NULL);

  comm = tree->config.parallel_config.communicator;

  sd.cptr = 0;
  MPI_Comm_size (comm, &sd.sz);

  vsg_prtree2@t@_distribute_nodes (tree,
                                   (VsgPRTree2@t@DistributionFunc)
                                   scatter_dist,
                                   &sd);
}

static gint concentrate_dist (VsgPRTree2@t@NodeInfo *node_info, gint *dst)
{
  return *dst;
}

/**
 * vsg_prtree2@t@_distribute_concentrate:
 * @tree: a #VsgPRTree2@t@.
 * @dst: destination processor number in the @tree's communicator.
 *
 * Performs parallel distribution of the nodes in @tree so as all
 * nodes will be sent to the processor @dst.
 */
void vsg_prtree2@t@_distribute_concentrate (VsgPRTree2@t@ *tree, gint dst)
{
  MPI_Comm comm;
  gint sz;

  g_return_if_fail (tree != NULL);

  comm = tree->config.parallel_config.communicator;
  MPI_Comm_size (comm, &sz);

  vsg_prtree2@t@_distribute_nodes (tree,
                                   (VsgPRTree2@t@DistributionFunc)
                                   concentrate_dist,
                                   &dst);

}
