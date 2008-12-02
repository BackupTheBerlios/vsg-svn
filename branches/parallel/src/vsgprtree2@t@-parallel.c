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

      tree->node = vsg_prtree2@t@node_alloc (&bounds[0], &bounds[1],
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

  vsg_prtree2@t@_traverse_custom_internal (tree, G_PRE_ORDER, NULL, NULL, NULL,
                                           (VsgPRTree2@t@InternalFunc) _migrate_traverse_point_send,
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

      vsg_prtree2@t@_traverse_custom_internal (tree, G_PRE_ORDER, NULL, NULL, NULL,
                                               (VsgPRTree2@t@InternalFunc) _migrate_traverse_region_send,
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

static void _destroy_children (VsgPRTree2@t@Node *node,
                               const VsgPRTree2@t@Config *config)
{
  if (PRTREE2@T@NODE_ISINT (node))
    {
      gint i;

      /* destroy children */
      for (i=0; i<4; i++)
        {
          vsg_prtree2@t@node_free (PRTREE2@T@NODE_CHILD (node, i), config);
          PRTREE2@T@NODE_CHILD (node, i) = NULL;
        }
    }

}

static void _node_remove_regions (VsgPRTree2@t@Node *node,
                                  const VsgPRTree2@t@Config *config)
{
  g_slist_foreach (node->region_list, 
                   (GFunc) config->parallel_config.region.destroy,
                   config->parallel_config.region.destroy_data);
  g_slist_free (node->region_list);
  node->region_list = NULL;
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
      _destroy_children (node, dd->config);
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
              _node_remove_regions (node, config);
              _destroy_children (node, config);
            }

          node->parallel_status.storage = storage;
          node->parallel_status.proc = dst;
        }

      return;
    }

  if (storage == VSG_PARALLEL_REMOTE)
    {
      _node_remove_regions (node, config);
      _destroy_children (node, config);
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

static void _traverse_flatten_remote (VsgPRTree2@t@Node *node,
				      VsgPRTree2@t@NodeInfo *node_info,
				      const VsgPRTree2@t@Config *config)
{
  if (PRTREE2@T@NODE_IS_REMOTE (node))
    {
      /* destroy remaining children */
      _destroy_children (node, config);

      /* remote nodes aren't aware of points and regions stored on another
       * processor */
      _node_remove_regions (node, config);

      node->point_count = 0;
      node->region_count = 0;
    }
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
  vsg_prtree2@t@_traverse_custom_internal (tree, G_POST_ORDER, NULL, NULL, NULL,
					   (VsgPRTree2@t@InternalFunc)
					   _traverse_flatten_remote,
					   &tree->config);

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


static gboolean _compute_one_visit (VsgPRTree2@t@ *tree,
                                    VsgNFConfig2@t@ *nfc);
static gboolean _propose_backward_send (VsgPRTree2@t@ *tree,
                                        VsgNFConfig2@t@ *nfc);

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

  node->parallel_status.storage = VSG_PARALLEL_REMOTE;
  node->parallel_status.proc = src;

  return wv;
}

static void _waiting_visitor_free (WaitingVisitor *wv)
{
  g_slice_free (WaitingVisitor, wv);
}

/*
 * Holds a message and the request associated with a particular transfer (send
 * or receive).
 */
typedef struct _VsgNFProcMsg VsgNFProcMsg;
struct _VsgNFProcMsg {
  VsgPackedMsg send_pm;

  MPI_Request request;
};

static void vsg_nf_proc_msg_init (VsgNFProcMsg *nfpm, MPI_Comm comm)
{
  vsg_packed_msg_init (&nfpm->send_pm, comm);
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

   nfc->waiting_visitors = NULL;
   nfc->done_visitors = NULL;

   nfc->pending_end_forward = nfc->sz-1;
   nfc->pending_backward_msgs = 0;

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

  g_slist_free (nfc->waiting_visitors);
  g_slist_free (nfc->done_visitors);
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
                                              VsgPRTreeKey2@t@ *id)
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
                                VsgPackedMsg *recv,
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
  vsg_packed_msg_recv_read (recv, &datapresent, 1, MPI_BYTE);
  vsg_packed_msg_recv_read (recv, &point_count, 1, MPI_INT);
  vsg_packed_msg_recv_read (recv, &region_count, 1, MPI_INT);

  /* unpack user_data */
  if (datapresent)
    {
      /* node_data should be already allocated */
      node_data->unpack (node->user_data, recv, node_data->unpack_data);
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

              point->unpack (pt, recv, point->unpack_data);

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

              point->unpack (pt, recv, point->unpack_data);

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

              region->unpack (pt, recv, region->unpack_data);

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

              region->unpack (pt, recv, region->unpack_data);

              regions = g_slist_next (regions);
            }
        }
    }
}

/*
 * Operates the send operations needed to dispatch a node on a list of
 * different processors.
 * Waits that a processors former send operation if completed before issuing
 * this message. While waiting, receive or send backward (or even visitor
 * computations) operations can be performed.
 */
static void _propose_node_forward (VsgPRTree2@t@ *tree,
                                   VsgNFConfig2@t@ *nfc,
                                   gint nprocs, gint *procs,
                                   VsgPRTree2@t@Node *node,
                                   VsgPRTreeKey2@t@ *id)
{
  gint index;
  gint i, flag;
  MPI_Request requests[nprocs];
  VsgNFProcMsg *nfpms[nprocs];

  for (i=0; i<nprocs; i ++)
    {
      nfpms[i] = vsg_nf_config2@t@_proc_msgs_lookup (nfc, procs[i]);
      requests[i] = nfpms[i]->request;
    }

/*   g_printerr ("%d : propose_node_forward begin (nbdests=%d)\n", nfc->rk, nprocs); */

  while (nprocs > 0)
    {
      MPI_Testany (nprocs, requests, &index, &flag, MPI_STATUS_IGNORE);

      if (flag)
        {
          VsgPackedMsg *msg;

          /* avoid problem when all requests are MPI_REQUEST_NULL */
          if (index == MPI_UNDEFINED) index = 0;

          msg = &nfpms[index]->send_pm;

/*           g_printerr ("%d : sending fw to %d - ", nfc->rk, procs[index]); */
/*           vsg_prtree_key2@t@_write (id, stderr); */
/*           g_printerr ("\n"); */
/*           fflush (stderr); */

          msg->position = 0;

          vsg_packed_msg_send_append (msg, id, 1, VSG_MPI_TYPE_PRTREE_KEY2@T@);
          _visit_pack_node (node, msg, &tree->config, DIRECTION_FORWARD, FALSE);

          vsg_packed_msg_isend (msg, procs[index], VISIT_FORWARD_TAG,
                                &nfpms[index]->request);

          nfc->pending_backward_msgs ++;

          if (index < nprocs-1)
            {
              procs[index] = procs[nprocs-1];
              requests[index] = requests[nprocs-1];
              nfpms[index] = nfpms[nprocs-1];
            }

          nprocs --;
        }
      else
        {
          if (!vsg_prtree2@t@_nf_check_receive (tree, nfc, MPI_ANY_TAG, FALSE))
            {
              if (_propose_backward_send (tree, nfc))
                {
                  /*
                   * since we have sent a backward message, we need to
                   * refresh the list of actual requests.
                   */
                  for (i=0; i<nprocs; i ++)
                    requests[i] = nfpms[i]->request;
                }
              else
                {
                  if (! _compute_one_visit (tree, nfc))
                    g_usleep (1);
                }
            }
        }
    }

/*   g_printerr ("%d : propose_node_forward end\n", nfc->rk); */
}

/*
 * Stores visitor's node info and near/far functions for a
 * vsg_prtree2@t@_traverse operation.
 */
typedef struct _NIAndFuncs NIAndFuncs;
struct _NIAndFuncs {
  VsgPRTree2@t@NodeInfo ref_info;
  VsgPRTree2@t@FarInteractionFunc far_func;
  VsgPRTree2@t@InteractionFunc near_func;
  gpointer user_data;
  gint8 done_flag;
};

/*
 * computes near/far interactions for a visiting node during a tree traversal.
 */
static void _traverse_visiting_nf (VsgPRTree2@t@Node *node,
                                   VsgPRTree2@t@NodeInfo *node_info,
                                   NIAndFuncs *niaf)
{
  if (PRTREE2@T@NODE_IS_LOCAL (node))
    {
      VsgPRTree2@t@NodeInfo *ref_info = &niaf->ref_info;
/*       gint rk; */
      guint8 nf;
      gboolean fardone;

/*       MPI_Comm_rank (MPI_COMM_WORLD, &rk); */

      /* interact only with interior local nodes of depth greater than visiting
       * node.
       */
      if (PRTREE2@T@NODE_ISINT (node) &&
          (node_info->id.depth < ref_info->id.depth))
        return;

      /* avoid interactions with empty local nodes. */
      if (PRTREE2@T@NODE_ISLEAF (node) && node->point_count == 0) return;

      nf = vsg_prtree_key2@t@_compare_near_far (&ref_info->id,
                                                &node_info->id);

      if (nf > 1 && node_info->id.depth > ref_info->id.depth) return;

      switch (nf) {
      case (1):
        if (! node_info->isleaf) return;
        if (ref_info->point_count == 0) return;

/*         g_printerr ("%d : near interaction [", rk); */
/*         vsg_prtree_key2@t@_write (&ref_info->id, stderr); */
/*         g_printerr ("] ["); */
/*         vsg_prtree_key2@t@_write (&node_info->id, stderr); */
/*         g_printerr ("]\n"); */

        niaf->near_func (ref_info, node_info, niaf->user_data);
        niaf->done_flag |= 1;
        break;
      case (2):
        if (node_info->id.depth != ref_info->id.depth) return;

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
        break;
      default:
        break;
      };

    }

}

/*
 * operates a traversal of near/far interactions for any visiting node left
 * in @nfc.
 */
static gboolean _compute_one_visit (VsgPRTree2@t@ *tree,
                                    VsgNFConfig2@t@ *nfc)
{
  GSList *first = nfc->waiting_visitors;
  WaitingVisitor *wv;
  NIAndFuncs niaf;

  if (first == NULL) return FALSE;

  niaf.far_func = nfc->far_func;
  niaf.near_func = nfc->near_func;
  niaf.user_data = nfc->user_data;

  nfc->waiting_visitors = g_slist_next (nfc->waiting_visitors);
  wv = (WaitingVisitor *) first->data;

  _vsg_prtree2@t@node_get_info (wv->node, &niaf.ref_info, NULL, 0);
  memcpy (&niaf.ref_info.id, &wv->id, sizeof (VsgPRTreeKey2@t@));

  niaf.done_flag = 0;

  first->next = NULL;

  vsg_prtree2@t@_traverse_custom_internal (tree, G_POST_ORDER, NULL,
                                           NULL, NULL,
                                           (VsgPRTree2@t@InternalFunc)
                                           _traverse_visiting_nf,
                                           &niaf);

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

  nfc->done_visitors = g_slist_concat (nfc->done_visitors, first);

  return TRUE;
}

/*
 * Inits a backward send operation for the first pending visitor node which
 * processor is available (ie. all send operations are currently complete for
 * this processor).
 */
static gboolean _propose_backward_send (VsgPRTree2@t@ *tree,
                                        VsgNFConfig2@t@ *nfc)
{
  GSList *first = nfc->done_visitors;
  GSList *current;
  WaitingVisitor *wv;
  VsgNFProcMsg *nfpm;
  gboolean sent = FALSE;
  gint flag;

  if (first == NULL) return FALSE;

  

  current = first;

  while (current != NULL)
    {
      wv = (WaitingVisitor *) current->data;

      nfpm = vsg_nf_config2@t@_proc_msgs_lookup (nfc, wv->src);

      MPI_Test (&nfpm->request, &flag, MPI_STATUS_IGNORE);

      if (flag)
        {
          VsgPackedMsg *msg = &nfpm->send_pm;

/*           g_printerr ("%d : sending bw to %d - ", nfc->rk, wv->src); */
/*           vsg_prtree_key2@t@_write (&wv->id, stderr); */
/*           g_printerr ("\n"); */
/*           fflush (stderr); */

          msg->position = 0;

          vsg_packed_msg_send_append (msg, &wv->id, 1,
                                      VSG_MPI_TYPE_PRTREE_KEY2@T@);
          _visit_pack_node (wv->node, msg, &tree->config,
                            DIRECTION_BACKWARD, FALSE);

/*           g_printerr ("%d : propose_backward_send isend begin\n", nfc->rk); */

          vsg_packed_msg_isend (msg, wv->src, VISIT_BACKWARD_TAG,
                                &nfpm->request);

/*           g_printerr ("%d : propose_backward_send isend end\n", nfc->rk); */

          first = g_slist_remove_link (first, current);

          _waiting_visitor_free (wv);

          g_slist_free_1 (current);

          sent = TRUE;

          /* stop once we have send _one_ backward visitor to avoid
           * entanglement (because we can be called by a function that waits
           * for sending a message).
           */
          current = NULL;
        }
      else
        {
          current = g_slist_next (current);
        }
    }

  nfc->done_visitors = first;

  return sent;
}

/* #include <sys/types.h> */
/* #include <unistd.h> */

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
          /* try to send some visitor backward */
          if (! _propose_backward_send (tree, nfc))
            {
              /* try to wait with some visitor computation */
              if (! _compute_one_visit (tree, nfc))
                {
                  /* fallback asleep for just a moment before rechecking */
                  g_usleep (1);
                }
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

        node = _new_visiting_node (tree, &id);

        wv = _waiting_visitor_new (node, &id, status.MPI_SOURCE);

        _visit_unpack_node (config, node, &nfc->recv, DIRECTION_FORWARD);

        nfc->waiting_visitors = g_slist_append (nfc->waiting_visitors, wv);
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

        _visit_unpack_node (config, node, &nfc->recv, DIRECTION_BACKWARD);

        nfc->pending_backward_msgs --;
/*             g_printerr ("bw recv(%d)\n", nfc->rk); */
/*         fflush (stderr); */
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
/*                   g_slist_length (nfc->waiting_visitors), */
/*                   g_slist_length (nfc->done_visitors)); */

    }
  return flag != FALSE;
}

typedef struct _NodeRemoteData NodeRemoteData;
struct _NodeRemoteData
{
  VsgPRTree2@t@NodeInfo *ref_info;
  gboolean *procs;
  GArray *procnums;
};

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
      guint8 nf;

      if (data->procs[proc]) return;

      nf = vsg_prtree_key2@t@_compare_near_far (&data->ref_info->id,
                                                &node_info->id);

      if (nf < 3 && nf > 0)
        {
          g_array_append_val (data->procnums, proc);

          data->procs[proc] = TRUE;
        }
    }
}

/*
 * checks wether some specified node is to be sent to distant processors in
 * order to complete a near/far traversal.
 */
void
vsg_prtree2@t@_node_check_parallel_near_far (VsgPRTree2@t@ *tree,
                                             VsgNFConfig2@t@ *nfc,
                                             VsgPRTree2@t@Node *node,
                                             VsgPRTree2@t@NodeInfo *info)
{
  gint rk;

  if (tree->config.parallel_config.communicator == MPI_COMM_NULL)
    return;

  MPI_Comm_rank (tree->config.parallel_config.communicator, &rk);

  vsg_prtree2@t@_nf_check_receive (tree, nfc, MPI_ANY_TAG, FALSE);

  if (VSG_PRTREE2@T@_NODE_INFO_IS_LOCAL (info))
    {
      NodeRemoteData nrd;
      gint nb;
      gint *procs;

      nrd.ref_info = info;
      nrd.procs = g_alloca (nfc->sz * sizeof (gboolean));
      memset (nrd.procs, 0, nfc->sz * sizeof (gboolean));
      nrd.procnums = g_array_new (FALSE, FALSE, sizeof (gint));

      vsg_prtree2@t@_traverse_custom_internal (tree, G_PRE_ORDER, NULL, NULL,
                                               NULL,
                                               (VsgPRTree2@t@InternalFunc)
                                               _traverse_check_remote_neighbours,
                                               &nrd);

      nb = nrd.procnums->len;
      procs = (gint *) g_array_free (nrd.procnums, FALSE);

      _propose_node_forward (tree, nfc, nb, procs, node, &info->id);

      g_free (procs);
    }
}

typedef struct _ConfigAndMsg ConfigAndMsg;

struct _ConfigAndMsg {
  VsgPRTree2@t@Config *config;
  VsgPackedMsg *msg;
};

/*
 * used in a traversal to unpack all shared nodes of a tree from a VsgPackedMsg
 * in a backward direction (that is: supposedly accumulating/reducing with the
 * output part of the near/far interaction).
 */
static void _bw_unpack_shared (VsgPRTree2@t@Node *node,
                               VsgPRTree2@t@NodeInfo *node_info,
                               ConfigAndMsg *cam)
{
  if (PRTREE2@T@NODE_IS_SHARED (node))
    {
/*       gint rk; */
/*       MPI_Comm_rank (MPI_COMM_WORLD, &rk); */

/*       g_printerr ("%d : unpacking shared ", rk); */
/*       vsg_prtree_key2@t@_write (&node_info->id, stderr); */
/*       g_printerr ("\n"); */

      _visit_unpack_node (cam->config, node, cam->msg, DIRECTION_BACKWARD); 
    }
}

/*
 * used in a traversal to pack all shared nodes of a tree into a VsgPackedMsg
 * in a backward direction (that is: supposedly accumulating/reducing with the
 * output part of the near/far interaction).
 */
static void _bw_pack_shared (VsgPRTree2@t@Node *node,
                             VsgPRTree2@t@NodeInfo *node_info,
                             ConfigAndMsg *cam)
{
  if (PRTREE2@T@NODE_IS_SHARED (node))
    {
/*       gint rk; */
/*       MPI_Comm_rank (MPI_COMM_WORLD, &rk); */

/*       g_printerr ("%d : packing shared ", rk); */
/*       vsg_prtree_key2@t@_write (&node_info->id, stderr); */
/*       g_printerr ("\n"); */

      _visit_pack_node (node, cam->msg, cam->config, DIRECTION_BACKWARD, TRUE); 
    }
}

/*
 * copies the semantic of the MPI_Allreduce operation but for all shared nodes
 * of a tree involved in a near/far interaction.
 */
static void _nf_allreduce_shared (VsgPRTree2@t@ *tree,
                                  VsgNFConfig2@t@ *nfc)
{
  VsgPackedMsg msg =
    VSG_PACKED_MSG_STATIC_INIT (tree->config.parallel_config.communicator);
  ConfigAndMsg cam_send = {&tree->config, &msg};
  ConfigAndMsg cam_recv = {&tree->config, &nfc->recv};
  MPI_Request requests[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
  gint step = 0;
  gint maxoffset, offset;
  gint quo;
  gint mod, div, src, dst;

  while ((1<<step) < nfc->sz)
    step ++;

  maxoffset = 1<<(step-1);

/*   g_printerr ("%d : allreduce steps=%d maxoffset=%d\n", nfc->rk, step, */
/*               maxoffset); */

  while (step > 0)
    {
      step --;
      offset = 1 << step;
      quo = offset << 1;
      mod = (nfc->rk + offset) % quo;
      div = nfc->rk / quo;

      dst = mod + div * quo;
      src = dst;
      if (src >= nfc->sz) src = src % maxoffset;

      if (dst < nfc->sz && dst != nfc->rk)
        {
          /* reinit msg for packing a new message */
          msg.position = 0;

          /* fill cam_send with shared nodes */
          vsg_prtree2@t@_traverse_custom_internal (tree, G_PRE_ORDER, NULL,
                                                   NULL, NULL,
                                                   (VsgPRTree2@t@InternalFunc)
                                                   _bw_pack_shared,
                                                   &cam_send);

          /* send to first destination */
/*           g_printerr ("%d : allreduce sending to %d\n", nfc->rk, dst); */
          vsg_packed_msg_isend (&msg, dst, VISIT_SHARED_TAG, &requests[0]);

          requests[1] = MPI_REQUEST_NULL;
          /* try alternate destination */
          dst += maxoffset;
          if (dst < nfc->sz && dst != nfc->rk)
            {
/*               g_printerr ("%d : allreduce sending2 to %d\n", nfc->rk, dst); */
              vsg_packed_msg_isend (&msg, dst, VISIT_SHARED_TAG, &requests[1]);
            }
        }

      if (src != nfc->rk)
        {
          /* receive from source */
/*           g_printerr ("%d : allreduce recv from %d\n", nfc->rk, src); */
          vsg_packed_msg_recv (&nfc->recv, src, VISIT_SHARED_TAG);

          /* add results to shared nodes */
          vsg_prtree2@t@_traverse_custom_internal (tree, G_PRE_ORDER, NULL,
                                                   NULL, NULL,
                                                   (VsgPRTree2@t@InternalFunc)
                                                   _bw_unpack_shared,
                                                   &cam_recv);
        }
/*       else */
/*           g_printerr ("%d : allreduce no recv\n", nfc->rk); */

      MPI_Waitall (2, requests, MPI_STATUSES_IGNORE);
    }

  vsg_packed_msg_drop_buffer (&msg);
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


/*   g_printerr ("%d(%d) : parallel_end begin\n", nfc->rk, getpid ()); */

  for (i=1; i<nfc->sz; i++)
    {
      dst = (nfc->rk+i) % nfc->sz;
      MPI_Send (&msg, 0, MPI_INT, dst, END_FW_TAG, comm);
    }

/*   g_printerr ("%d(%d) : end fw sent\n", nfc->rk, getpid ()); */

  /* check all remaining messages */
  while (vsg_prtree2@t@_nf_check_receive (tree, nfc, MPI_ANY_TAG, FALSE) ||
         nfc->pending_end_forward > 0)
/*     g_printerr ("%d : check %d\n", nfc->rk,  nfc->end_forward_received); */
    /* pass */;
/*   g_printerr ("%d(%d) : end fw received\n", nfc->rk, getpid ()); */

  /* now, no forward visitor should be left incoming */

  /* do all remaining stuff */
  while (nfc->pending_backward_msgs > 0)
    {
/*       g_printerr ("%d(%d) : pending bw %d\n", */
/*                   nfc->rk, getpid (), nfc->pending_backward_msgs); */
      vsg_prtree2@t@_nf_check_receive (tree, nfc, MPI_ANY_TAG, TRUE);
    }

/*   g_printerr ("%d(%d) : pending bw recv ok\n", nfc->rk, getpid ()); */

  /* all receives have occured, we now only have pending sends */
  while (nfc->done_visitors != NULL ||
         nfc->waiting_visitors != NULL)
    {
/*       g_printerr ("%d(%d) : remaining computes=%d send=%d\n", nfc->rk, getpid (), */
/*                   g_slist_length (nfc->done_visitors), */
/*                   g_slist_length (nfc->waiting_visitors)); */
      _compute_one_visit (tree, nfc);
      _propose_backward_send (tree, nfc);
    }

  g_hash_table_foreach (nfc->procs_msgs, (GHFunc) _wait_procs_msgs, NULL);

  MPI_Barrier (comm);
/*   g_printerr ("%d : allreduce begin\n", nfc->rk); */


    _nf_allreduce_shared (tree, nfc);

  g_printerr ("%d : parallel_end done (%f seconds)\n", nfc->rk,
              g_timer_elapsed (timer, NULL));
  g_timer_destroy (timer);
}
