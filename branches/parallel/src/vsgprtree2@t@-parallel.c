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
  vm->vtable->migrate_pack (pt, vm->msg, vm->vtable->migrate_pack_data);
  vm->vtable->destroy (pt, TRUE, vm->vtable->destroy_data);
}

static void _send_region (VsgRegion2 rg, VTableAndMsg *vm)
{
  vm->vtable->migrate_pack (rg, vm->msg, vm->vtable->migrate_pack_data);
  vm->vtable->destroy (rg, TRUE, vm->vtable->destroy_data);
}

static void _send_shared_region (VsgRegion2 rg, VTableAndMsg *vm)
{
  vm->vtable->migrate_pack (rg, vm->msg, vm->vtable->migrate_pack_data);
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
            sizeof (VsgPRTreeParallelConfig));;

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

              pt_vtable->migrate_unpack (pt, &pt_recv,
                                         pt_vtable->migrate_unpack_data);
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

              rg_vtable->migrate_unpack (rg, &rg_recv,
                                         rg_vtable->migrate_unpack_data);
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

  vsg_prtree2@t@_traverse_custom_internal (tree, G_PRE_ORDER, NULL,
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

          pt_vtable->migrate_unpack (pt, recv, pt_vtable->migrate_unpack_data);
          vsg_prtree2@t@_insert_point (tree, pt);
          
        }

      vsg_comm_buffer_drop_recv_buffer (cb, src);
    }

  if (rg_vtable->migrate_pack != NULL)
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

      vsg_prtree2@t@_traverse_custom_internal (tree, G_PRE_ORDER, NULL,
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

              rg_vtable->migrate_unpack (rg, recv,
                                         rg_vtable->migrate_unpack_data);
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

static void _migrate_pack_node_header (VsgPRTree2@t@Node *node,
                                       VsgPackedMsg *msg, gint dst,
                                       VsgPRTree2@t@Config *config)
{
  vsg_packed_msg_send_append (msg, &node->center, 1, VSG_MPI_TYPE_VECTOR2@T@);
  vsg_packed_msg_send_append (msg, &dst, 1, MPI_INT);
}

static void _migrate_pack_node (VsgPRTree2@t@Node *node, VsgPackedMsg *msg,
                                VsgPRTree2@t@Config *config,
                                gboolean shared)
{
  VsgPRTreeParallelConfig *pc = &config->parallel_config;
  gint datapresent, point_count, region_count;

  datapresent = pc->node_data.migrate_pack != NULL && node->user_data != NULL;
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
      pc->node_data.migrate_pack (node->user_data, msg,
                                  pc->node_data.migrate_pack_data);

      if (pc->node_data.destroy)
        pc->node_data.destroy (node->user_data, TRUE,
                               pc->node_data.destroy_data);
      else
        g_boxed_free (config->user_data_type, node->user_data);

      node->user_data = NULL;
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

  g_assert (PRTREE2@T@NODE_ISINT (node));

  if (PRTREE2@T@NODE_IS_SHARED (PRTREE2@T@NODE_INT (node).children[0]))
    children_proc = -1;
  else
    children_proc =
      PRTREE2@T@NODE_PROC (PRTREE2@T@NODE_INT (node).children[0]);

  for (i=1; i<4; i++)
    {
      gint proc =
        PRTREE2@T@NODE_PROC (PRTREE2@T@NODE_INT (node).children[i]);

      if (PRTREE2@T@NODE_IS_SHARED (PRTREE2@T@NODE_INT (node).children[i]) ||
          proc != children_proc)
        children_proc = -1;
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
          _migrate_pack_node_header (node, dd->bcast, dst, dd->config);

          /* choose between ptp communication or bcast */
          if (new_storage == VSG_PARALLEL_REMOTE)
            msg = vsg_comm_buffer_get_send_buffer (dd->cb, dst);
          else
            msg = dd->bcast;

          /* write node to the corresponding message */
          _migrate_pack_node (node, msg, dd->config, dst < 0);
        }

    }

  if (old_storage == VSG_PARALLEL_SHARED && new_storage == VSG_PARALLEL_REMOTE)
    {

      /* free all shared regions */
      g_slist_foreach (node->region_list, 
                       (GFunc) dd->config->parallel_config.region.destroy,
                       dd->config->parallel_config.region.destroy_data);
      g_slist_free (node->region_list);
      node->region_list = NULL;

      _destroy_children (node, dd->config);
    }

  if (old_storage == VSG_PARALLEL_REMOTE) return; /* check migrations later */

  /* update node's parallel status */
  node->parallel_status.storage = new_storage;
  node->parallel_status.proc = dst;
}

void vsg_prtree2@t@node_insert_child (VsgPRTree2@t@Node *node,
                                      const VsgPRTree2@t@Config *config,
                                      VsgParallelStorage storage,
                                      gint dst,
                                      VsgVector2@t@ *center,
                                      gpointer user_data,
                                      GSList *points,
                                      GSList *regions)
{
  const VsgPRTreeParallelConfig *pc = &config->parallel_config;

  if (vsg_vector2@t@_dist (&node->center, center) > @epsilon@)
    {
      vsgloc2 child = CALL_POINT2@T@_LOC (config, center, &node->center);
      gint children_proc;

      if (PRTREE2@T@NODE_ISLEAF (node))
        {
          vsg_prtree2@t@node_make_int (node, config);
        }

      vsg_prtree2@t@node_insert_child (PRTREE2@T@NODE_CHILD (node, child),
                                       config, storage, dst, center, user_data,
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

          /* ? */
          if (storage == VSG_PARALLEL_REMOTE)
            {
/*               gint rk; */

/*               MPI_Comm_rank (pc->communicator, &rk); */
/*               g_printerr ("%d: destroy children of ", rk); */
/*               vsg_vector2@t@_write (&node->center, stderr); */
/*               g_printerr ("\n"); */

              _destroy_children (node, config);
            }

          node->parallel_status.storage = storage;
          node->parallel_status.proc = dst;
        }

      return;
    }

  /* ? */
  if (storage == VSG_PARALLEL_REMOTE) _destroy_children (node, config);

  if (user_data != NULL)
    {
      if (pc->node_data.destroy)
        pc->node_data.destroy (node->user_data, TRUE,
                               pc->node_data.destroy_data);
      else
        g_boxed_free (config->user_data_type, node->user_data);

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
  vsg_prtree2@t@_traverse_custom_internal
    (tree, G_POST_ORDER, NULL,
     (VsgPRTree2@t@InternalFunc) _traverse_distribute_nodes, &dd);

/*   g_printerr ("%d: after gather\n", rk); */

  /* send/receive pending messages */
  vsg_comm_buffer_set_bcast (bcastcb, &bcast);

/*   g_printerr ("%d: before share\n", rk); */

  vsg_comm_buffer_share (bcastcb);
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
          VsgVector2@t@ center;
          gint dst;

          /* unpack the header */
          vsg_packed_msg_recv_read (hdrmsg, &center, 1,
                                    VSG_MPI_TYPE_VECTOR2@T@);
          vsg_packed_msg_recv_read (hdrmsg, &dst, 1, MPI_INT);

/*           g_printerr ("%d: unpacking src=%d center=", rk, src); */
/*           vsg_vector2@t@_write (&center, stderr);  */
/*           g_printerr (" dst=%d\n", dst); */

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

                  pc->node_data.migrate_unpack (data, unpack,
                                                pc->node_data.migrate_unpack_data);
                }

              g_assert (point_count == 0 || storage == VSG_PARALLEL_LOCAL);

              for (i=0; i<point_count; i++)
                {
                  VsgPoint2 pt =
                    pc->point.alloc (TRUE, pc->point.alloc_data);

                  pc->point.migrate_unpack (pt, unpack,
                                            pc->point.migrate_unpack_data);

                  points = g_slist_prepend (points, pt);
                }

              for (i=0; i<region_count; i++)
                {
                  VsgRegion2 pt =
                    pc->region.alloc (TRUE, pc->region.alloc_data);

                  pc->region.migrate_unpack (pt, unpack,
                                             pc->region.migrate_unpack_data);

                  regions = g_slist_prepend (regions, pt);
                }

              /* insert the node in the tree */
              vsg_prtree2@t@node_insert_child (tree->node, &tree->config,
                                               storage, dst, &center, data,
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
                                               &center, NULL, NULL, NULL);
              
            }

/*           g_printerr ("%d: unpacking done\n", rk); */
        }
    }

  vsg_comm_buffer_free (bcastcb);
  vsg_comm_buffer_free (cb);
}
