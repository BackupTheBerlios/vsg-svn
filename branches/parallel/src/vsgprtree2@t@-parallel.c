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
      VsgPackedMsg pt_send = VSG_PACKED_MSG_NULL;
      VsgPackedMsg rg_send = VSG_PACKED_MSG_NULL;
      VTableAndMsg pt_vm = {&tree->config.parallel_config.point, &pt_send};
      VTableAndMsg rg_vm = {&tree->config.parallel_config.region, &rg_send};

      vsg_packed_msg_init (&pt_send, comm);
      vsg_packed_msg_init (&rg_send, comm);

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
      VsgPackedMsg pt_recv = VSG_PACKED_MSG_NULL;
      VsgPackedMsg rg_recv = VSG_PACKED_MSG_NULL;

      vsg_packed_msg_init (&pt_recv, comm);
      vsg_packed_msg_init (&rg_recv, comm);

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

static void _print_region (gpointer rg, FILE *file)
{
  fprintf (file, "%p ", rg);
}

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

