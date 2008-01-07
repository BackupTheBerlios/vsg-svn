#include "vsgprtree2@t@-parallel.h"
#include "vsgprtree2@t@-private.h"

#include "string.h"

typedef struct _VTableAndMsg VTableAndMsg;

struct _VTableAndMsg
{
  VsgParallelVTable *vtable;
  VsgPackedMsg *msg;
};

static void _set_parallel_send_point (VsgPoint2 pt, VTableAndMsg *vm)
{
  vm->vtable->migrate_pack (pt, vm->msg, vm->vtable->migrate_pack_data);
  vm->vtable->destroy (pt, TRUE, vm->vtable->destroy_data);
}

static void _set_parallel_send_region (VsgRegion2 rg, VTableAndMsg *vm)
{
  vm->vtable->migrate_pack (rg, vm->msg, vm->vtable->migrate_pack_data);
  vm->vtable->destroy (rg, TRUE, vm->vtable->destroy_data);
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
      VsgPackedMsg pt_send = {comm, NULL, 0, 0};
      VsgPackedMsg rg_send = {comm, NULL, 0, 0};
      VTableAndMsg pt_vm = {&tree->config.parallel_config.point, &pt_send};
      VTableAndMsg rg_vm = {&tree->config.parallel_config.region, &rg_send};

      /* send points to 0 */
      cnt = vsg_prtree2@t@_point_count (tree);
      vsg_packed_msg_send_append (&pt_send, &cnt, 1, MPI_INT);
      vsg_prtree2@t@_foreach_point (tree, (GFunc) _set_parallel_send_point,
                                    &pt_vm);

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
      vsg_prtree2@t@_foreach_region (tree, (GFunc) _set_parallel_send_region,
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
      VsgPackedMsg pt_recv = {comm, NULL, 0, 0};
      VsgPackedMsg rg_recv = {comm, NULL, 0, 0};

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

