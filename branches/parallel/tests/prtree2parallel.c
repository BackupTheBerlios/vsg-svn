/* basic point insertion, removal and traversal on a cloned VsgPRTree2d */

#include "vsg-config.h"

#include "vsg/vsgd.h"

#include <math.h>
#include <glib.h>

static gint mytid, numtasks;

static GPtrArray *points = NULL;

static gpointer pt_alloc (gboolean resident, gpointer user_data)
{
  gpointer ret;
  ret = g_malloc (sizeof (VsgVector2d));
  g_ptr_array_add (points, ret);
/*   g_printerr ("%d: alloc 1 VsgVector2d (%p)\n", mytid, ret); */
  return ret;
}

static void pt_destroy (gpointer data, gboolean resident,
                        gpointer user_data)
{
/*   g_printerr ("%d: destroy 1 VsgVector2d (%p)\n", mytid, data); */
  g_ptr_array_remove (points, data);
  g_free (data);
}

static void pt_migrate_pack (VsgVector2d *pt, VsgPackedMsg *pm,
                             gpointer user_data)
{
/*   g_printerr ("%d: pack ", mytid); */
/*   vsg_vector2d_write (pt, stderr); */
/*   g_printerr ("\n"); */
  vsg_packed_msg_send_append (pm, pt, 1, VSG_MPI_TYPE_VECTOR2D);
}

static void pt_migrate_unpack (VsgVector2d *pt, VsgPackedMsg *pm,
                               gpointer user_data)
{
  vsg_packed_msg_recv_read (pm, pt, 1, VSG_MPI_TYPE_VECTOR2D);

/*   g_printerr ("%d: unpack ", mytid); */
/*   vsg_vector2d_write (pt, stderr); */
/*   g_printerr ("\n"); */
}

static VsgPRTreeParallelConfig pconfig = {
  MPI_COMM_WORLD,
  {pt_alloc, NULL,
   pt_destroy, NULL,
   (VsgMigrablePackDataFunc) pt_migrate_pack, NULL,
   (VsgMigrablePackDataFunc) pt_migrate_unpack, NULL,
  },
};

static void empty_array (gpointer var, gpointer data)
{
/*   g_printerr ("%d: static destroy 1 VsgVector2d (%p)\n", mytid, var); */
  g_free (var);
}

gint main (gint argc, gchar ** argv)
{
  gint ret = 0;

  VsgPRTree2d *tree;
  gint i;

  VsgVector2d lb;
  VsgVector2d ub;


  MPI_Init (&argc, &argv);

  MPI_Comm_size (MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank (MPI_COMM_WORLD, &mytid);

  if (argc > 1 && g_strncasecmp (argv[1], "--version", 9) == 0)
    {
      if (mytid == 0)
        g_print ("%s\n", PACKAGE_VERSION);
      return 0;
    }

  vsg_init_gdouble ();

  points = g_ptr_array_new ();

  if (mytid == 0)
    {
      VsgVector2d *pt;
      lb.x = -1.; lb.y = -1.;
      ub.x = 0.; ub.y = 0.;

      pt = pt_alloc (TRUE, NULL);
      pt->x = -0.5; pt->y = -0.5;
    }
  else
    {
      VsgVector2d *pt;

      lb.x = 0.; lb.y = 0.;
      ub.x = 1.*mytid; ub.y = 1.*mytid;

      pt = pt_alloc (TRUE, NULL);
      pt->x = 0.5*mytid; pt->y = 0.5*mytid;

      pt = pt_alloc (TRUE, NULL);
      pt->x = 0.65*mytid; pt->y = 0.65*mytid;

      pt = pt_alloc (TRUE, NULL);
      pt->x = 0.75*mytid; pt->y = 0.75*mytid;
    }

  /* create the tree */
  tree =
    vsg_prtree2d_new_full (&lb, &ub,
                           (VsgPoint2dLocFunc) vsg_vector2d_vector2d_locfunc,
                           (VsgPoint2dDistFunc) vsg_vector2d_dist,
                           NULL, 2);

  /* insert some points */
  for (i=0; i<points->len; i++)
    {
      vsg_prtree2d_insert_point (tree, g_ptr_array_index (points, i));
    }

  MPI_Barrier (MPI_COMM_WORLD);

/*   g_printerr ("%d: set_parallel begin\n", mytid); */

  vsg_prtree2d_set_parallel (tree, &pconfig);

/*   MPI_Barrier (MPI_COMM_WORLD); */
/*   g_printerr ("%d: set_parallel ok\n", mytid); */

/*   if (mytid == 0) */
/*     vsg_prtree2d_write (tree, stderr); */

  /* destroy the points */
  g_ptr_array_foreach (points, empty_array, NULL);
  g_ptr_array_free (points, TRUE);

  /* destroy the cloned tree */
  vsg_prtree2d_free (tree);

  MPI_Finalize ();

  return ret;
}
