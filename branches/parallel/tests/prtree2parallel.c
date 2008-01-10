/* basic point insertion, removal and traversal on a cloned VsgPRTree2d */

#include "vsg-config.h"

#include "vsg/vsgd.h"

#include "glib/gprintf.h"

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

static void empty_array (gpointer var, gpointer data)
{
/*   g_printerr ("%d: static destroy 1 VsgVector2d (%p)\n", mytid, var); */
  g_free (var);
}

static GPtrArray *regions = NULL;

typedef struct _Circle Circle;

struct _Circle {
  VsgVector2d center;
  gdouble radius;
};

static vsgrloc2 _circle_loc2 (Circle *circle, VsgVector2d *center)
{
  VsgVector2d tmp;
  gdouble dist;
  vsgloc2 centerpos;
  vsgrloc2 ret;

  vsg_vector2d_sub (center, &circle->center, &tmp);

  dist = vsg_vector2d_norm (&tmp);

  if (dist <= circle->radius) return VSG_RLOC2_MASK;

  ret = 0;

  centerpos = vsg_vector2d_vector2d_locfunc (&circle->center, center);

  ret |= VSG_RLOC2_COMP (centerpos);

  if (fabs (tmp.x) <= circle->radius)
    {
      vsgloc2 itmp = centerpos & ~ (VSG_LOC2_X & centerpos);
      ret |= VSG_RLOC2_COMP (itmp);
    }

  if (fabs (tmp.y) <= circle->radius)
    {
      vsgloc2 itmp = centerpos & ~ (VSG_LOC2_Y & centerpos);
      ret |= VSG_RLOC2_COMP (itmp);
    }

  return ret;
}

static gpointer rg_alloc (gboolean resident, gpointer user_data)
{
  gpointer ret;
  ret = g_malloc (sizeof (Circle));
  g_ptr_array_add (regions, ret);
/*   g_printerr ("%d: alloc 1 Circle (%p)\n", mytid, ret); */
  return ret;
}

static void rg_destroy (gpointer data, gboolean resident,
                        gpointer user_data)
{
/*   g_printerr ("%d: destroy 1 Circle (%p)\n", mytid, data); */
  g_ptr_array_remove (regions, data);
  g_free (data);
}

static void rg_migrate_pack (Circle *rg, VsgPackedMsg *pm,
                             gpointer user_data)
{
/*   g_printerr ("%d: pack circle {c=", mytid); */
/*   vsg_vector2d_write (&rg->center, stderr); */
/*   g_printerr (" r=%lf}\n", rg->radius); */
  vsg_packed_msg_send_append (pm, &rg->center, 1, VSG_MPI_TYPE_VECTOR2D);
  vsg_packed_msg_send_append (pm, &rg->radius, 1, MPI_DOUBLE);
}

static void rg_migrate_unpack (Circle *rg, VsgPackedMsg *pm,
                               gpointer user_data)
{
  vsg_packed_msg_recv_read (pm, &rg->center, 1, VSG_MPI_TYPE_VECTOR2D);
  vsg_packed_msg_recv_read (pm, &rg->radius, 1, MPI_DOUBLE);

/*   g_printerr ("%d: unpack circle {c=", mytid); */
/*   vsg_vector2d_write (&rg->center, stderr); */
/*   g_printerr (" r=%lf}\n", rg->radius); */
}

static VsgPRTreeParallelConfig pconfig = {
  MPI_COMM_WORLD,
  {pt_alloc, NULL,
   pt_destroy, NULL,
   (VsgMigrablePackDataFunc) pt_migrate_pack, NULL,
   (VsgMigrablePackDataFunc) pt_migrate_unpack, NULL,
  },
  {rg_alloc, NULL,
   rg_destroy, NULL,
   (VsgMigrablePackDataFunc) rg_migrate_pack, NULL,
   (VsgMigrablePackDataFunc) rg_migrate_unpack, NULL,
  },
};


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
  regions = g_ptr_array_new ();

  if (mytid == 0)
    {
      VsgVector2d *pt;
      Circle *c;

      lb.x = -1.; lb.y = -1.;
      ub.x = 0.; ub.y = 0.;

      pt = pt_alloc (TRUE, NULL);
      pt->x = -0.5; pt->y = -0.5;

      c = rg_alloc (TRUE, NULL);
      c->center.x = -0.6; c->center.y = -0.6;
      c->radius = 0.1;
    }
  else
    {
      VsgVector2d *pt;
      Circle *c;

      lb.x = 0.; lb.y = 0.;
      ub.x = 1.*mytid; ub.y = 1.*mytid;

      pt = pt_alloc (TRUE, NULL);
      pt->x = 0.5*mytid; pt->y = 0.5*mytid;

      pt = pt_alloc (TRUE, NULL);
      pt->x = 0.65*mytid; pt->y = 0.65*mytid;

      pt = pt_alloc (TRUE, NULL);
      pt->x = 0.75*mytid; pt->y = 0.75*mytid;

      c = rg_alloc (TRUE, NULL);
      c->center.x = 0.6*mytid; c->center.y = 0.6*mytid;
      c->radius = 0.1;
    }

  /* create the tree */
  tree =
    vsg_prtree2d_new_full (&lb, &ub,
                           (VsgPoint2dLocFunc) vsg_vector2d_vector2d_locfunc,
                           (VsgPoint2dDistFunc) vsg_vector2d_dist,
                           (VsgRegion2dLocFunc) _circle_loc2, 2);

  /* insert the points */
  for (i=0; i<points->len; i++)
    {
      vsg_prtree2d_insert_point (tree, g_ptr_array_index (points, i));
    }

  /* insert the regions */
  for (i=0; i<regions->len; i++)
    {
      vsg_prtree2d_insert_region (tree, g_ptr_array_index (regions, i));
    }

/*   MPI_Barrier (MPI_COMM_WORLD); */
/*   g_printerr ("%d: set_parallel begin\n", mytid); */

  vsg_prtree2d_set_parallel (tree, &pconfig);

/*   MPI_Barrier (MPI_COMM_WORLD); */
/*   g_printerr ("%d: set_parallel ok\n", mytid); */

  {
    VsgVector2d *pt;
    Circle *c;

    pt = pt_alloc (TRUE, NULL);
    pt->x = 0.5*mytid; pt->y = 0.75*mytid;
    vsg_prtree2d_insert_point (tree, pt);

    c = rg_alloc (TRUE, NULL);
    c->center.x = 0.; c->center.y = 0.6*mytid;
    c->radius = 0.1;
    vsg_prtree2d_insert_region (tree, c);
  }

/*   MPI_Barrier (MPI_COMM_WORLD); */
/*   g_printerr ("%d: migrate_flush begin\n", mytid); */

  vsg_prtree2d_migrate_flush (tree);

/*   MPI_Barrier (MPI_COMM_WORLD); */
/*   g_printerr ("%d: migrate_flush ok\n", mytid); */

  MPI_Barrier (MPI_COMM_WORLD);
  {
    gchar fn[1024];
    FILE *f;

    g_sprintf (fn, "prtree2parallel%03d.txt", mytid);
    f = fopen (fn, "w");
    vsg_prtree2d_write (tree, f);
    fclose (f);
  }

  /* destroy the points */
  g_ptr_array_foreach (points, empty_array, NULL);
  g_ptr_array_free (points, TRUE);

  /* destroy the circles */
  g_ptr_array_foreach (regions, empty_array, NULL);
  g_ptr_array_free (regions, TRUE);

  /* destroy the tree */
  vsg_prtree2d_free (tree);

  MPI_Finalize ();

  return ret;
}
