/* basic point insertion, removal and traversal on a cloned VsgPRTree2d */

#include "vsg-config.h"

#include "vsg/vsgd.h"

#include "glib/gprintf.h"

#include <stdlib.h>
#include <math.h>
#include <glib.h>


static gboolean _do_write = FALSE;

static gint rk, sz;

static GPtrArray *points = NULL;

static gint _np = 12;
static guint32 _random_seed = 0;
static gint _flush_interval = 10;

static gpointer pt_alloc (gboolean resident, gpointer user_data)
{
  gpointer ret;
  ret = g_malloc (sizeof (VsgVector2d));
  g_ptr_array_add (points, ret);
/*   g_printerr ("%d: alloc 1 VsgVector2d (%p)\n", rk, ret); */
  return ret;
}

static void pt_destroy (gpointer data, gboolean resident,
                        gpointer user_data)
{
/*   g_printerr ("%d: destroy 1 VsgVector2d (%p)\n", rk, data); */
  g_ptr_array_remove (points, data);
  g_free (data);
}

static void pt_migrate_pack (VsgVector2d *pt, VsgPackedMsg *pm,
                             gpointer user_data)
{
/*   g_printerr ("%d: pack ", rk); */
/*   vsg_vector2d_write (pt, stderr); */
/*   g_printerr ("\n"); */
  vsg_packed_msg_send_append (pm, pt, 1, VSG_MPI_TYPE_VECTOR2D);
}

static void pt_migrate_unpack (VsgVector2d *pt, VsgPackedMsg *pm,
                               gpointer user_data)
{
  vsg_packed_msg_recv_read (pm, pt, 1, VSG_MPI_TYPE_VECTOR2D);

/*   g_printerr ("%d: unpack ", rk); */
/*   vsg_vector2d_write (pt, stderr); */
/*   g_printerr ("\n"); */
}

static void empty_array (gpointer var, gpointer data)
{
/*   g_printerr ("%d: static destroy 1 VsgVector2d (%p)\n", rk, var); */
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

/*   g_printerr ("%d: {c=", rk); */
/*   vsg_vector2d_write (&circle->center, stderr); */
/*   g_printerr (" r=%g} ref=", circle->radius); */
/*   vsg_vector2d_write (center, stderr); */

  ret |= VSG_RLOC2_COMP (centerpos);

/*   g_printerr (" / 0x%X", ret); */

  if (fabs (tmp.x) <= circle->radius)
    {
      vsgloc2 itmp = centerpos ^ VSG_LOC2_X;
      ret |= VSG_RLOC2_COMP (itmp);
/*       g_printerr (" / x<rad 0x%X (%x %x)", itmp, centerpos, VSG_LOC2_X); */
    }

  if (fabs (tmp.y) <= circle->radius)
    {
      vsgloc2 itmp = centerpos ^ VSG_LOC2_Y;
      ret |= VSG_RLOC2_COMP (itmp);
/*       g_printerr (" / y<rad 0x%X (%x %x)", itmp, centerpos, VSG_LOC2_Y); */
    }

/*   g_printerr (" => 0x%X\n", ret); */

  return ret;
}

static gpointer rg_alloc (gboolean resident, gpointer user_data)
{
  gpointer ret;
  ret = g_malloc (sizeof (Circle));
  g_ptr_array_add (regions, ret);

/*   g_printerr ("%d: alloc 1 Circle (%p)\n", rk, ret); */

  return ret;
}

static void rg_destroy (gpointer data, gboolean resident,
                        gpointer user_data)
{
/*   g_printerr ("%d: destroy 1 Circle (%p)\n", rk, data); */
/*   contiguous_distribute_nodes (tree); */

  g_ptr_array_remove (regions, data);
  g_free (data);
}

static void rg_migrate_pack (Circle *rg, VsgPackedMsg *pm,
                             gpointer user_data)
{
/*   g_printerr ("%d: pack circle {c=", rk); */
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

/*   g_printerr ("%d: unpack circle {c=", rk); */
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


static gint concentrate_dist (VsgPRTree2dNodeInfo *node_info, gint *dst)
{
/*   g_printerr ("%d: dist node center ", rk); */
/*   vsg_vector2d_write (&node_info->center, stderr); */
/*   g_printerr (" dst=%d\n", *dst); */

  return *dst;
}

static gint scatter_dist (VsgPRTree2dNodeInfo *node_info, gint *cptr)
{
  if (VSG_PRTREE2D_NODE_INFO_IS_LOCAL (node_info))
    {
      gint dst = *cptr;

/*       g_printerr ("%d: dist node center ", rk); */
/*       vsg_vector2d_write (&node_info->center, stderr);  */
/*       g_printerr (" dst=%d\n", dst); */

      (*cptr) ++;
      (*cptr) = (*cptr) % sz;

      return dst;
    }

  return -1;
}

static void _pt_write (VsgVector2d *pt, FILE *file)
{
  fprintf (file, "<circle cx=\"%g\" cy=\"%g\" r=\"%g\" " \
           "style=\"stroke:#000000;fill:#00ff00;\"/>\n",
           pt->x, -pt->y, 0.02);
}

static void _rg_write (Circle *c, FILE *file)
{
  fprintf (file, "<circle cx=\"%g\" cy=\"%g\" r=\"%g\" " \
           "style=\"stroke:#000000;fill:none;\"/>\n",
           c->center.x, -c->center.y, c->radius);

}

static void _traverse_bg_write (VsgPRTree2dNodeInfo *node_info, FILE *file)
{
  gdouble x = node_info->lbound.x;
  gdouble y = -node_info->ubound.y;
  gdouble w = node_info->ubound.x - x;
  gdouble h = -node_info->lbound.y - y;
  gchar *fill = "none";

  if (!node_info->isleaf) return;

  if (VSG_PRTREE2D_NODE_INFO_IS_REMOTE (node_info))
    fill = "#ff0000";

  fprintf (file, "<rect x=\"%g\" y=\"%g\" width=\"%g\" height=\"%g\" " \
           "rx=\"0\" style=\"stroke:#000000; " \
           "stroke-linejoin:miter; stroke-linecap:butt; fill:%s;\"/>\n",
           x, y, w, h, fill);
}

static void _traverse_fg_write (VsgPRTree2dNodeInfo *node_info, FILE *file)
{
  g_slist_foreach (node_info->region_list, (GFunc) _rg_write, file);
  g_slist_foreach (node_info->point_list, (GFunc) _pt_write, file);
}

static void _tree_write (VsgPRTree2d *tree)
{
  gchar fn[1024];
  FILE *f;
  VsgVector2d lb, ub;
  gdouble x;
  gdouble y;
  gdouble w;
  gdouble h;

  vsg_prtree2d_get_bounds (tree, &lb, &ub);

  x = lb.x;
  y = -ub.y;
  w = ub.x - x;
  h = -lb.y - y;

  g_sprintf (fn, "prtree2parallel-%03d.svg", rk);
  f = fopen (fn, "w");

  fprintf (f, "<?xml version=\"1.0\" standalone=\"no\"?>\n" \
           "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n" \
           "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n" \
           "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"10cm\" " \
           "height=\"10cm\" viewBox=\"%g %g %g %g\">\n" \
           "<g style=\"stroke-width:0.01; stroke:black; fill:none\">\n",
           x, y, w, h);

  vsg_prtree2d_traverse (tree, G_PRE_ORDER,
                         (VsgPRTree2dFunc) _traverse_bg_write,
                         f);
  vsg_prtree2d_traverse (tree, G_PRE_ORDER,
                         (VsgPRTree2dFunc) _traverse_fg_write,
                         f);

  fprintf (f, "</g>\n</svg>\n");

  fclose (f);

}

static gint reference_total_points = 0;

static void init_total_points_count ()
{
  gint local_points = points->len;

  MPI_Allreduce (&local_points, &reference_total_points, 1, MPI_INT, MPI_SUM,
                 MPI_COMM_WORLD);
}

static void _local_points_count (VsgPRTree2dNodeInfo *node_info, gint *cnt)
{
  *cnt += g_slist_length (node_info->point_list);
}

static gint check_points_number (VsgPRTree2d *tree)
{
  gint ret = 0;
  gint local_points = 0;
  gint total_points = 0;

  vsg_prtree2d_traverse (tree, G_PRE_ORDER,
                         (VsgPRTree2dFunc) _local_points_count, &local_points);

  MPI_Allreduce (&local_points, &total_points, 1, MPI_INT, MPI_SUM,
                 MPI_COMM_WORLD);

  if (total_points != reference_total_points)
    {
      g_printerr ("%d: total points mismatch : ref=%d verif=%d\n",
                  rk, reference_total_points, total_points);

      ret ++;
    }

  return ret;
}

static gint reference_total_regions = 0;

static void init_total_regions_count ()
{
  gint local_regions = regions->len;

  MPI_Allreduce (&local_regions, &reference_total_regions, 1, MPI_INT, MPI_SUM,
                 MPI_COMM_WORLD);
}

static void _local_regions_count (VsgPRTree2dNodeInfo *node_info, gint *cnt)
{
  if (VSG_PRTREE2D_NODE_INFO_IS_LOCAL (node_info) || rk == 0)
    *cnt += g_slist_length (node_info->region_list);
}

static gint check_regions_number (VsgPRTree2d *tree)
{
  gint ret = 0;
  gint local_regions = 0;
  gint total_regions = 0;

  vsg_prtree2d_traverse (tree, G_PRE_ORDER,
                         (VsgPRTree2dFunc) _local_regions_count,
                         &local_regions);

  MPI_Allreduce (&local_regions, &total_regions, 1, MPI_INT, MPI_SUM,
                 MPI_COMM_WORLD);

  if (total_regions != reference_total_regions)
    {
      g_printerr ("%d: total regions mismatch : ref=%d verif=%d\n",
                  rk, reference_total_regions, total_regions);

      ret ++;
    }

  return ret;
}

static void _local_leaves_count (VsgPRTree2dNodeInfo *node_info, GArray *array)
{
  if (! VSG_PRTREE2D_NODE_INFO_IS_SHARED (node_info))
    {
      /* detect root of local (remote) subtrees */
      if (node_info->father_info == NULL ||
          VSG_PRTREE2D_NODE_INFO_IS_SHARED (node_info->father_info))
        {
          gint i = 0;
          g_array_append_val (array, i);
        }

      /* increment local leaves */
      if (node_info->isleaf && VSG_PRTREE2D_NODE_INFO_IS_LOCAL (node_info))
        g_array_index (array, gint, array->len-1) ++;
    }
}

typedef struct _ContiguousDistData ContiguousDistData;
struct _ContiguousDistData {
  GArray *array;        /* array of number of leaves on each local subtree */
  gint total_leaves;   /* total number of leaves in the tree (sum of array) */
  gint q, r, m;        /* distribution variables */
  gint current_index;  /* number of already checked array entries */
  gint current_lcount; /* number of already checked leaves */
};

static gint contiguous_dist (VsgPRTree2dNodeInfo *node_info,
                             ContiguousDistData *cda)
{
  if (VSG_PRTREE2D_NODE_INFO_IS_LOCAL (node_info))
    {
      gint ret = 0;

      if (cda->current_lcount >= cda->m)
        ret = (cda->current_lcount - cda->r) / cda->q;
      else
        ret = cda->current_lcount / (cda->q + 1);

      cda->current_lcount ++;

      return ret;
    }

  g_assert (VSG_PRTREE2D_NODE_INFO_IS_REMOTE (node_info));

  cda->current_lcount += g_array_index (cda->array, gint, cda->current_index);
  cda->current_index ++;

  return -1;
}

static void contiguous_distribute_nodes (VsgPRTree2d *tree)
{
  ContiguousDistData cda;
  GArray *array = g_array_sized_new (FALSE, FALSE, sizeof (gint), 1024);
  GArray *reduced;
  gint i;

  /* accumulate number of local leaves in the local subtrees' roots */
  vsg_prtree2d_traverse (tree, G_PRE_ORDER,
                         (VsgPRTree2dFunc) _local_leaves_count,
                         array);

  /* prepare reduced for storing the results of the Allreduce */
  reduced = g_array_sized_new (FALSE, TRUE, sizeof (gint), array->len);
  reduced = g_array_set_size (reduced, array->len);

  MPI_Allreduce (array->data, reduced->data, array->len, MPI_INT, MPI_MAX,
                 MPI_COMM_WORLD);

  g_array_free (array, TRUE);

  cda.array = reduced;
  cda.total_leaves = 0;
  for (i=0; i<reduced->len; i++)
    cda.total_leaves += g_array_index (reduced, gint, i);

  cda.q = cda.total_leaves / sz;
  cda.r = cda.total_leaves % sz;
  cda.m = (cda.q+1) * cda.r;
  cda.current_lcount = 0;
  cda.current_index = 0;

  vsg_prtree2d_distribute_nodes (tree,
                                 (VsgPRTree2dDistributionFunc) contiguous_dist,
                                 &cda);

  g_array_free (reduced, TRUE);

}

static void random_fill (VsgPRTree2d *tree, guint np);

static void (*_fill) (VsgPRTree2d *tree, guint np) = random_fill;

static void random_fill (VsgPRTree2d *tree, guint np)
{
  gint i;
  VsgVector2d *pt;
  Circle *c;
  VsgVector2d lb, ub;
  GRand *rand = g_rand_new_with_seed (_random_seed);

  vsg_prtree2d_get_bounds (tree, &lb, &ub);

  for (i=0; i< np; i++)
    {
      gdouble x1, x2, y1, y2, r;

      x1 = g_rand_double_range (rand, lb.x, ub.x);
      y1 = g_rand_double_range (rand, lb.y, ub.y);
      x2 = g_rand_double_range (rand, lb.x, ub.x);
      y2 = g_rand_double_range (rand, lb.y, ub.y);
      r = g_rand_double_range (rand, 0., 0.1);

      if (i%_flush_interval == 0) vsg_prtree2d_migrate_flush (tree);

      if (i%sz != rk) continue;

      pt = pt_alloc (TRUE, NULL);
      pt->x = x1;
      pt->y = y1;
      vsg_prtree2d_insert_point (tree, pt);

      c = rg_alloc (TRUE, NULL);
      c->center.x = x2;
      c->center.y = y2;
      c->radius = r;
      vsg_prtree2d_insert_region (tree, c);

/*       g_printerr ("%d: %d\n", rk, i); */
    }

  vsg_prtree2d_migrate_flush (tree);

  g_rand_free (rand);
}

static
void parse_args (int argc, char **argv)
{
  int iarg = 1;
  char *arg;

  while (iarg < argc)
    {
      arg = argv[iarg];

      if (g_ascii_strcasecmp (arg, "-np") == 0)
	{
	  guint tmp = 0;
	  iarg ++;

	  arg = (iarg<argc) ? argv[iarg] : NULL;

	  if (sscanf (arg, "%u", &tmp) == 1)
            _np = tmp;
	  else
	    g_printerr ("Invalid particles number (-np %s)\n", arg);
	}
      else if (g_strncasecmp (arg, "--write", 7) == 0)
        {
          _do_write = TRUE;
        }
      else if (g_ascii_strcasecmp (arg, "--version") == 0)
	{
	  g_printerr ("%s version %s\n", argv[0], PACKAGE_VERSION);
	  exit (0);
	}
      else
	{
	  g_printerr ("Invalid argument \"%s\"\n", arg);
	}

      iarg ++;
    }
}


gint main (gint argc, gchar ** argv)
{
  gint ret = 0;
  gint i;

  VsgPRTree2d *tree;

  VsgVector2d lb;
  VsgVector2d ub;


  MPI_Init (&argc, &argv);

  MPI_Comm_size (MPI_COMM_WORLD, &sz);
  MPI_Comm_rank (MPI_COMM_WORLD, &rk);

  vsg_init_gdouble ();

  parse_args (argc, argv);

  points = g_ptr_array_new ();
  regions = g_ptr_array_new ();

  lb.x = -1.; lb.y = -1.;
  ub.x = 1.; ub.y = 1.;

  /* create the tree */
  tree =
    vsg_prtree2d_new_full (&lb, &ub,
                           (VsgPoint2dLocFunc) vsg_vector2d_vector2d_locfunc,
                           (VsgPoint2dDistFunc) vsg_vector2d_dist,
                           (VsgRegion2dLocFunc) _circle_loc2, 2);

/*   MPI_Barrier (MPI_COMM_WORLD); */
/*   g_printerr ("%d: set_parallel begin\n", rk); */

  vsg_prtree2d_set_parallel (tree, &pconfig);

/*   MPI_Barrier (MPI_COMM_WORLD); */
/*   g_printerr ("%d: set_parallel ok\n", rk); */

  _fill (tree, _np);

  /* update total points and regions count */
  init_total_points_count ();
  init_total_regions_count ();

/*   MPI_Barrier (MPI_COMM_WORLD); */
/*   g_printerr ("%d: distribute_nodes begin\n", rk); */

/*   i = 0; */
/*   vsg_prtree2d_distribute_nodes (tree, */
/*                                  (VsgPRTree2dDistributionFunc) */
/*                                  scatter_dist, */
/*                                  &i); */

  contiguous_distribute_nodes (tree);

  ret += check_points_number (tree);
  ret += check_regions_number (tree);

/*   MPI_Barrier (MPI_COMM_WORLD); */
/*   g_printerr ("%d: distribute_nodes ok\n", rk); */

  if (_do_write)
    {
      MPI_Barrier (MPI_COMM_WORLD);
      _tree_write (tree);
    }

  if (_do_write)
    {
      gchar fn[1024];
      FILE *f;

      g_sprintf (fn, "prtree2parallel-%03d.txt", rk);
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