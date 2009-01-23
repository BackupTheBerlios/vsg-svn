/* basic point insertion, removal and traversal on a cloned VsgPRTree2d */

#include "vsg-config.h"

#include "vsg/vsgd.h"

#include <string.h>

typedef struct _TestData TestData;

struct _TestData {
  VsgPRTreeKey2d k1;
  VsgPRTreeKey2d k2;
  guint8 ldi_ref;
  VsgPRTreeKey2d dpa_ref;
  vsgloc2 loc1_ref;
  vsgloc2 loc2_ref;
  guint64 dist_ref;
};

static const TestData _test_data_sentinel = {
  {0x0, 0x0, 0}, {0x0, 0x0, 0},  0, {0x0, 0x0, 0}, 0, 0, 0x0
};

static TestData _test_data[] = {
  {{0x132, 0x0, 10}, {0x132, 0x0, 10}, 0, {0x132, 0x0, 10}, 0, 0, 0x0},
  {{0x132, 0x0, 10}, {0x133, 0x0, 10}, 1, {0x99, 0x0, 9}, 0, 1, 0x1},
  {{0x132, 0x0, 10}, {0x72, 0x0, 10}, 9, {0x0, 0x0, 1}, 1, 0, 0xc0},
  {{0x72, 0x0, 10}, {0x132, 0x0, 10}, 9, {0x0, 0x0, 1}, 0, 1, 0xc0},
  {{0x0, 0x132, 10}, {0x0, 0x132, 10}, 0, {0x0, 0x132, 10}, 0, 0, 0x0},
  {{0x0, 0x132, 10}, {0x0, 0x72, 10}, 9, {0x0, 0x0, 1}, 2, 0, 0xc0},
  {{0x0aab5, 0x1579b, 17}, {0x555aa, 0xabcde, 20},  3, {0xaab5, 0x1579b, 17}, 0, 3, },

  {{0x0, 0x0, 2}, {0x1, 0x2, 2},  2, {0x0, 0x0, 0}, 0, 3, 0x2},
  {{0x0, 0x0, 2}, {0x1, 0x0, 2},  1, {0x0, 0x0, 1}, 0, 1, 0x1},

  {{0xff, 0xef, 8}, {0x70, 0xef, 8},  8, {0x0, 0x0, 0}, 1, 0, 0x8f},
  {{0xffaa, 0xefaa, 16}, {0x70aa, 0xefaa, 16},  16, {0x0, 0x0, 0}, 1, 0, 0x8f00},

  {{0x1, 0x0, 1}, {0x2, 0x0, 2},  0, {0x2, 0x0, 2}, 0, 0, 0x0},
  {{0x0, 0x1, 1}, {0x0, 0x2, 2},  0, {0x0, 0x2, 2}, 0, 0, 0x0},
  {{0x1, 0x0, 1}, {0x4, 0x0, 3},  0, {0x4, 0x0, 3}, 0, 0, 0x0},
  {{0x0, 0x1, 1}, {0x0, 0x4, 3},  0, {0x0, 0x4, 3}, 0, 0, 0x0},

  {{0x2, 0x0, 2}, {0x1, 0x0, 1},  0, {0x2, 0x0, 2}, 0, 0, 0x0},
  {{0x0, 0x2, 2}, {0x0, 0x1, 1},  0, {0x0, 0x2, 2}, 0, 0, 0x0},
  {{0x4, 0x0, 3}, {0x1, 0x0, 1},  0, {0x4, 0x0, 3}, 0, 0, 0x0},
  {{0x0, 0x4, 3}, {0x0, 0x1, 1},  0, {0x0, 0x4, 3}, 0, 0, 0x0},
  //{{0x, 0x, }, {0x, 0x, },  , {0x, 0x, }, , },
  {{0x0, 0x0, 0}, {0x0, 0x0, 0},  0, {0x0, 0x0, 0}, 0, 0, 0x0}, // sentinel
};

static gboolean _check_test_data (TestData *data)
{
  gint8 ldi;
  VsgPRTreeKey2d dpa = {0};
  vsgloc2 loc1, loc2;
  guint64 dist;

  ldi = vsg_prtree_key2d_first_different_index (&data->k1, &data->k2);

  if (ldi != data->ldi_ref)
    {
      g_printerr ("Error on first_different_index:\n");
      g_printerr ("k1: ");
      vsg_prtree_key2d_write (&data->k1, stderr);
      g_printerr ("\nk2: ");
      vsg_prtree_key2d_write (&data->k2, stderr);

      g_printerr ("\ngot %d, should be %d\n\n", ldi, data->ldi_ref);
    }

  vsg_prtree_key2d_deepest_common_ancestor(&data->k1, &data->k2, &dpa);

  if (memcmp (&data->dpa_ref, &dpa, sizeof (VsgPRTreeKey2d)) != 0)
    {
      g_printerr ("Error on deepest_common_ancestor:\n");
      g_printerr ("k1: ");
      vsg_prtree_key2d_write (&data->k1, stderr);
      g_printerr ("\nk2: ");
      vsg_prtree_key2d_write (&data->k2, stderr);

      g_printerr ("\ngot (");
      vsg_prtree_key2d_write (&dpa, stderr);
      g_printerr (")\nshould be (");
      vsg_prtree_key2d_write (&data->dpa_ref, stderr);
      g_printerr (")\n\n");
    }

  loc1 = vsg_prtree_key2d_loc2 (&data->k1, &data->k2);
  if (loc1 != data->loc1_ref)
    {
      g_printerr ("Error on loc2 (k1, k2):\n");
      g_printerr ("k1: ");
      vsg_prtree_key2d_write (&data->k1, stderr);
      g_printerr ("\nk2: ");
      vsg_prtree_key2d_write (&data->k2, stderr);

      g_printerr ("\ngot %d, should be %d\n\n", loc1, data->loc1_ref);
    }

  loc2 = vsg_prtree_key2d_loc2 (&data->k2, &data->k1);
  if (loc2 != data->loc2_ref)
    {
      g_printerr ("Error on loc2 (k2, k1):\n");
      g_printerr ("k1: ");
      vsg_prtree_key2d_write (&data->k1, stderr);
      g_printerr ("\nk2: ");
      vsg_prtree_key2d_write (&data->k2, stderr);

      g_printerr ("\ngot %d, should be %d\n\n", loc2, data->loc2_ref);
    }

  dist = vsg_prtree_key2d_distance (&data->k1, &data->k2);
  if (dist != data->dist_ref)
    {
      g_printerr ("Error on dist (k1, k2):\n");
      g_printerr ("k1: ");
      vsg_prtree_key2d_write (&data->k1, stderr);
      g_printerr ("\nk2: ");
      vsg_prtree_key2d_write (&data->k2, stderr);

      g_printerr ("\ngot %#"G_GINT64_MODIFIER"x, should be %#"G_GINT64_MODIFIER"x\n\n", dist, data->dist_ref);

    }

  return memcmp (data, &_test_data_sentinel, sizeof (TestData)) == 0;
}

typedef struct _NFTestData NFTestData;

struct _NFTestData {
  VsgPRTreeKey2d k1;
  VsgPRTreeKey2d k2;
  guint8 mindepth;
  guint8 cmp_mindepth_ref;
};

static NFTestData _nf_test_data[] = {
  {{0x04, 0x0, 3}, {0x19, 0x0, 5}, 4, 3},
  {{0x1, 0x1, 2}, {0x0, 0x1, 1}, 1, 1},
  {{0x1, 0x1, 2}, {0x0, 0x1, 1}, 2, 1},
  {{0x1, 0x0, 2}, {0x0, 0x1, 1}, 2, 2},
  {{0x1, 0x1, 3}, {0x0, 0x2, 2}, 3, 3},
  {{0x358, 0x356, 10}, {0xd66, 0xd58, 12}, 12, 3},
  {{0x0, 0x0, 1}, {0x6, 0x0, 3}, 3, 3},
  {{0x1, 0x2, 2}, {0x4, 0x5, 3}, 2, 1},
  {{0x165, 0x13a, 11}, {0x2cc, 0x274, 12}, 11, 1},
  {{0x0, 0x0, 0}, {0x0, 0x0, 0}, 0, 0},
};

static NFTestData _nf_test_data_sentinel = {
  {0x0, 0x0, 0}, {0x0, 0x0, 0}, 0, 0
};

gboolean _check_nf_mindepth (NFTestData *data)
{
  guint8 cmp_mindepth;

  if (memcmp (data, &_nf_test_data_sentinel, sizeof (NFTestData)) == 0)
    return TRUE;

  cmp_mindepth = vsg_prtree_key2d_compare_near_far_mindepth (&data->k1,
                                                             &data->k2,
                                                             data->mindepth);
  if (cmp_mindepth != data->cmp_mindepth_ref)
    {
      g_printerr ("Error on compare_near_far_mindepth (k1, k2, mindepth):\n");
      g_printerr ("k1: ");
      vsg_prtree_key2d_write (&data->k1, stderr);
      g_printerr ("\nk2: ");
      vsg_prtree_key2d_write (&data->k2, stderr);

      g_printerr ("\nmindepth: %u\ngot %u, should be %u\n\n", data->mindepth,
                  cmp_mindepth,
                  data->cmp_mindepth_ref);
    }

  return FALSE;
}

void _nf_error (VsgPRTreeKey2d *k1, VsgPRTreeKey2d *k2, gint8 got,
                gint8 shouldbe)
{
  g_printerr ("Error on compare_near_far (k1, k2):\n");
  g_printerr ("k1: ");
  vsg_prtree_key2d_write (k1, stderr);
  g_printerr ("\nk2: ");
  vsg_prtree_key2d_write (k2, stderr);

  g_printerr ("\ngot %d should be %d\n\n", got, shouldbe);

}

void _check_nf (VsgPRTreeKey2d *k1, guint8 lvl, gint8 window)
{
  /* check k1 against neighbours */
  if (k1->depth > 1)
    {
      VsgPRTreeKey2d k2;
      guint64 i, j;
      guint64 n = (1 << k1->depth) - 1;

      for (i=MAX (0, k1->x-window); i<MIN (n, k1->x+window); i++)
        for (j=MAX (0, k1->y-window); j<MIN (n, k1->y+window); j++)
          {
            gint8 nf;
            gint64 dx, dy, diff, fdx, fdy, fdiff;

            k2.x = i;
            k2.y = j;
            k2.depth = k1->depth;

            nf = vsg_prtree_key2d_compare_near_far (k1, &k2);

            dx = k2.x-k1->x;
            dy = k2.y-k1->y;
            diff = MAX (ABS (dx), ABS (dy));

            fdx = (k2.x>>1)-(k1->x>>1);
            fdy = (k2.y>>1)-(k1->y>>1);
            fdiff = MAX (ABS (fdx), ABS (fdy));

            if (nf < 0) 
              {
                /* ancestor order */
              }
            else if (diff == 0)
              {
                if (nf != 0)
                  _nf_error (k1, &k2, nf, 0);
              }
            else if (diff == 1)
              {
                if (nf != 1)
                  _nf_error (k1, &k2, nf, 1);
              }
            else if (diff == 2)
              {
                if (nf != 2)
                  _nf_error (k1, &k2, nf, 2);
              }
            else if (diff == 3 && fdiff <= 1)
              {
                if (nf != 2)
                  _nf_error (k1, &k2, nf, 2);
              }
            else
              {
                if (nf != 3)
                  _nf_error (k1, &k2, nf, 3);
              }
          }
      
    }

  /* step down in the tree */
  if (lvl > 0)
    {
      vsgloc2 i;

      for (i=0; i<4; i++)
        {
          VsgPRTreeKey2d child;

          vsg_prtree_key2d_build_child (k1, i, &child);

          _check_nf (&child, lvl-1, window);
        }
    }
}

gint main (gint argc, gchar ** argv)
{
  gint ret = 0;
  gint i = 0;
  vsg_init_gdouble ();

  if (argc > 1 && g_strncasecmp (argv[1], "--version", 9) == 0)
    {
      g_print ("%s\n", PACKAGE_VERSION);
      return 0;
    }

  while (! _check_test_data (&_test_data[i]))
    i ++;

  {
    VsgPRTreeKey2d k1 = vsg_prtree_key2d_root;
    _check_nf (&k1, 8, 5);
  }

  i = 0;
  while (! _check_nf_mindepth (&_nf_test_data[i]))
    i ++;

  return ret;
}
