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
};

static const TestData _test_data_sentinel = {
  {0x0, 0x0, 0}, {0x0, 0x0, 0},  0, {0x0, 0x0, 0}, 0, 0
};

static TestData _test_data[] = {
  {{0x132, 0x0, 10}, {0x132, 0x0, 10}, 15, {0x132, 0x0, 10}, 0, 0},
  {{0x132, 0x0, 10}, {0x72, 0x0, 10}, 6, {0x32, 0x0, 6}, 0, 1},
  {{0x72, 0x0, 10}, {0x132, 0x0, 10}, 6, {0x32, 0x0, 6}, 1, 0},
  {{0x0, 0x132, 10}, {0x0, 0x132, 10}, 15, {0x0, 0x132, 10}, 0, 0},
  {{0x0, 0x132, 10}, {0x0, 0x72, 10}, 6, {0x0, 0x32, 6}, 0, 2},
  {{0x05557, 0xabcde, 17}, {0xaaa57, 0xabcde, 20},  8, {0x57, 0xde, 8}, 1, 0},
  {{0x0, 0x0, 2}, {0x1, 0x2, 2},  0, {0x0, 0x0, 0}, 0, 3},
  {{0x0, 0x0, 2}, {0x1, 0x0, 2},  0, {0x0, 0x0, 0}, 0, 1},
  //{{0x, 0x, }, {0x, 0x, },  , {0x, 0x, }, , },
  {{0x0, 0x0, 0}, {0x0, 0x0, 0},  0, {0x0, 0x0, 0}, 0, 0}, // sentinel
};

static gboolean _check_test_data (TestData *data)
{
  gint8 ldi;
  VsgPRTreeKey2d dpa;
  vsgloc2 loc1, loc2;

  ldi = vsg_prtree_key2d_first_different_index (&data->k1, &data->k2);

  if (ldi != data->ldi_ref)
    {
      g_printerr ("Error on first_different_index:\n");
      g_printerr ("k1: ");
      vsg_prtreekey2d_write (&data->k1, stderr);
      g_printerr ("\nk2: ");
      vsg_prtreekey2d_write (&data->k2, stderr);

      g_printerr ("\ngot %d, should be %d\n\n", ldi, data->ldi_ref);
    }

  vsg_prtree_key2d_deepest_common_ancestor(&data->k1, &data->k2, &dpa);

  if (memcmp (&data->dpa_ref, &dpa, sizeof (VsgPRTreeKey2d)) != 0)
    {
      g_printerr ("Error on deepest_common_ancestor:\n");
      g_printerr ("k1: ");
      vsg_prtreekey2d_write (&data->k1, stderr);
      g_printerr ("\nk2: ");
      vsg_prtreekey2d_write (&data->k2, stderr);

      g_printerr ("\ngot %#"G_GINT64_MODIFIER"x %#"G_GINT64_MODIFIER"x d=%d\n",
                  dpa.x, dpa.y, dpa.depth);
      g_printerr ("should be %#"G_GINT64_MODIFIER"x %#"G_GINT64_MODIFIER"x "
                  "d=%d\n\n",
                  data->dpa_ref.x, data->dpa_ref.y, data->dpa_ref.depth);
    }

  loc1 = vsg_prtree_key2d_loc2 (&data->k1, &data->k2);
  if (loc1 != data->loc1_ref)
    {
      g_printerr ("Error on loc2 (k1, k2):\n");
      g_printerr ("k1: ");
      vsg_prtreekey2d_write (&data->k1, stderr);
      g_printerr ("\nk2: ");
      vsg_prtreekey2d_write (&data->k2, stderr);

      g_printerr ("\ngot %d, should be %d\n\n", loc1, data->loc1_ref);
    }

  loc2 = vsg_prtree_key2d_loc2 (&data->k2, &data->k1);
  if (loc2 != data->loc2_ref)
    {
      g_printerr ("Error on loc2 (k2, k1):\n");
      g_printerr ("k1: ");
      vsg_prtreekey2d_write (&data->k1, stderr);
      g_printerr ("\nk2: ");
      vsg_prtreekey2d_write (&data->k2, stderr);

      g_printerr ("\ngot %d, should be %d\n\n", loc2, data->loc2_ref);
    }

  return memcmp (data, &_test_data_sentinel, sizeof (TestData)) == 0;
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

  return ret;
}
