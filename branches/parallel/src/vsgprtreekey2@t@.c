#include "vsgprtreekey2@t@.h"

#include <string.h>

#include <glib/gprintf.h>

/* number of bytes in a key */
#define KEY_SIZE (sizeof (@key_type@))
/* number of bits in a key */
#define KEY_BITS (8*KEY_SIZE)

#define INDEX_MASK(index) ( \
(1<<((index)+1)) - 1 \
)

typedef struct _BitMaskData BitMaskData;

struct _BitMaskData {
  @key_type@ mask;
  guint8 base;
};

/* a list of binary sieves used to analyze keys */
static guint8 _bitmasks_number = 0;
static BitMaskData *_bitmasks = NULL;

/* list of number of sieves for any  */
static guint8 _number_of_sieves[KEY_BITS+1];

static inline void _set_bitmasks ()
{
  if (G_UNLIKELY (_bitmasks_number == 0))
    {
      guint8 ks = KEY_BITS;
      gint8 i, j;
      @key_type@ mask = ~0;

      while (ks>>(_bitmasks_number+1) != 0) _bitmasks_number++;

      _bitmasks = g_malloc (_bitmasks_number * sizeof (BitMaskData));

      for (i=_bitmasks_number-1; i>=0; i--)
        {
          ks >>= 1; /* ks = 1 << i */
          mask ^= (mask << ks);
          _bitmasks[i].mask = ~ mask;
          _bitmasks[i].base = 1 << (i);
/*           g_printerr ("num=%u offset=%u mask mask=%#@kmod@x\n", */
/*                       i, _bitmasks[i].base, _bitmasks[i].mask); */
        }

      _number_of_sieves[0] = 0;
      _number_of_sieves[1] = 0;
      for (i=0; i<_bitmasks_number; i++)
        {
          for (j=_bitmasks[i].base/2+1; j<=_bitmasks[i].base; j ++)
            {
              _number_of_sieves[j] = i+1;
/*               g_printerr ("powers %d %d\n", j,  _number_of_sieves[j]); */
            }
        }
   }
}

void vsg_prtreekey2@t@_write (VsgPRTreeKey2@t@ *key, FILE *file)
{
  g_fprintf (file, "x=%#@kmod@x, y=%#@kmod@x, d=%d",
             key->x, key->y, key->depth);
}

static void _key_scale_up (VsgPRTreeKey2@t@ *key, guint8 offset,
                        VsgPRTreeKey2@t@ *result)
{
  result->x = key->x << offset;
  result->y = key->y << offset;
  result->depth = key->depth + offset;
}

static guint8 _single_key_first_true_bit (@key_type@ key, guint8 maxdepth)
{
  gint8 sieves_number;
  gint8 i, ret;

  _set_bitmasks ();

  if (!key) return 0;

  sieves_number = _number_of_sieves[maxdepth]-1;
  ret = INDEX_MASK (sieves_number);

  for (i=sieves_number; i >= 0; i--)
    {
      @key_type@ masked = key & _bitmasks[i].mask;
      if (!masked) ret -= _bitmasks[i].base;
      else key = masked;
/*       g_printerr ("%d, %#@kmod@x, %#@kmod@x, %#@kmod@x, %u\n", i, _bitmasks[i].mask, key, masked, ret); */
/*       g_printerr ("%d, %u %u\n", i, _bitmasks[i].base, maxdepth); */
    }

  return ret + 1;
}

void vsg_prtree_key2@t@_xor (VsgPRTreeKey2@t@ *one,
                             VsgPRTreeKey2@t@ *other,
                             VsgPRTreeKey2@t@ *result)
{

  guint8 d;

  if (one->depth >= other->depth)
    {
      d = one->depth - other->depth;
      result->x = one->x ^ (other->x << d);
      result->y = one->y ^ (other->y << d);
      result->depth = one->depth;
    }
  else
    {
      d = other->depth - one->depth;
      result->x = (one->x << d) ^ other->x;
      result->y = (one->y << d) ^ other->y;
      result->depth = other->depth;
    }
}

static guint8
vsg_prtree_key2@t@_first_different_index_internal (VsgPRTreeKey2@t@ *one,
                                                   VsgPRTreeKey2@t@ *other,
                                                   VsgPRTreeKey2@t@ *xor)
{
  guint8 x, y;

  vsg_prtree_key2@t@_xor (one, other, xor);

  x = _single_key_first_true_bit (xor->x, xor->depth);
  y = _single_key_first_true_bit (xor->y, xor->depth);

  return MAX (x, y);
}

guint8 vsg_prtree_key2@t@_first_different_index (VsgPRTreeKey2@t@ *one,
                                                 VsgPRTreeKey2@t@ *other)
{
  VsgPRTreeKey2@t@ xor;

  return vsg_prtree_key2@t@_first_different_index_internal (one, other, &xor);
}

void vsg_prtree_key2@t@_deepest_common_ancestor (VsgPRTreeKey2@t@ *one,
                                                 VsgPRTreeKey2@t@ *other,
                                                 VsgPRTreeKey2@t@ *result)
{
  VsgPRTreeKey2@t@ *max = one;
  guint8 index;

  if (one->depth < other->depth) max = other;

  index = vsg_prtree_key2@t@_first_different_index (one, other);

  result->x = max->x >> index;
  result->y = max->y >> index;

  result->depth = max->depth - index;
}

void vsg_prtree_key2@t@_build_child (VsgPRTreeKey2@t@ *father,
                                     vsgloc2 child_num,
                                     VsgPRTreeKey2@t@ *result)
{
  if (father == NULL)
    {
      result->x = 0;
      result->y = 0;
      result->depth = 0;
    }
  else
    {
      result->x = (father->x << 1) | (child_num & VSG_LOC2_X);
      result->y = (father->y << 1) | ((child_num & VSG_LOC2_Y)>>1);

      result->depth = father->depth+1;
    }
}

vsgloc2 vsg_prtree_key2@t@_loc2 (VsgPRTreeKey2@t@ *key,
                                 VsgPRTreeKey2@t@ *center)
{
  VsgPRTreeKey2@t@ tmp;
  VsgPRTreeKey2@t@ *k = key;
  VsgPRTreeKey2@t@ *c = center;
  vsgloc2 loc = 0;

  if (k->depth > c->depth)
    {
      _key_scale_up (c, k->depth - c->depth, &tmp);
      c = &tmp;
    }
  else if (k->depth < c->depth)
    {
      _key_scale_up (k, c->depth - k->depth, &tmp);
      k = &tmp;
    }

  if (k->x > c->x) loc |= VSG_LOC2_X;
  if (k->y > c->y) loc |= VSG_LOC2_Y;

  return loc;
}
