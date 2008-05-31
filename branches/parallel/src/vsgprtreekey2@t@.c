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

/* list of immediately greater than a bits number power of 2 */
static guint8 _greater_power_of_2[KEY_BITS+1];

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
          _bitmasks[i].mask = mask;
          _bitmasks[i].base = ks;
/*           g_printerr ("offset=%u mask num=%u mask=%#@kmod@x\n", */
/*                       ks, i, mask); */
        }

      _greater_power_of_2[0] = 0;
      _greater_power_of_2[1] = 0;
      for (i=0; i<_bitmasks_number; i++)
        {
          for (j=_bitmasks[i].base+1; j<=2*_bitmasks[i].base; j ++)
            {
              _greater_power_of_2[j] = i+1;
/*               g_printerr ("powers %d %d\n", j,  _greater_power_of_2[j]); */
            }
        }
   }
}

void vsg_prtreekey2@t@_write (VsgPRTreeKey2@t@ *key, FILE *file)
{
  g_fprintf (file, "x=%#@kmod@x, y=%#@kmod@x, d=%d",
             key->x, key->y, key->depth);
}

static guint8 _single_key_first_true_bit (@key_type@ one, guint8 maxdepth)
{
  gint8 i, ret = 0;

  _set_bitmasks ();

/*   for (i=0; i<_bitmasks_number && _bitmasks[i].base<maxdepth; i++) */
/*   for (i=_bitmasks_number-1; i >= 0; i--) */
    for (i=_greater_power_of_2[maxdepth]-1; i >= 0; i--)
    {
      @key_type@ masked = one & _bitmasks[i].mask;
      if (!masked) ret += _bitmasks[i].base;
      else one = masked;
/*       g_printerr ("%d, %#@kmod@x, %#@kmod@x, %#@kmod@x, %u\n", i, _bitmasks[i].mask, one, masked, ret); */
/*       g_printerr ("%d, %u %u\n", i, _bitmasks[i].base, maxdepth); */
    }

/*   g_printerr ("size:%d\n", sizeof (@key_type@)); */

  return ret;
}

void vsg_prtree_key2@t@_xor (VsgPRTreeKey2@t@ *one,
                             VsgPRTreeKey2@t@ *other,
                             VsgPRTreeKey2@t@ *result)
{
  result->x = one->x ^ other->x;
  result->y = one->y ^ other->y;
  result->depth = MAX (one->depth, other->depth);
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

  return MIN (x, y);
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
  guint8 index = vsg_prtree_key2@t@_first_different_index (one, other);

  if (index > MIN (one->depth, other->depth))
    {
      memcpy (result, one, sizeof (VsgPRTreeKey2@t@));

      result->depth = MIN (one->depth, other->depth);
    }
  else
    {
      @key_type@ mask = INDEX_MASK (index-1);

      result->x = one->x & mask;
      result->y = one->y & mask;

      result->depth = index;
    }
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
      result->x = father->x |
        ((child_num & VSG_LOC2_X)<<father->depth);
      result->y = father->y |
        (((child_num & VSG_LOC2_Y)>>1)<<father->depth);
      result->depth = father->depth+1;
    }
}

vsgloc2 vsg_prtree_key2@t@_loc2 (VsgPRTreeKey2@t@ *key,
                                 VsgPRTreeKey2@t@ *center)
{
  vsgloc2 loc = 0;
  VsgPRTreeKey2@t@ xor;
  @key_type@ x, y;
  guint8 xi, yi;

  vsg_prtree_key2@t@_xor (key, center, &xor);

  xi = _single_key_first_true_bit (xor.x, xor.depth);
  yi = _single_key_first_true_bit (xor.y, xor.depth);


  x = key->x >> (xi);
  y = key->y >> (yi);

/*   g_printerr ("%#@kmod@x, %#@kmod@x\n", x, y); */

  if (x & 1) loc |= VSG_LOC2_X;
  if (y & 1) loc |= VSG_LOC2_Y;

  return loc;
}
