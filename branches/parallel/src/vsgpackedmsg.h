/* LIBVSG - Visaurin Geometric Library
 * Copyright (C) 2006-2007 Pierre Gay
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __VSG_PACKEDMSG_H__ /* __VSG_PACKEDMSG_H__ */
#define __VSG_PACKEDMSG_H__

#include <glib.h>
#include <glib-object.h>

#include <vsg/vsgmpi.h>

G_BEGIN_DECLS;

typedef struct _VsgPackedMsg VsgPackedMsg;

struct _VsgPackedMsg
{
  MPI_Comm communicator;

  /* <private> */

  gchar *buffer;
  gint position;
  gint size;
  gboolean own_buffer;
};

static const VsgPackedMsg VSG_PACKED_MSG_NULL = {
  MPI_COMM_NULL,
  NULL,
  0,
  0,
  TRUE
};

VsgPackedMsg *vsg_packed_msg_new (MPI_Comm comm);

void vsg_packed_msg_init (VsgPackedMsg *pm, MPI_Comm comm);

void vsg_packed_msg_set_reference (VsgPackedMsg *pm, VsgPackedMsg *model);

void vsg_packed_msg_send_append (VsgPackedMsg *pm, gpointer buf,
                                 gint count, MPI_Datatype type);


void vsg_packed_msg_recv_read (VsgPackedMsg *pm, gpointer buf,
                               gint count, MPI_Datatype type);

void vsg_packed_msg_send (VsgPackedMsg *pm, gint dst, gint tag);

void vsg_packed_msg_recv (VsgPackedMsg *pm, gint src, gint tag);

VsgPackedMsg * vsg_packed_msg_recv_new (MPI_Comm comm, gint src, gint tag);

void vsg_packed_msg_drop_buffer (VsgPackedMsg *pm);

void vsg_packed_msg_free (VsgPackedMsg *pm);

G_END_DECLS;

#endif /* __VSG_PACKEDMSG_H__ */
