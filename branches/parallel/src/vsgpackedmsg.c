#include "vsgpackedmsg.h"

/**
 * VsgPackedMsg:
 * @communicator: the MPI comunicator on which the #VsgPackedMsg is to be used.
 *
 * Stores heterogenous MPI message data.
 */

/**
 * vsg_packed_msg_new:
 * @comm: a MPI communicator.
 *
 * Creates a new packed message buffer for use inside @comm.
 *
 * Returns: newly allocated #VsgPackedMsg.
 */
VsgPackedMsg *vsg_packed_msg_new (MPI_Comm comm)
{
  VsgPackedMsg *ret;

  ret = g_malloc (sizeof (VsgPackedMsg));

  vsg_packed_msg_init (ret, comm);

  return ret;
}

/**
 * vsg_packed_msg_init:
 * @pm: a #VsgPackedMsg.
 * @comm: a MPI communicator.
 *
 * Initializes @pm for holding @comm messages. Will free all previously stored
 * message data.
 */
void vsg_packed_msg_init (VsgPackedMsg *pm, MPI_Comm comm)
{
  g_return_if_fail (pm != NULL);

  pm->communicator = comm;

  vsg_packed_msg_drop_buffer (pm);
}

/**
 * vsg_packed_msg_send_append:
 * @pm: a #VsgPackedMsg.
 * @buf: pointer to the beginning of data to be stored.
 * @count: number of @type data to store.
 * @type: type of the data to be stored.
 *
 * Appends @count instances of @type data to the message buffer.
 */
void vsg_packed_msg_send_append (VsgPackedMsg *pm, gpointer buf,
                                 gint count, MPI_Datatype type)
{
  gint pos, size, addsize;
  gint ierr;

  g_return_if_fail (pm != NULL);

  pos = pm->position;
  size = pm->size;

  /* compute size of this new message part */
  MPI_Pack_size (count, type, pm->communicator, &addsize);

  /* allocate enough memory in message msg to store this message */
  if ((addsize + pos) > size)
    {
      size = MAX (size + 1024, addsize+pos);

      pm->buffer = g_realloc (pm->buffer, size * sizeof (char));

      pm->size = size;
    }

  ierr = MPI_Pack (buf, count, type, pm->buffer, size,
		   &pm->position, pm->communicator);

  if (ierr != MPI_SUCCESS) vsg_mpi_error_output (ierr);

}

/**
 * vsg_packed_msg_recv_read:
 * @pm: a #VsgPackedMsg.
 * @buf: pointer to the beginning of data to be read.
 * @count: number of @type data to read.
 * @type: type of the data to be read.
 *
 * Reads @count instance of @type data from the current position in the buffer.
 * The position then is updated to the first byte after read data.
 */
void vsg_packed_msg_recv_read (VsgPackedMsg *pm, gpointer buf,
                               gint count, MPI_Datatype type)
{
  gint size;
  gint ierr;

  g_return_if_fail (pm != NULL);

  size = pm->size;

/*   { */
/*     gint mytid; */
/*     MPI_Comm_rank (pm->communicator, &mytid); */
/*     g_printerr ("%d: unpacking %d data at %d\n", mytid, count, pm->position); */
/*   } */

  ierr = MPI_Unpack (pm->buffer, size, &pm->position, 
		     buf, count, type,  pm->communicator);
  
  if (ierr != MPI_SUCCESS) vsg_mpi_error_output (ierr);

}

/**
 * vsg_packed_msg_send:
 * @pm: a #VsgPackedMsg.
 * @dst: the destination task id.
 * @tag: an integer message tag.
 *
 * Sends stored message to the specified destination with the specified tag.
 */
void vsg_packed_msg_send (VsgPackedMsg *pm, gint dst, gint tag)
{
  gint ierr;

  ierr = MPI_Send (pm->buffer, pm->position, MPI_PACKED, dst, tag,
                   pm->communicator);

  if (ierr != MPI_SUCCESS) vsg_mpi_error_output (ierr);
}

/**
 * vsg_packed_msg_recv:
 * @pm: a #VsgPackedMsg.
 * @src: the source task id. Can be %MPI_ANY_SOURCE.
 * @tag: an integer message tag. Can be %MPI_ANY_TAG.
 *
 * Receives a message from source @src with @tag message tag and stores it in
 * @pm. any previously stored data will be lost.
 */
void vsg_packed_msg_recv (VsgPackedMsg *pm, gint src, gint tag)
{
  MPI_Status status;
  gint ierr;
  gint rsize = 0;

  MPI_Probe (src, tag, pm->communicator, &status);

  MPI_Get_count (&status, MPI_PACKED, &rsize);

  pm->buffer = g_realloc (pm->buffer, rsize);
  pm->size = rsize;

/*   { */
/*     gint mytid; */
/*     MPI_Comm_rank (pm->communicator, &mytid); */
/*     g_printerr ("%d: received %d bytes from %d\n", mytid, rsize, src); */
/*   } */

  ierr = MPI_Recv (pm->buffer, rsize, MPI_PACKED, src, tag,
                   pm->communicator, &status);

  pm->position = 0;

  if (ierr != MPI_SUCCESS) vsg_mpi_error_output (ierr);
}

/**
 * vsg_packed_msg_recv_new:
 * @comm: a MPI communicator
 * @src: the source task id. Can be %MPI_ANY_SOURCE.
 * @tag: an integer message tag. Can be %MPI_ANY_TAG.
 *
 * Receives a message from source @src with @tag message tag and stores it in
 * a newly allocated #VsgPackedMsg.
 *
 * Returns: the received message.
 */
VsgPackedMsg * vsg_packed_msg_recv_new (MPI_Comm comm, gint src, gint tag)
{
  VsgPackedMsg *ret = vsg_packed_msg_new (comm);

  vsg_packed_msg_recv (ret, src, tag);

  return ret;
}

/**
 * vsg_packed_msg_drop_buffer:
 * @pm: a #VsgPackedMsg.
 *
 * Drops data stored in @pm buffer.
 */
void vsg_packed_msg_drop_buffer (VsgPackedMsg *pm)
{
  g_return_if_fail (pm != NULL);

  if (pm->buffer != NULL)
    g_free (pm->buffer);

  pm->buffer = NULL;
  pm->position = 0;
  pm->size = 0;
}

/**
 * vsg_packed_msg_free:
 * @pm: a #VsgPackedMsg.
 *
 * deallocates @pm and all data associated with it.
 */
void vsg_packed_msg_free (VsgPackedMsg *pm)
{
  g_return_if_fail (pm != NULL);

  vsg_packed_msg_drop_buffer (pm);

  g_free (pm);
}
