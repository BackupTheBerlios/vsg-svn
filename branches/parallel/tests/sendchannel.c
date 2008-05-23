#include "vsg-config.h"

#include "vsg/vsg.h"

#include "vsg/vsgpackedmsg.h"

gint main (gint argc, gchar ** argv)
{
  gint ret = 0;
  gint mytid, numtasks;

  vsg_init_gdouble ();

  MPI_Init (&argc, &argv);

  MPI_Comm_size (MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank (MPI_COMM_WORLD, &mytid);

  if (argc > 1 && g_strncasecmp (argv[1], "--version", 9) == 0)
    {
      if (mytid == 0)
	g_print ("%s\n", PACKAGE_VERSION);
      return 0;
    }
 
/*   g_printerr ("%d %d\n", mytid, numtasks); */

  if (mytid == 0)
    {
      gint i;
      VsgSendChannel *sc;

      sc = vsg_send_channel_new (MPI_COMM_WORLD, 1, 123, 2*sizeof (gint));

      i = 1;
      vsg_send_channel_append (sc, &i, 1, MPI_INT);
      i = 2;
      vsg_send_channel_append (sc, &i, 1, MPI_INT);

      if (sc->pm.position != 0) g_error ("%d: message not sent after second "
                                         "integer as it should be\n", mytid);

      i = 3;
      vsg_send_channel_append (sc, &i, 1, MPI_INT);

      vsg_send_channel_flush (sc);

      vsg_send_channel_free (sc);

    }

  if (mytid == 1)
    {
      gint i;
      VsgPackedMsg pm = VSG_PACKED_MSG_STATIC_INIT (MPI_COMM_WORLD);

      vsg_packed_msg_recv (&pm, 0, 123);

      vsg_packed_msg_recv_read (&pm, &i, 1, MPI_INT);
      if (i != 1) g_error ("%d: first int should be %d (unpacked %d)\n", mytid,
                           1, i);

      vsg_packed_msg_recv_read (&pm, &i, 1, MPI_INT);
      if (i != 2) g_error ("%d: second int should be %d (unpacked %d)\n",
                           mytid, 2, i);

      vsg_packed_msg_recv (&pm, 0, 123);

      vsg_packed_msg_recv_read (&pm, &i, 1, MPI_INT);
      if (i != 3) g_error ("%d: third int should be %d (unpacked %d)\n", mytid,
                           3, i);

      vsg_packed_msg_drop_buffer (&pm);

    }

  MPI_Finalize ();

  return ret;
}
