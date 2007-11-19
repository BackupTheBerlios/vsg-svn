#include "vsg-config.h"

#include "vsg/vsgcommbuffer.h"

gint main (gint argc, gchar ** argv)
{
  gint ret = 0;
  VsgCommBuffer *cb;
  gint mytid, numtasks;
  gint i;

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

  cb = vsg_comm_buffer_new (MPI_COMM_WORLD);

  for (i=0; i<numtasks; i++)
    {
      gchar *a="a";

      /* proc #2 won't receive any data: empty message */
      if (i == 2) continue;

      vsg_comm_buffer_send_append (cb, i, &i, 1, MPI_INT);
      vsg_comm_buffer_send_append (cb, i, a, 2, MPI_CHAR);
    }

  vsg_comm_buffer_share (cb);

  if (mytid != 2)
    for (i=0; i<numtasks; i++)
      {
	gchar ra[2] = "";
	gint ri = 0;

	vsg_comm_buffer_recv_read (cb, i, &ri, 1, MPI_INT);
	vsg_comm_buffer_recv_read (cb, i, ra, 2, MPI_CHAR);

	if (ri != mytid)
	  {
	    g_printerr ("%d: integer msg (from %d) error %d (should be %d).\n",
			mytid, i, ri, mytid);
	  }

	if (g_strncasecmp (ra, "a", 1) != 0)
	  {
	    g_printerr ("%d: character msg (from %d) error %s (should be %s).\n",
			mytid, i, ra, "a");
	  }
      }

  vsg_comm_buffer_free (cb);

  MPI_Finalize ();

  return ret;
}
