/* Import a serial Exodus II mesh (exo) and attach real value attributes
   (att) read from an auxilliary file. Export a new Exodus II mesh
   file or partitioned Nemesis I files with the new attributes attached */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#ifdef _MSTK_HAVE_MPI
#include <mpi.h>
#endif

#include "MSTK.h" 

char progname[16] = "exoatt";
void import_attributes(Mesh_ptr mesh, char *attfname);
void import_sets(Mesh_ptr mesh, char *elsetfname);


int main(int argc, char *argv[]) {
  char infname[256], outfname[256], attfname[256], setfname[256];
  Mesh_ptr mesh;
  int len, ok=0;
  int build_classfn=1, partition=0, weave=0, use_geometry=0, parallel_check=0;
  int num_ghost_layers=0, partmethod=0, attfile_given=0, setfile_given=0;
  MshFmt inmeshformat, outmeshformat;

  if (argc < 3) {
    fprintf(stderr,"\n");
    fprintf(stderr,"usage: mpirun -n NP %s <--partition=y|1|n|0> <--partition-method=0|1|2> <--parallel-check=y|1|n|0> <--attfile=attfilename> <--setfile=setfilename> infilename outfilename\n\n",progname);
    fprintf(stderr,"partition-method = 0, METIS\n");
    fprintf(stderr,"                 = 1, ZOLTAN with GRAPH partioning\n");
    fprintf(stderr,"                 = 2, ZOLTAN with RCB partitioning\n");
    fprintf(stderr,"Choose 2 if you want to avoid partitioning models\n");
    fprintf(stderr,"with high aspect ratio along the short directions\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"infilename   = Input Exodus II File\n");
    fprintf(stderr,"outfilename  = Output Exodus II file name\n");
    fprintf(stderr,"attfilename  = Auxiliary file with real valued attributes\n");
    fprintf(stderr,"setfilename  = Auxiliary file with entity set specifications\n");
    fprintf(stderr,"\n");
    exit(-1);
  }

#ifdef MSTK_HAVE_MPI
  MPI_Init(&argc,&argv);
#endif


  MSTK_Init();


  int rank=0, numprocs=1;
#ifdef MSTK_HAVE_MPI

  MSTK_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&numprocs);

#else
  MSTK_Comm comm = NULL;
#endif
  

  if (argc > 3) {
    int i;
    for (i = 1; i < argc-2; i++) {
      if (strncmp(argv[i],"--partition=",12) == 0) {
        if (strncmp(argv[i]+12,"y",1) == 0 ||
            strncmp(argv[i]+12,"1",1) == 0)
          partition = 1;
        else if (strncmp(argv[i]+12,"n",1) == 0 ||
                 strncmp(argv[i]+12,"0",1) == 0) 
          partition = 0;
        else
          MSTK_Report(progname,
                      "--partition option should be y, n, 1 or 0",
                      MSTK_FATAL);          
      }
      else if (strncmp(argv[i],"--partition-method",18) == 0 ||
               strncmp(argv[i],"--partition_method",18) == 0) {
        sscanf(argv[i]+19,"%d",&partmethod);
        partition = 1;
      }
      else if (strncmp(argv[i],"--parallel-check",16) == 0 ||
               strncmp(argv[i],"--parallel_check",16) == 0) {
        if (strncmp(argv[i]+17,"y",1) == 0 ||
            strncmp(argv[i]+17,"1",1) == 1)
          parallel_check = 1;
        else if (strncmp(argv[i]+17,"n",13) == 0 ||
                 strncmp(argv[i]+17,"0",13) == 0) 
          parallel_check = 0;
      }
      else if (strncmp(argv[i],"--attfile",9) == 0) {
        sscanf(argv[i]+10,"%s",attfname);
        attfile_given=1;
      }
      else if (strncmp(argv[i],"--setfile",9) == 0) {
        sscanf(argv[i]+10,"%s",setfname);
        setfile_given=1;
      }
      else
        fprintf(stderr,"Unrecognized option...Ignoring\n");
    }

    /* If running on multiple processors, but partition is 0, assume
       that you are partitioning a serial mesh */

    if (numprocs > 1 && partition == 0) {
      MSTK_Report(progname,
                  "Running on multiple processors but partition option not specified",
                  MSTK_WARN);
      MSTK_Report(progname,
                  "Output mesh will be partitioned into as many parts as processors",
                  MSTK_WARN);
      partition = 1; 
    }
  }

  strcpy(infname,argv[argc-2]);
  strcpy(outfname,argv[argc-1]);


  len = strlen(infname);
  if (len > 4 && strncmp(&(infname[len-4]),".exo",4) == 0)
    inmeshformat = EXODUSII;
  else
    MSTK_Report(progname,"Input file must have .exo extension", MSTK_FATAL);

  len = strlen(outfname);
  if (len > 4 && strncmp(&(outfname[len-4]),".exo",4) == 0)
    outmeshformat = EXODUSII;
  else if (len > 4 && strncmp(&(outfname[len-4]),".par",4) == 0)
    outmeshformat = EXODUSII;
  else
    MSTK_Report(progname,"Output file must have .exo or .par extension", MSTK_FATAL);



  /* Import the Exodus II mesh */

  mesh = MESH_New(F1);

#ifdef MSTK_HAVE_MPI
  if (rank == 0) {
#endif

    int opts[5]={0,0,0,0,0};
  
    ok = 0;
    fprintf(stderr,"Importing mesh from ExodusII file...");
    opts[0] = 0; /* don't partition while importing - do it later */
    opts[1] = 0;
    opts[2] = 0; /* no ghost layers */  
    opts[3] = 0;
    ok = MESH_ImportFromFile(mesh,infname,"exodusii",opts,comm);
    if (ok)
      fprintf(stderr,"done\n\n");
    else {
      MSTK_Report(progname,"Failed\n",MSTK_FATAL);
    }
    

    /* Read in the attributes only on rank 0 and attach it to the mesh */
    if (attfile_given)
      import_attributes(mesh, attfname);


    /* Import element set definitions and create in the mesh */
    if (setfile_given)
      import_sets(mesh, setfname);

#ifdef MSTK_HAVE_MPI
  }
#endif


#ifdef MSTK_HAVE_MPI

  if (partition) {

    /* Need to make sure that we have a mesh only on processor 0 */
    /* For now we cannot repartition                             */

    int prepartitioned = 0, mesh_present = 0;

    int dim = 3;
    if (rank > 0) {
      if (mesh && MESH_Num_Vertices(mesh) != 0)
        mesh_present = 1;
    }
    
    MPI_Reduce(&mesh_present,&prepartitioned,1,MPI_INT,MPI_MAX,0,comm);
    MPI_Bcast(&prepartitioned,1,MPI_INT,0,comm);

    if (!prepartitioned) {
    
      Mesh_ptr mesh0=NULL;
      if (rank == 0) {
        fprintf(stderr,"Partitioning mesh into %d parts...",numprocs);
        
        if (MESH_Num_Regions(mesh))
          dim = 3;
        else if (MESH_Num_Faces(mesh))
          dim = 2;
        else {
          fprintf(stderr,"Partitioning possible only for 2D or 3D meshes\n");
          exit(-1);
        }
        
        mesh0 = mesh;
        mesh = NULL;
      }
      
      int ring = 0; /* No ghost ring of elements */
      int with_attr = 1; /* Do allow exchange of attributes */
      int del_inmesh = 1; /* Delete input mesh after partitioning */

      ok = MSTK_Mesh_Distribute(mesh0, &mesh, &dim, ring, with_attr,
                                    partmethod, del_inmesh, comm);
      
      if (rank == 0) {
        if (ok)
          fprintf(stderr,"done\n");
        else {
          fprintf(stderr,"failed\n");
          exit(-1);
        }
      }

      /*      if (rank == 0)
	      MESH_Delete(mesh0);
      */

    }    
  }

#endif /* MSTK_HAVE_MPI */


#ifdef MSTK_HAVE_MPI
  if (numprocs > 1 && parallel_check == 1) {

    /* Do a parallel consistency check too */

    ok = ok && MESH_Parallel_Check(mesh,comm);
    
  }
#endif

  fprintf(stderr,"\nExporting mesh to ExodusII format...");
  ok = MESH_ExportToFile(mesh,outfname,"exodusii",-1,NULL,NULL,comm);

  if (ok)
    fprintf(stderr,"Done\n");
  else {
    fprintf(stderr,"Failed\n");
    exit(-1);
  }


  MESH_Delete(mesh);
 
  // while (1);

#ifdef MSTK_HAVE_MPI
  MPI_Finalize();
#endif

  return 1;
}







void import_attributes(Mesh_ptr mesh, char *attfname) {
  FILE *fp;
  int i, j, idx, status, entid, natt, nent, celldim;
  char mesg[256];
  MRegion_ptr mr;
  MFace_ptr mf;
  MVertex_ptr mv;
  char attname[256], typestr[256];
  double attval;
  MAttrib_ptr att;
  

  fp = fopen(attfname,"r");
  if (!fp) {
    sprintf(mesg,"Could not open attribute file %s",attfname);
    MSTK_Report(progname,mesg,MSTK_FATAL);
  }

  status = fscanf(fp,"%d",&natt);
  if (status == EOF)
    MSTK_Report(progname,"Attributes file ended prematurely",MSTK_FATAL);

  celldim = MESH_Num_Regions(mesh) ? 3 : 2;
  
  for (i = 0; i < natt; ++i) {

    status = fscanf(fp,"%s %s %d\n",attname,typestr,&nent);

    fprintf(stderr,"\nCreating %s attribute %s\n",typestr,attname);

    if (strcasecmp(typestr,"NODE") == 0) {
      att = MAttrib_New(mesh,attname,DOUBLE,MVERTEX);

      if (nent == -1) { /* data specified on the whole mesh */
        idx = 0;
        while ((mv = MESH_Next_Vertex(mesh,&idx))) {
          status = fscanf(fp,"%lf",&attval);
          if (status == EOF)
            MSTK_Report(progname,"Attributes file ended prematurely",MSTK_FATAL);
          
          MEnt_Set_AttVal(mv,att,0,attval,NULL);
        }
      }
      else {
        for (j = 0; j < nent; ++j) {
          status = fscanf(fp,"%d %lf",&entid,&attval);
          if (status == EOF)
            MSTK_Report(progname,"Attributes file ended prematurely",MSTK_FATAL);
          
          mv = MESH_VertexFromID(mesh,entid);
          if (!mv) {
            sprintf(mesg,"Cannot find mesh element with ID = %-d",entid);
            MSTK_Report(progname,mesg,MSTK_FATAL);
          }
            
          MEnt_Set_AttVal(mv,att,0,attval,NULL);
        }
      }
    }
    else if (strcasecmp(typestr,"SIDESET") == 0) {
      MSTK_Report(progname,"SIDESET attribute import not implemented",MSTK_FATAL);
    }
    else if (strcasecmp(typestr,"CELL") == 0) {
      if (celldim == 3) {

        att = MAttrib_New(mesh,attname,DOUBLE,MREGION);

        if (nent == -1) { /* data specified on the whole mesh */
          idx = 0;
          while ((mr = MESH_Next_Region(mesh,&idx))) {
            status = fscanf(fp,"%lf",&attval);
            if (status == EOF)
              MSTK_Report(progname,"Attributes file ended prematurely",MSTK_FATAL);

            MEnt_Set_AttVal(mr,att,0,attval,NULL);
          }
        }
        else {
          for (j = 0; j < nent; ++j) {
            status = fscanf(fp,"%d %lf",&entid,&attval);
            if (status == EOF)
              MSTK_Report(progname,"Attributes file ended prematurely",MSTK_FATAL);

            mr = MESH_RegionFromID(mesh,entid);
            if (!mr) {
              sprintf(mesg,"Cannot find mesh element with ID = %-d",entid);
              MSTK_Report(progname,mesg,MSTK_FATAL);
            }
            
            MEnt_Set_AttVal(mr,att,0,attval,NULL);
          }
        }
      }
      else {

        att = MAttrib_New(mesh,attname,DOUBLE,MFACE);

        if (nent == -1) { /* data specified on the whole mesh */
          idx = 0;
          while ((mf = MESH_Next_Face(mesh,&idx))) {
            status = fscanf(fp,"%lf",&attval);
            if (status == EOF)
              MSTK_Report(progname,"Attributes file ended prematurely",MSTK_FATAL);

            MEnt_Set_AttVal(mf,att,0,attval,NULL);
          }
        }
        else {
          for (j = 0; j < nent; ++j) {
            status = fscanf(fp,"%d %lf",&entid,&attval);
            if (status == EOF)
              MSTK_Report(progname,"Attributes file ended prematurely",MSTK_FATAL);

            mf = MESH_FaceFromID(mesh,entid);
            if (!mf) {
              sprintf(mesg,"Cannot find mesh element with ID = %-d",entid);
              MSTK_Report(progname,mesg,MSTK_FATAL);
            }
            
            MEnt_Set_AttVal(mf,att,0,attval,NULL);
          }
        }        
      }
    }
    else
      MSTK_Report(progname,
                  "Unknown entity type for attribute: Valid types are NODE, CELL or SIDESET (case does not matter)",
                  MSTK_FATAL);
  }

}


// Import entity set specifications from auxiliary files and create
// them in the MSTK mesh so that the output Exodus II mesh file will
// contain these sets

void import_sets(Mesh_ptr mesh, char *setfname) {
  FILE *fp;
  int i, j, k, l, m;
  int idx, status, setid, entid, nset, nent, celldim;
  char mesg[256];
  MRegion_ptr mr;
  MFace_ptr mf;
  MEdge_ptr me;
  MVertex_ptr mv;
  char setname[256], typestr[256];
  MSet_ptr mset;
  
  fp = fopen(setfname,"r");
  if (!fp) {
    sprintf(mesg,"Could not open entity set file %s",setfname);
    MSTK_Report(progname,mesg,MSTK_FATAL);
  }

  status = fscanf(fp,"%d",&nset);
  if (status == EOF)
    MSTK_Report(progname,"Entity sets file ended prematurely",MSTK_FATAL);

  celldim = MESH_Num_Regions(mesh) ? 3 : 2;
  
  for (i = 0; i < nset; ++i) {

    status = fscanf(fp,"%d %s %d\n",&setid,typestr,&nent);
    
    if (strcasecmp(typestr,"NODE") == 0) {

      sprintf(setname,"nodeset_%-d",setid);
      fprintf(stderr,"Creating %s set %s\n",typestr,setname);

      mset = MSet_New(mesh,setname,MVERTEX);

      if (nent <= 0)
        MSTK_Report(progname,"Number of entities in set must be greater than zero",MSTK_FATAL);
      for (j = 0; j < nent; ++j) {
        status = fscanf(fp,"%d",&entid);
        if (status == EOF)
          MSTK_Report(progname,"Entity sets file ended prematurely",MSTK_FATAL);
        
        mv = MESH_VertexFromID(mesh,entid);
        if (!mv) {
          sprintf(mesg,"Cannot find mesh vertex with ID = %-d",entid);
          MSTK_Report(progname,mesg,MSTK_FATAL);
        }
            
        MSet_Add(mset,mv);
      }
    }
    else if (strcasecmp(typestr,"FACE") == 0) {
      int lfnum;

      sprintf(setname,"sideset_%-d",setid);
      fprintf(stderr,"Creating %s set %s\n",typestr,setname);

      if (celldim == 3) {
        mset = MSet_New(mesh,setname,MFACE);
        
        if (nent <= 0)
          MSTK_Report(progname,"Number of entities in set must be greater than zero",MSTK_FATAL);

        for (j = 0; j < nent; ++j) {

          int num_face_verts;
          int face_verts[MAXPV2];

          status = fscanf(fp,"%d",&num_face_verts);
          if (status == EOF)
            MSTK_Report(progname,"Entity sets file ended prematurely",MSTK_FATAL);

          for (k = 0; k < num_face_verts; k++) 
            status = fscanf(fp,"%d",&(face_verts[k]));
          
          MVertex_ptr fverts[MAXPV2];
          for (k = 0; k < num_face_verts; k++)
            fverts[k] = MESH_VertexFromID(mesh,face_verts[0]);

          mf = MVs_CommonFace(num_face_verts,fverts);
          if (mf)
            MSet_Add(mset,mf);
          else {
            MSTK_Report(progname,"Sideset face not found",MSTK_WARN);
            fprintf(stderr,"Did not find common face for verts ");
            for (k = 0; k < num_face_verts; k++)
              fprintf(stderr,"%-d",face_verts[k]);
          }
        }
      }
      else { /* assume celldim == 2 */
        mset = MSet_New(mesh,setname,MEDGE);
        
        if (nent <= 0)
          MSTK_Report(progname,"Number of entities in set must be non-zero",MSTK_FATAL);

        int edge_verts[2];
        for (j = 0; j < nent; ++j) {
          status = fscanf(fp,"%d %d",&(edge_verts[0]),&(edge_verts[1]));
          if (status == EOF)
            MSTK_Report(progname,"Entity sets file ended prematurely",MSTK_FATAL);

          MVertex_ptr ev[2];
          ev[0] = MESH_VertexFromID(mesh,edge_verts[0]);
          ev[1] = MESH_VertexFromID(mesh,edge_verts[1]);          
          me = MVs_CommonEdge(ev[0],ev[1]);
          if (me) 
            MSet_Add(mset,me);
          else {
            MSTK_Report(progname,"Sideset edge not found",MSTK_WARN);
            fprintf(stderr,"Did not find common edge between vertices %-d %-d",
                    edge_verts[0],edge_verts[1]);
          }
        }
      }
    }
    else if (strcasecmp(typestr,"CELL") == 0) {

      sprintf(setname,"elemset_%-d",setid);
      fprintf(stderr,"Creating %s set %s\n",typestr,setname);

      if (celldim == 3) {
        mset = MSet_New(mesh,setname,MREGION);

        if (nent <= 0)
          MSTK_Report(progname,"Number of entities in set must be non-zero",MSTK_FATAL);
        for (j = 0; j < nent; ++j) {
          status = fscanf(fp,"%d",&entid);
          if (status == EOF)
            MSTK_Report(progname,"Entity sets file ended prematurely",MSTK_FATAL);
          
          mr = MESH_RegionFromID(mesh,entid);
          if (!mr) {
            sprintf(mesg,"Cannot find mesh region with ID = %-d",entid);
            MSTK_Report(progname,mesg,MSTK_FATAL);
          }
            
          MSet_Add(mset,mr);
        }
      }
      else { /* assume celldim == 2 */
        mset = MSet_New(mesh,setname,MFACE);

        if (nent <= 0)
          MSTK_Report(progname,"Number of entities in set must be non-zero",MSTK_FATAL);
        for (j = 0; j < nent; ++j) {
          status = fscanf(fp,"%d",&entid);
          if (status == EOF)
            MSTK_Report(progname,"Entity sets file ended prematurely",MSTK_FATAL);
          
          mf = MESH_FaceFromID(mesh,entid);
          if (!mf) {
            sprintf(mesg,"Cannot find mesh face with ID = %-d",entid);
            MSTK_Report(progname,mesg,MSTK_FATAL);
          }
            
          MSet_Add(mset,mf);
        }
      }
    }
    else
      MSTK_Report(progname,
                  "Unknown entity type for set: Valid types are NODE, FACE or CELL",
                  MSTK_FATAL);
  }

  fclose(fp);

} /* import_sets */


