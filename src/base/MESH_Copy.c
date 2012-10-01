
#include "string.h"
#include "MSTK.h"
#include "MSTK_private.h"


#ifdef __cplusplus
extern "C" {
#endif

  /* Duplicate a mesh optionally copying the attributes (with_attr) and 
     mesh sets (with_sets) */


  int MESH_Copy(Mesh_ptr srcmesh, Mesh_ptr dstmesh, int with_attr, int with_sets) {
    MVertex_ptr mv, mv_copy;
    MEdge_ptr me, me_copy;
    MFace_ptr mf, mf_copy;
    MRegion_ptr mr, mr_copy;
    MEntity_ptr ment, ment_copy;
    int i, idx, idx2, ival;
    double rval;
    void *pval;
    char funcname[256] = "MESH_Copy";

    if (srcmesh == NULL)
      MSTK_Report(funcname,"Source mesh is NULL",MSTK_FATAL);
    if (dstmesh == NULL)
      MSTK_Report(funcname,"Destination mesh is NULL",MSTK_FATAL);


    MAttrib_ptr copyatt = MAttrib_New(srcmesh,"copyatt",POINTER,MALLTYPE);

    idx = 0;
    while ((mv = MESH_Next_Vertex(srcmesh,&idx))) {
      mv_copy = MV_New(dstmesh);

      double xyz[3];
      MV_Coords(mv,xyz);
      MV_Set_Coords(mv_copy,xyz);

      MV_Set_GEntDim(mv_copy,MV_GEntDim(mv));
      MV_Set_GEntID(mv_copy,MV_GEntID(mv));
      MV_Set_ID(mv_copy,MV_ID(mv));

#ifdef MSTK_HAVE_MPI
      MV_Set_PType(mv_copy,MV_PType(mv));
      MV_Set_MasterParID(mv_copy,MV_MasterParID(mv));
      MV_Set_GlobalID(mv_copy,MV_GlobalID(mv));
#endif

      MEnt_Set_AttVal(mv,copyatt,0,0.0,mv_copy);
    }


    if (MESH_RepType(srcmesh) < R1) { /* only full representations have edges */

      idx = 0;
      while ((me = MESH_Next_Edge(srcmesh,&idx))) {
        me_copy = ME_New(dstmesh);

        MVertex_ptr ev, ev_copy;

        ev = ME_Vertex(me,0);        
        MEnt_Get_AttVal(ev,copyatt,&ival,&rval,&ev_copy);
        ME_Set_Vertex(me_copy,0,ev_copy);

        ev = ME_Vertex(me,1);
        MEnt_Get_AttVal(ev,copyatt,&ival,&rval,&ev_copy);
        ME_Set_Vertex(me_copy,1,ev_copy);

        ME_Set_GEntDim(me_copy,ME_GEntDim(me));
        ME_Set_GEntID(me_copy,ME_GEntID(me));
        ME_Set_ID(me_copy,ME_ID(me));

#ifdef MSTK_HAVE_MPI
        ME_Set_PType(me_copy,ME_PType(me));
        ME_Set_MasterParID(me_copy,ME_MasterParID(me));
        ME_Set_GlobalID(me_copy,ME_GlobalID(me));
#endif

        MEnt_Set_AttVal(me,copyatt,0,0.0,me_copy);
      }

    }


    if (MESH_RepType(srcmesh) < R1) { /* Faces in terms of edges */

      idx = 0;
      while ((mf = MESH_Next_Face(srcmesh,&idx))) {
        mf_copy = MF_New(dstmesh);
                
        List_ptr fedges = MF_Edges(mf,1,0);
        int nfe = List_Num_Entries(fedges);
        MEdge_ptr fe_copy[MAXPV2];
        int fedirs[MAXPV2];

        for (i = 0; i < nfe; i++) {
          MEdge_ptr fe = List_Entry(fedges,i);
          MEnt_Get_AttVal(fe,copyatt,&ival,&rval,&(fe_copy[i]));
          fedirs[i] = MF_EdgeDir_i(mf,i);
        }
        
        List_Delete(fedges);

        MF_Set_Edges(mf_copy,nfe,fe_copy,fedirs);

        MF_Set_GEntDim(mf_copy,MF_GEntDim(mf));
        MF_Set_GEntID(mf_copy,MF_GEntID(mf));
        MF_Set_ID(mf_copy,MF_ID(mf));

#ifdef MSTK_HAVE_MPI
        MF_Set_PType(mf_copy,MF_PType(mf));
        MF_Set_MasterParID(mf_copy,MF_MasterParID(mf));
        MF_Set_GlobalID(mf_copy,MF_GlobalID(mf));
#endif

        MEnt_Set_AttVal(mf,copyatt,0,0.0,mf_copy);
      }
        
    }
    else { /* Faces in terms of vertices */

      idx = 0;
      while ((mf = MESH_Next_Face(srcmesh,&idx))) {
        mf_copy = MF_New(dstmesh);
                
        int nfv = MF_Num_Vertices(mf);
        MVertex_ptr fv_copy[MAXPV2];
        List_ptr fverts = MF_Vertices(mf,1,0);

        for (i = 0; i < nfv; i++) {
          MVertex_ptr fv = List_Entry(fverts,i);
          MEnt_Get_AttVal(fv,copyatt,&ival,&rval,&(fv_copy[i]));
        }
        List_Delete(fverts);

        MF_Set_Vertices(mf_copy,nfv,fv_copy);

        MF_Set_GEntDim(mf_copy,MF_GEntDim(mf));
        MF_Set_GEntID(mf_copy,MF_GEntID(mf));
        MF_Set_ID(mf_copy,MF_ID(mf));

#ifdef MSTK_HAVE_MPI
        MF_Set_PType(mf_copy,MF_PType(mf));
        MF_Set_MasterParID(mf_copy,MF_MasterParID(mf));
        MF_Set_GlobalID(mf_copy,MF_GlobalID(mf));
#endif

        MEnt_Set_AttVal(mf,copyatt,0,0.0,mf_copy);
      }
        
     
      
    }

    
    if (MESH_RepType(srcmesh) != R1 || MESH_RepType(srcmesh) != R2) { /* Regions --> Faces */

      idx = 0;
      while ((mr = MESH_Next_Region(srcmesh,&idx))) {
        mr_copy = MR_New(dstmesh);

        List_ptr rfaces = MR_Faces(mr);
        int nrf = List_Num_Entries(rfaces);
        MFace_ptr rf_copy[MAXPF3];
        int rfdirs[MAXPF3];
        
        for (i = 0; i < nrf; i++) {
          MFace_ptr rf = List_Entry(rfaces,i);
          MEnt_Get_AttVal(rf,copyatt,&ival,&rval,&(rf_copy[i]));
          rfdirs[i] = MR_FaceDir_i(mr,i);
        }

        List_Delete(rfaces);

        MR_Set_Faces(mr_copy,nrf,rf_copy,rfdirs);

        MR_Set_GEntDim(mr_copy,MR_GEntDim(mr));
        MR_Set_GEntID(mr_copy,MR_GEntID(mr));
        MR_Set_ID(mr_copy,MR_ID(mr));

#ifdef MSTK_HAVE_MPI
        MR_Set_PType(mr_copy,MR_PType(mr));
        MR_Set_MasterParID(mr_copy,MR_MasterParID(mr));
        MR_Set_GlobalID(mr_copy,MR_GlobalID(mr));
#endif

        MEnt_Set_AttVal(mr,copyatt,0,0.0,mr_copy);
      }

    }
    else { /* Regions --> Vertices */

      idx = 0;
      while ((mr = MESH_Next_Region(srcmesh,&idx))) {
        mr_copy = MR_New(dstmesh);

        List_ptr rverts = MR_Vertices(mr);
        int nrv = List_Num_Entries(rverts);
        MFace_ptr rv_copy[MAXPV3];
        
        for (i = 0; i < nrv; i++) {
          MVertex_ptr rv = List_Entry(rverts,i);
          MEnt_Get_AttVal(rv,copyatt,&ival,&rval,&(rv_copy[i]));
        }
        List_Delete(rverts);

        MR_Set_Vertices(mr_copy,nrv,rv_copy,0,NULL); 

        MR_Set_GEntDim(mr_copy,MR_GEntDim(mr));
        MR_Set_GEntID(mr_copy,MR_GEntID(mr));
        MR_Set_ID(mr_copy,MR_ID(mr));

#ifdef MSTK_HAVE_MPI
        MR_Set_PType(mr_copy,MR_PType(mr));
        MR_Set_MasterParID(mr_copy,MR_MasterParID(mr));
        MR_Set_GlobalID(mr_copy,MR_GlobalID(mr));
#endif

        MEnt_Set_AttVal(mr,copyatt,0,0.0,mr_copy);
      }

    }

    

    if (with_attr) {

      MAttrib_ptr att, att_copy;

      idx = 0;
      while ((att = MESH_Next_Attrib(srcmesh,&idx))) {
        int numcomps = 0;

        char attname[256];
        MAttrib_Get_Name(att,attname);
        if (strcmp(attname,"copyatt") == 0) continue;

        MAttType atttype = MAttrib_Get_Type(att);
        MType attentdim = MAttrib_Get_EntDim(att);

        MAttrib_ptr att_copy;
        if (atttype == VECTOR || atttype == TENSOR) {
          numcomps = MAttrib_Get_NumComps(att);
          att_copy = MAttrib_New(dstmesh,attname,atttype,attentdim,numcomps);
        }
        else
          att_copy = MAttrib_New(dstmesh,attname,atttype,attentdim);

        if (attentdim == MVERTEX || attentdim == MALLTYPE) {
          idx2 = 0;
          while ((mv = MESH_Next_Vertex(srcmesh,&idx2))) {
            MEnt_Get_AttVal(mv,copyatt,&ival,&rval,&mv_copy);
            
            MEnt_Get_AttVal(mv,att,&ival,&rval,&pval);
            
            if (atttype == INT && ival) 
              MEnt_Set_AttVal(mv_copy,att_copy,ival,rval,pval);
            else if (atttype == DOUBLE && rval)
              MEnt_Set_AttVal(mv_copy,att_copy,ival,rval,pval);
            else if ((atttype == VECTOR || atttype == TENSOR) && pval) {
              double *vec = (double *) MSTK_malloc(numcomps*sizeof(double));
              memcpy(vec,pval,numcomps*sizeof(double));
              MEnt_Set_AttVal(mv_copy,att_copy,ival,rval,vec);
            }
            else {
              MSTK_Report(funcname,"Does not make sense to copy POINTER attributes from one mesh to another",MSTK_WARN);
            }
          }
        }
      
        if (attentdim == MEDGE || attentdim == MALLTYPE) {
          idx2 = 0;
          while ((me = MESH_Next_Edge(srcmesh,&idx2))) {
            MEnt_Get_AttVal(me,copyatt,&ival,&rval,&me_copy);
            
            MEnt_Get_AttVal(me,att,&ival,&rval,&pval);
            
            if (atttype == INT && ival)
              MEnt_Set_AttVal(me_copy,att_copy,ival,rval,pval);
            else if (atttype == DOUBLE && rval)
              MEnt_Set_AttVal(me_copy,att_copy,ival,rval,pval);
            else if ((atttype == VECTOR || atttype == TENSOR) && pval) {
              double *vec = (double *) MSTK_malloc(numcomps*sizeof(double));
              memcpy(vec,pval,numcomps*sizeof(double));
              MEnt_Set_AttVal(me_copy,att_copy,ival,rval,vec);
            }
            else {
              MSTK_Report(funcname,"Does not make sense to copy POINTER attributes from one mesh to another",MSTK_WARN);
            }
          }
        }
      
        if (attentdim == MFACE || attentdim == MALLTYPE) {
          idx2 = 0;
          while ((mf = MESH_Next_Face(srcmesh,&idx2))) {
            MEnt_Get_AttVal(mf,copyatt,&ival,&rval,&mf_copy);
            
            MEnt_Get_AttVal(mf,att,&ival,&rval,&pval);
            
            if (atttype == INT && ival)
              MEnt_Set_AttVal(mf_copy,att_copy,ival,rval,pval);
            else if (atttype == DOUBLE && rval)
              MEnt_Set_AttVal(mf_copy,att_copy,ival,rval,pval);
            else if ((atttype == VECTOR || atttype == TENSOR) && pval) {
              double *vec = (double *) MSTK_malloc(numcomps*sizeof(double));
              memcpy(vec,pval,numcomps*sizeof(double));
              MEnt_Set_AttVal(mf_copy,att_copy,ival,rval,vec);
            }
            else {
              MSTK_Report(funcname,"Does not make sense to copy POINTER attributes from one mesh to another",MSTK_WARN);
            }
          }
        }
      
        if (attentdim == MREGION || attentdim == MALLTYPE) {
          idx2 = 0;
          while ((mr = MESH_Next_Region(srcmesh,&idx2))) {
            MEnt_Get_AttVal(mr,copyatt,&ival,&rval,&mr_copy);
            
            MEnt_Get_AttVal(mr,att,&ival,&rval,&pval);
            
            if (atttype == INT && ival)
              MEnt_Set_AttVal(mr_copy,att_copy,ival,rval,pval);
            else if (atttype == DOUBLE && rval)
              MEnt_Set_AttVal(mr_copy,att_copy,ival,rval,pval);
            else if ((atttype == VECTOR || atttype == TENSOR) && pval) {
              double *vec = (double *) MSTK_malloc(numcomps*sizeof(double));
              memcpy(vec,pval,numcomps*sizeof(double));
              MEnt_Set_AttVal(mr_copy,att_copy,ival,rval,vec);
            }
            else {
              MSTK_Report(funcname,"Does not make sense to copy POINTER attributes from one mesh to another",MSTK_WARN);
            }
          }
        }

      } /* while (att = MESH_Next_Attrib...) */

    }




    if (with_sets) {

      MSet_ptr set;
      idx = 0;
      while ((set = MESH_Next_MSet(srcmesh,&idx))) {

        char setname[256];
        MSet_Name(set,setname);

        MType setentdim = MSet_EntDim(set);

        MSet_ptr set_copy = MSet_New(dstmesh,setname,setentdim);

        idx2 = 0;
        while ((ment = MSet_Next_Entry(set,&idx2))) {
          MEnt_Get_AttVal(ment,copyatt,&ival,&rval,&ment_copy);
          MSet_Add(set_copy,ment_copy);
        }

      } /* while (set = MESH_Next_MSet...) */

    }




    idx = 0;
    while ((mv = MESH_Next_Vertex(srcmesh,&idx)))
      MEnt_Rem_AttVal(mv,copyatt);
      
    idx = 0;
    while ((me = MESH_Next_Edge(srcmesh,&idx)))
      MEnt_Rem_AttVal(me,copyatt);

    idx = 0;
    while ((mf = MESH_Next_Face(srcmesh,&idx)))
      MEnt_Rem_AttVal(mf,copyatt);
      
    idx = 0;
    while ((mr = MESH_Next_Region(srcmesh,&idx)))
      MEnt_Rem_AttVal(mr,copyatt);

    MAttrib_Delete(copyatt);

    return 1;
  }

#ifdef __cplusplus
}
#endif
