#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MSTK_private.h"
#include "MSTK.h"



void Graph_Build_LevelStructure(int nnodes, int adj[], int nadj[], 
                                int offset[], int start_node, 
                                int nodemap[], int level[], 
                                int *depth, int *maxwidth);
void Graph_Find_ExtremalNodes(int nnodes, int adj[], int nadj[], int offset[], 
                              int *n_beg, int *n_end);
void Graph_Assign_Levels_MinWidth(int nnodes, int adj[], int nadj[], 
                                  int offset[], int n_1, int n_k, 
                                  int level[], int *depth, int *maxwidth);
void Graph_Assign_Numbers(int nnodes, int adj[], int nadj[], int offset[], 
                          int n_beg, int n_end, int level[], int nodemap[]);


/* 
   Renumber nodes of a graph according to the Gibbs-Pool-Stockmeyer algorithm 
   with George-Liu modifications to find the pseudo-peripheral nodes

   nnodes   (in)      -- Number of nodes in the graph
   nstart   (in)      -- Start node (optional - can be 0 if not known)
   nadj     (in)      -- number of neighbors for each node
   adj      (in)      -- array containing the neighbors of each node
   nodemap  (out)     -- new numbers for the nodes (starting from 0)
   depth    (out)     -- Number of levels of nodes
   maxwidth (out)     -- maxwidth of any level
*/

int Graph_Renumber_GPS(int nnodes, int nstart, int *nadj, int *adj, 
                       int *nodemap, int *depth, int *maxwidth) {
  int i, j, k, d, l, n_1, n_k, nnew, ncset;
  int curlev, curnode;
  int curdepth, curmaxwidth, minmaxwidth, newmaxwidth, newdepth;
  int nlastlev, mindeg, maxdeg, ntmp;
  int *offset;
  int *level;

  /* Compute the min and max degree of the mesh */

  mindeg = 1e+9;
  maxdeg = 0;
  for (i = 0; i < nnodes; i++) {
    mindeg = (nadj[i] < mindeg) ? nadj[i] : mindeg;
    maxdeg = (nadj[i] > maxdeg) ? nadj[i] : maxdeg;
  }

  /* Compute the offsets for each node into the adjacency array based
     on the number of neighbors of each node */

  offset = (int *) malloc(nnodes*sizeof(int));
  offset[0] = 0;
  for (i = 1; i < nnodes; i++)
    offset[i] = offset[i-1] + nadj[i-1];

  n_1 = nstart;




  /* Step I: Finding two extremal nodes representing the pseudo-diameter of
     the problem */

  Graph_Find_ExtremalNodes(nnodes, adj, nadj, offset, &n_1, &n_k);

  /*
  fprintf(stderr,"Extremal nodes found %-d  %-d\n",n_1,n_k);
  */


  /* Step II: assigning final level numbers to minimize level width */

  level =  (int *) malloc(nnodes*sizeof(int));
  Graph_Assign_Levels_MinWidth(nnodes, adj, nadj, offset, n_1, n_k, level, 
                               depth, maxwidth);
   
  /* Step III: Numbering */

  
  Graph_Assign_Numbers(nnodes, adj, nadj, offset, n_1, n_k, level, nodemap);


  free(offset);
  free(level);

  return 1;
}



/* Given a graph and a starting node, build a level structure and
   return some information about it.

   nnodes   (in)    --- Number of nodes
   adj      (in)    --- Array listing adjacent nodes of nodes
   nadj     (in)    --- Number of adjacent nodes for each node
   offset   (in)    --- Starting location or offset in 'adj' array where the 
                        adjacency information for each node starts. 
   n_beg    (in)    --- Starting node for building the level structure
   nodemap  (out)   --- Mapping of the nodes to new numbers (or in other words, 
                        the order they are added into the level structure)
   level    (out)   --- Level of each node
   depth    (out)   --- Depth of the level structure (number of levels)
   maxwidth (out)   --- Maximum number of nodes in any one level
*/

void Graph_Build_LevelStructure(int nnodes, int adj[], int nadj[], 
                                int offset[], int start_node, 
                                int nodemap[], int level[], 
                                int *depth, int *maxwidth) {
  int i, j;
  int curlev, curnode, nnew, lastnode;

  for (i = 0; i < nnodes; i++) level[i] = 0;

  /* Start building a level structure */
      
  curlev = 1;
  curnode = start_node;
  nodemap[0] = curnode;
  level[curnode] = curlev;
  nnew = 1;
    
  for (i = 0; i < nnodes; i++) {

    /* Visit the next node */
        
    curnode = nodemap[i];
    curlev = level[curnode];

    /* Visit its neighbors */

    for (j = 0; j < nadj[curnode]; j++) {
      int off = offset[curnode];
      int adjnode = adj[off+j];
      
      if (level[adjnode] != 0) continue; /* already visited */
      
      /* If the neighbor has not been visited and tagged with a level
         already, tag it with the next level and add it to the 
         new list */
      
      level[adjnode] = curlev+1;
      nodemap[nnew++] = adjnode;
    }
    
  }

  if (nnew != nnodes) {
    fprintf(stderr,"Error in algorithm - did not traverse all nodes (nnodes %-d, nnew %-d)\n",nnodes,nnew);
    exit(-1);
  }

  lastnode = nodemap[nnodes-1];
  *depth = level[lastnode];
  
  *maxwidth = 0;
  for (i = 1; i <= *depth; i++) {
    int width = 0;
    for (j = 0; j < nnodes; j++)
      width = (level[j] == i) ? width+1 : width;
    *maxwidth = (width > *maxwidth) ? width : *maxwidth;
  }
    
}
 

/* Find the extremal nodes of a pseudo diameter of a graph - These
   nodes will result in a maximum depth level structure. This routine
   uses a variation of George and Liu's suggestion for shrinking the
   last level from which to do the peripheral search. This cuts down
   the search time and gives good enough results (for a structured
   mesh, the full GPS algorithm does not give the results we want
   anyway) 

   nnodes (in)     --- Number of nodes
   adj    (in)     --- Array listing adjacent nodes of nodes
   nadj   (in)     --- Number of adjacent nodes for each node
   offset (in)     --- Starting location or offset in 'adj' array where the 
                       adjacency information for each node starts. 
   n_beg  (in/out) --- Guess and final result for first extremal node of pseudo-diameter
   n_end  (out)    --- Second extremal node of pseudo-diameter
*/


void Graph_Find_ExtremalNodes(int nnodes, int adj[], int nadj[], int offset[], 
                              int *n_beg, int *n_end) {
  int i, j, l, outerdone, nlastlev, ntmp, restart;
  int curdepth, curmaxwidth, newdepth, newmaxwidth, minmaxwidth;
  int mindeg_last, maxdeg_last;
  int *level, *nodemap, *lastlevnodes;
  

  level = (int *) malloc(nnodes*sizeof(int));
  nodemap = (int *) malloc(nnodes*sizeof(int));
  

  /* Step I: Finding two extremal nodes representing the pseudo-diameter of
     the problem */

  outerdone = 0;
  while (!outerdone) {

    Graph_Build_LevelStructure(nnodes, adj, nadj, offset, *n_beg, 
                               nodemap, level, &curdepth, &curmaxwidth);

    /*
    fprintf(stderr,"Found level structure with starting at node %-d \n",*n_beg);
    fprintf(stderr,"Level structure has %-d levels and maximum width of %-d\n",curdepth,curmaxwidth);
    */


    /* Generate level structures from the nodes of the deepest level
       sorted according to increasing degree */

    /* Count the number of nodes on the last level */

    nlastlev = 0;
    for (i = nnodes-1; i >= 0; i--) {
      int node = nodemap[i];
      if (level[node] == curdepth) 
        nlastlev++;
      else
        break;
    }

    /* collect the last level nodes */

    lastlevnodes = (int *) malloc(nlastlev*sizeof(int));
    ntmp = 0;
    for (i = 0; i < nlastlev; i++) {
      int j = nnodes-1-i;
      int node = nodemap[j];
      lastlevnodes[ntmp++] = node;
    }

    /* sort them according to degree */

    mstk_quicksort(lastlevnodes, nadj, 0, nlastlev-1);

    /* Keep only the ones with the minimal degree */

    for (i = 1; i < nlastlev; i++) {
      if (nadj[lastlevnodes[i]] > nadj[lastlevnodes[0]]) {
        nlastlev = i+1; /* truncate */
        break;
      }
    }
    
    /* Build a level structure using each of these minimal degree last
       level nodes as a starting node */

    minmaxwidth = 1e+9;
    restart = 0;
    for (l = 0; l < nlastlev; l++) {

      Graph_Build_LevelStructure(nnodes, adj, nadj, offset, lastlevnodes[l],
                                 nodemap, level, &newdepth, &newmaxwidth);

      if (newmaxwidth < minmaxwidth) {
        minmaxwidth = newmaxwidth;
        *n_end = lastlevnodes[l];
      }

      if (newdepth > curdepth) {
        *n_beg = lastlevnodes[l];
        restart = 1;
        break;
      }      
    }

    free(lastlevnodes);
          
    if (restart) 
      continue;
    else 
      outerdone = 1;

  } /* while (!outerdone) */

  free(level);
  free(nodemap);
}


/* Given a graph and the two extremal nodes of a pseudo-diameter, build a level structure 
   of minimum width 

   nnodes   (in)    --- Number of nodes
   adj      (in)    --- Array listing adjacent nodes of nodes
   nadj     (in)    --- Number of adjacent nodes for each node
   offset   (in)    --- Starting location or offset in 'adj' array where the 
                        adjacency information for each node starts. 
   n_beg    (in)    --- First extremal node of pseudo-diameter for building the level structure
   n_end    (in)    --- Second extremal node of pseudo-diameter for building the level structure
   nodemap  (out)   --- Mapping of the nodes to new numbers (or in other words, 
                        the order they are added into the level structure)
   level    (out)   --- Assigned level of each node
   depth    (out)   --- Depth of the level structure (number of levels)
   maxwidth (out)   --- Maximum number of nodes in any one level
*/

void Graph_Assign_Levels_MinWidth(int nnodes, int adj[], int nadj[], int offset[], 
                                  int n_1, int n_k, int level[], int *depth, int *maxwidth) {

  int i, j, k, done;
  int *nodemap, *zeros, (*levelpairs)[2];
  int *baselevelwidth, *deltawidth1, *deltawidth2, *cset, *tag, *maxcset;
  int ncset, nmaxcset, maxwidth1, maxwidth2, curdepth, curmaxwidth;

  zeros = (int *) calloc(nnodes,sizeof(int));

  nodemap = (int *) malloc(nnodes*sizeof(int));

  levelpairs = (int (*)[2]) malloc(nnodes*sizeof(int [2]));

  /* Build levels starting from node n_1 - they form the first number in the level pairs */

  Graph_Build_LevelStructure(nnodes, adj, nadj, offset, n_1,
                             nodemap, level, &curdepth, &curmaxwidth);

  for (i = 0; i < nnodes; i++)
    levelpairs[i][0] = level[i];

  /* Build levels starting from node n_2 - they form the second number in the level pairs */

  Graph_Build_LevelStructure(nnodes, adj, nadj, offset, n_k,
                             nodemap, level, &curdepth, &curmaxwidth);

  for (i = 0; i < nnodes; i++)
    levelpairs[i][1] = curdepth+1-level[i];

  free(nodemap);


  /* reset the info in the array level1 */
  memcpy(level,zeros,nnodes*sizeof(int));
  
  /* Initiate a set of level widths and also width increments */
  baselevelwidth = (int *) malloc(curdepth*sizeof(int));
  deltawidth1 = (int *) malloc(curdepth*sizeof(int));
  deltawidth2 = (int *) malloc(curdepth*sizeof(int));


  /* For those nodes whose levels match from both directions, the
     level assignment is trivial */

  memcpy(baselevelwidth,zeros,curdepth*sizeof(int));
  for (i = 0; i < nnodes; i++)
    if (levelpairs[i][0] == levelpairs[i][1]) {
      int lev = levelpairs[i][0];
      level[i] = lev;
      baselevelwidth[lev-1]++;
    }


  /* Connected sets of nodes whose level pairs have a discrepancy are
     assigned the first or send level based on which will minimize the
     maximum width of the level structure (which in turn determines
     bandwidth) */

  /* Allocate some space for connected sets of nodes with unassigned levels */

  cset = (int *) malloc(nnodes*sizeof(int)); /* worst case allocation */
  maxcset = (int *) malloc(nnodes*sizeof(int));
  tag = (int *) malloc(nnodes*sizeof(int));


  /* Find sets of nodes that are connected and don't have a level
     assigned to them (level1[i] = 0) because the levels computed from
     them in two different directions are not the same. Process them
     in increasing order of size until none are left */

  done = 0;
  while (!done) { 
    nmaxcset = 0;
    for (i = 0; i < nnodes; i++) {

      if (level[i] == 0) { 
      
        /* level is not set for this node which means that there was a
           conflict between the two level predictions */
        /* Find the connected set of nodes which have not been assigned
           a level */

        memcpy(tag,zeros,nnodes*sizeof(int));
      
        cset[0] = i;
        ncset = 1;
        tag[i] = 1;
      
        for (j = 0; j < ncset; j++) {
          int node = cset[j];
          int off = offset[node];
        
          for (k = 0; k < nadj[node]; k++) {
            int adjnode = adj[off+k];
            if (level[adjnode] == 0 && tag[adjnode] == 0) {
              cset[ncset++] = adjnode;
              tag[adjnode] = 1;
            }
          }
        }

        if (ncset > nmaxcset) {
          nmaxcset = ncset;
          memcpy(maxcset,cset,ncset*sizeof(int));
        }

      }
    
    }

    if (nmaxcset > 0) {

      /* Find what the overall widths will be if the first or
         the second candidate for level numbers were chosen for
         these nodes */
    
      memcpy(deltawidth1,zeros,curdepth*sizeof(int));
      memcpy(deltawidth2,zeros,curdepth*sizeof(int));
    
      for (j = 0; j < nmaxcset; j++) {
        int node = maxcset[j];
        int lev;
        
        /* increment the width of the corresponding level if the first
           candidate level is chosen for this node */
      
        lev = levelpairs[node][0];
        deltawidth1[lev-1]++;
      
        /* increment the width of the corresponding level if the
           second candidate level is chosen for this node */
      
        lev = levelpairs[node][1];
        deltawidth2[lev-1]++;
      }
      
      /* Using the base levelwidths and width increments we can tell
         what the maximum width will be for each choice */
    
      maxwidth1 = 0;
      maxwidth2 = 0;
      for (i = 0; i < curdepth; i++) {
        int newwidth;
      
        newwidth = baselevelwidth[i]+deltawidth1[i];
        maxwidth1 = (newwidth > maxwidth1) ? newwidth : maxwidth1;
      
        newwidth = baselevelwidth[i]+deltawidth2[i];
        maxwidth2 = (newwidth > maxwidth2) ? newwidth : maxwidth2;
      }
    
      if (maxwidth1 < maxwidth2) { /* use the first candidate levels */
        for (j = 0; j < nmaxcset; j++) {
          int node = maxcset[j];
          int lev = levelpairs[node][0];
          level[node] = lev;
          baselevelwidth[lev-1]++;
        }
      }
      else {  /* use the second candidate levels */
        for (j = 0; j < nmaxcset; j++) {
          int node = maxcset[j];
          int lev = levelpairs[node][1];
          level[node] = lev;
          baselevelwidth[lev-1]++;
        }
      }
    }
    else 
      done = 1;
  } 

  *maxwidth = 0;
  for (i = 0; i < curdepth; i++)
    *maxwidth = baselevelwidth[i] > (*maxwidth) ? baselevelwidth[i] : (*maxwidth);


  free(baselevelwidth);
  free(deltawidth1);
  free(deltawidth2);
  free(cset);
  free(tag);
  free(maxcset);
  free(levelpairs);
  free(zeros);
}


/* Given a graph and a starting node, build a level structure and
   return some information about it.

   nnodes   (in)    --- Number of nodes
   adj      (in)    --- Array listing adjacent nodes of nodes
   nadj     (in)    --- Number of adjacent nodes for each node
   offset   (in)    --- Starting location or offset in 'adj' array where the 
                        adjacency information for each node starts. 
   n_beg    (in)    --- First extremal node of pseudo-diameter for building the level structure
   n_end    (in)    --- Second extremal node of pseudo-diameter for building the level structure   
   level    (out)   --- Level of each node
   nodemap  (out)   --- Mapping of the nodes to new numbers to minimize bandwidth
*/

void Graph_Assign_Numbers(int nnodes, int adj[], int nadj[], int offset[], 
                          int n_beg, int n_end, int level[], int nodemap[]) {

  int i, j, nnew, ncurlevnodes, lev, nnbrs_curlev, nnbrs_nxtlev;
  int *node0lev, *curlevnodes, *nbrs_curlev, *nbrs_nxtlev, *zeros, *tag;
  int depth, curlev, maxdeg;
  
  depth = 0;
  for (i = 0; i < nnodes; i++)
    if (level[i] > depth) depth = level[i];

  maxdeg = 0;
  for (i = 0; i < nnodes; i++)
    if (nadj[i] > maxdeg) maxdeg = nadj[i];

  /* If degree of n_end is greater than degree of n_beg, swap the two and
     reverse the levels */

  if (nadj[n_end] < nadj[n_beg]) {
    int tmp = n_beg;
    n_beg = n_end;
    n_end = tmp;

    for (i = 0; i < nnodes; i++)
      level[i] = depth-level[i]+1;
  }

  tag = (int *) calloc(nnodes,sizeof(int));

  for (i = 0; i < nnodes; i++) nodemap[i] = -1;
  nodemap[n_beg] = 0; 
  tag[n_beg] = 1;
  nnew = 1;


  /* March through level by level and number the nodes */

  curlevnodes = (int *) malloc(nnodes*sizeof(int)); /* assume the worst */
  nbrs_curlev = (int *) malloc(maxdeg*sizeof(int));
  nbrs_nxtlev = (int *) malloc(maxdeg*sizeof(int));

  for (lev = 1; lev <= depth; lev++) {
    
    /* collect all nodes of level 'lev' */

    ncurlevnodes = 0;
    for (i = 0; i < nnodes; i++) {
      if (level[i] == lev && nodemap[i] != -1) {
        curlevnodes[ncurlevnodes] = i;
        ncurlevnodes++;
      }
    }

    /* Sort the nodes in the level by their new node numbers */

    mstk_quicksort(curlevnodes, nodemap, 0, ncurlevnodes-1);


    /* For each of these nodes number the adjacent nodes that are in
       the SAME level by increasing degree (valence) */

    i = 0;
    while (i < ncurlevnodes) {
      int off, node;

      node = curlevnodes[i];
      if (nodemap[node] == -1) continue; /* unnumbered node */

      off = offset[node];
      
      /* collect adjacent nodes in the same level */

      nnbrs_curlev = 0;
      for (j = 0; j < nadj[node]; j++) {
        int adjnode = adj[off+j];
        if (level[adjnode] == lev && nodemap[adjnode] == -1) { 
          nbrs_curlev[nnbrs_curlev++] = adjnode;
          tag[adjnode] = 1;
        }
      }

      /* sort them according to degree */

      mstk_quicksort(nbrs_curlev, nadj, 0, nnbrs_curlev-1);
    
      /* Number any unnumbered neighbors in this level and add them
       to the nodes in this level */

      for (j = 0; j < nnbrs_curlev; j++) {
        int adjnode = nbrs_curlev[j];

        nodemap[adjnode] = nnew++;

        curlevnodes[ncurlevnodes] = adjnode;
        ncurlevnodes++;
      }

      i++;
    }

    /* Next, for each of these nodes number the adjacent nodes that are in
       the NEXT level by increasing degree (valence) */

    for (i = 0; i < ncurlevnodes; i++) {
      int node = curlevnodes[i];
      int off = offset[node];
      
      /* collect adjacent nodes in the next level */

      nnbrs_nxtlev = 0;
      for (j = 0; j < nadj[node]; j++) {
        int adjnode = adj[off+j];
        if (level[adjnode] == lev+1 && tag[adjnode] == 0) { 
          nbrs_nxtlev[nnbrs_nxtlev++] = adjnode;
          tag[adjnode] = 1;
        }
      }

      /* sort them according to degree */

      mstk_quicksort(nbrs_nxtlev, nadj, 0, nnbrs_nxtlev-1);
    
      /* Number any unnumbered ones in this level */

      for (j = 0; j < nnbrs_nxtlev; j++) {
        int adjnode = nbrs_nxtlev[j];
        nodemap[adjnode] = nnew++;
      }

    }
  } /* for each level */

  free(nbrs_curlev);
  free(nbrs_nxtlev);
  free(curlevnodes);
  free(tag);
}
