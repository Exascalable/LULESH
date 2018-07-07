#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#if USE_MPI
#include <mpi.h>
#endif
#include "lulesh.h"

/* Helper function for converting strings to ints, with error checking */
int StrToInt(const char *token, int *retVal)
{
   const char *c ;
   char *endptr ;
   const int decimal_base = 10 ;

   if (token == NULL)
      return 0 ;
   
   c = token ;
   *retVal = (int)strtol(c, &endptr, decimal_base) ;
   if((endptr != c) && ((*endptr == ' ') || (*endptr == '\0')))
      return 1 ;
   else
      return 0 ;
}

static void PrintCommandLineOptions(char *execname, int myRank)
{
   if (myRank == 0) {

      printf("Usage: %s [opts]\n", execname);
      printf(" where [opts] is one or more of:\n");
      printf(" -q              : quiet mode - suppress all stdout\n");
      printf(" -i <iterations> : number of cycles to run\n");
      printf(" -s <size>       : length of cube mesh along side\n");
      printf(" -r <numregions> : Number of distinct regions (def: 11)\n");
      printf(" -b <balance>    : Load balance between regions of a domain (def: 1)\n");
      printf(" -c <cost>       : Extra cost of more expensive regions (def: 1)\n");
      printf(" -f <numfiles>   : Number of files to split viz dump into (def: (np+10)/9)\n");
      printf(" -p              : Print out progress\n");
      printf(" -v              : Output viz file (requires compiling with -DVIZ_MESH\n");
      printf(" -h              : This message\n");
      printf(" -cp             : Checkpoint frequency\n");
      printf("\n\n");
   }
}

static void ParseError(const char *message, int myRank)
{
   if (myRank == 0) {
      printf("%s\n", message);
#if USE_MPI      
      MPI_Abort(MPI_COMM_WORLD, -1);
#else
      exit(-1);
#endif
   }
}

void ParseCommandLineOptions(int argc, char *argv[],
                             int myRank, struct cmdLineOpts *opts)
{
   if(argc > 1) {
      int i = 1;

      while(i < argc) {
         int ok;
         /* -i <iterations> */
         if(strcmp(argv[i], "-i") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -i", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->its));
            if(!ok) {
               ParseError("Parse Error on option -i integer value required after argument\n", myRank);
            }
            i+=2;
         }
         /* -s <size, sidelength> */
         else if(strcmp(argv[i], "-s") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -s\n", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->nx));
            if(!ok) {
               ParseError("Parse Error on option -s integer value required after argument\n", myRank);
            }
            i+=2;
         }
	 /* -r <numregions> */
         else if (strcmp(argv[i], "-r") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -r\n", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->numReg));
            if (!ok) {
               ParseError("Parse Error on option -r integer value required after argument\n", myRank);
            }
            i+=2;
         }
	 /* -f <numfilepieces> */
         else if (strcmp(argv[i], "-f") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -f\n", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->numFiles));
            if (!ok) {
               ParseError("Parse Error on option -f integer value required after argument\n", myRank);
            }
            i+=2;
         }
         /* -p */
         else if (strcmp(argv[i], "-p") == 0) {
            opts->showProg = 1;
            i++;
         }
         /* -q */
         else if (strcmp(argv[i], "-q") == 0) {
            opts->quiet = 1;
            i++;
         }
         else if (strcmp(argv[i], "-b") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -b\n", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->balance));
            if (!ok) {
               ParseError("Parse Error on option -b integer value required after argument\n", myRank);
            }
            i+=2;
         }
         else if (strcmp(argv[i], "-c") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -c\n", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->cost));
            if (!ok) {
               ParseError("Parse Error on option -c integer value required after argument\n", myRank);
            }
            i+=2;
         }
         /* -v */
         else if (strcmp(argv[i], "-v") == 0) {
#if VIZ_MESH            
            opts->viz = 1;
#else
            ParseError("Use of -v requires compiling with -DVIZ_MESH\n", myRank);
#endif
            i++;
         }
         /* -h */
         else if (strcmp(argv[i], "-h") == 0) {
            PrintCommandLineOptions(argv[0], myRank);
#if USE_MPI            
            MPI_Abort(MPI_COMM_WORLD, 0);
#else
            exit(0);
#endif
         }
         else if(strcmp(argv[i], "-cp") == 0) {
        	 if (i+1 >= argc) {
        		 ParseError("Missing integer argument to -cp", myRank);
        	 }
        	 ok = StrToInt(argv[i+1], &(opts->cpInterval));
        	 if(!ok) {
        		 ParseError("Parse Error on option -cp integer value required after argument\n", myRank);
        	 }
        	 i+=2;
         }
         else {
            char msg[80];
            PrintCommandLineOptions(argv[0], myRank);
            sprintf(msg, "ERROR: Unknown command line argument: %s\n", argv[i]);
            ParseError(msg, myRank);
         }
      }
   }
}

/////////////////////////////////////////////////////////////////////

void VerifyAndWriteFinalOutput(Real_t elapsed_time,
                               Domain& locDom,
                               Int_t nx,
                               Int_t numRanks)
{
   // GrindTime1 only takes a single domain into account, and is thus a good way to measure
   // processor speed indepdendent of MPI parallelism.
   // GrindTime2 takes into account speedups from MPI parallelism 
   Real_t grindTime1 = ((elapsed_time*1e6)/locDom.cycle())/(nx*nx*nx);
   Real_t grindTime2 = ((elapsed_time*1e6)/locDom.cycle())/(nx*nx*nx*numRanks);

   Index_t ElemId = 0;
   printf("Run completed:  \n");
   printf("   Problem size        =  %i \n",    nx);
   printf("   MPI tasks           =  %i \n",    numRanks);
   printf("   Iteration count     =  %i \n",    locDom.cycle());
   printf("   Final Origin Energy = %12.6e \n", locDom.e(ElemId));

   Real_t   MaxAbsDiff = Real_t(0.0);
   Real_t TotalAbsDiff = Real_t(0.0);
   Real_t   MaxRelDiff = Real_t(0.0);

   for (Index_t j=0; j<nx; ++j) {
      for (Index_t k=j+1; k<nx; ++k) {
         Real_t AbsDiff = FABS(locDom.e(j*nx+k)-locDom.e(k*nx+j));
         TotalAbsDiff  += AbsDiff;

         if (MaxAbsDiff <AbsDiff) MaxAbsDiff = AbsDiff;

         Real_t RelDiff = AbsDiff / locDom.e(k*nx+j);

         if (MaxRelDiff <RelDiff)  MaxRelDiff = RelDiff;
      }
   }

   // Quick symmetry check
   printf("   Testing Plane 0 of Energy Array on rank 0:\n");
   printf("        MaxAbsDiff   = %12.6e\n",   MaxAbsDiff   );
   printf("        TotalAbsDiff = %12.6e\n",   TotalAbsDiff );
   printf("        MaxRelDiff   = %12.6e\n\n", MaxRelDiff   );

   // Timing information
   printf("\nElapsed time         = %10.2f (s)\n", elapsed_time);
   printf("Grind time (us/z/c)  = %10.8g (per dom)  (%10.8g overall)\n", grindTime1, grindTime2);
   printf("FOM                  = %10.8g (z/s)\n\n", 1000.0/grindTime2); // zones per second

   return ;
}

/////////////////////////////////////////////////////////
// Write application state in checkpoint file

void ApplicationCheckpointWrite(int rank, Domain& locDom, struct cmdLineOpts &opts, double start) {
	std::ofstream ofs;
	char filename[20];

	sprintf(filename, "tmp_%d", rank);
	ofs.open(filename, std::ios::out | std::ios::app | std::ios::binary);
	ofs.write(reinterpret_cast<char *> (&locDom.m_cycle), sizeof(Int_t));
	ofs.close();

	sprintf(filename, "check_%d_%d", rank, locDom.m_cycle);
	ofs.open(filename, std::ios::out | std::ios::binary);

	int size;
		size = locDom.m_x.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_x[i]), sizeof(Real_t));
		}

		size = locDom.m_y.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_y[i]), sizeof(Real_t));
		}

		size = locDom.m_z.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_z[i]), sizeof(Real_t));
		}

		size = locDom.m_xd.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_xd[i]), sizeof(Real_t));
		}

		size = locDom.m_yd.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_yd[i]), sizeof(Real_t));
		}

		size = locDom.m_zd.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_zd[i]), sizeof(Real_t));
		}

		size = locDom.m_xdd.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_xdd[i]), sizeof(Real_t));
		}

		size = locDom.m_ydd.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_ydd[i]), sizeof(Real_t));
		}

		size = locDom.m_zdd.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_zdd[i]), sizeof(Real_t));
		}

		size = locDom.m_fx.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_fx[i]), sizeof(Real_t));
		}

		size = locDom.m_fy.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_fy[i]), sizeof(Real_t));
		}

		size = locDom.m_fz.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_fz[i]), sizeof(Real_t));
		}

		size = locDom.m_nodalMass.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_nodalMass[i]), sizeof(Real_t));
		}

		size = locDom.m_symmX.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_symmX[i]), sizeof(Index_t));
		}

		size = locDom.m_symmY.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_symmY[i]), sizeof(Index_t));
		}

		size = locDom.m_symmZ.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_symmZ[i]), sizeof(Index_t));
		}

		ofs.write(reinterpret_cast<char *>(&locDom.m_numReg), sizeof(Int_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_cost), sizeof(Int_t));
		ofs.write(reinterpret_cast<char *>(locDom.m_regElemSize), sizeof(Index_t) * locDom.m_numReg);
		//int  nElem = locDom.m_numElem;
		ofs.write(reinterpret_cast<char *>(&locDom.m_numElem), sizeof(Index_t));
		ofs.write(reinterpret_cast<char *>(locDom.m_regNumList), sizeof(Index_t) * locDom.m_numElem);

		for(int i = 0; i < locDom.m_numReg; i++)
		{
			ofs.write(reinterpret_cast<char *>(locDom.m_regElemlist[i]), sizeof(Index_t) * locDom.regElemSize(i));
		}

		size = locDom.m_nodelist.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_nodelist[i]), sizeof(Index_t));
		}

		size = locDom.m_lxim.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_lxim[i]), sizeof(Index_t));
		}

		size = locDom.m_lxip.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_lxip[i]), sizeof(Index_t));
		}

		size = locDom.m_letam.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_letam[i]), sizeof(Index_t));
		}

		size = locDom.m_letap.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_letap[i]), sizeof(Index_t));
		}

		size = locDom.m_lzetam.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_lzetam[i]), sizeof(Index_t));
		}

		size = locDom.m_lzetap.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_lzetap[i]), sizeof(Index_t));
		}

		size = locDom.m_elemBC.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_elemBC[i]), sizeof(Int_t));
		}

		size = locDom.m_dxx.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_dxx[i]), sizeof(Real_t));
		}

		size = locDom.m_dyy.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_dyy[i]), sizeof(Real_t));
		}

		size = locDom.m_dzz.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_dzz[i]), sizeof(Real_t));
		}

		size = locDom.m_delv_xi.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_delv_xi[i]), sizeof(Real_t));
		}

		size = locDom.m_delv_eta.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_delv_eta[i]), sizeof(Real_t));
		}

		size = locDom.m_delv_zeta.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_delv_zeta[i]), sizeof(Real_t));
		}

		size = locDom.m_delx_xi.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_delx_xi[i]), sizeof(Real_t));
		}

		size = locDom.m_delx_eta.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_delx_eta[i]), sizeof(Real_t));
		}

		size = locDom.m_delx_zeta.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_delx_zeta[i]), sizeof(Real_t));
		}

		size = locDom.m_e.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_e[i]), sizeof(Real_t));
		}

		size = locDom.m_p.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_p[i]), sizeof(Real_t));
		}

		size = locDom.m_q.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_q[i]), sizeof(Real_t));
		}

		size = locDom.m_ql.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_q[i]), sizeof(Real_t));
		}

		size = locDom.m_qq.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_qq[i]), sizeof(Real_t));
		}

		size = locDom.m_v.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_v[i]), sizeof(Real_t));
		}

		size = locDom.m_volo.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_volo[i]), sizeof(Real_t));
		}

		size = locDom.m_vnew.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_vnew[i]), sizeof(Real_t));
		}

		size = locDom.m_delv.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_delv[i]), sizeof(Real_t));
		}

		size = locDom.m_vdov.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_vdov[i]), sizeof(Real_t));
		}

		size = locDom.m_arealg.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_arealg[i]), sizeof(Real_t));
		}

		size = locDom.m_ss.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_ss[i]), sizeof(Real_t));
		}

		size = locDom.m_elemMass.size();
		ofs.write(reinterpret_cast<char *>(&size), sizeof(int));
		for(int i = 0; i < size; i++) {
			ofs.write(reinterpret_cast<char *>(&locDom.m_elemMass[i]), sizeof(Real_t));
		}

		ofs.write(reinterpret_cast<char *>(&locDom.m_dtcourant), sizeof(Real_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_dthydro), sizeof(Real_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_cycle), sizeof(Int_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_dtfixed), sizeof(Real_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_time), sizeof(Real_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_deltatime), sizeof(Real_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_deltatimemultlb), sizeof(Real_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_deltatimemultub), sizeof(Real_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_dtmax), sizeof(Real_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_stoptime), sizeof(Real_t));

		ofs.write(reinterpret_cast<char *>(&locDom.m_numRanks), sizeof(Int_t));

		ofs.write(reinterpret_cast<char *>(&locDom.m_colLoc), sizeof(Index_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_rowLoc), sizeof(Index_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_planeLoc), sizeof(Index_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_tp), sizeof(Index_t));

		ofs.write(reinterpret_cast<char *>(&locDom.m_sizeX), sizeof(Index_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_sizeY), sizeof(Index_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_sizeZ), sizeof(Index_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_numElem), sizeof(Index_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_numNode), sizeof(Index_t));

		ofs.write(reinterpret_cast<char *>(&locDom.m_maxPlaneSize), sizeof(Index_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_maxEdgeSize), sizeof(Index_t));

#if _OPENMP
		ofs.write(reinterpret_cast<char *>(locDom.m_nodeElemStart), sizeof(Index_t) * (locDom.m_numNode + 1));
		int elem = locDom.m_nodeElemStart[locDom.numNode()];
		ofs.write(reinterpret_cast<char *>(locDom.m_nodeElemCornerList), sizeof(Index_t) * elem);
#endif

		ofs.write(reinterpret_cast<char *>(&locDom.m_rowMin), sizeof(Index_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_rowMax), sizeof(Index_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_colMin), sizeof(Index_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_colMax), sizeof(Index_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_planeMin), sizeof(Index_t));
		ofs.write(reinterpret_cast<char *>(&locDom.m_planeMax), sizeof(Index_t));

		//Size of comm send/recv buffer
		ofs.write(reinterpret_cast<char *>(&locDom.commBufSize), sizeof(int));
		ofs.write(reinterpret_cast<char *>(locDom.commDataSend), sizeof(Real_t) * locDom.commBufSize);
		ofs.write(reinterpret_cast<char *>(locDom.commDataRecv), sizeof(Real_t) * locDom.commBufSize);

		//struct
		ofs.write(reinterpret_cast<char *>(&opts), sizeof(opts));

	//Time
#if USE_MPI
		ofs.write(reinterpret_cast<char *>(&start), sizeof(double));
#endif
		ofs.close();
}

/////////////////////////////////////////////////////////
// Read application state from checkpoint file

void ApplicationCheckpointRead(int rank, Domain& locDom, struct cmdLineOpts &opts, double &start) {
	std::ifstream ifs;
	char filename[20];

	sprintf(filename, "tmp_%d", rank);
	ifs.open(filename, std::ios::in | std::ios::binary);
	int second_to_last, last;
	std::vector<int> cp_count;
	while (true) {
		ifs.read(reinterpret_cast<char *> (&last), sizeof(Int_t));
		if (ifs.eof()) break;
		cp_count.push_back(last);
	}
	ifs.close();
	second_to_last = cp_count[cp_count.size() - 2];
	sprintf(filename, "check_%d_%d", rank, second_to_last);
	ifs.open(filename, std::ios::in | std::ios::binary);

	int sz;
	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_x.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_x[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_y.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_y[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_z.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_z[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_xd.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_xd[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_yd.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_yd[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_zd.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_zd[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_xdd.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_xdd[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_ydd.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_ydd[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_zdd.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_zdd[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_fx.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_fx[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_fy.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_fy[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_fz.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_fz[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_nodalMass.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_nodalMass[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_symmX.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_symmX[i]), sizeof(Index_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_symmY.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_symmY[i]), sizeof(Index_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_symmZ.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_symmZ[i]), sizeof(Index_t));
	}

	ifs.read(reinterpret_cast<char *>(&locDom.m_numReg), sizeof(Int_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_cost), sizeof(Int_t));

	locDom.m_regElemSize = new Index_t[locDom.m_numReg];
	ifs.read(reinterpret_cast<char *>(locDom.m_regElemSize), sizeof(Index_t) * locDom.m_numReg);

	Index_t nElem;
	ifs.read(reinterpret_cast<char *>(&nElem), sizeof(Index_t));
	locDom.m_regNumList = new Index_t[nElem];
	ifs.read(reinterpret_cast<char *>(locDom.m_regNumList), sizeof(Index_t) * nElem);

	locDom.m_regElemlist = new Index_t*[locDom.m_numReg];
	for (int i = 0; i < locDom.m_numReg; i++) {
		locDom.m_regElemlist[i] = new Index_t[locDom.regElemSize(i)];
	}
	for(int i = 0; i < locDom.m_numReg; i++)
	{
		ifs.read(reinterpret_cast<char *>(locDom.m_regElemlist[i]), sizeof(Index_t) * locDom.regElemSize(i));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_nodelist.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_nodelist[i]), sizeof(Index_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_lxim.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_lxim[i]), sizeof(Index_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_lxip.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_lxip[i]), sizeof(Index_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_letam.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_letam[i]), sizeof(Index_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_letap.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_letap[i]), sizeof(Index_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_lzetam.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_lzetam[i]), sizeof(Index_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_lzetap.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_lzetap[i]), sizeof(Index_t));
	}


	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_elemBC.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_elemBC[i]), sizeof(Int_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_dxx.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_dxx[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_dyy.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_dyy[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_dzz.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_dzz[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_delv_xi.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_delv_xi[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_delv_eta.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_delv_eta[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_delv_zeta.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_delv_zeta[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_delx_xi.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_delx_xi[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_delx_eta.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_delx_eta[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_delx_zeta.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_delx_zeta[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_e.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_e[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_p.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_p[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_q.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_q[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_ql.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_ql[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_qq.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_qq[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_v.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_v[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_volo.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_volo[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_vnew.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_vnew[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_delv.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_delv[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_vdov.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_vdov[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_arealg.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_arealg[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_ss.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_ss[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&sz), sizeof(int));
	locDom.m_elemMass.resize(sz);
	for (int i = 0; i < sz; i++) {
		ifs.read(reinterpret_cast<char *>(&locDom.m_elemMass[i]), sizeof(Real_t));
	}

	ifs.read(reinterpret_cast<char *>(&locDom.m_dtcourant), sizeof(Real_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_dthydro), sizeof(Real_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_cycle), sizeof(Int_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_dtfixed), sizeof(Real_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_time), sizeof(Real_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_deltatime), sizeof(Real_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_deltatimemultlb), sizeof(Real_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_deltatimemultub), sizeof(Real_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_dtmax), sizeof(Real_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_stoptime), sizeof(Real_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_numRanks), sizeof(Int_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_colLoc), sizeof(Index_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_rowLoc), sizeof(Index_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_planeLoc), sizeof(Index_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_tp), sizeof(Index_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_sizeX), sizeof(Index_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_sizeY), sizeof(Index_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_sizeZ), sizeof(Index_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_numElem), sizeof(Index_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_numNode), sizeof(Index_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_maxPlaneSize), sizeof(Index_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_maxEdgeSize), sizeof(Index_t));

#if _OPENMP
	locDom.m_nodeElemStart = new Index_t[locDom.m_numNode + 1];
	ifs.read(reinterpret_cast<char *>(locDom.m_nodeElemStart), sizeof(Index_t) * (locDom.m_numNode + 1));
	int elem = locDom.m_nodeElemStart[locDom.numNode()];
	locDom.m_nodeElemCornerList = new Index_t[elem];
	ifs.read(reinterpret_cast<char *>(locDom.m_nodeElemCornerList), sizeof(Index_t) * elem);
#endif

	ifs.read(reinterpret_cast<char *>(&locDom.m_rowMin), sizeof(Index_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_rowMax), sizeof(Index_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_colMin), sizeof(Index_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_colMax), sizeof(Index_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_planeMin), sizeof(Index_t));
	ifs.read(reinterpret_cast<char *>(&locDom.m_planeMax), sizeof(Index_t));
	ifs.read(reinterpret_cast<char *>(&locDom.commBufSize), sizeof(int));
	locDom.commDataSend = new Real_t[locDom.commBufSize];
	ifs.read(reinterpret_cast<char *>(locDom.commDataSend), sizeof(Real_t) * locDom.commBufSize);
	locDom.commDataRecv = new Real_t[locDom.commBufSize];
	ifs.read(reinterpret_cast<char *>(locDom.commDataRecv), sizeof(Real_t) * locDom.commBufSize);

	//struct
	ifs.read(reinterpret_cast<char *>(&opts), sizeof(opts));

	// time
#if USE_MPI
	ifs.read(reinterpret_cast<char *>(&start), sizeof(double));
#endif

	ifs.close();
}
