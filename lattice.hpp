#ifndef LATTICE_HPP
#define LATTICE_HPP


class Lattice{
/* 3-D lattice with seperate anion/cation heterogeneity, 
   absorbing and creating boundary conditions for cations, 
   reflecting boundary conditions for anions
*/
  public:
    Lattice(int,int,float,float,float,float);
    void  initialize();
		void  get_drpair0();
		void  get_drpair1(int);
    float get_energy();
	
		void  moveIons(int); // perform a monte carlo iteration
		void  outputAnalysis(int,int);
 
  private:

    // todo: change variables to use camelCase??

   	int   SideLength; // number of lattice sites on a side
		int   SLZ;        // number of lattice sites between electrodes
		float rhop;       // fraction of lattice sites inaccessible to positive ions
		float rhon;       // fraction of lattice sites inaccessible to negative ions
		float icon;       // fraction of lattice sites populated with pos/neg ions
    float gam;        // parameter that is inversely related to strength of coulomb force
    int   Ntot;       // = SideLength*Sidelength*SLZ (total number of lattice sites)
    int   Nion;       // = Ntot*icon*2 (total number of ions)
    int   Npos;       // = Ntot*rhop (Number of positive inacessible sites)
    int   Nneg;       // = Ntot*rhon (Number of negative inacessible sites)
    int   NIpairs;    // = Nion! (Total number of ion pairs)
    float Eng;        // energy of the lattice syste

    int MP; // number of positive ion moves accepted since last output
    int MN; // number of negative ion moves accepted since last output
    int CP; // number of ions absorbed since last output

    const static int numNeighbors = 6; //dimensionality of lattice
    long int idum;        //for random number generator: ran2(&idum)

		std::vector< std::vector<int> >                sites;    // sites[siteIndex][x/y/z] = position in x/y/z
		std::vector< std::vector<int> >                neighb;   // neighb[siteIndex][0-5] = index of nearest neighbor 0-6
    std::vector< std::vector< std::vector<int> > > ssind;    // ssind[x][y][z] = siteIndex
    std::vector< std::vector<int> >                pos;      // pos[ionIndex][x/y/z] = position in x/y/z of ion
	  std::vector<float>                             charge;   // charge[ionIndex] = ion's charge
    std::vector< std::vector< std::vector<int> > > occ;      // occ[x][y][z] = 0/1 depending on if site is occupied by an ion
		std::vector<float>                             drpair;   // drpair[ionPairIndex] = 3-D distance between ions in ion pair
		std::vector< std::vector<int> >                pairind;  // pairind[ionPairIndex][0/1] = index of 1st/2nd ion in pair
		std::vector< std::vector<int> >                pairind2; // pairind2[ionIndex 1][ionIndex 2] = ionPairIndex
		std::vector< std::vector< std::vector<int> > > netp;     // netp[x][y][z] = 0/1 depending on if site is blocked to positive ions
		std::vector< std::vector< std::vector<int> > > netn;     // netn[x][y][z] = 0/1 depending on if site is blocked to negative ions
    std::vector< std::vector<int> >                density;  // density[i][0]: number of times a timestep ended with a positive ion in that slice
                                                             // denisty[i][1]: number of times a timestep ended with a negative ion in that slice
    
      
};

#endif
