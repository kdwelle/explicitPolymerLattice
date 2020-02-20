#ifndef POLYMER_HPP
#define POLYMER_HPP


class Polymer{
/* 3-D lattice with seperate anion/cation heterogeneity, 
   absorbing and creating boundary conditions for cations, 
   reflecting boundary conditions for anions
*/
  public:
    Polymer(int,int,int,int); //sidelength, SLZ, numChains, sizePolymer, ionConc
    void  initialize();
    float get_energy();
	
	void  movePolymer(int); // ts
    void  moveIons(int);    //ts
	void  outputAnalysis(int,int);
    
 
  private:

    // todo: change variables to use camelCase??
   	int   SideLength;  // number of lattice sites on a side
	int   SLZ;         // number of lattice sites between electrodes
    int   numChains;   // number of polymer chains in the lattice
    int   sizePolymer; // number of connected lattice sites in a polymer chain
    float ionConc;
    float gamma;

    int   Ntot;        // = SideLength*Sidelength*SLZ (total number of lattice sites)
    float Eng;         // energy of the lattice system

    int   numBeads;   // numPolymer*sizePolymer
    int   numIons;     
    int   numPairs;   //(numBead+numIons)*(numBead+numIons-1)/2

    const static int numNeighbors = 6; //dimensionality of lattice
    long int idum;        //for random number generator: ran2(&idum)    

	std::vector< std::vector<int> >                sites;    // sites[siteIndex][x/y/z] = position in x/y/z
	std::vector< std::vector<int> >                neighb;   // neighb[siteIndex][0-5] = index of nearest neighbor 0-6
    std::vector< std::vector< std::vector<int> > > ssind;    // ssind[x][y][z] = siteIndex
    std::vector< std::vector<int> >                pos;      // pos[beadIndex][x/y/z] = position in x/y/z of bead
    std::vector< std::vector<int> >                ionPos;   // pos[ionIndex][x/y/z] = position in x/y/z of the ion
    std::vector< std::vector< std::vector<int> > > occ;      // occ[x][y][z] = 0/1 depending on if site is occupied by a polymer
    std::vector<float>                             charge;   // charge[ionIndex] = ion's charge
    std::vector< std::vector< std::vector<int> > > getPair;  // getpair[0/1/2][bead/ionIndex][bead/ionIndex] = pairIndex
    std::vector< std::vector <int> >               pairInd;  // pairInd[pairIndex][0/1] = ion/beadIndex
    std::vector<int>                               pairType; // pairType[pairIndex] = 0/1/2
    std::vector<float>                             pairDist; // distance between pair
    std::vector<float>                             pairEng;  // Energy contribution from pair


    
    bool              endMove(int);
    bool              kinkJump(int);
    bool              crankshaft(int);
    bool              slitherSnek(int);
    bool              isEnd(int);         //checks if bead at index is an end bead
    std::array<int,3> get2ndFromEnd(int); //gets the coordinates of a connected bead site of an end bead

    std::array<int,3> getVector(int, int); //get the vector of a bond between two monomers (by index)
    bool              areParallel(std::array<int,3>, std::array<int,3>);
    std::array<int,3> crossVector(std::array<int,3>, std::array<int,3>);
    std::array<int,3> addVector(int, std::array<int,3>); // add a vector to the position of the bead with index int
    
    void              getAllDist();            //refresh distance and energy vectors for all pairs
    void              getOneDist(int,int);     //refresh distance and energy vectors for pairs that incude the bead/ion with index
      
};

#endif
