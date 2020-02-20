#include "polymer.cpp"
int main(int argc, char *argv[]){
    
  if (argc != 5){
    // SideLength: number of lattice sites on an edge (must be less than 100)
    // SLZ: number of lattice sites between electrodes
    // numPolymer: 
    // sizePolymer: 
    cerr << "Usage Error: SideLength, SLZ, numPolymer, sizePolymer" << endl;
    cerr << "you had " << argc << "arguments" << endl;
    exit(1); 
  }
  
  // SideLength = atoi( argv[1] );
  // SLZ = atoi( argv[2] );
  // numChains = atof( argv[3] );
  // sizePolymer = atof( argv[4] );

  

  // Make a lattice object:
  Polymer polymer(atoi( argv[1] ),atoi( argv[2] ),atof( argv[3] ),atof( argv[4] ));
  polymer.outputAnalysis(0,100);
  
}