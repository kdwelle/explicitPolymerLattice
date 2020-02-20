#include "lattice.cpp"
int main(int argc, char *argv[]){
    
  if (argc != 8){
    // SideLength: number of lattice sites on an edge (must be less than 100)
    // SLZ: number of lattice sites between electrodes
    // rho(pos): fraction of sites which are inaccessible to positive ions (between 0 and 1)
    // rho(neg): fraction of sites which are inaccessible to negative ions (between 0 and 1)
    // concentration: fraction of sites occupied by (each positive and negative) ions
    // gamma: parameter which is inversely related to strength of coulombic forces
    // outind: determines type of output
    //   0: "timestep energy" to stderr
    //      .xyz file to stdout
    //   1: timestep | CP | MP | MN | Energy | "Con"
    //      timestep | i | density[i][0] | density[i][1] | "Den"
    //   2: timestep | ind1 | ind2
    //   3: timestep | MP | "MP"
    //      timestep | i | pos[i][0] | pos[i][1] | pos[i][2] | charge[i] | "Pos"
    cerr << "Usage Error: SideLength, SLZ, rho(pos), rho(neg), concentration, gamma, outind" << endl;
    exit(1); 
  }
  
  // SideLength = atoi( argv[1] );
  // SLZ = atoi( argv[2] );
  // rhop = atof( argv[3] );
  // rhon = atof( argv[4] );
  // icon = atof( argv[5] );
  // gam = atof( argv[6] );
  
  int outind = atoi( argv[7] );
  
  if(atof( argv[5] )>.1){
    cerr << "concentration is too high, icon<.1" << endl;
    exit(1);
  }

  // Make a lattice object:
  Lattice lattice(atoi( argv[1] ),atoi( argv[2] ),atof( argv[3] ),atof( argv[4] ),atof( argv[5] ),atof( argv[6] ));

  // equilibrate for 500 Monte Carlo sweeps and then generate a 20000 MC sweep trajectory
  // todo: make number of timesteps not a magic number?
  for(int ts = -500 ; ts <= 20000 ; ts++){
    lattice.moveIons(ts); // Monte Carlo sweep
    
    if (outind ==4){
      if (ts == 20000){
        lattice.outputAnalysis(ts,outind);
      }
    }
    else if (outind==1){
      if (ts%1000==0){
        lattice.outputAnalysis(ts,outind);
      }
    }
    else if (ts >= 0){
      lattice.outputAnalysis(ts,outind);
    }
    
  }
}