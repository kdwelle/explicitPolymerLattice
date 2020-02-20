#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "lattice.hpp"
#include "randomGen.cpp"
using namespace std;

Lattice::Lattice(int SideLength, int SLZ, float rhop, float rhon, float icon, float gam): SideLength(SideLength), SLZ(SLZ), rhop(rhop), rhon(rhon), icon(icon), gam(gam)
{
  idum=-time(NULL); // todo: random number seed shouldn't be regenenerated for every new object
  Ntot = SideLength*SideLength*SLZ;
  this->initialize();
  this->get_drpair0();
}

void Lattice::initialize(){
	/*initialization for the lattice model
    set sizes of arrays, zero them out, add ions etc.
  */

  // Build sites, ssind, netp, netn, occ 
  sites.reserve(Ntot); // allocate memory for the arrays now
  ssind.reserve(SideLength);
  netp.reserve(SideLength);
  netn.reserve(SideLength);
  occ.reserve(SideLength);
  neighb.reserve(Ntot);
  // Some temperary arrays for building the matricies
  vector<int> sitesRow; //can reuse on each loop
  vector< vector <int> > ssindY;
  vector<int> ssindZ (SLZ); //initialize with SLZ elements set to default value (0)
  vector<int> neighbRow (numNeighbors);
  
  int count = 0;
  for(int i = 0 ; i < SideLength ; i++){
    ssindY.clear(); //use same ssindY vector on each run so memory only has to be allocated once
    for(int j = 0 ; j < SideLength ; j++){
      for(int k = 0 ; k < SLZ ; k++){
        sitesRow.clear(); 
        sitesRow.push_back(i);
        sitesRow.push_back(j);
        sitesRow.push_back(k);
        sites.push_back(sitesRow);

        ssindZ[k] = count;
        count++;

        neighb.push_back(neighbRow);
      }
      ssindY.push_back(ssindZ); // add row to the slice (ssindY)
    }
    ssind.push_back(ssindY); // add slice to the ssind matrix
  }

  // Build arrays that start as all zeroes
  vector< vector <int> > tempY;
  tempY.reserve(SideLength);
  vector <int> tempZ (SLZ,0); //initialize with zeroes
  for(int j = 0 ; j < SideLength ; j++){
    tempY.push_back(tempZ); // build a slice of all zeroes
  }
  for(int i = 0 ; i < SideLength ; i++){
    netp.push_back(tempY); //push_back slice into the matricies
    netn.push_back(tempY);
    occ.push_back(tempY);
  }

  // Build pos
 
  int ind1,ind2;
  //generate the neighbor lists (periodic boundary conditions)
  for(int i = 0 ; i < SideLength ; i++){
    for(int j = 0 ; j < SideLength ; j++){
      for(int k = 0 ; k < SLZ ; k++){
      	count=0;
      	ind2=ssind[i][j][k];
      	ind1=i-1;
      	if(ind1<0) {ind1=SideLength-1;}
      	neighb[ind2][count]=ssind[ind1][j][k];
      	count++;
      	ind1=i+1;
      	if(ind1==SideLength) {ind1=0;}
      	neighb[ind2][count]=ssind[ind1][j][k];
      	count++;
      	ind1=j-1;
      	if(ind1<0) {ind1=SideLength-1;}
      	neighb[ind2][count]=ssind[i][ind1][k];
      	count++;
      	ind1=j+1;
      	if(ind1==SideLength) {ind1=0;}
      	neighb[ind2][count]=ssind[i][ind1][k];
      	count++;
      	ind1=k-1;
      	if(ind1<0) {ind1=SLZ-1;}
      	neighb[ind2][count]=ssind[i][j][ind1];
      	count++;
      	ind1=k+1;
      	if(ind1==SLZ) {ind1=0;}
      	neighb[ind2][count]=ssind[i][j][ind1];
      }
    }
  }  
  
  int itemp[4];
  bool found;
  
  Npos=Ntot*rhop;  
  for(int i = 0 ; i < Npos ; i++){ // add blocked sites for positive ions randomly
    found=false;
    while(!found){
      itemp[0]=int(ran2(&idum)*SideLength);
      itemp[1]=int(ran2(&idum)*SideLength);
      itemp[2]=int(ran2(&idum)*SLZ);
      if(netp[itemp[0]][itemp[1]][itemp[2]]==0){
      	netp[itemp[0]][itemp[1]][itemp[2]]=1;
      	found=true;
      }
    }
  }
  
  Nneg=Ntot*rhon;
  for(int i = 0 ; i < Nneg ; i++){ // add blocked sites for negative ions randomly
    found=false;
    while(!found){
      itemp[0]=int(ran2(&idum)*SideLength);
      itemp[1]=int(ran2(&idum)*SideLength);
      itemp[2]=int(ran2(&idum)*SLZ);
      if(netn[itemp[0]][itemp[1]][itemp[2]]==0){
      	netn[itemp[0]][itemp[1]][itemp[2]]=1;
      	found=true;
      }
    }
  }
  
  itemp[3]=Ntot*icon;
  Nion=itemp[3]*2; 
  if(Nion/2 >= Ntot-(Npos+Nneg)/2){ // todo: fix this
    cerr << "ERROR: too many ions" << endl;
    exit(1);
  }

  charge.resize(Nion); // size charge matrix appropiately
  drpair.resize(Nion*(Nion-1)); // size drpair matrix

  pos.reserve(Nion);   // build pos matrix
  vector<int> posTemp(3);
  for (int i = 0; i < Nion; ++i){
    pos.push_back(posTemp);
  }

  count=0;
  for(int i = 0 ; i < itemp[3] ; i++){ // add positive ions randomly
    found=false;
    while(!found){
      itemp[0]=int(ran2(&idum)*SideLength);
      itemp[1]=int(ran2(&idum)*SideLength);
      itemp[2]=int(ran2(&idum)*SLZ);
      if(netp[itemp[0]][itemp[1]][itemp[2]]==0 && occ[itemp[0]][itemp[1]][itemp[2]]==0){
      	occ[itemp[0]][itemp[1]][itemp[2]]=1;
      	pos[count][0]=itemp[0];
      	pos[count][1]=itemp[1];
      	pos[count][2]=itemp[2];
      	charge[count]=1;
      	found=true;
      	count++;
      }
    }
  }
  
  for(int i = 0 ; i < itemp[3] ; i++){ // add negative ions randomly
    found=false;
    while(!found){
      itemp[0]=int(ran2(&idum)*SideLength);
      itemp[1]=int(ran2(&idum)*SideLength);
      itemp[2]=int(ran2(&idum)*SLZ);
      if(netn[itemp[0]][itemp[1]][itemp[2]]==0 && occ[itemp[0]][itemp[1]][itemp[2]]==0){
      	occ[itemp[0]][itemp[1]][itemp[2]]=1;
      	pos[count][0]=itemp[0];
      	pos[count][1]=itemp[1];
      	pos[count][2]=itemp[2];
      	charge[count]=-1;
      	found=true;
      	count++;
      }
    }
  }  

  NIpairs=0;
  vector<int> pairindTemp;
  vector<int> pairind2Temp(Nion);
  pairind2.reserve(Nion);
  for(int i = 0 ; i < Nion-1 ; i++){ // Build list of ion pairs
    pairind2.push_back(pairind2Temp);
    for(int j = i+1 ; j < Nion ; j++){
      pairindTemp.clear();
      pairindTemp.push_back(i);
      pairindTemp.push_back(j);
      pairind.push_back(pairindTemp);
      
      pairind2[i][j] = NIpairs;
      NIpairs++;
    }
  }
  pairind2.push_back(pairind2Temp);

  for(int i = 0 ; i < Nion-1 ; i++){ // mirror pairind2
    for(int j = i+1 ; j < Nion ; j++){
      pairind2[j][i] = pairind2[i][j];
    }
  }

  drpair.resize(NIpairs); //initialize drpairs
  this->get_drpair0(); //needs to be called before get_energy
  Eng = this->get_energy();
 
  density.reserve(SLZ);
  vector<int> temp(2); //vector [0,0]
  for(int i = 0 ; i < SLZ ; i++){ //initialize density
    density.push_back(temp);
  }

}

void Lattice::get_drpair0(){

  int dr[3];
  int ind1,ind2;
  //#pragma omp parallel for private(ind1,ind2,dr)
  for(int i = 0 ; i < NIpairs ; i++){
    ind1=pairind[i][0];
    ind2=pairind[i][1];
    for(int k = 0 ; k < 3 ; k++){
      dr[k]=pos[ind1][k]-pos[ind2][k];
      if(k<2){
	if(dr[k]>SideLength/2.) dr[k]-=SideLength;
	if(dr[k]<-SideLength/2.) dr[k]+=SideLength;
      }
    }
    drpair[i]=sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
  } 
}

void Lattice::get_drpair1(int ind1){

  int dr[3];
  int ind0;
  //#pragma omp parallel for private(ind0,dr)
  for(int ind2 = 0 ; ind2 < Nion ; ind2++){
    if(ind2!=ind1){
      ind0=pairind2[ind1][ind2];
      for(int k = 0 ; k < 3 ; k++){
      	dr[k]=pos[ind1][k]-pos[ind2][k];
      	if(k<2){
      	  if(dr[k]>SideLength/2.) dr[k]-=SideLength;
      	  if(dr[k]<-SideLength/2.) dr[k]+=SideLength;
      	}
      }
      drpair[ind0]=sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
    }
  }
}

float Lattice::get_energy(){

  float eng=0.0;
  float ftemp;
  int ind1,ind2;
  
  for(int i = 0 ; i < NIpairs ; i++){
    ind1=pairind[i][0];
    ind2=pairind[i][1];
    ftemp=charge[ind1]*charge[ind2]/(drpair[i]*gam);
    eng+=ftemp;
  }
  this->Eng = eng;
  return eng;
}

void Lattice::moveIons(int ts){ //todo: should this really be passed as a parameter?
  bool acc; // declare some variables for use inside the loop
  bool found;
  int iter;
  int ind1,ind2,ind3,ind4;
  int x0,y0,z0;
  int x1,y1,z1;
  float Eng2,deltaE;


  if(ts==0){ // todo: this seems inefficient
      CP=0;
      MN=0;
      MP=0;
    }

  //perform a MC sweep (make N=Nion trial moves)
  for(int n = 0 ; n < Nion ; n++){
    //start off assuming we will accept the trial move
    acc=true;

    //choose random ion to move
    ind1=int(Nion*ran2(&idum));

    //x0,y0,z0 are the ion's current position
    x0=pos[ind1][0];
    y0=pos[ind1][1];
    z0=pos[ind1][2];

    //choose a random nearest neighbor site to hop to
    ind2=int(6*ran2(&idum));
    ind4=ssind[x0][y0][z0];
    ind3=neighb[ind4][ind2];

    //x1,y1,z1 are the positions of the site the ion is trying to hop to
    x1=sites[ind3][0];
    y1=sites[ind3][1];
    z1=sites[ind3][2];
    
    //dynamic constraints
    if(ts>=0){
      
      //for anions
      if(charge[ind1]<0){
      
        //(z0-z1)**2 is only greater than 1 if the ion is hopping through the periodic boundaries in the z-direction. The electrodes prevent this for anions.
        if((z0-z1)*(z0-z1)>1){
          //reject move
          acc=false;
        }
      }
      
      //for cations
      if(charge[ind1]>0){
        //cations can't move back into left-hand electrode
        if(z0==0 && z1==(SLZ-1)){
          acc=false;
        }
        if(z0==(SLZ-1) && z1==0){
          //if cation is trying to hop into the right-hand electrode we need to find an empty site to create a new cation at the left-hand electrode. This type of move is automatically accepted witout evaluating agains metropolis condition
          acc=false;
          found=false;
          iter=0;
          while(!found){
            //choose random site in the x-y plane
            x1=int(SideLength*ran2(&idum));
            y1=int(SideLength*ran2(&idum));
            
            //check to see if site is empty
            if(occ[x1][y1][z1]==0){
              //put ion at new position in left-hand electrode
              pos[ind1][0]=x1;
              pos[ind1][1]=y1;
              pos[ind1][2]=z1;
              
              //find distance between ion and all other ions in the system
              get_drpair1(ind1);
              
              //compute energy of system
              Eng=get_energy();
              
              occ[x0][y0][z0]=0;
              occ[x1][y1][z1]=1;
              CP++;
              found=true;
            }
            iter++;
            //fail if there are no empty places to put new ion
            if(iter>1000){
              cerr << "ERROR: Could not find landing spot" << endl;
              exit(1);
            }
          }
        }
      }
    }
    
    if(acc){
      if(charge[ind1]>0){
        //make sure the trial move site is not blocked by cation spatial contraint or by another ion
        if(netp[x1][y1][z1]==0 && occ[x1][y1][z1]==0){
          //move particle
          pos[ind1][0]=x1;
          pos[ind1][1]=y1;
          pos[ind1][2]=z1;
          get_drpair1(ind1);
          
          //compute new energy
          Eng2=get_energy();

          //change in energy associate with move
          deltaE=Eng2-Eng;

          //metropolis conditions is below
          if(deltaE<0){
            occ[x0][y0][z0]=0;
            occ[x1][y1][z1]=1;
            Eng=Eng2;
            MP++;
          }
          else{
            if(ran2(&idum)<=exp(-deltaE)){
              occ[x0][y0][z0]=0;
              occ[x1][y1][z1]=1;
              Eng=Eng2;
              MP++;
            }else{
              pos[ind1][0]=x0;
              pos[ind1][1]=y0;
              pos[ind1][2]=z0;
              get_drpair1(ind1);
            }
          }
        }
      }
      
      //same thing for anion
      if(charge[ind1]<0){
        if(netn[x1][y1][z1]==0 && occ[x1][y1][z1]==0){
          pos[ind1][0]=x1;
          pos[ind1][1]=y1;
          pos[ind1][2]=z1;
          get_drpair1(ind1);
          Eng2=get_energy();
          deltaE=Eng2-Eng;
          if(deltaE<0){
            MN++;
            occ[x0][y0][z0]=0;
            occ[x1][y1][z1]=1;
            Eng=Eng2;
          }
          else{
            if(ran2(&idum)<=exp(-deltaE)){
              MN++;
              occ[x0][y0][z0]=0;
              occ[x1][y1][z1]=1;
              Eng=Eng2;
            }else{
              pos[ind1][0]=x0;
              pos[ind1][1]=y0;
              pos[ind1][2]=z0;
              get_drpair1(ind1);
            }
          }
        }
      }
    }
    //gather some statistics
    for(int i = 0 ; i < Nion ; i++){
      if(charge[i]>0){
        density[pos[i][2]][0]++;
      }
      if(charge[i]<0){
        density[pos[i][2]][1]++;
      }
    }
  }

}

void Lattice::outputAnalysis(int ts, int outind){
  // outind: determines type of output
  //   0: "timestep energy" to stderr
  //      .xyz file to stdout
  //   1: timestep | CP | MP | MN | Energy | "Con"
  //      timestep | i | density[i][0] | density[i][1] | "Den"
  //   2: timestep | ind1 | ind2
  //   3: timestep | MP | "MP"
  //      timestep | i | pos[i][0] | pos[i][1] | pos[i][2] | charge[i] | "Pos"
  //   4: CP


  //create xyz file of trajectory
  // nN: site unavailible to negative ion
  // nP: site unavailible to positive ion
  // iN: negative ion site
  // iP: positive ion site
  if(outind==0){   
    cerr << ts << "  " << Eng << endl;
    cout << Nion+Npos+Nneg << endl << endl;
    for(int i = 0 ; i < Nion ; i++){
      if(charge[i]>0){
        cout << "iP  " << pos[i][0] << "  " << pos[i][1] << "  " << pos[i][2] << endl;
      }
      if(charge[i]<0){
        cout << "iN  " << pos[i][0] << "  " << pos[i][1] << "  " << pos[i][2] << endl;
      }
    }
    for(int i = 0 ; i < SideLength ; i++){
      for(int j = 0 ; j < SideLength ; j++){
        for(int k = 0 ; k < SLZ ; k++){
          if(netp[i][j][k]==1){
            cout << "nP  " << i << "  " << j << "  " << k << endl;
          }
          if(netn[i][j][k]==1){
            cout << "nN  " << i << "  " << j << "  " << k << endl;
          }
        }
      }
    } 
  }

  // timestep | CP | MP | MN | Energy | "Con"
  // timestep | i | density[i][0] | density[i][1] | "Den"
  //
  // CP: # of ions absorbed since last output
  // MP: number of positive ion moves accepted since last output
  // MN: number of negative ion moves accepted since last output
  // i: slice position. Goes from 0 to SLZ (distance between electrodes)
  // density[i][0]: number of times a timestep ended with a positive ion in that slice
  // denistry[i][1]: number of times a timestep ended with a negative ion in that slice
  if(outind==1){
    cout << ts << "  " << CP << "  " << MP << "  " << MN << "  " << Eng << "  Con" << endl;
    CP=0;
    MP=0;
    MN=0;
    for(int i = 0 ; i < SLZ ; i++){
      cout << ts << "  " << i << "  " << density[i][0] << "  " << density[i][1] << "  Den" << endl;
      density[i][0]=0;
      density[i][1]=0;
    }
  }
  
  
  // timestep | ind1 | ind2
  //
  // ind1: number of positive ions
  // ind2: number of negative ions
  if(outind==2){
    for(int i = 0 ; i < SLZ ; i++){
      density[i][0]=0;
      density[i][1]=0;
    }
    int ind1=0;
    int ind2=0;
    for(int i = 0 ; i < Nion ; i++){
      if(charge[i]>0){
        density[pos[i][2]][0]++;
        ind1++;
      }
      if(charge[i]<0){
        density[pos[i][2]][1]++;
        ind2++;
      }
    }
    cout << ts << "  " << ind1 << "  " << ind2 << endl;
    /*
    for(int i = 0 ; i < SideLength ; i++){
    cout << ts << "  " << i << "  " << density[i][0]/float(ind1) << "  " << density[i][1]/float(ind2) << "  " << ind1 << "  " << ind2 << endl;
    }
    */
  }
  
  // timestep | MP | "MP"
  // timestep | i | pos[i][0] | pos[i][1] | pos[i][2] | charge[i] | "Pos"
  //
  // MP: number of positive ion moves accepted in that timestep
  // i: index of ion
  // pos[i][2]: position along inter-electrode distance
  if(outind==3){
    cout << ts << "  " << MP << "  MP" << endl;
    MP=0;
    for(int i = 0 ; i < Nion ; i++){
     cout << ts << "  " << i << "  " << pos[i][0] << "  " << pos[i][1] << "  " << pos[i][2] << "  " << charge[i] <<  "  Pos" << endl;
    }
  }

  // output total number of ions absorbed since last output
  if(outind==4){
    cout << CP << endl;
    CP = 0;
  }

}


