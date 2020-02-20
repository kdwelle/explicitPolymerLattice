#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>

#include "polymer.hpp"
#include "randomGen.cpp"
using namespace std;

Polymer::Polymer(int SideLength, int SLZ, int numChains, int sizePolymer): SideLength(SideLength), SLZ(SLZ), numChains(numChains), sizePolymer(sizePolymer)
{
  idum=-time(NULL); // todo: random number seed shouldn't be regenenerated for every new object
  gamma = 1.0;
  Ntot = SideLength*SideLength*SLZ;
  numBeads = numChains*sizePolymer;
  // numIons = int(ionCon*Ntot)*2;
  ionConc = 0.00001;
  numIons = Ntot*ionConc*2;
  numPairs = (numBeads+numIons)*(numBeads+numIons-1)/2;
  this->initialize();
  this->getAllDist();
  Eng = get_energy();
}

void Polymer::initialize(){
	/*initialization for the lattice model
    set sizes of arrays, zero them out, add ions etc.
  */

  // Build sites, ssind, occ 
  sites.reserve(Ntot); // allocate memory for the arrays now
  ssind.reserve(SideLength);
  occ.reserve(SideLength);
  neighb.reserve(Ntot);

  charge.resize(numIons);
  pairDist.resize(numPairs);
  pairEng.resize(numPairs);
  pairType.resize(numPairs);
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
  vector< vector< vector<int> > > occTemp;      
  vector< vector <int> > tempY;
  tempY.reserve(SideLength);
  vector <int> tempZ (SLZ,0); //initialize with zeroes
  for(int j = 0 ; j < SideLength ; j++){
    tempY.push_back(tempZ); // build a slice of all zeroes
  }
  for(int i = 0 ; i < SideLength ; i++){
    occ.push_back(tempY);
    occTemp.push_back(tempY);
  }

 
  int ind1,ind2,ind3,ind4;
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
      	if(ind1==SideLength) {ind1=0;} //pbc
      	neighb[ind2][count]=ssind[i][ind1][k];
      	count++;
      	ind1=k-1;
      	if(ind1<0) {ind1=k+1;} // reflecting conditions
      	neighb[ind2][count]=ssind[i][j][ind1];
      	count++;
      	ind1=k+1;
      	if(ind1==SLZ) {ind1=k-1;} // reflecting conditions
      	neighb[ind2][count]=ssind[i][j][ind1];
      }
    }
  }  
  

  bool done,added;
  int counter;
  int x0,y0,z0,x1,y1,z1;

  pos.reserve(numBeads);   // build pos matrix
  vector<int> posTemp(3);
  for (int i = 0; i < numBeads; ++i){
    pos.push_back(posTemp);
  }

  pos.reserve(numIons);   // and for the ion positions
  for (int i=0; i < numIons; ++i){
    ionPos.push_back(posTemp);
  }

  //Pair Indecies
  int pairIndex = 0;
  vector<int> tempBeads(numBeads); 
  vector<int> tempIons(numIons);
  vector<int> tempPairs(2);
  vector< vector<int> > beadsBeads;
  for (int i=0; i < numBeads-1; ++i){
    beadsBeads.push_back(tempBeads);
    for (int j=i+1; j < numBeads; ++j){
      beadsBeads[i][j] = pairIndex;
      pairInd.push_back(tempPairs);
      pairInd[pairIndex][0] = i;
      pairInd[pairIndex][1] = j;
      pairType[pairIndex] = 0;
      ++pairIndex;
    }
  }
  vector< vector<int> > ionsIons;
  for (int i=0; i < numIons-1; ++i){
    ionsIons.push_back(tempIons);
    for (int j=i+1; j < numIons; ++j){
      ionsIons[i][j] = pairIndex;
      pairInd.push_back(tempPairs);
      pairInd[pairIndex][0] = i;
      pairInd[pairIndex][1] = j;
      pairType[pairIndex] = 1;
      ++pairIndex;
    }
  }
  vector< vector<int> > beadsIons;
  for (int i=0; i < numBeads; ++i){
    beadsIons.push_back(tempIons);
    for (int j=0; j < numIons; ++j){
      beadsIons[i][j] = pairIndex;
      pairInd.push_back(tempPairs);
      pairInd[pairIndex][0] = i;
      pairInd[pairIndex][1] = j;
      pairType[pairIndex] = 2;      
      ++pairIndex;
    }
  }
  // make sure mirror images exist in getPair
  beadsBeads.push_back(tempBeads);
  for (int i=0; i < numBeads-1; ++i){
    for (int j=0i+1; j < numBeads; ++j){
      beadsBeads[j][i] = beadsBeads[i][j];
    }
    //beads-ions is not hermetian
  }
  ionsIons.push_back(tempIons);
  for (int i=0; i < numIons-1; ++i){
    for (int j=i+1; j < numIons; ++j){
      ionsIons[j][i] = ionsIons[i][j];
    }
  }
  getPair.push_back(beadsBeads);
  getPair.push_back(ionsIons);
  getPair.push_back(beadsIons);


  // build polymer
  vector< vector< vector<int> > > occTemp2 = occTemp;
  for (int i=0; i<numChains; ++i){
    done=false;
    added=true;
    counter=0;
    while(!done){
      ++counter;
      added=true;
      occTemp = occ;
      // pick a random starting site
      x0=int(ran2(&idum)*SideLength);
      y0=int(ran2(&idum)*SideLength);
      z0=int(ran2(&idum)*SLZ);

      pos[i*sizePolymer][0] =x0;
      pos[i*sizePolymer][1] =y0;
      pos[i*sizePolymer][2] =z0;

      if (occTemp[x0][y0][z0] == 0){ //site availible
        occTemp[x0][y0][z0] = 1; //mark as occupied
        for(int j = 1 ; j < sizePolymer ; ++j){ // add monomers
          //x0,y0,z0 are the end bead's current position
          x0=pos[i*sizePolymer+j-1][0];
          y0=pos[i*sizePolymer+j-1][1];
          z0=pos[i*sizePolymer+j-1][2];
          //choose a random nearest neighbor site to add on
          ind2=int(6*ran2(&idum));
          ind4=ssind[x0][y0][z0];
          ind3=neighb[ind4][ind2];

          //x1,y1,z1 are the positions of the site the ion is trying to hop to
          x1=sites[ind3][0];
          y1=sites[ind3][1];
          z1=sites[ind3][2];

          if (occTemp[x1][y1][z1] == 0){ // if the site is availible
            occTemp[x1][y1][z1] = 1;
            pos[i*sizePolymer+j][0] = x1;
            pos[i*sizePolymer+j][1] = y1;
            pos[i*sizePolymer+j][2] = z1;
          }else{
            added=false;
            break;
          }
        }
        if (added){ //only do if whole polymer is added
          cerr << "Chain # " << i << " added!" << endl;
          occ=occTemp;
          occTemp=occTemp2; //zero out occTemp
          done=true;
        }
      }
    }
  }
  bool found;
  count=0;
  vector<int> itemp(3);
  // Add positive ions randomly
  for (int i=0; i<numIons/2; ++i){
    found=false;
    while(!found){
      itemp[0]=int(ran2(&idum)*SideLength);
      itemp[1]=int(ran2(&idum)*SideLength);
      itemp[2]=int(ran2(&idum)*SLZ);
      if(occ[itemp[0]][itemp[1]][itemp[2]]==0){
        occ[itemp[0]][itemp[1]][itemp[2]]=1;
        ionPos[count][0]=itemp[0];
        ionPos[count][1]=itemp[1];
        ionPos[count][2]=itemp[2];
        charge[count]=1;
        found=true;
        count++;
      }
    }
  }
  // Add negative ions randomly
  for (int i=0; i<numIons/2; ++i){
    found=false;
    while(!found){
      itemp[0]=int(ran2(&idum)*SideLength);
      itemp[1]=int(ran2(&idum)*SideLength);
      itemp[2]=int(ran2(&idum)*SLZ);
      if(occ[itemp[0]][itemp[1]][itemp[2]]==0){
        occ[itemp[0]][itemp[1]][itemp[2]]=1;
        ionPos[count][0]=itemp[0];
        ionPos[count][1]=itemp[1];
        ionPos[count][2]=itemp[2];
        charge[count]=-1;
        found=true;
        count++;
      }
    }
  }

}


float Polymer::get_energy(){
  float eng=0.0;  
  for (int i=0; i<numPairs; ++i){
    eng += pairEng[i];
  }
  return eng;
}

void Polymer::movePolymer(int ts){ //todo: should this really be passed as a parameter?
  int ind1;
  if(ts<0){
    return;
  }
  //perform a MC sweep (make N=numMoves trial moves)
  for(int n = 0 ; n < numBeads ; ++n){
    //choose random bead to move
    ind1=int(numBeads*ran2(&idum));
    
    if (isEnd(ind1)){ //end bead
      
      if (ran2(&idum)<0.5){
        endMove(ind1);

      }else{ //slithering snake
        slitherSnek(ind1);
      }

    }else{ //internal bead

      if (ran2(&idum)<0.5){ //kink jump
        kinkJump(ind1);

      }else{ //crankshaft move
        crankshaft(ind1);
      }
    }
  }
}

void Polymer::moveIons(int ts){ //todo: same question about passing ts as a parameter
  int ind1,ind2,ind3,ind4;
  int x0,y0,z0;
  int x1,y1,z1;
  float Eng2,deltaE;
  if(ts<0){
    return;
  }
  //perform a MC sweep (make N=numIons trial moves)
  for(int n = 0; n < numIons; ++n){
    //choose random bead to move
    ind1=int(numIons*ran2(&idum));

    //x0,y0,z0 are the bead's current position
    x0=ionPos[ind1][0];
    y0=ionPos[ind1][1];
    z0=ionPos[ind1][2];

    //choose a random nearest neighbor site
    ind4 = ssind[x0][y0][z0];
    ind2=int(6*ran2(&idum));
    ind3=neighb[ind4][ind2]; //neighb list already includes pbc
    x1=sites[ind3][0];
    y1=sites[ind3][1];
    z1=sites[ind3][2];

    if (!(occ[x1][y1][z1])){ //if site is unoccupied 
      //make move
      ionPos[ind1][0]=x1;
      ionPos[ind1][1]=y1;
      ionPos[ind1][2]=z1;
      //compute new energy
      getOneDist(ind1,0); //0 means it's an ion
      Eng2=get_energy();
      deltaE=Eng2-Eng;
    
      // metropolis
      // note: || does lazy evaluation
      if(deltaE<0 || ran2(&idum)<=exp(-deltaE)){ //accept
        occ[x0][y0][z0]=0;
        occ[x1][y1][z1]=1;
        Eng=Eng2;
      }else{ //reject
        ionPos[ind1][0]=x0;
        ionPos[ind1][1]=y0;
        ionPos[ind1][2]=z0;
        getOneDist(ind1,0); //0 means it's an ion
      }
    }
  }
}

bool Polymer::endMove(int ind1){ //ind1 is the index of the bead the move is attempted on
  int ind2,ind3,ind4;
  int x0,y0,z0;
  int x1,y1,z1;
  float Eng2,deltaE;
  //x0,y0,z0 are the bead's current position
  x0=pos[ind1][0];
  y0=pos[ind1][1];
  z0=pos[ind1][2];

  //choose a random nearest neighbor site of the next most inward bead
  array<int,3> beta = get2ndFromEnd(ind1);
  ind4 = ssind[beta[0]][beta[1]][beta[2]];
  ind2=int(6*ran2(&idum));
  ind3=neighb[ind4][ind2]; //neighb list already includes pbc
  x1=sites[ind3][0];
  y1=sites[ind3][1];
  z1=sites[ind3][2];

  if (!(occ[x1][y1][z1])){ //if site is unoccupied 
    //make move
    pos[ind1][0]=x1;
    pos[ind1][1]=y1;
    pos[ind1][2]=z1;
    //compute new energy
    getOneDist(ind1,1); //1 means it's a polymer bead
    Eng2=get_energy();
    deltaE=Eng2-Eng;
  
    // metropolis
    // note: || does lazy evaluation
    if(deltaE<0 || ran2(&idum)<=exp(-deltaE)){ //accept
      occ[x0][y0][z0]=0;
      occ[x1][y1][z1]=1;
      Eng=Eng2;
      return true;
    }else{ //reject
      pos[ind1][0]=x0;
      pos[ind1][1]=y0;
      pos[ind1][2]=z0;
      getOneDist(ind1,1); //1 means it's a polymer bead
    }
  }
  return false;
}

bool Polymer::kinkJump(int ind1){
  int x0,y0,z0;
  int x1,y1,z1;
  float Eng2,deltaE;

  //x0,y0,z0 are the bead's current position
  x0=pos[ind1][0];
  y0=pos[ind1][1];
  z0=pos[ind1][2];

  //get two nearest neighbors
  array<int,3> vec1 = getVector(ind1,ind1+1);
  array<int,3> vec2 = getVector(ind1,ind1-1);

  if (!areParallel(vec1,vec2)){
    array<int,3> newPosition = addVector(ind1+1,vec2);
    x1 = newPosition[0];
    y1 = newPosition[1];
    z1 = newPosition[2];

    if (!(occ[x1][y1][z1])){ //if site is unoccupied
      //make move
      pos[ind1][0]=x1;
      pos[ind1][1]=y1;
      pos[ind1][2]=z1;
      //compute new energy
      getOneDist(ind1,1); //1 means it's a polymer bead
      Eng2=get_energy();
      deltaE=Eng2-Eng;

      //metropolis 
      if(deltaE<0 || ran2(&idum)<=exp(-deltaE)){ //accept
        occ[x0][y0][z0]=0;
        occ[x1][y1][z1]=1;
        Eng=Eng2;
        return true;
      }else{
        pos[ind1][0]=x0;
        pos[ind1][1]=y0;
        pos[ind1][2]=z0;
        getOneDist(ind1,1); //1 means it's a polymer bead
      }
    }
  }
  return false;
}

bool Polymer::crankshaft(int ind1){
  if(isEnd(ind1+1)){
    return false;
  } 

  int x0,y0,z0;
  int x01,y01,z01;
  int x1,y1,z1,x2,y2,z2;
  float Eng2,deltaE;

  x0=pos[ind1][0]; //current position
  y0=pos[ind1][1];
  z0=pos[ind1][2];

  array<int,3> vec1 = getVector(ind1,ind1+1);
  array<int,3> vec2 = getVector(ind1,ind1-1);                             
  array<int,3> vec3 = getVector(ind1+1, ind1+2);

  if (vec2 == vec3 && !areParallel(vec1,vec2)){ //can do a crankshaft

    array<int,3> crossProduct = crossVector(vec1,vec2);
    array<int,3> new1 = addVector(ind1-1,crossProduct); //new position for ind1
    array<int,3> new2 = addVector(ind1+2,crossProduct); //new position for ind1+1

    x01=pos[ind1+1][0];  //initial positions of ind1+1
    y01=pos[ind1+1][1];
    z01=pos[ind1+1][2]; 

    x1 = new1[0];
    y1 = new1[1];
    z1 = new1[2];

    x2 = new2[0];
    y2 = new2[1];
    z2 = new2[2];
    if(!occ[new1[0]][new1[1]][new1[2]] && !occ[new2[0]][new2[1]][new2[2]]){
      //make move
      pos[ind1][0] = x1;
      pos[ind1][1] = y1;
      pos[ind1][2] = z1;
      pos[ind1+1][0] = x2;
      pos[ind1+1][1] = y2;
      pos[ind1+1][2] = z2;

      //compute new energy
      getOneDist(ind1,1); //1 means it's a polymer bead
      getOneDist(ind1+1,1); //1 means it's a polymer bead
      Eng2=get_energy();
      deltaE=Eng2-Eng;
    
      //metropolis 
      if(deltaE<0 || ran2(&idum)<=exp(-deltaE)){ //accept
        occ[x0][y0][z0]=0;
        occ[x01][y01][z01]=0;
        occ[x1][y1][z1]=1;
        occ[x2][y2][z2]=1;
        Eng=Eng2;
        return true;
      }else{
        //reject
        pos[ind1][0] = x0;
        pos[ind1][1] = y0;
        pos[ind1][2] = z0;
        pos[ind1+1][0] = x01;
        pos[ind1+1][1] = y01;
        pos[ind1+1][2] = z01;
        getOneDist(ind1,1); //1 means it's a polymer bead
        getOneDist(ind1+1,1); //1 means it's a polymer bead
      }
    }
  }
  return false;
}

bool Polymer::slitherSnek(int ind1){
  int ind2,ind3,ind4;
  int x0,y0,z0;
  int x1,y1,z1;
  int x01,y01,z01;
  float Eng2,deltaE;

  x0=pos[ind1][0]; //current position
  y0=pos[ind1][1];
  z0=pos[ind1][2];

  //find random nearest neighbor of end bead
  ind4 = ssind[x0][y0][z0];
  ind2=int(6*ran2(&idum));
  ind3=neighb[ind4][ind2];

  x1=sites[ind3][0];
  y1=sites[ind3][1];
  z1=sites[ind3][2];

  //find other end of chain
  int otherInd = (ind1%sizePolymer == 0) ? ind1+sizePolymer-1 : ind1-(sizePolymer-1);
  //this is the site that would vanish if we accept the move, so set its coordinates as x0 y0 z0
  x0=pos[otherInd][0];
  y0=pos[otherInd][1];
  z0=pos[otherInd][2];

  if(!occ[x1][y1][z1]){
    //make move
    pos[otherInd][0] = x1;
    pos[otherInd][1] = y1;
    pos[otherInd][2] = z1;
    
    //compute new energy
    getOneDist(otherInd,1); //1 means it's a polymer bead
    Eng2=get_energy();
    deltaE=Eng2-Eng;
  
    //metropolis 
    if(deltaE<0 || ran2(&idum)<=exp(-deltaE)){ //accept
      occ[x0][y0][z0]=0;
      occ[x1][y1][z1]=1;
      Eng=Eng2;
      
      //rearrange order of monomer indicies to avoid confusion
      if (otherInd < ind1){ // 1->2, 2->3, etc.
        x01=x1;
        y01=y1;
        z01=z1;
        for(int i=ind1; i>otherInd; --i){
          x1=pos[i][0];
          y1=pos[i][1];
          z1=pos[i][2]; //save old coordinates
          pos[i][0]=x01;
          pos[i][1]=y01;
          pos[i][2]=z01; //move index to previous index's corrdinates
          x01=x1;
          y01=y1;
          z01=z1; //move saved coordinates over
        }
        pos[otherInd][0]=x01;
        pos[otherInd][1]=y01;
        pos[otherInd][2]=z01;

      }else{ // 2->1, 3->2, etc.
        x01=x1;
        y01=y1;
        z01=z1; //store coordinates of otherInd (new added site)
        for(int i=ind1; i<otherInd; ++i){
          x1=pos[i][0];
          y1=pos[i][1];
          z1=pos[i][2]; //save old coordinates
          pos[i][0]=x01;
          pos[i][1]=y01;
          pos[i][2]=z01; //move index to previous index's corrdinates
          x01=x1;
          y01=y1;
          z01=z1; //move saved coordinates over
        }
        pos[otherInd][0]=x01;
        pos[otherInd][1]=y01;
        pos[otherInd][2]=z01;
      }
      getAllDist();  // todo: see if this is actually nescessary
      return true;
    }else{
      pos[otherInd][0] = x0;
      pos[otherInd][1] = y0;
      pos[otherInd][2] = z0;
      getOneDist(otherInd,1); //1 means it's a polymer bead
    }
  }
  return false;
}


void Polymer::outputAnalysis(int ts, int outind){
  if (outind == 100){
    cout << numBeads+numIons << endl;
    cout << endl;
    cerr << Eng << endl;
    for(int i=0; i < numBeads; ++i){
      cout << "POL  " << pos[i][0] << "  " << pos[i][1] << "  " << pos[i][2] << endl;
    }
    for(int i=0; i < numIons; ++i){
      string ionType = (charge[i] > 0) ? "POS" : "NEG";
      cout << ionType << "  " << ionPos[i][0] << "  " << ionPos[i][1] << "  " << ionPos[i][2] << endl;
    }
  }else{
    cout << ts << endl;
  }
}

bool Polymer::isEnd(int beadIndex){ //checks if bead at index is an end bead
  if (beadIndex % sizePolymer == 0 || beadIndex % sizePolymer == sizePolymer-1){
    return true;
  }
  return false;
} 

array<int,3> Polymer::get2ndFromEnd(int beadIndex){
  array<int,3> returnCoordinates;
  int newIndex;
  if (beadIndex % sizePolymer == 0){
    newIndex = beadIndex + 1;
    
  }else if (beadIndex % sizePolymer == sizePolymer-1){
    newIndex = beadIndex - 1;

  }else{
    throw invalid_argument( "index to get2ndFromEnd must be the index of an end bead.");
  }

  returnCoordinates.at(0) = pos[newIndex][0];
  returnCoordinates.at(1) = pos[newIndex][1];
  returnCoordinates.at(2) = pos[newIndex][2];
  return returnCoordinates;

}

array<int,3> Polymer::getVector(int bead1, int bead2){ //only intended for adjacent beads
  array<int,3> solution;
  int x = pos[bead2][0] - pos[bead1][0];
  int y = pos[bead2][1] - pos[bead1][1];
  int z = pos[bead2][2] - pos[bead1][2];
  if (abs(x) >= (SideLength-1)) { x = (x>0) ? x - SideLength : x + SideLength; } //periodic boundary conditions
  if (abs(y) >= (SideLength-1)) { y = (y>0) ? y - SideLength : y + SideLength; }
  if (abs(z) >= (SLZ-1)) { z = (z>0) ? z - SLZ : z + SLZ; }
  solution.at(0) = x;
  solution.at(1) = y;
  solution.at(2) = z;
  
  return solution;
}

array<int,3> Polymer::crossVector(array<int,3> vec1, array<int,3> vec2){
  array<int,3> cross;
  cross.at(0) = int( (vec1[1]*vec2[2]) - (vec1[2]*vec2[1]) );
  cross.at(1) = int(  -(vec1[0]*vec2[2]) + (vec1[2]*vec2[0]) );
  cross.at(2) = int( (vec1[0]*vec2[1]) - (vec1[1]*vec2[0]) ); 

  return cross;
}

bool Polymer::areParallel(array<int,3> vec1, array<int,3> vec2){
  array<int,3> zero = {{0,0,0}}; 
  if (crossVector(vec1,vec2) == zero){
    return true;
  }
  return false;
}

array<int,3> Polymer::addVector(int beadIndex, array<int,3> vect){ //deals with periodic bondary conditions as well
  
  array<int,3> newPosition;

  int x = pos[beadIndex][0]; // get coordinates of bead
  int y = pos[beadIndex][1];
  int z = pos[beadIndex][2];

  int x1 = x+vect[0];
  int y1 = y+vect[1];
  int z1 = z+vect[2];

  if (z1<0){ z1 = -1*z1; }

  newPosition.at(0) = (x1 < 0) ? SideLength+x1 : x1%SideLength; //pbc
  newPosition.at(1) = (y1 < 0) ? SideLength+y1 : y1%SideLength;
  newPosition.at(2) = (z1 >= SideLength) ? (SideLength-1)-(z1%SideLength) : z1; //reflecting condition
  

  return newPosition;
}

void Polymer::getAllDist(){

  int dr[3];
  int ind1,ind2;

  for(int i=0; i<numPairs; ++i){
    ind1=pairInd[i][0];
    ind2=pairInd[i][1];
    int type = pairType[i];
    
    if (type == 0){ //bead-bead
      for(int k = 0 ; k < 3 ; k++){
        dr[k]=pos[ind1][k]-pos[ind2][k];
        if(k<2){
          if(dr[k] > SideLength/2.) dr[k]-=SideLength;
          if(dr[k] < -SideLength/2.) dr[k]+=SideLength;
        }
      }
      pairDist[i]=sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
      pairEng[i] = 1/(pairDist[i]*pairDist[i]);
    
    }else if (type == 1){ //ion-ion
      for(int k = 0 ; k < 3 ; k++){
        dr[k]=ionPos[ind1][k]-ionPos[ind2][k];
        if(k<2){
          if(dr[k] > SideLength/2.) dr[k]-=SideLength;
          if(dr[k] < -SideLength/2.) dr[k]+=SideLength;
        }
      }
      pairDist[i]=sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
      pairEng[i] = charge[ind1]*charge[ind2]/(pairDist[i]*gamma);

    }else if (type == 2){ //bead-ion
      for(int k = 0 ; k < 3 ; k++){
        dr[k]=pos[ind1][k]-ionPos[ind2][k];
        if(k<2){
          if(dr[k] > SideLength/2.) dr[k]-=SideLength;
          if(dr[k] < -SideLength/2.) dr[k]+=SideLength;
        }
      }
      pairDist[i]=sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
      pairEng[i] = charge[ind2]/(pairDist[i]*pairDist[i]);
    }   
  }
}
void Polymer::getOneDist(int index, int polymer){
  int pairIndex;
  int dr[3];
  if (polymer==1){
    int pairType = 0; //bead-bead
    for (int i=0; i<numBeads; ++i){
      if (i != index){
        for(int k = 0 ; k < 3 ; k++){
          dr[k]=pos[i][k]-pos[index][k];
          if(k<2){
            if(dr[k] > SideLength/2.) dr[k]-=SideLength;
            if(dr[k] < -SideLength/2.) dr[k]+=SideLength;
          }
        }
        pairIndex = getPair[pairType][index][i];
        pairDist[pairIndex] = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
        pairEng[pairIndex] = 1/(pairDist[pairIndex]*pairDist[pairIndex]);
      }
    }
    pairType = 2; //bead-ion
    for (int i=0; i<numIons; ++i){
      for(int k = 0 ; k < 3 ; k++){
        dr[k]=ionPos[i][k]-pos[index][k];
        if(k<2){
          if(dr[k] > SideLength/2.) dr[k]-=SideLength;
          if(dr[k] < -SideLength/2.) dr[k]+=SideLength;
        }
      }
      pairIndex = getPair[pairType][index][i];
      pairDist[pairIndex] = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
      pairEng[pairIndex] = charge[i]/(pairDist[pairIndex]*pairDist[pairIndex]);
    }

  }else{ //it's an ion
    int pairType = 2; //ion-bead
    for (int i=0; i<numBeads; ++i){
      for(int k = 0 ; k < 3 ; k++){
        dr[k]=pos[i][k]-ionPos[index][k];
        if(k<2){
          if(dr[k] > SideLength/2.) dr[k]-=SideLength;
          if(dr[k] < -SideLength/2.) dr[k]+=SideLength;
        }
      }
      pairIndex = getPair[pairType][i][index];
      pairDist[pairIndex] = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
      pairEng[pairIndex] = charge[index]/(pairDist[pairIndex]*pairDist[pairIndex]);
    }
    
    pairType = 1; //ion-ion
    for (int i=0; i<numIons; ++i){
      if (i != index){
        for(int k = 0 ; k < 3 ; k++){
          dr[k]=ionPos[i][k]-ionPos[index][k];
          if(k<2){
            if(dr[k] > SideLength/2.) dr[k]-=SideLength;
            if(dr[k] < -SideLength/2.) dr[k]+=SideLength;
          }
        }
        pairIndex = getPair[pairType][index][i];
        pairDist[pairIndex] = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
        pairEng[pairIndex] = charge[i]*charge[index]/(pairDist[pairIndex]*gamma);
      }
    } 
  }
}





