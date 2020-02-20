#include "polymer.cpp"
#include <iostream>
#include "randomGen.cpp"


int main(){
  // long int idum=-time(NULL);

  // Make a polymer object:
  Polymer polymer(100, 100, 3, 60); //100x100x100, 2 chains, 30 beads/chain
  polymer.outputAnalysis(0,100);

  for (int i=0; i<100000; ++i){
    polymer.moveIons(i);
    if(i%100 == 0){
      polymer.movePolymer(i);
      polymer.outputAnalysis(i,100);
    }
  }
    
  

  // bool go = false;
  
  // cerr << "end wiggle test" << endl;
  // cout << "end wiggle test" << endl;
  // while(!go){
  //   go = polymer.endMove(0);
  // }
  // polymer.outputAnalysis(1,100);

  // cerr << "slither snek test 1" << endl;
  // cout << "slither snek test 1" << endl;
  // go=false;
  // while(!go){
  //   go = polymer.slitherSnek(0);
  // }
  // polymer.outputAnalysis(2,100);

  // cerr << "slither snek test 2" << endl;
  // cout << "slither snek test 2" << endl;
  // go=false;
  // while(!go){
  //   go = polymer.slitherSnek(29);
  // }
  // polymer.outputAnalysis(3,100);

  // cerr << "slither snek test 3" << endl;
  // cout << "slither snek test 3" << endl;
  // go=false;
  // while(!go){
  //   go = polymer.slitherSnek(30);
  // }
  // polymer.outputAnalysis(4,100);

  // cerr << "crankshaft test" << endl;
  // cout << "crankshaft test" << endl;
  // go=false;
  // while(!go){
  //   int ind1=int(90*ran2(&idum));
  //   if (ind1 % 30 != 0 && ind1 % 30 != 29){
  //     go = polymer.crankshaft(ind1);
  //   }
  // }
  // polymer.outputAnalysis(5,100);

  // cerr << "kink jump test" << endl;
  // cout << "kink jump test" << endl;
  // go=false;
  // while(!go){
  //   int ind1=int(90*ran2(&idum));
  //   if (ind1 % 30 != 0 && ind1 % 30 != 29){
  //     go = polymer.kinkJump(ind1);
  //   }
  // }
  // polymer.outputAnalysis(6,100);

  
}