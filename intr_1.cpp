/*
This code has been designed by Mitradip Das, NISER, Bhubaneswar
For documentation, please look into doc.html
*/

#include "analysis.cpp"
using namespace std;

int main(int argc, char const *argv[]) {
  string pdb_name="Data/sample_1000_frames.pdb";
  long int st=line_skip(pdb_name,0,1);
  trajectory tt;
  long int pp=tt.read_frames(pdb_name,st,1000,50.0);
  std::vector<array<long,2>> v;
  tt.intr_cnt_3("C","H","O",0.9,1.5,1.5,1.8,0,180,v);
  for(int i=0;i<v.size();++i) {
    std::cout << (v[i][0]*5) << " " << v[i][1] << '\n';
  }
  return 0;
}
