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
  std::vector<array<double,2>> v;
  string pr[]={"C","CH"};
  tt.pair_rdf_fn(0,26.0,0.1,pr,v);
  for(int i=0;i<v.size();++i) {
    std::cout << v[i][0] << " " << v[i][1] << '\n';
  }
  return 0;
}
