/*
This code has been designed by Mitradip Das, NISER, Bhubaneswar
For documentation, please look into doc.html
*/

#include "cpp_include/libraries.h"
using namespace std;

class atom {
public:
  long int atm_id;
  string atm_name;
  long int res_id;
  float coord[3];

  void show() {
    std::cout << atm_id <<"\t"<<atm_name <<"\t"<<res_id;
    for(auto i=0;i<3;i++) {
      std::cout << "\t" << coord[i];
    }
    std::cout  << '\n';
  }
};

class frame {
public:
  std::vector<atom> atmx;
  double box_len;

  long int read_file(string file_name,long int start_line) {
    long int ctr=0;
    fstream pdb_file;
    pdb_file.open(file_name,ios::in);
    string line;
    atom this_atom;
    pdb_file.seekg(start_line,ios::beg);
    while(true) {
      line="";
      getline(pdb_file,line);
      if(trim(line).compare("END")==0)
        break;
      this_atom.atm_id=stoi(line.substr(6,5));
      this_atom.atm_name=trim(line.substr(12,4));
      this_atom.res_id=stoi(line.substr(22,4));
      this_atom.coord[0]=stof(line.substr(30,8));
      this_atom.coord[1]=stof(line.substr(38,8));
      this_atom.coord[2]=stof(line.substr(46,8));
      atmx.push_back(this_atom);
    }
    ctr=pdb_file.tellg();
    pdb_file.close();
    std::sort(atmx.begin(), atmx.end(), [](const atom& lhs, const atom& rhs){ return lhs.atm_id < rhs.atm_id; });
    return ctr;
  }

  void clear_all() {
    atmx.clear();
    box_len=0;
  }
  void set_box(double len) {
    box_len=len;
  }

  void show_frame() {
    for(auto i=0;i<atmx.size();++i)
      atmx[i].show();
  }


  double calc_dist_pbc(long int id1,long int id2,bool check_diff_res) {
    double dist=-1.0,diff,dst;
    if(atmx[id1].res_id!=atmx[id2].res_id or check_diff_res==false){
      dst=0;
      for(int j=0;j<3;++j){
        diff=atmx[id2].coord[j]-atmx[id1].coord[j];
        if(abs(diff)>0.5*box_len) {
          if(diff<0) diff+=box_len;
          else diff-=box_len;
        }
        dst+=pow(diff,2);
      }
      dist=sqrt(dst);
    }
    return dist;
  }

  void get_idx_of(string src,std::vector<int> &v) {
    for(int i=0;i<atmx.size();++i)
      if (src.compare(atmx[i].atm_name)==0) v.push_back(i);
  }

  void pair_rdf_fn(double start,double stop,double interval,const string pairs[],std::vector<array<double,2>> &rdf){
    int nrdf=(int)((stop-start)/interval);
    std::array <double,2> rdf_comp={0,0};
    int count=0;
    double d;
    for(int i=0;i<=nrdf;++i){
      rdf_comp={start+i*interval+interval/2.0,0.0};
      rdf.push_back({start+i*interval+interval/2.0,0.0});
    }
    std::vector<int> id1;
    std::vector<int> id2;
    for(int i=0;i<atmx.size();++i){
      if (pairs[0].compare(atmx[i].atm_name)==0) id1.push_back(i);
      if (pairs[1].compare(atmx[i].atm_name)==0) id2.push_back(i);
    }
    for(int i=0;i<id1.size();++i){
      for(int j=0;j<id2.size();++j){
	if(i!=j){
          d=calc_dist_pbc(id1[i],id2[j],true);
          count+=1;
          if(d>=start && d<=stop && d<=0.5*box_len) {
            rdf[(int)((d-start)/interval)][1]+=1.0;
	 }
        }
      }
    }
    if(count!=0) {
      for(int i=0;i<rdf.size();++i) {
        rdf[i][1]/=(M_PI*4*pow(rdf[i][0],2)*interval);
        rdf[i][1]/=count;
        rdf[i][1]*=pow(box_len,3);
      }
    }

  }

  bool interact_2(long id1,long id2,double d12l,double d12h) {
    double dist1=0,dst=0;
    for(int j=0;j<3;++j)
      dst+=pow(atmx[id2].coord[j]-atmx[id1].coord[j],2);
    dist1=sqrt(dst);
    if(dist1>=d12l and dist1<=d12h)
      return true;
    else
      return false;
  }

  bool interact_3(long id1,long id2,long id3,double d12l,double d12h,double d23l,double d23h,double a123l,double a123h){
    double dotp=0,dist1=0,dist2=0,dst1=0,dst2=0;
    double v1=0,v2=0;
    for(int j=0;j<3;++j) {
      v1=atmx[id2].coord[j]-atmx[id1].coord[j];
      v2=atmx[id2].coord[j]-atmx[id3].coord[j];
      dotp+=v1*v2;
      dst1+=pow(v1,2);
      dst2+=pow(v2,2);
    }
    dist1=sqrt(dst1);
    dist2=sqrt(dst2);
    double ang=acos(dotp/(dist1*dist2))*180/M_PI;
    int qty=0;
    if(dist1>=d12l and dist1<=d12h){
      qty+=1;
    }
    if(dist2>=d23l and dist2<=d23h){
      qty+=1;
    }
    if(ang>=a123l and ang<=a123h){
      qty+=1;
    }
    if(qty==3)
      return true;
    else
      return false;
  }

  long cnt_intr_2(const string atm1,const string atm2,double d12l,double d12h) {
    std::vector<int> id1,id2;

    get_idx_of(atm1,id1);
    get_idx_of(atm2,id2);
    long ctr=0;
    for(auto i=0;i<id1.size();++i) {
      for(auto j=0;j<id2.size();++j) {
        if(interact_2(id1[i],id2[j],d12l,d12h)) {
          ++ctr;
        }
      }
    }

    return ctr;
  }

  long cnt_intr_3(const string atm1,const string atm2,const string atm3,double d12l,double d12h,double d23l,double d23h,double a123l,double a123h) {
    std::vector<int> id1,id2,id3;

    get_idx_of(atm1,id1);
    get_idx_of(atm2,id2);
    get_idx_of(atm3,id3);
    long ctr=0;
    for(auto i=0;i<id1.size();++i) {
      for(auto j=0;j<id2.size();++j) {
        if(interact_2(id1[i],id2[j],d12l,d12h)) {
          for(auto k=0;k<id3.size();++k) {
            if(interact_3(id1[i],id2[j],id3[k],d12l,d12h,d23l,d23h,a123l,a123h))
              ++ctr;
          }
        }
      }
    }
    return ctr;
  }
};

class trajectory {
public:
  std::vector<frame> framx;


  long int read_frames(string file_name,long int st_pos,long int nframes,double boxl){
    frame my_frame;
    long int tsps=st_pos;
    for(long int i=0;i<nframes;++i) {
      tsps=my_frame.read_file(file_name,tsps);
      my_frame.set_box(boxl);
      framx.push_back(my_frame);
      my_frame.clear_all();
    }
    return tsps;
  }

  void clear_all() {
    framx.clear();
  }

  void pair_rdf_fn(double start,double stop,double interval,const string pairs[],std::vector<array<double,2>> &rdf){
    int nrdf=(int)((stop-start)/interval);
    std::array <double,2> rdf_comp={0,0};
    int count=0;
    double d;
    std::vector<array<double,2>> v;
    for(int i=0;i<=nrdf;++i){
      rdf_comp={start+i*interval+interval/2.0,0.0};
      rdf.push_back({start+i*interval+interval/2.0,0.0});
    }
    for(int i=0;i<framx.size();++i) {

      framx[i].pair_rdf_fn(start,stop,interval,pairs,v);
      for(int j=0;j<rdf.size();++j){
        rdf[j][1]+=(v[j][1]/framx.size());
      }
      v.clear();
    }
  }

  double idx_msd_fn(int dt,const int idx) {
    double msd=0;
    int ctr=0;
    for(int i=0;i<framx.size()-dt;++i) {
      ++ctr;
      for(int j=0;j<3;++j)
        msd+=(pow(framx[i].atmx[idx].coord[j]-framx[i+dt].atmx[idx].coord[j],2)/(framx.size()-dt));
    }
    return msd;
  }

  double type_msd_fn(int dt,const string src) {
    double msd=0;
    std::vector<int> idx;
    framx[0].get_idx_of(src,idx);
    for(int i=0;i<idx.size();++i)
      msd+=(idx_msd_fn(dt,idx[i])/idx.size());
    return msd;
  }

  void type_msd_fn_all_dt(const string src,std::vector<array<double,2>> &v) {
    array<double,2> ar;
    for(int dt=1;dt<framx.size();++dt) {
      ar[0]=dt;
      ar[1]=type_msd_fn(dt,src);
      v.push_back(ar);
    }
  }

  void intr_cnt_3(const string atm1,const string atm2,const string atm3,double d12l,double d12h,double d23l,double d23h,double a123l,double a123h,std::vector<array<long,2>> &v) {
    array<long,2> indv_info;
    for(auto t=0;t<framx.size();++t){
      
      indv_info[0]=t;
      indv_info[1]=framx[t].cnt_intr_3(atm1,atm2,atm3,d12l,d12h,d23l,d23h,a123l,a123h);
      v.push_back(indv_info);
    }
  }
  void intr_cnt_2(const string atm1,const string atm2,double d12l,double d12h,std::vector<array<long,2>> &v) {
    array<long,2> indv_info;
    for(auto t=0;t<framx.size();++t){
      indv_info[0]=t;
      indv_info[1]=framx[t].cnt_intr_2(atm1,atm2,d12l,d12h);
      v.push_back(indv_info);
    }
  }
};
