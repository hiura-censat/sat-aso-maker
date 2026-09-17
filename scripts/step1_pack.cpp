#include <zlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdint>
using namespace std;
uint32_t encode(const string&s){if(s.size()!=16)throw runtime_error("k != 16");uint32_t k=0;for(char c:s){k<<=2;switch(c){case 'A':break;case 'C':k|=1;break;case 'G':k|=2;break;case 'T':k|=3;break;default:throw runtime_error("non ACGT");}}return k;}
void put(gzFile f,uint32_t k,uint64_t c){if(gzwrite(f,&k,4)!=4||gzwrite(f,&c,8)!=8)throw runtime_error("write failed");}
int main(int argc,char**argv){try{
 string mode=argc>1?argv[1]:"";
 if(mode=="pack"&&argc==3){gzFile out=gzopen(argv[2],"wb1");if(!out)throw runtime_error("open failed");string s;uint64_t n,total=0,distinct=0;while(cin>>s>>n){if(n){put(out,encode(s),n);total+=n;distinct++;}}if(!cin.eof())throw runtime_error("bad dump");if(gzclose(out)!=Z_OK)throw runtime_error("close failed");cout<<distinct<<"\t"<<total<<endl;}
 else if(mode=="pool"&&argc==4){
  struct C{uint64_t sum=0;uint64_t maximum=0;};unordered_map<uint32_t,C> all;ifstream list(argv[2]);string path;
  while(getline(list,path)){gzFile f=gzopen(path.c_str(),"rb");if(!f)throw runtime_error("open failed "+path);uint32_t k;uint64_t n;int got;
   while((got=gzread(f,&k,4))==4){if(gzread(f,&n,8)!=8)throw runtime_error("truncated count");auto &c=all[k];c.sum+=n;c.maximum=max(c.maximum,n);}
   if(got!=0||gzclose(f)!=Z_OK)throw runtime_error("read failed");
  }
  vector<uint32_t> keys;for(auto &p:all)if(p.second.sum>=100||p.second.maximum>=10)keys.push_back(p.first);sort(keys.begin(),keys.end());
  ofstream out(argv[3],ios::binary);for(auto k:keys){out.write((char*)&k,4);}if(!out)throw runtime_error("write failed");
  cout<<"{\"distinct_target_kmers\":"<<all.size()<<",\"preliminary_candidates\":"<<keys.size()<<"}"<<endl;
 }else throw runtime_error("usage: pack out.bin.gz OR pool input-list output.u32");
}catch(exception&e){cerr<<e.what()<<endl;return 1;}}
