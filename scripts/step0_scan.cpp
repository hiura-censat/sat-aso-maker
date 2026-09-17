#include <zlib.h>
#include <emmintrin.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <stdexcept>
using namespace std;
using I=pair<size_t,size_t>;
int main(int argc,char**argv){try{
 if(argc!=5) throw runtime_error("usage: scan FASTA region-dir hap output");
 string dir=argv[2],hap=argv[3];
 map<string,map<string,vector<I>>> regions;
 for(string family:{"HSat2","HSat3","HSat23","ambiguous","background"}){
  ifstream f(dir+"/"+family+".bed"); if(!f) throw runtime_error("missing bed");
  string c,label;size_t s,e;while(f>>c>>s>>e>>label)regions[c][family].push_back({s,e});
 }
 map<string,size_t> lengths; ifstream fai(string(argv[1])+".fai"); string line;
 while(getline(fai,line)){istringstream f(line);string c;size_t n;f>>c>>n;if(regions.count(c)) lengths[c]=n;}
 ofstream out(argv[4]);out<<"hap\tchromosome\tregion\tinterval_count\tbp\tACGT_bp\tnon_ACGT_bp\tvalid_16mer_starts\n";
 string name,seq;size_t seen=0;
 auto emit=[&](){
  if(!regions.count(name))return;
  if(seq.size()!=lengths.at(name))throw runtime_error("length mismatch "+name);
  ++seen;
  for(string fam:{"HSat2","HSat3","HSat23","ambiguous","background"}){
   size_t bp=0,acgt=0,valid=0;auto &iv=regions[name][fam];
   for(auto x:iv){size_t run=0;bp+=x.second-x.first;
    for(size_t i=x.first;i<x.second;i++){
     if(i+16<=x.second){
      __m128i v=_mm_and_si128(_mm_loadu_si128((const __m128i*)(seq.data()+i)),_mm_set1_epi8(0xDF));
      __m128i ok16=_mm_or_si128(_mm_or_si128(_mm_cmpeq_epi8(v,_mm_set1_epi8('A')),_mm_cmpeq_epi8(v,_mm_set1_epi8('C'))),_mm_or_si128(_mm_cmpeq_epi8(v,_mm_set1_epi8('G')),_mm_cmpeq_epi8(v,_mm_set1_epi8('T'))));
      if(_mm_movemask_epi8(ok16)==65535){acgt+=16;valid+=(run>=15?16:run+1);run+=16;i+=15;continue;}
     }
     char c=seq[i];bool ok=c=='A'||c=='C'||c=='G'||c=='T'||c=='a'||c=='c'||c=='g'||c=='t';if(ok){acgt++;if(++run>=16)valid++;}else run=0;}
   }
   out<<hap<<'\t'<<name<<'\t'<<fam<<'\t'<<iv.size()<<'\t'<<bp<<'\t'<<acgt<<'\t'<<bp-acgt<<'\t'<<valid<<'\n';
  }
 };
 gzFile gz=gzopen(argv[1],"rb");if(!gz)throw runtime_error("open fasta");gzbuffer(gz,1<<20);
 vector<char> buf(1<<20);bool header=false,bol=true;string h;int n;
 while((n=gzread(gz,buf.data(),buf.size()))>0){for(int i=0;i<n;i++){char c=buf[i];
  if(bol&&c=='>'){if(!name.empty())emit();seq.clear();h.clear();header=true;bol=false;continue;}
  if(c=='\n'){if(header){name=h.substr(0,h.find_first_of(" \t\r"));header=false;}bol=true;continue;}
  if(header)h+=c;else if(c!='\r')seq+=c;bol=false;
 }}
 if(n<0)throw runtime_error("gzip read error");if(!name.empty())emit();if(gzclose(gz)!=Z_OK)throw runtime_error("gzip close error");
 if(seen!=lengths.size())throw runtime_error("missing chromosome"); if(!out)throw runtime_error("output write error");
}catch(exception&e){cerr<<e.what()<<endl;return 1;} }
