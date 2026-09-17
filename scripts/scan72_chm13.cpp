// One-pass Aho-Corasick exact 16-mer scanner for gzipped CHM13 FASTA.
#include <zlib.h>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
struct Hit {int q; char strand;};
struct Node {array<int,4> next{{-1,-1,-1,-1}};int fail=0;vector<Hit> out;};
static int base(char c){switch(c){case 'A':case 'a':return 0;case 'C':case 'c':return 1;case 'G':case 'g':return 2;case 'T':case 't':return 3;default:return -1;}}
static void writegz(gzFile g,const string& s){if(gzwrite(g,s.data(),s.size())!=(int)s.size()){cerr<<"gzip write failed\n";exit(4);}}
int main(int argc,char**argv){
 if(argc!=4){cerr<<"usage: scan72_chm13 queries.tsv fasta.gz hits.tsv.gz\n";return 2;}
 ifstream in(argv[1]);if(!in){cerr<<"queries unavailable\n";return 2;}
 string line;getline(in,line);vector<string> ids;vector<Node> trie(1);
 auto insert=[&](const string& motif,int qi,char strand){int n=0;for(char c:motif){int b=base(c);if(b<0){cerr<<"non-ACGT motif\n";exit(3);}if(trie[n].next[b]<0){trie[n].next[b]=trie.size();trie.emplace_back();}n=trie[n].next[b];}trie[n].out.push_back({qi,strand});};
 while(getline(in,line)){
  if(line.empty())continue;vector<string> fields;stringstream ss(line);string f;while(getline(ss,f,'\t'))fields.push_back(f);
  if(fields.size()<4 || fields[2].size()!=16 || fields[3].size()!=16){cerr<<"bad query line\n";return 3;}
  int qi=ids.size();ids.push_back(fields[1]);insert(fields[2],qi,'+');if(fields[2]!=fields[3])insert(fields[3],qi,'-');
 }
 queue<int> todo;for(int b=0;b<4;b++){int v=trie[0].next[b];if(v<0)trie[0].next[b]=0;else{trie[v].fail=0;todo.push(v);}}
 while(!todo.empty()){int u=todo.front();todo.pop();for(int b=0;b<4;b++){int v=trie[u].next[b];if(v<0)trie[u].next[b]=trie[trie[u].fail].next[b];else{trie[v].fail=trie[trie[u].fail].next[b];auto &o=trie[trie[v].fail].out;trie[v].out.insert(trie[v].out.end(),o.begin(),o.end());todo.push(v);}}}
 gzFile fa=gzopen(argv[2],"rb");gzFile out=gzopen(argv[3],"wb6");if(!fa||!out){cerr<<"cannot open FASTA/output\n";return 2;}
 writegz(out,"kmer_id\tchromosome\tstart_0based\tend_0based\tstrand_of_canonical_target\n");
 char buf[1<<20];string chrom,header;long long pos=0,total=0;int state=0;bool line_start=true,in_header=false;
 while(true){int got=gzread(fa,buf,sizeof(buf));if(got<0){cerr<<"gzip read failed\n";return 4;}if(!got)break;
  for(int i=0;i<got;i++){
   char c=buf[i];if(line_start&&c=='>'){in_header=true;header.clear();state=0;pos=0;line_start=false;continue;}
   if(c=='\n'||c=='\r'){
    if(in_header&&c=='\n'){size_t end=header.find_first_of(" \t");chrom=header.substr(0,end);cerr<<"scan "<<chrom<<"\n";in_header=false;}
    if(c=='\n')line_start=true;continue;
   }
   line_start=false;if(in_header){header.push_back(c);continue;}
   int b=base(c);state=(b<0)?0:trie[state].next[b];
   if(b>=0&&!trie[state].out.empty())for(const Hit& h:trie[state].out){long long start=pos-15;if(start>=0){string row=ids[h.q]+"\t"+chrom+"\t"+to_string(start)+"\t"+to_string(start+16)+"\t"+h.strand+"\n";writegz(out,row);total++;}}
   pos++;
  }
 }
 gzclose(fa);gzclose(out);cerr<<"DONE queries="<<ids.size()<<" hits="<<total<<"\n";return 0;
}
