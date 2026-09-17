#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
// Exact rolling 16-mer counter. Output: uint64 [candidate, chromosome, strand].
int main(int argc,char**argv){try{
 std::ios::sync_with_stdio(false);std::cin.tie(nullptr);
 if(argc!=5)throw std::runtime_error("usage: step3_count codes.u32 chromosomes.txt counts.u64 stats.tsv < regions.fa");
 std::ifstream input(argv[1],std::ios::binary);if(!input)throw std::runtime_error("codes open failed");
 std::unordered_map<uint32_t,uint32_t> index;std::vector<uint32_t> codes;uint32_t v;
 while(input.read(reinterpret_cast<char*>(&v),4))codes.push_back(v);
 if(!input.eof()||codes.empty())throw std::runtime_error("invalid codes");
 index.reserve(codes.size()*2);for(size_t i=0;i<codes.size();++i)if(!index.emplace(codes[i],i).second)throw std::runtime_error("duplicate code");
 std::ifstream cf(argv[2]);std::vector<std::string> names;std::unordered_map<std::string,size_t> chroms;std::string line;
 while(std::getline(cf,line))if(!line.empty()){chroms.emplace(line,names.size());names.push_back(line);}
 if(names.empty()||chroms.size()!=names.size())throw std::runtime_error("invalid chromosomes");
 std::vector<uint64_t> prefix(1<<18,0);for(auto code:codes){uint32_t p=code>>8;prefix[p>>6]|=uint64_t(1)<<(p&63);}
 size_t nc=names.size();std::vector<uint64_t> counts(codes.size()*nc*2,0),valid(nc,0),bases(nc,0),records(nc,0);
 int chr=-1,run=0;uint32_t fw=0,rc=0;
 while(std::getline(std::cin,line)){
  if(!line.empty()&&line[0]=='>'){
   std::string name=line.substr(1,line.find_first_of(" \t")-1);auto colon=name.rfind(':');if(colon!=std::string::npos)name.resize(colon);
   auto it=chroms.find(name);if(it==chroms.end())throw std::runtime_error("unknown FASTA chromosome "+name);
   chr=it->second;records[chr]++;run=0;fw=rc=0;continue;
  }
  if(chr<0&&!line.empty())throw std::runtime_error("sequence without header");
  for(unsigned char ch:line){if(ch=='\r')continue;int b=-1;switch(ch){case 'A':case 'a':b=0;break;case 'C':case 'c':b=1;break;case 'G':case 'g':b=2;break;case 'T':case 't':b=3;break;}
   bases[chr]++;if(b<0){run=0;fw=rc=0;continue;}fw=(fw<<2)|b;rc=(rc>>2)|(uint32_t(b^3)<<30);if(run<16)run++;if(run<16)continue;
   valid[chr]++;uint32_t key=std::min(fw,rc),p=key>>8;if(!(prefix[p>>6]&(uint64_t(1)<<(p&63))))continue;auto it=index.find(key);if(it!=index.end())counts[(size_t(it->second)*nc+chr)*2+(fw>rc)]++;
  }
 }
 if(!std::cin.eof())throw std::runtime_error("FASTA read failed");
 std::ofstream out(argv[3],std::ios::binary);out.write(reinterpret_cast<char*>(counts.data()),counts.size()*8);if(!out)throw std::runtime_error("counts write failed");
 std::ofstream stats(argv[4]);stats<<"chromosome\tbp\tvalid_16mer_starts\trecords\n";
 for(size_t c=0;c<nc;++c){stats<<names[c]<<'\t'<<bases[c]<<'\t'<<valid[c]<<'\t'<<records[c]<<'\n';}if(!stats)throw std::runtime_error("stats write failed");
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
