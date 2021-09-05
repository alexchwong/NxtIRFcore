#include "FastaReader.h"

void FastaReader::SetInputHandle(std::istream *in_stream) {
  IN = in_stream;
  FirstSeq = true;
}

void FastaReader::Profile() {
  // Profiles the FASTQ file to know the chr_names and chr_lens
  chr_names.clear();
  chr_lens.clear();
  IN->seekg (0, std::ios_base::beg);
  
  while(!IN->eof() && !IN->fail()) {
    ReadSeq();
    chr_names.push_back(seqname);
    chr_lens.push_back((int32_t)sequence.length());
  }
  IN->clear();
  IN->seekg (0, std::ios_base::beg);
  FirstSeq = true;
}

bool FastaReader::ReadSeq() {
  std::string myLine;
  std::string sequence_raw;
  std::string line;
  std::string subline;
  std::string subline2;
  
  sequence.clear();
  if(FirstSeq) {
    std::getline(*IN, myLine, '>');
    FirstSeq = false;
  }
  std::getline(*IN, myLine, '\n');
  // ensure tabs and spaces are removed
  std::istringstream iss;
  iss.str(myLine);
  std::getline(iss, myLine, '\r');  // chromosome name is tab-terminated
  std::istringstream iss2;
  iss2.str(myLine);
  std::getline(iss2, myLine, '\t');  // chromosome name is tab-terminated
  std::istringstream iss3;
  iss3.str(myLine);
  std::getline(iss3, seqname, ' ');  // chromosome name is tab-terminated
  
  std::getline(*IN, sequence_raw, '>');

  sequence_raw.erase( std::remove(sequence_raw.begin(), sequence_raw.end(), ' '), sequence_raw.end() );
  sequence_raw.erase( std::remove(sequence_raw.begin(), sequence_raw.end(), '\r'), sequence_raw.end() );
  sequence_raw.erase( std::remove(sequence_raw.begin(), sequence_raw.end(), '\n'), sequence_raw.end() );
  sequence.append(sequence_raw);
  return(true);
}