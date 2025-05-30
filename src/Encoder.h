#ifndef ENCODER_H
#define ENCODER_H

#include <string>
#include <vector>
#include "MatchFinder.h" 

class Encoder {
   public:
      void encode(
        const std::string& reference,
        const std::string& target,
        const std::vector<Match>& matches,
        int kmer_size,
        const std::string& out_filename = "encoded.txt"
      );
};

#endif