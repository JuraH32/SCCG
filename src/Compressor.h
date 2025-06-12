#ifndef SCCG_COMPRESSOR_H
#define SCCG_COMPRESSOR_H

#include <string>
#include <vector>
#include <cstddef> // For size_t
#include <functional>
#include "Position.h"
#include "MatchFinder.h"


class Compressor {
public:
    Compressor(); // Default constructor

    void compress(const std::string& ref_fasta_path,
                  const std::string& target_fasta_path,
                  const std::string& output_intermediate_path,
                  int kmer_size,
                  int k0,
                  size_t segment_length,
                  int search_range_limit,
                  double mismatch_threshold_T1, // Changed from int to double
                  int mismatch_threshold_T2,
                  bool only_global = false);

private:
    std::string meta_data;
    bool saveToFile(const std::string& filename, const std::string& content);
    std::string readSequenceFromFile(const std::string& filename);
    [[nodiscard]] std::vector<Position> getPositions(
        const std::string& sequence,
        std::function<bool(char)> is_of_interest
    ) const;
    void savePositionsToFile(const std::string& filename, const std::vector<Position>& positions, const std::string& label_prefix); // Added label_prefix
    void saveAlignmentSegmentsToFile(
        const std::string& filename, 
        const std::vector<AlignmentSegment>& segments,
        int ref_segment_offset // Added for local segment processing
    );
    void postprocess(std::string& temp_file, const std::string& output_path);
};

#endif // SCCG_COMPRESSOR_H
