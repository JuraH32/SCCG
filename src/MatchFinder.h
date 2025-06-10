#ifndef MATCH_FINDER_H
#define MATCH_FINDER_H

#include <string>
#include <vector>
#include <unordered_map>
#include "KmerIndex.h"
#include "Position.h" // Include the Position struct

// Defines the result of a match or a segment of mismatch
struct AlignmentSegment {
    Position region;             // Coordinates of the segment.
                                 // For matches: all 4 coords are relevant.
                                 // For mismatches: only target coords are relevant, ref coords are -1.
    std::string mismatched_sequence; // Stores the sequence of the mismatch if is_match is false.
    bool is_match;               // True if this entry represents a match, false if it's a mismatch block.

    // Constructor for a successful match
    AlignmentSegment(int ref_start_pos, int tar_start_pos, int length)
        : mismatched_sequence(""), is_match(true) {
        if (length <= 0) {
            is_match = true; // Assuming length > 0 based on call site logic
        }
        region.setStartInReference(ref_start_pos);
        region.setEndInReference(ref_start_pos + length - 1);
        region.setStartInTarget(tar_start_pos);
        region.setEndInTarget(tar_start_pos + length - 1);
    }

    AlignmentSegment()
        : mismatched_sequence(""), is_match(true) {
        // Default constructor initializes to a match with empty sequence
        region.setStartInReference(-1);
        region.setEndInReference(-1);
        region.setStartInTarget(-1);
        region.setEndInTarget(-1);
    }

    // Constructor for a mismatch
    AlignmentSegment(const std::string& mismatch_str, int tar_start_pos)
        : mismatched_sequence(mismatch_str), is_match(false) {
        region.setStartInTarget(tar_start_pos);
        if (!mismatch_str.empty()) {
            region.setEndInTarget(tar_start_pos + static_cast<int>(mismatch_str.length()) - 1);
        } else {
            region.setEndInTarget(tar_start_pos - 1); // Or tar_start_pos if 0-length mismatch means single point
        }
        // Reference positions in 'region' will remain -1 (default from Position constructor)
    }
};

class MatchFinder {
public:
    MatchFinder(const std::string& reference, const std::string& target, 
                size_t kmer_size, size_t search_range = 100); // search_range default was in .cpp, better in .h
    
    std::vector<AlignmentSegment> findMatches(bool global);

private:
    const std::string& reference;
    const std::string& target;
    size_t kmer_size;
    size_t search_range;
    KmerIndex kmerIndex; // Was std::unordered_multimap<std::string, size_t> kmer_hash;

    void buildKmerIndex(const std::string& sequence); // Renamed from buildKmerHash
    size_t extendMatch(size_t ref_pos, size_t target_pos);

    static const size_t DEFAULT_MIN_MATCH_LEN = 3; // Example, if needed
};

#endif // MATCH_FINDER_H
