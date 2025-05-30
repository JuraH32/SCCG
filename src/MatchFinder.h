#ifndef MATCH_FINDER_H
#define MATCH_FINDER_H

#include <string>
#include <vector>
#include "KmerIndex.h" // Assumed to be in the same directory or include path

// Structure to represent a match between reference and target sequences
struct Match {
    size_t ref_pos;    // Starting position of the match in the reference sequence
    size_t target_pos; // Starting position of the match in the target sequence
    size_t length;     // Length of the matched segment

    // Default constructor
    Match(size_t r_pos = 0, size_t t_pos = 0, size_t len = 0)
            : ref_pos(r_pos), target_pos(t_pos), length(len) {}

    // Optional: For debugging or more complex logic, an equality operator
    bool operator==(const Match& other) const {
        return ref_pos == other.ref_pos &&
               target_pos == other.target_pos &&
               length == other.length;
    }
};

class MatchFinder {
private:
    const KmerIndex& ref_kmer_index;  // Reference to the k-mer index of the reference genome
    const std::string& reference_seq; // Reference to the reference genome sequence string
    const std::string& target_seq;    // Reference to the target genome sequence string
    int k_mer_size;                   // The k-mer size used for indexing and initial seeding of matches

    // Minimum length for a match to be considered significant.
    // Matches shorter than this (or k_mer_size, whichever is larger) might be ignored.
    static const size_t DEFAULT_MIN_MATCH_LEN = 5;

public:
    /**
     * @brief Constructs a MatchFinder object.
     * @param kmer_index_param A const reference to a pre-built KmerIndex of the reference sequence.
     * @param reference_sequence_param A const reference to the reference genome sequence.
     * @param target_sequence_param A const reference to the target genome sequence.
     * @param k_val_param The k-mer size used for the KmerIndex and for seeding matches. Must be positive.
     */
    MatchFinder(const KmerIndex& kmer_index_param,
                const std::string& reference_sequence_param,
                const std::string& target_sequence_param,
                int k_val_param);

    /**
     * @brief Finds significant, non-overlapping matches between the target and reference sequences.
     *
     * This method employs a greedy strategy:
     * 1. Iterate through the target sequence.
     * 2. At each position, attempt to find the longest possible match with the reference sequence,
     * seeded by k-mers.
     * 3. If a match meeting the minimum length criteria is found, it's added to the results,
     * and the search in the target sequence jumps past this match.
     * 4. If no significant match is found, the search advances by one character in the target.
     *
     * @return A vector of Match objects, ordered by their appearance in the target sequence.
     * The regions in the target sequence not covered by these matches are considered mismatches
     * by the subsequent encoding steps.
     */
    std::vector<Match> findSignificantMatches();
};

#endif // MATCH_FINDER_H
