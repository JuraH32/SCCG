#include "MatchFinder.h"
#include <iostream>
#include <algorithm>

MatchFinder::MatchFinder(const KmerIndex& kmer_index_param,
                         const std::string& reference_sequence_param,
                         const std::string& target_sequence_param,
                         int k_val_param)
        : ref_kmer_index(kmer_index_param),
          reference_seq(reference_sequence_param),
          target_seq(target_sequence_param),
          k_mer_size(k_val_param) {
    // Validate k_mer_size
    if (this->k_mer_size <= 0) {
        // Using this->k_mer_size as k_mer_size is now a member
        std::cerr << "Error: K-mer size must be positive. Using default of 12 instead of " << k_val_param << "." << std::endl;
        this->k_mer_size = 12;
    }
}

std::vector<Match> MatchFinder::findSignificantMatches() {
    std::vector<Match> all_matches;
    size_t current_target_idx = 0;

    // Ensure k_mer_size (as size_t) is valid and get effective minimum match length
    const size_t current_k_mer_size_t = static_cast<size_t>(this->k_mer_size);
    const size_t effective_min_match_len = std::max(DEFAULT_MIN_MATCH_LEN, current_k_mer_size_t);

    // Basic sanity checks for sequence lengths
    // Target or reference must be long enough for at least one k-mer and one minimal match
    if (target_seq.length() < effective_min_match_len || reference_seq.length() < effective_min_match_len ||
        target_seq.length() < current_k_mer_size_t || reference_seq.length() < current_k_mer_size_t) {
        // Not enough characters in one or both sequences to find meaningful matches
        return all_matches;
    }

    // Loop while there's enough space in the target sequence for a potential minimal match
    while (current_target_idx + effective_min_match_len <= target_seq.length()) {
        // Stores the best match found starting at current_target_idx
        Match best_match_this_iteration(0, current_target_idx, 0);

        // Ensure we can extract a k-mer from the current position in the target sequence
        if (current_target_idx + current_k_mer_size_t > target_seq.length()) {
            break;
        }

        std::string current_kmer = target_seq.substr(current_target_idx, current_k_mer_size_t);
        const std::vector<size_t>* kmer_hits_in_ref = ref_kmer_index.findKmer(current_kmer);

        if (kmer_hits_in_ref) {
            // Iterate through all positions in the reference where this k-mer was found
            for (size_t ref_start_pos : *kmer_hits_in_ref) {
                // Extend this k-mer match as far as possible
                size_t current_match_len = current_k_mer_size_t; // Start with k-mer length

                // Continue extending while characters match and we are within bounds
                while (ref_start_pos + current_match_len < reference_seq.length() &&
                       current_target_idx + current_match_len < target_seq.length() &&
                       reference_seq[ref_start_pos + current_match_len] == target_seq[current_target_idx + current_match_len]) {
                    current_match_len++;
                }

                // If this extended match is better (longer) than what we've found so far for current_target_idx
                if (current_match_len > best_match_this_iteration.length) {
                    best_match_this_iteration.ref_pos = ref_start_pos;
                    // best_match_this_iteration.target_pos is already current_target_idx
                    best_match_this_iteration.length = current_match_len;
                }
            }
        }

        // After checking all k-mer hits, if a significant match was found:
        if (best_match_this_iteration.length >= effective_min_match_len) {
            all_matches.push_back(best_match_this_iteration);
            current_target_idx += best_match_this_iteration.length; // Advance past the found match
        } else {
            // No significant match found starting at current_target_idx.
            // Advance by 1 to check the next position in the target sequence.
            // The character at current_target_idx will be part of a mismatch region.
            current_target_idx++;
        }
    }
    return all_matches;
}
