#include "MatchFinder.h"
#include <algorithm>
#include <vector> // Required for std::vector

const size_t MatchFinder::DEFAULT_MIN_MATCH_LEN;

MatchFinder::MatchFinder(const std::string& ref, const std::string& tgt, 
                        size_t k_size, size_t range)
    : reference(ref), target(tgt), kmer_size(k_size), search_range(range), kmerIndex(k_size) { // Initialize kmerIndex here
    buildKmerIndex(reference);
}

void MatchFinder::buildKmerIndex(const std::string& sequence) {
    // KmerIndex is already initialized with kmer_size in MatchFinder's constructor initializer list.
    // kmerIndex = KmerIndex(kmer_size); // This would re-assign, not ideal. Better to init in MIL.
    kmerIndex.buildIndex(sequence);
}

size_t MatchFinder::extendMatch(size_t ref_pos, size_t target_pos) {
    size_t length = 0;
    size_t ref_length = reference.length();
    size_t target_length = target.length();

    // Extend the match as long as characters match and within bounds
    while (ref_pos + length < ref_length && target_pos + length < target_length &&
           reference[ref_pos + length] == target[target_pos + length]) {
        length++;
    }

    return length;
}

std::vector<AlignmentSegment> MatchFinder::findMatches(bool global) {
    std::vector<AlignmentSegment> segments; // Changed from std::vector<Match> matches
    size_t L = target.length();
    size_t index = 0;
    size_t last_match_end_ref = 0;  // For global matching, tracks end in reference

    while (index < L) { // Loop until end of target
        if (index + kmer_size > L) { // Not enough characters left for a k-mer
            if (index < L) { // Remaining characters are mismatches
                segments.push_back(AlignmentSegment(
                    target.substr(index, L - index),
                    static_cast<int>(index)
                ));
            }
            break; // End of target processing
        }

        size_t current_target_kmer_start = index;
        std::string current_kmer = target.substr(current_target_kmer_start, kmer_size);
        const std::vector<size_t>* kmer_positions_ptr = kmerIndex.findKmer(current_kmer);
        
        // If no k-mer match found in hash table or kmer_positions_ptr is null
        if (!kmer_positions_ptr || kmer_positions_ptr->empty()) {
            segments.push_back(AlignmentSegment(
                target.substr(current_target_kmer_start, 1),
                static_cast<int>(current_target_kmer_start)
            ));
            index++;
            continue;
        }
        const std::vector<size_t>& kmer_positions = *kmer_positions_ptr;


        size_t lmax1 = 0, lmax2 = 0;
        size_t pn1 = 0, pn2 = 0;
        // ln1 and ln2 will store the actual match lengths found
        size_t best_ln1 = 0, best_ln2 = 0; 


        // Check all matching k-mer positions
        for (size_t ref_kmer_start_pos : kmer_positions) {
            // Ensure ref_kmer_start_pos is valid before calling extendMatch
            if (ref_kmer_start_pos + kmer_size > reference.length()) continue;

            size_t current_match_len = extendMatch(ref_kmer_start_pos, current_target_kmer_start);

            if (global && last_match_end_ref > 0) {
                // Using ref_kmer_start_pos for distance calculation
                int64_t distance = static_cast<int64_t>(ref_kmer_start_pos) - static_cast<int64_t>(last_match_end_ref);
                if (std::abs(distance) <= static_cast<int64_t>(search_range)) {
                    if (current_match_len == lmax2) { // lmax2 stores best length for global constrained
                        if (pn2 == 0 || // if pn2 is not set yet
                            std::abs(static_cast<int64_t>(ref_kmer_start_pos) - static_cast<int64_t>(last_match_end_ref)) <
                            std::abs(static_cast<int64_t>(pn2) - static_cast<int64_t>(last_match_end_ref))) {
                            pn2 = ref_kmer_start_pos;
                            // best_ln2 = current_match_len; // lmax2 already holds this length
                        }
                    } else if (current_match_len > lmax2) {
                        lmax2 = current_match_len;
                        pn2 = ref_kmer_start_pos;
                        // best_ln2 = current_match_len;
                    }
                }
            }

            // General best match (lmax1)
            if (current_match_len == lmax1) { // lmax1 stores best length overall
                 if (pn1 == 0 || // if pn1 is not set yet
                    (global && last_match_end_ref > 0 && std::abs(static_cast<int64_t>(ref_kmer_start_pos) - static_cast<int64_t>(last_match_end_ref)) <
                     std::abs(static_cast<int64_t>(pn1) - static_cast<int64_t>(last_match_end_ref))) ||
                    (!global && ref_kmer_start_pos < pn1) // Simple tie-break for non-global: smallest ref_pos
                 ) {
                    pn1 = ref_kmer_start_pos;
                    // best_ln1 = current_match_len; // lmax1 already holds this length
                }
            } else if (current_match_len > lmax1) {
                lmax1 = current_match_len;
                pn1 = ref_kmer_start_pos;
                // best_ln1 = current_match_len;
            }
        }
        
        // Determine final pn and ln based on global flag and findings
        size_t final_pn = 0;
        size_t final_ln = 0;

        if (global && pn2 != 0 && lmax2 > 0) { // Check lmax2 > 0 as well
            final_pn = pn2;
            final_ln = lmax2;
        } else if (pn1 != 0 && lmax1 > 0) { // Check lmax1 > 0
            final_pn = pn1;
            final_ln = lmax1;
        }
        // If neither condition met, final_ln remains 0

        if (final_ln > 0) {
            segments.push_back(AlignmentSegment(
                static_cast<int>(final_pn),
                static_cast<int>(current_target_kmer_start),
                static_cast<int>(final_ln)
            ));
            last_match_end_ref = final_pn + final_ln; // Update for next global search
            index = current_target_kmer_start + final_ln; // Advance index past the match
        } else {
            // No effective match found (ln is 0), treat as mismatch of one char
            segments.push_back(AlignmentSegment(
                target.substr(current_target_kmer_start, 1),
                static_cast<int>(current_target_kmer_start)
            ));
            index = current_target_kmer_start + 1; // Advance index by 1
        }
    }
    return segments;
}
