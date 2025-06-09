#include "MatchFinder.h"
#include <algorithm>
#include <vector>
#include <iostream>

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
    size_t L = target.length();
    std::vector<AlignmentSegment> segments;
    segments.reserve(L / 2);
    size_t index = 0;
    int64_t last_match_end_ref = -1;

    while (index < L - kmer_size + 1) {
        if (index + kmer_size > L) {
            if (index < L) {
                segments.emplace_back(target.substr(index, L - index), static_cast<int>(index));
            }
            break;
        }

        size_t tgt_kmer_start = index;
        std::string kmer = target.substr(tgt_kmer_start, kmer_size);
        std::vector<size_t> ref_kmer_positions_vec = kmerIndex.findKmer(kmer);

        if (ref_kmer_positions_vec.empty()) {
            // No k-mer match at all
            segments.emplace_back(target.substr(tgt_kmer_start, 1), static_cast<int>(tgt_kmer_start));
            index++;
            continue;
        }

        size_t best_ref_pos = 0, best_len = 0;
        bool found_within_range = false;

        if (global && last_match_end_ref >= 0) {
            for (size_t ref_pos : ref_kmer_positions_vec) {
                if (ref_pos + kmer_size > reference.length()) continue;
                int64_t dist = static_cast<int64_t>(ref_pos) - last_match_end_ref;

                if (std::abs(dist) <= static_cast<int64_t>(search_range)) {
                    size_t len = extendMatch(ref_pos, tgt_kmer_start);
                    if (len > best_len || (len == best_len && std::abs(dist) < std::abs(static_cast<int64_t>(best_ref_pos) - last_match_end_ref))) {
                        best_ref_pos = ref_pos;
                        best_len = len;
                        found_within_range = true;
                    }
                }
            }
        }

        if (!found_within_range) {
            for (size_t ref_pos : ref_kmer_positions_vec) {
                if (ref_pos + kmer_size > reference.length()) continue;
                size_t len = extendMatch(ref_pos, tgt_kmer_start);
                if (len > best_len || (len == best_len && ref_pos < best_ref_pos)) {
                    best_ref_pos = ref_pos;
                    best_len = len;
                }
            }
        }

        if (best_len > 0) {
            segments.emplace_back(static_cast<int>(best_ref_pos), static_cast<int>(tgt_kmer_start), static_cast<int>(best_len));
            last_match_end_ref = best_ref_pos + best_len - 1;
            index = tgt_kmer_start + best_len;
        } else {
            segments.emplace_back(target.substr(tgt_kmer_start, 1), static_cast<int>(tgt_kmer_start));
            index++;
        }
    }
    return segments;
}