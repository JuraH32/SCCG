#include "KmerIndex.h"
#include <iostream> // For potential error messages (though not used heavily here)

// Constructor for KmerIndex.
// Initializes the k-mer size. If k is not positive, it defaults to 12.
KmerIndex::KmerIndex(size_t k) : k_size(k), p_sequence_(nullptr) { // Initialize p_sequence_
    if (this->k_size <= 0) {
        std::cerr << "Warning: K-mer size must be positive. Using default of 12 instead of " << k << "." << std::endl;
        this->k_size = 12; // Fallback to a common default k-mer size.
    }
}

// Implements FNV-1a hashing algorithm for strings
// This function is static as per the header and its logic.
size_t KmerIndex::hashKmer(const std::string& kmer) {
    // FNV-1a hash algorithm constants for 64-bit
    const size_t FNV_PRIME = 1099511628211ULL;
    const size_t FNV_OFFSET_BASIS = 14695981039346656037ULL;

    size_t hash = FNV_OFFSET_BASIS;

    for (char c : kmer) {
        hash ^= static_cast<size_t>(c);
        hash *= FNV_PRIME;
    }

    return hash;
}

// Builds the k-mer index from the provided sequence using a custom hash table.
void KmerIndex::buildIndex(const std::string& sequence) {
    p_sequence_ = &sequence; // Store pointer to the sequence

    // Clear previous index data
    kmer_location_heads_.assign(MAX_HASH_TABLE_SIZE, -1);
    if (sequence.empty()) {
        next_kmer_positions_.clear();
        return;
    }
    next_kmer_positions_.assign(sequence.length(), -1);


    // Check if the sequence is long enough to extract any k-mers.
    if (sequence.length() < k_size) {
        return;
    }

    // Iterate through the sequence to extract k-mers.
    for (size_t i = 0; i <= sequence.length() - k_size; ++i) {
        // Extract the k-mer starting at position i with length k_size.
        std::string kmer_str = sequence.substr(i, k_size);

        size_t hash_val = KmerIndex::hashKmer(kmer_str);
        size_t scaled_hash = hash_val % MAX_HASH_TABLE_SIZE;

        // Add k-mer position to the custom hash table (linked list style)
        // Cast i to int because kmer_location_heads_ and next_kmer_positions_ store int
        next_kmer_positions_[i] = kmer_location_heads_[scaled_hash];
        kmer_location_heads_[scaled_hash] = static_cast<int>(i);
    }
}

// Finds a k-mer in the index and returns a vector of its occurrences.
// Returns an empty vector if the k-mer is not found or if arguments are invalid.
std::vector<size_t> KmerIndex::findKmer(const std::string& kmer) const {
    std::vector<size_t> positions;
    if (!p_sequence_ || kmer.length() != k_size || k_size == 0) {
        return positions; // Index not built, kmer length mismatch, or invalid k_size
    }

    size_t hash_val = KmerIndex::hashKmer(kmer);
    size_t scaled_hash = hash_val % MAX_HASH_TABLE_SIZE;

    int current_pos_idx = kmer_location_heads_[scaled_hash];

    while (current_pos_idx != -1) {
        // Boundary check for substr
        if (static_cast<size_t>(current_pos_idx) + k_size <= p_sequence_->length()) {
            // Verify actual k-mer to handle collisions
            if (p_sequence_->substr(static_cast<size_t>(current_pos_idx), k_size) == kmer) {
                positions.push_back(static_cast<size_t>(current_pos_idx));
            }
        }
        current_pos_idx = next_kmer_positions_[static_cast<size_t>(current_pos_idx)];
    }
    // The positions will be in reverse order of their appearance in the sequence for a given hash bucket.
    return positions;
}

// Getter for the k-mer size.
size_t KmerIndex::getKmerSize() const {
    return k_size;
}