#include "KmerIndex.h"
#include <iostream> // For potential error messages (though not used heavily here)

// Constructor for KmerIndex.
// Initializes the k-mer size. If k is not positive, it defaults to 12.
KmerIndex::KmerIndex(size_t k) : k_size(k) {
    if (this->k_size <= 0) {
        std::cerr << "Warning: K-mer size must be positive. Using default of 12 instead of " << k << "." << std::endl;
        this->k_size = 12; // Fallback to a common default k-mer size.
    }
}

// Builds the k-mer index from the provided sequence.
void KmerIndex::buildIndex(const std::string& sequence) {
    // Clear any existing index data to allow re-indexing with a new sequence.
    index.clear();

    // Check if the sequence is long enough to extract any k-mers.
    if (sequence.length() < static_cast<size_t>(k_size)) {
        return;
    }

    // Iterate through the sequence to extract k-mers.
    for (size_t i = 0; i <= sequence.length() - static_cast<size_t>(k_size); ++i) {
        // Extract the k-mer starting at position i with length k_size.
        std::string kmer = sequence.substr(i, static_cast<size_t>(k_size));
        // Add the starting position 'i' to the list of occurrences for this k-mer.
        // If the k-mer is not yet in the map, it will be added.
        index[kmer].push_back(i);
    }
}

// Finds a k-mer in the index and returns a pointer to its list of occurrences.
// Returns nullptr if the k-mer is not found.
const std::vector<size_t>* KmerIndex::findKmer(const std::string& kmer) const {
    // Use 'find' to search for the k-mer in the unordered_map.
    auto it = index.find(kmer);
    // Check if the k-mer was found.
    if (it != index.end()) {
        // If found, return a pointer to the vector of positions.
        return &(it->second);
    }
    // If not found, return nullptr.
    return nullptr;
}

// Getter for the k-mer size.
size_t KmerIndex::getKmerSize() const {
    return k_size;
}
