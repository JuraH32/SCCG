#ifndef KMER_INDEX_H
#define KMER_INDEX_H

#include <string>
#include <vector>
#include <unordered_map>

// KmerIndex class: Builds and provides access to an index of k-mers from a sequence.
// This index maps each k-mer (a string of length k) to a list of starting positions
// where that k-mer occurs in the sequence.
class KmerIndex {
private:
    size_t k_size;
    // Pointer to the sequence being indexed. KmerIndex does not own this memory.
    const std::string* p_sequence_ = nullptr;
    // Stores the head of the linked list for each hash bucket. Indexed by scaled hash.
    std::vector<int> kmer_location_heads_;
    // Stores the next k-mer start position in a collision chain. Indexed by k-mer start position in sequence.
    std::vector<int> next_kmer_positions_;

    // A prime number for the hash table size, helps in distributing hash values.
    static const size_t MAX_HASH_TABLE_SIZE = 10000019;
    
    /**
     * @brief Hashes a k-mer string using FNV-1a algorithm
     * @param kmer The k-mer string to hash
     * @return A hash value for the k-mer
     */
    static size_t hashKmer(const std::string& kmer) ;

public:
    /**
     * @brief Constructs a KmerIndex object.
     * @param k The size of the k-mers to index. Defaults to 12 if not specified or if an invalid value is given.
     * It's recommended to choose k based on sequence characteristics and memory constraints.
     */
    KmerIndex(size_t k = 12);

    /**
     * @brief Builds the k-mer index from a given DNA/protein sequence.
     *
     * Iterates through the sequence, extracts all possible k-mers of the specified k_size,
     * and stores their starting positions in the index. If the sequence is shorter than k_size,
     * the index will remain empty.
     *
     * @param sequence The biological sequence (e.g., DNA or protein string) to index.
     */
    void buildIndex(const std::string& sequence);

    /**
     * @brief Finds all occurrences of a given k-mer in the indexed sequence.
     * @param kmer The k-mer string to search for. Its length should match k_size.
     * @return A vector of starting positions (size_t) where the k-mer is found.
     * Returns an empty vector if the k-mer is not found, if the index is not built,
     * or if kmer length does not match k_size.
     */
    std::vector<size_t> findKmer(const std::string& kmer) const;

    /**
     * @brief Gets the k-mer size used for this index.
     * @return The k-mer size (int).
     */
    size_t getKmerSize() const;
};

#endif // KMER_INDEX_H