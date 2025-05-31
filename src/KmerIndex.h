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
    // The index itself: a hash map where keys are k-mers (strings)
    // and values are vectors of starting positions (size_t) of that k-mer.
    std::unordered_map<std::string, std::vector<size_t>> index;

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
     * @return A const pointer to a vector of starting positions (size_t) where the k-mer is found.
     * Returns nullptr if the k-mer is not found in the index or if its length
     * does not match k_size (though this implementation doesn't explicitly check kmer length here).
     */
    const std::vector<size_t>* findKmer(const std::string& kmer) const;

    /**
     * @brief Gets the k-mer size used for this index.
     * @return The k-mer size (int).
     */
    size_t getKmerSize() const;
};

#endif // KMER_INDEX_H
