#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <algorithm>
#include <limits>
#include <cctype>
#include <sys/stat.h>
#include <cstdio>
#include <cstdlib>

// Kmer class used to store k-mer information
class newkmer {
private:
    int kmerstart; // Starting position of the k-mer in the sequence
    std::string kmer; // The k-mer string itself

public:
    void setkmerstart(int start) { kmerstart = start; }

    int getkmerstart() const { return kmerstart; }

    void setkmer(const std::string &k) { kmer = k; }

    std::string getkmer() const { return kmer; }
};

// Position class used to store start and end positions of matches in target and reference sequences
class Position {
private:
    int startinTar, endinTar, startinRef, endinRef; // Coordinates for matches

public:
    void setstartinTar(int start) { startinTar = start; }

    int getstartinTar() const { return startinTar; }

    void setendinTar(int end) { endinTar = end; }

    int getendinTar() const { return endinTar; }

    void setstartinRef(int start) { startinRef = start; }

    int getstartinRef() const { return startinRef; }

    void setendinRef(int end) { endinRef = end; }

    int getendinRef() const { return endinRef; }
};

// Enum to define the scale of processing (local or global)
enum ScaleType {
    LOCAL = 0,  // Processing is done in segments/locally
    GLOBAL = 1 // Processing is done on the entire sequence
};

// Enum to define the type of sequence (reference or target)
enum SequenceType {
    REFERENCE = 0, // The sequence is a reference genome
    TARGET = 1    // The sequence is a target genome to be compressed
};


class SequenceProcessor {
public:
    // Hashmap to store k-mers for local matching; maps hash value to a list of k-mers
    static std::unordered_map<int, std::vector<newkmer> > hashmap;
    // Vector to store index of the next k-mer with the same hash in the reference sequence (for global matching)
    static std::vector<int> next_kmer;
    // Location of the head (start index) of each k-mer in the reference sequence (for global matching)
    static std::vector<int> kmer_location;
    static const int maxchar = 268435456; // 2^28
    static const int maxseq = 268435456 * 2; // 2^29

    // Calculates a hash code for a given string
    static int stringHashCode(const std::string &str) {
        int hash = 0;
        for (char c: str) {
            hash = 31 * hash + static_cast<int>(c); // Common string hashing algorithm
        }
        return hash;
    }

    // Builds a hash table for k-mers from the input sequence (for local matching)
    static void buildLhashtable(const std::string &read,
                                int kmer_length_param) {
        int current_length = read.length(), i = kmer_length_param;
        std::string Nkmer = "";
        while (i > 0) {
            Nkmer += "N";
            i--;
        }

        int Nkey = stringHashCode(Nkmer); // Hash code for a k-mer of all 'N's
        i = 0;

        // Iterate through the sequence to extract and hash k-mers
        while (i < current_length - kmer_length_param + 1) {
            std::string kmer_str = read.substr(i, kmer_length_param);
            newkmer newKMer;
            newKMer.setkmerstart(i);
            newKMer.setkmer(kmer_str);
            int key = stringHashCode(kmer_str);

            // Add the k-mer to the hashmap
            if (hashmap.find(key) != hashmap.end()) { // If key exists, append to list
                std::vector<newkmer> &list = hashmap[key];
                list.push_back(newKMer);
            } else { // If key doesn't exist, create new list
                std::vector<newkmer> list;
                list.push_back(newKMer);
                hashmap[key] = list;
            }
            i++;
            // Skip consecutive 'N' or 'n' characters efficiently
            if (key == Nkey) {
                while (i < current_length - kmer_length_param + 1 &&
                       (read.substr(i, 1) == "n" || read.substr(i, 1) == "N")) {
                    i++;
                }
            }
        }
    }

    // Builds a hash table for k-mers from the input sequence (for global matching)
    static void buildGhashtable(const std::string &read, int kmer_length_param) {
        int iteration = read.length() - kmer_length_param + 1; // Total number of k-mers

        // Initialize kmer_location array with -1 (indicating no k-mer at that hash yet)
        for (int i = 0; i < maxseq; i++) {
            kmer_location[i] = -1;
        }
        for (int i = 0; i < iteration; i++) {
            std::string kmer_str = read.substr(i, kmer_length_param); // Extract k-mer
            long key = std::abs(static_cast<long>(stringHashCode(kmer_str))); // Calculate hash and take absolute value

            if (key == -2147483648LL) { // Handle Integer.MIN_VALUE case for abs
                key = 0;
            }

            // Ensure key is within the bounds of kmer_location array
            while (key > maxseq - 1) {
                key = key / 2;
            }

            // Store k-mer's starting position 'i' using a linked-list style in arrays
            // next_kmer[i] points to the previous k-mer that had the same hash
            next_kmer[i] = kmer_location[static_cast<int>(key)];
            // kmer_location[key] stores the latest k-mer's start position for this hash
            kmer_location[static_cast<int>(key)] = i;
        }
    }

    // Identifies and returns regions of lowercase characters in the sequence
    static std::vector<Position> lowercase_position(const std::string &sequence) {
        std::vector<Position> list; // List to store positions of lowercase segments
        bool successive = false;    // Flag to track if currently in a lowercase segment
        int start = 0, current_end = 0; // start: start of current segment, current_end: end of segment or current position

        for (int i = 0; i < sequence.length(); i++) {
            if (std::islower(sequence[i])) {
                if (successive) { // If already in a lowercase segment
                    current_end += 1; // Extend the segment
                } else { // New lowercase segment starts
                    start = i;
                    current_end += 1;
                    successive = true;
                }
            } else { // If character is not lowercase
                if (successive) { // If a lowercase segment just ended
                    Position position;
                    position.setstartinTar(start);
                    position.setendinTar(current_end - 1);
                    list.push_back(position);
                }
                successive = false;
                start = 0;
                current_end = i + 1;
            }
        }

        // After loop, if the sequence ends with a lowercase segment, store it
        if (successive) {
            Position position;
            position.setstartinTar(start);
            position.setendinTar(current_end);
            list.push_back(position);
        }
        return list;
    }
};

// Initialize static members of SequenceProcessor
std::unordered_map<int, std::vector<newkmer> > SequenceProcessor::hashmap;
std::vector<int> SequenceProcessor::next_kmer;
std::vector<int> SequenceProcessor::kmer_location;

class Matcher {
public:
    // Performs local matching between reference and target sequences
    static std::vector<Position>
    Lmatch(const std::string &ref_seq, const std::string &tar_seq, int kmer_len) {
        std::vector<Position> list;
        int index = 0, startIndex;
        int increment, most_incre, key_val;
        int kmerstart_val, endinRef_val, endinTar_val, Refendindex, Tarendindex;
        std::string kmer_str;

        SequenceProcessor::buildLhashtable(ref_seq, kmer_len); // Build k-mer hash table for the reference

        while (true) { // Iterate through the target sequence
            increment = 0;
            most_incre = 0; // Longest extension found for the current k-mer
            if (index + kmer_len > tar_seq.length()) { // Stop if remaining target is shorter than k-mer length
                break;
            }
            kmer_str = tar_seq.substr(index, kmer_len);
            key_val = SequenceProcessor::stringHashCode(kmer_str);

            // If k-mer hash not found in reference's hash table, move to next position in target
            if (SequenceProcessor::hashmap.find(key_val) == SequenceProcessor::hashmap.end()) {
                startIndex = std::numeric_limits<int>::max(); // Reset startIndex
                index = index + 1;
                if (index + kmer_len > tar_seq.length()) {
                    break;
                }
                continue;
            }

            // K-mer hash found, retrieve list of reference k-mers with this hash
            std::vector<newkmer> &klist = SequenceProcessor::hashmap[key_val];
            startIndex = std::numeric_limits<int>::max(); // Initialize with a large value
            most_incre = 0;
            for (const newkmer &newKMer: klist) { // Iterate through all reference k-mers with the same hash
                if (newKMer.getkmer() == kmer_str) { // Verify exact k-mer match (collision check)
                    kmerstart_val = newKMer.getkmerstart(); // Start of k-mer in reference
                    endinRef_val = kmerstart_val + kmer_len - 1; // End of k-mer in reference
                    endinTar_val = index + kmer_len - 1;      // End of k-mer in target
                    Refendindex = ref_seq.length() - 1;     // Last index of reference sequence
                    Tarendindex = tar_seq.length() - 1;     // Last index of target sequence

                    // Calculate how many additional characters match beyond the k-mer
                    increment = get_incre(ref_seq, tar_seq, endinRef_val, endinTar_val, Refendindex,
                                          Tarendindex);

                    if (klist.size() > 1) { // If multiple ref k-mers have the same hash and string
                        if (increment == most_incre) { // If this match has the same extension length as current best
                            // Prioritize the match closest to the end of the previous match in reference
                            if (list.size() > 1) {
                                int lastEIR = list.back().getendinRef();
                                if (std::abs(kmerstart_val - lastEIR) <
                                    std::abs(startIndex - lastEIR))
                                    startIndex = kmerstart_val;
                            } else if (list.empty()) { // If it's the first match
                                startIndex = kmerstart_val;
                            }
                        } else if (increment > most_incre) { // If this match has a longer extension
                            most_incre = increment;
                            startIndex = kmerstart_val;
                        }
                    } else { // Only one ref k-mer matched
                        most_incre = increment;
                        startIndex = kmerstart_val;
                        break; // No need to check other (non-existent) k-mers in klist
                    }
                }
            }

            // If no exact k-mer match was found (e.g. hash collision but different k-mer string) or no valid start
            if (startIndex == std::numeric_limits<int>::max()) {
                index = index + 1; // Move to next position in target
                if (index + kmer_len > tar_seq.length()) {
                    break;
                }
                continue;
            }

            // A match is found, store its position
            Position position_obj;
            position_obj.setstartinTar(index);
            position_obj.setendinTar(index + kmer_len + most_incre - 1);
            position_obj.setstartinRef(startIndex);
            position_obj.setendinRef(startIndex + kmer_len + most_incre - 1);
            list.push_back(position_obj);

            // Move index in target past the current match
            index = index + kmer_len + most_incre + 1;
            if (index + kmer_len > tar_seq.length()) { // Check bounds after advancing index
                break;
            }
        }
        return list;
    }

    // Performs global matching between reference and target sequences
    static std::vector<Position>
    Gmatch(const std::string &ref_seq, const std::string &tar_seq, int kmer_len, int limit_param) {
        std::vector<Position> list; // List to store found matches
        int index = 0; // Current position in the target sequence
        int startIndex; // Starting position of the best match in the reference sequence
        int lastEIR = 0; // End position in reference of the previously found match
        std::string kmer_str;

        SequenceProcessor::buildGhashtable(ref_seq, kmer_len); // Build k-mer hash table for the reference

        while (true) { // Iterate through the target sequence
            int increment = 0, most_incre = 0; // Variables for match extension length
            if (index + kmer_len > tar_seq.length()) { // Stop if remaining target is shorter than k-mer
                break;
            }
            kmer_str = tar_seq.substr(index, kmer_len); // Get current k-mer from target
            int key_val = std::abs(SequenceProcessor::stringHashCode(kmer_str)); // Hash the target k-mer

            if (key_val == -2147483648) { // Handle Integer.MIN_VALUE for abs
                key_val = 0;
            }
            // Ensure key is within bounds of kmer_location array
            while (key_val > SequenceProcessor::maxseq - 1) {
                key_val = key_val / 2;
            }

            // If k-mer hash not found in reference's hash table, move to next position in target
            if (SequenceProcessor::kmer_location[key_val] == -1) {
                startIndex = std::numeric_limits<int>::max();
                index = index + 1;
                if (index + kmer_len > tar_seq.length()) {
                    break;
                }
                continue;
            }

            startIndex = std::numeric_limits<int>::max(); // Initialize with a large value
            most_incre = 0;
            bool match_found = false;

            // First pass: iterate through reference k-mers with same hash, considering the distance limit
            for (int k = SequenceProcessor::kmer_location[key_val]; k != -1; k = SequenceProcessor::next_kmer[k]) {
                increment = 0;
                if (k + kmer_len > ref_seq.length()) continue; // Boundary check for ref_seq
                std::string Rkmer = ref_seq.substr(k, kmer_len); // Get k-mer from reference
                if (kmer_str != Rkmer) { // Verify exact k-mer match (collision check)
                    continue;
                }

                // Get end of last reference match to check distance limit
                try { // Safely access last element
                    if (!list.empty()) {
                        lastEIR = list.back().getendinRef();
                    } else {
                        lastEIR = 0; // No previous match
                    }
                } catch (...) { // Catch any potential exception during access
                    lastEIR = 0;
                }

                // If the current ref k-mer start 'k' is too far from the end of the last ref match
                if (std::abs(k - lastEIR) > limit_param) {
                    continue; // Skip this k-mer as it's outside the allowed proximity
                }

                match_found = true;
                int ref_idx = k + kmer_len;
                int tar_idx = index + kmer_len;
                while (ref_idx < ref_seq.length() && tar_idx < tar_seq.length()
                       && (ref_seq.substr(ref_idx, 1) == tar_seq.substr(tar_idx, 1))) {
                    ref_idx++;
                    tar_idx++;
                    increment++;
                }
                // Update best match if current extension is better or equally good but closer
                if (increment == most_incre) {
                    if (!list.empty()) { // Consider proximity if not the first match
                        if (std::abs(k - lastEIR) < std::abs(startIndex - lastEIR))
                            startIndex = k;
                    } else if (list.empty() && startIndex == std::numeric_limits<int>::max()) { // First match consideration
                        startIndex = k;
                    }
                } else if (increment > most_incre) {
                    most_incre = increment;
                    startIndex = k;
                }
            }
            // Second pass (fallback): if no match was found within the distance limit, search without the limit
            if (!match_found) {
                for (int k = SequenceProcessor::kmer_location[key_val]; k != -1; k = SequenceProcessor::next_kmer[k]) {
                    increment = 0;
                    if (k + kmer_len > ref_seq.length()) continue; // Boundary check
                    std::string Rkmer = ref_seq.substr(k, kmer_len);
                    if (kmer_str != Rkmer) { // Verify exact k-mer match
                        continue;
                    }
                    // No distance limit check in this fallback loop

                    // Extend the match
                    int ref_idx = k + kmer_len;
                    int tar_idx = index + kmer_len;

                    while (ref_idx < ref_seq.length() && tar_idx < tar_seq.length()
                           && (ref_seq.substr(ref_idx, 1) == tar_seq.substr(tar_idx, 1))) {
                        ref_idx++;
                        tar_idx++;
                        increment++;
                    }
                    // Update best match (same logic as above, but without the limit consideration active)
                    if (increment == most_incre) {
                        if (!list.empty()) {
                            if (std::abs(k - lastEIR) < std::abs(startIndex - lastEIR)) // Still prefer closer if extensions are equal
                                startIndex = k;
                        } else if (list.empty() && startIndex == std::numeric_limits<int>::max()) {
                            startIndex = k;
                        }
                    } else if (increment > most_incre) {
                        most_incre = increment;
                        startIndex = k;
                    }
                }
            }

            // If no valid start index was found after checking all candidates
            if (startIndex == std::numeric_limits<int>::max()) {
                index = index + 1; // Move to next position in target
                if (index + kmer_len > tar_seq.length()) {
                    break;
                }
                continue;
            }

            // A match is found, store its position
            Position position_obj;
            position_obj.setstartinTar(index);
            position_obj.setendinTar(index + kmer_len + most_incre - 1);
            position_obj.setstartinRef(startIndex);
            position_obj.setendinRef(startIndex + kmer_len + most_incre - 1);
            list.push_back(position_obj);

            // Move index in target past the current match
            index = index + kmer_len + most_incre + 1; // Advance index past the matched segment
            if (index + kmer_len > tar_seq.length()) { // Check bounds
                break;
            }
        }
        return list;
    }

private:
    // Calculates the number of additional characters that match beyond an initial k-mer match
    static int get_incre(const std::string& reference_seq, const std::string& target_seq,
                         int endinRef_param, int endinTar_param, int Refendindex, int Tarendindex) {
        int position = 0; // Number of additional matching characters
        int endIndex;

        // Determine the maximum possible extension based on the shorter remaining sequence length
        if (Refendindex - endinRef_param <= Tarendindex - endinTar_param) {
            endIndex = Refendindex - endinRef_param + 1; // Max extension limited by reference
        } else {
            endIndex = Tarendindex - endinTar_param + 1; // Max extension limited by target
        }

        // Iterate and count matching characters
        for (int i = 1; i < endIndex; i++) { // Start from 1 (next char after k-mer)
            // Boundary check (though max_possible_extension should prevent this, defensive check)
            if ((endinTar_param + i) >= target_seq.length() || (endinRef_param + i) >= reference_seq.length()) {
                // This error indicates an issue with logic or input parameters if reached
                std::cerr << "OUT OF BOUNDS in get_incre: endinTar_param + i = " << (endinTar_param + i)
                          << " target_seq.length() = " << target_seq.length()
                          << ", endinRef_param + i = " << (endinRef_param + i)
                          << " reference_seq.length() = " << reference_seq.length() << std::endl;
                break;
            }
            return position;
        }
    }
};

class FileUtils {
public:
    static std::string meta_data; // Stores metadata (e.g., header line from FASTA)
    static int length;           // Stores line length, typically for FASTA format consistency

    // Reads sequence from a file
    static std::string readSeq(const std::string &sequenceFileName, ScaleType scaleType = ScaleType::LOCAL, SequenceType sequenceType = SequenceType::TARGET) {
        std::ifstream file(sequenceFileName);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + sequenceFileName);
        }

        std::stringstream stringbuilder; // Efficiently builds the sequence string
        std::string line;
        int line_length;
        char temp_ch;
        bool base = false; // Flag to capture initial line length for local processing

        std::getline(file, line); // Read the first line

        // Metadata is stored if it's a target sequence in local processing mode
        if (sequenceType == SequenceType::TARGET && scaleType == ScaleType::LOCAL) {
            meta_data = line;
        }

        // Process remaining lines
        while (std::getline(file, line)) {
            line_length = line.length();
            if (sequenceType == SequenceType::REFERENCE && scaleType == ScaleType::GLOBAL) {
                // For global reference, convert to uppercase and exclude 'N'
                for (int i = 0; i < line_length; i++) {
                    temp_ch = std::toupper(line[i]);
                    if (temp_ch != 'N') {
                        stringbuilder << temp_ch;
                    }
                }
            } else {
                // For other cases, append line as is
                stringbuilder << line;
            }
            // For local processing, capture the length of the first sequence data line
            if (!base && scaleType == ScaleType::LOCAL && !line.empty()) { // Ensure line is not empty
                length = line_length;
                base = true;
            }
        }
        file.close();
        return stringbuilder.str();
    }

    // Reads target sequence for global processing, also extracts and writes N/L lists.
    static std::string GreadtarSeq(const std::string &sequenceFileName, const std::string &outFileName) {
        std::ifstream file(sequenceFileName);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + sequenceFileName);
        }

        // Open output file in append mode, as it will write metadata and N/L lists
        std::ofstream out(outFileName, std::ios::app);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot open output file: " + outFileName);
        }

        std::stringstream Llist;
        std::stringstream Nlist;

        std::string meta_data_local, line; // meta_data_local for this function's scope

        int line_length_val, totallength = 0, Llen = 0, Nlen = 0;
        int last_L_end = 0, last_N_end = 0; // Track end of last L/N region for delta calculation
        int start_L = 0, start_N = 0; // Start position of current L/N region
        char temp_ch;
        bool is_in_lowercase_region = false; // True if current char is part of a lowercase region
        bool is_in_N_region = false;      // True if current char is part of an 'N' region
        bool first_data_line = true;       // Flag to write line length once

        std::stringstream sequence_builder; // To build the actual sequence string (uppercase, N's removed)

        std::getline(file, meta_data_local); // Read metadata line
        out << meta_data_local << '\n';      // Write metadata to output

        while (std::getline(file, line)) { // Process each line of the sequence file
            line_length_val = line.length();

            if (first_data_line) { // Write original line length for the first data line
                out << std::to_string(line_length_val) << '\n';
                first_data_line = false;
            }

            for (int i = 0; i < line_length_val; i++) {
                temp_ch = line[i];

                if (std::islower(temp_ch)) { // Handle lowercase characters
                    if (is_in_lowercase_region) { // Continuing a lowercase region
                        Llen++;
                    } else { // Starting a new lowercase region
                        is_in_lowercase_region = true;
                        start_L = i + totallength;
                        Llist << (start_L - last_L_end) << ' '; // Delta from end of last L region
                        Llen++;
                    }
                    temp_ch = std::toupper(temp_ch); // Convert to uppercase for the main sequence
                } else { // Not lowercase
                    if (is_in_lowercase_region) { // End of a lowercase region
                        Llen--;
                        Llist << Llen << ' '; // Length of the just-ended L region
                        last_L_end = start_L + Llen;
                    }
                    is_in_lowercase_region = false;
                    Llen = 0;
                }

                if (temp_ch == 'N') { // Handle 'N' characters (already uppercased if it was 'n')
                    if (is_in_N_region) { // Continuing an 'N' region
                        Nlen++;
                    } else { // Starting a new 'N' region
                        is_in_N_region = true;
                        start_N = i + totallength;
                        Nlist << (start_N - last_N_end) << ' '; // Delta from end of last N region
                        Nlen++;
                    }
                    // 'N's are not added to sequence_builder
                } else { // Not 'N'
                    if (is_in_N_region) { // End of an 'N' region
                        Nlist << Nlen << ' '; // Length of the just-ended N region
                        last_N_end = start_N + Nlen;
                    }
                    is_in_N_region = false;
                    Nlen = 0;
                    sequence_builder << temp_ch; // Add valid character to sequence
                }
            }
            totallength += line_length_val; // Accumulate total length processed
        }

        // Finalize Llist string
        std::string LlistStr = Llist.str();

        if (LlistStr.length() > 0) {
            if (is_in_lowercase_region) { // If sequence ends in a lowercase region
                LlistStr += std::to_string(Llen);
            } else { // Remove trailing space if last was a length
                LlistStr.pop_back();
            }
        }

        // Finalize Nlist string
        std::string NlistStr = Nlist.str();
        if (NlistStr.length() > 0) {
            if (is_in_N_region) { // If sequence ends in an 'N' region
                NlistStr += std::to_string(Nlen);
            } else { // Remove trailing space
                NlistStr.pop_back();
            }
        }

        file.close();

        out << LlistStr << '\n'; // Write Llist data
        out << NlistStr << '\n'; // Write Nlist data
        out.flush();
        out.close();
        return sequence_builder.str(); // Return the processed sequence (uppercase, no 'N's)
    }

    // Writes text content to a file
    static void
    write(const std::string &filename, const std::string &text_content, bool append) {
        std::ofstream output;
        if (append) { // Open in append mode if true
            output.open(filename, std::ios::app);
        } else { // Otherwise, overwrite
            output.open(filename);
        }
        if (!output.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        output << text_content; // Write content
        output.flush();
        output.close();
    }

    // Writes a list of Position objects to a file in a specific format
    static void
    write(std::string filename, std::vector<Position> list_param, bool append, std::string auxiliary) {
        std::ofstream output;
        if (append) { // Open in append mode
            output.open(filename, std::ios::app);
        } else { // Overwrite mode
            output.open(filename);
        }

        if (!output.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        std::ostringstream text_stream; // Use ostringstream for efficient string building
        int temp_end = 0; // Tracks the end position of the previous segment in target
        output << auxiliary; // Write auxiliary information (like metadata) first

        // Format: delta_start length
        for (Position position_obj: list_param) {
            int start_val = position_obj.getstartinTar();
            int end_val = position_obj.getendinTar();
            text_stream << (start_val - temp_end) << " " << (end_val - start_val) << " ";
            temp_end = end_val;
        }
        output << text_stream.str(); // Write the formatted position data
        output << "\n"; // Add a newline at the end
        output.flush();
        output.close();
    }

    // Compresses a file using 7zip with PPMD method
    static void use7zip(std::string filename) {
        struct stat buffer;
        // Check if the file to be compressed exists
        if (stat(filename.c_str(), &buffer) != 0) {
            std::cerr << "Warning: File " << filename << " not found for 7zip." << std::endl;
            return;
        }

        // Construct the 7zip command
        std::string exec = "./7za a " + filename + ".7z " + filename + " -m0=PPMD";
        try {
            int result = std::system(exec.c_str()); // Execute the command
            if (result != 0) { // Check command execution status
                std::cerr << "Warning: 7zip command failed with result " << result << " for file " << filename << std::endl;
            }
        } catch (const std::exception &e) { // Catch any exceptions during system call
            std::cerr << "Exception in use7zip for " << filename << ": " << e.what() << std::endl;
        }
    }
};


class ORGC { // Main class for Off-Reference Genome Compressor
public:
    // Static members to hold current reference and target sequences/segments being processed
    static std::string reference, target;
    static std::string text; // Buffer for accumulating output text before writing
    // Parameters for k-mer matching and sequence processing
    static int kmer_length;  // Length of k-mers used for matching
    static int sub_length;   // Length of sub-sequences for local processing
    static int limit;        // Distance limit for global matching
    static double T1;        // Threshold for mismatch ratio in local processing
    static int T2;           // Threshold for number of mismatched segments in local processing
    // Start and end coordinates for current segments in target (sot, eot) and reference (sor, eor)
    static int sot, eot, sor, eor;
    static int mismatch;     // Counter for mismatched segments in local processing
    static int endref;       // End position of the last match in the reference segment
    static bool local;        // Flag indicating if local compression attempt was successful

    // Array of chromosome names (standard human genome)
    static const std::string chrName[];


    // Constructor
    ORGC() {
        text = "";
        SequenceProcessor::hashmap.clear(); // Clear hashmap for new ORGC object/segment processing
    }

    // Gets current CPU time in nanoseconds
    static long getCPUTime() {
        auto now = std::chrono::high_resolution_clock::now();
        auto duration = now.time_since_epoch();
        auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);
        return nanoseconds.count();
    }

    // Formats a list of matches (Position objects) into a string representation for local processing.
    // Also identifies and includes unaligned segments (mismatches) from the target.
    static Position format_matches(const std::vector<Position> &list_param) {
        int startinTar_val, startinRef_val, endinRef_val;
        int trouble = 0; // Accumulates length of mismatched segments within current block
        if (list_param.empty()) {
            return Position(); // Return an empty Position if no matches
        }

        for (int i = 0; i < list_param.size(); i++) {
            // For the first match in the list
            if (i == 0) {
                startinTar_val = list_param[i].getstartinTar();
                startinRef_val = list_param[i].getstartinRef();
                endinRef_val = list_param[i].getendinRef();
                if (endinRef_val >= endref) { // Update overall end of reference match
                    endref = endinRef_val;
                }

                // If the match doesn't start at the beginning of the target segment, write the preceding part as mismatch
                if (startinTar_val > 0) {
                    std::string preamble = target.substr(0, startinTar_val);
                    text += preamble + "\n"; // Append mismatch to output text
                    trouble += preamble.length();
                }
                // Append the match coordinates (adjusted by sor for global reference position)
                text += "" + std::to_string(startinRef_val + sor) + "," + std::to_string(endinRef_val + sor) +
                        "\n";
                continue; // Move to next match in list
            }

            // For subsequent matches
            startinTar_val = list_param[i].getstartinTar();
            startinRef_val = list_param[i].getstartinRef();
            endinRef_val = list_param[i].getendinRef();
            if (endinRef_val >= endref) { // Update overall end of reference match
                endref = endinRef_val;
            }

            int endinTarPrev = list_param[i - 1].getendinTar(); // End of previous match in target
            // If there's a gap between previous match and current match in target, write it as mismatch
            if (endinTarPrev + 1 < startinTar_val) {
                std::string mismatch_str = target.substr(endinTarPrev + 1, startinTar_val - (endinTarPrev +
                                                                                             1));
                if (mismatch_str.length() > 0) {
                    text += mismatch_str + "\n";
                    trouble += mismatch_str.length();
                }
            }
            // Append current match coordinates
            text += std::to_string(startinRef_val + sor) + "," + std::to_string(endinRef_val + sor) +
                    "\n";
        }
        // Check if accumulated mismatches exceed threshold T1 for the current sub-segment
        if (trouble > (sub_length * T1))
        {
            mismatch++; // Increment overall mismatch counter for the chromosome
        }

        return list_param.back(); // Return the last match processed
    }

    // Formats a list of matches and writes them to a file (used in global processing)
    static void format_matches(const std::vector<Position> &list_param, const std::string &fileName) {
        std::stringstream stringbuilder; // For building the output string

        int startinTar_val, startinRef_val, endinRef_val, endinTar_val = 0; // Initialize endinTar_val
        for (int i = 0; i < list_param.size(); i++) {
            if (i == 0) { // First match
                startinTar_val = list_param[i].getstartinTar();
                endinTar_val = list_param[i].getendinTar();
                startinRef_val = list_param[i].getstartinRef();
                endinRef_val = list_param[i].getendinRef();
                // If match doesn't start at beginning of target, write preamble
                if (startinTar_val > 0) {
                    std::string preamble = target.substr(0, startinTar_val);
                    stringbuilder << preamble << "\n";
                }
                // Write match coordinates (raw, not adjusted by sor as in local)
                stringbuilder << startinRef_val << "," << endinRef_val << "\n";
                continue;
            }

            // Subsequent matches
            startinTar_val = list_param[i].getstartinTar();
            startinRef_val = list_param[i].getstartinRef();
            endinRef_val = list_param[i].getendinRef();
            endinTar_val = list_param[i].getendinTar(); // Update endinTar_val for current match
            int endinTarPrev = list_param[i - 1].getendinTar(); // End of previous target match

            // Handle gap between matches (mismatch)
            if (endinTarPrev + 1 < startinTar_val) {
                std::string mismatch_str = target.substr(endinTarPrev + 1, startinTar_val - (endinTarPrev +
                                                                                             1));
                if (mismatch_str.length() > 0) {
                    stringbuilder << mismatch_str << "\n";
                }
            }

            stringbuilder << startinRef_val << "," << endinRef_val << "\n"; // Write current match
        }
        // After all matches, if there's a remaining part of the target sequence, write it
        if (endinTar_val < (target.length() - 1)) {
            // Boundary check for substr
            if ( (endinTar_val + 1) < target.length() ) { // Ensure there's something to write
                stringbuilder << target.substr(endinTar_val + 1, (target.length()) - (endinTar_val + 1));
            }
        }
        FileUtils::write(fileName, stringbuilder.str(), true); // Append formatted output to file
    }

    // Post-processes the formatted match file to merge consecutive reference segments and reformat coordinates
    static void postprocess(std::string filename, std::string final_file) {
        std::ifstream input(filename); // Input is the interim file with matches and mismatches
        if (!input.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        std::string line;
        std::ostringstream stringbuilder_merged; // String builder for merging phase
        std::vector<int> num_list; // Stores start and end of a potential sequence of reference coordinates

        // Phase 1: Merge consecutive reference segments
        while (std::getline(input, line)) {
            if (line.find(",") != std::string::npos) { // Line contains reference coordinates (e.g., "start,end")
                size_t comma_pos = line.find(",");
                int begin = std::stoi(line.substr(0, comma_pos));
                int end_val = std::stoi(line.substr(comma_pos + 1));

                if (num_list.size() > 0) { // If there's an existing segment being built
                    int prevEnd = num_list.back();
                    if (begin != prevEnd + 1) { // If current segment is not consecutive to previous
                        // Write out the previously accumulated segment
                        stringbuilder_merged << num_list[0] << "," << num_list.back() << "\n";
                        num_list.clear(); // Start a new segment
                    }
                }
                num_list.push_back(begin);
                num_list.push_back(end_val);
            } else if (!line.empty()) {

                if (!num_list.empty()) {
                    stringbuilder_merged << num_list[0] << "," << num_list.back() << "\n";
                }
                // Handle lines that might have special meaning (like "^") though not fully clear from context
                if (line.find("^") == std::string::npos) {
                    stringbuilder_merged << line << "\n";
                }
                num_list.clear();
            }
        }
        // After loop, write out any remaining reference segment
        if (!num_list.empty()) {
            stringbuilder_merged << num_list[0] << "," << num_list.back() << "\n";
        }
        input.close();

        // Phase 2: Reformat merged reference coordinates to delta encoding (start, length_delta)
        std::istringstream inputStream(stringbuilder_merged.str()); // Use the output of phase 1 as input
        stringbuilder_merged.str("");
        stringbuilder_merged.clear();
        std::ostringstream stringbuilder_final; // String builder for final output format
        num_list.clear(); // Re-use num_list for new format
        std::vector<std::string> stringList; // Store lines for processing

        while (std::getline(inputStream, line)) {
            stringList.push_back(line);
        }

        int prev_ref_end = 0; // Tracks the end of the previous reference segment for delta calculation
        bool successive_ref_segment = false; // Tracks if currently processing a sequence of reference segments
        for (size_t i = 0; i < stringList.size(); i++) {
            std::string str = stringList[i];

            if (str.find(",") != std::string::npos) { // If it's a reference coordinate line
                size_t comma_pos = str.find(",");
                int begin = std::stoi(str.substr(0, comma_pos));
                int end_val = std::stoi(str.substr(comma_pos + 1));
                if (!successive_ref_segment) {
                    num_list.push_back(begin);
                    num_list.push_back(end_val - begin);
                    prev_ref_end = end_val;

                    if (!num_list.empty()) stringbuilder_merged << num_list[0] << "," << num_list.back() << "\n";

                    successive_ref_segment = true;
                } else {
                    num_list.push_back(begin - prev_ref_end);
                    num_list.push_back(end_val - begin);
                    prev_ref_end = end_val;
                    if (!num_list.empty()) stringbuilder_merged << num_list[0] << "," << num_list.back() << "\n";
                }
                num_list.clear();
            } else if (str.length() > 0) {
                stringbuilder_merged << str << "\n";
            }
        }

        if (!num_list.empty()) {
            stringbuilder_merged << num_list[0] << "," << num_list.back() << "\n";
        }

        FileUtils::write(final_file, stringbuilder_merged.str(), true);
    }
private:
    // Processes a target genome locally against a reference genome segment by segment
    static void process_local(const std::string& ref_genome_path, const std::string& tar_genome_path, const std::string& final_file_path, const std::string& tempfile) {
        int controuble = 0; // Counter for consecutive unmatchable segments
        bool is_con = false; // Flag indicating if the previous segment was unmatchable

        std::cout << tar_genome_path << " is compressing (local attempt)...\n" << std::endl;

        struct stat buffer;
        // Remove existing final file to start fresh for this chromosome
        if (stat(final_file_path.c_str(), &buffer) == 0) {
            if (std::remove(final_file_path.c_str()) != 0) {
                std::cerr << "Warning: Failed to remove existing final file: " << final_file_path << std::endl;
            }
        }

        // Read reference and target sequences
        std::string reference_seq_content = FileUtils::readSeq(ref_genome_path, LOCAL, REFERENCE);
        std::string target_seq_content = FileUtils::readSeq(tar_genome_path, LOCAL, TARGET); // This also sets FileUtils::meta_data and FileUtils::length

        // Error handling for empty sequences
        if (reference_seq_content.empty()) {
            std::cerr << "Error: Reference sequence is empty for " << ref_genome_path << std::endl;
            ORGC::local = false; // Mark local processing as failed
            return;
        }
        if (target_seq_content.empty()) {
            std::cerr << "Error: Target sequence is empty for " << tar_genome_path << std::endl;
            ORGC::local = false; // Mark local processing as failed
            return;
        }

        int target_length_val = target_seq_content.length();

        // Adjust thresholds T1, T2 based on target length; if too short, local is unsuitable
        if (target_length_val < ORGC::sub_length * 5) { // If target is too short
            ORGC::local = false; // Switch to global
            return;
        } else if (target_length_val < ORGC::sub_length * 1333) { // For moderately sized targets
            ORGC::T1 = 0.1; // Lower mismatch ratio threshold
            ORGC::T2 = 0;   // Lower segment mismatch threshold (effectively 1 allowed)
        }

        // Prepare auxiliary information (metadata, line length) for the output file
        std::string auxiliary = FileUtils::meta_data + "\n" + std::to_string(FileUtils::length) + "\n";
        // Find and write positions of lowercase segments from target
        std::vector<Position> L_list = SequenceProcessor::lowercase_position(target_seq_content);
        FileUtils::write(final_file_path, L_list, false, auxiliary); // false: overwrite, write auxiliary info
        FileUtils::write(final_file_path, "\n", true); // Add a newline separator

        // Convert sequences to uppercase for matching
        std::transform(reference_seq_content.begin(), reference_seq_content.end(), reference_seq_content.begin(),
                       ::toupper);
        std::transform(target_seq_content.begin(), target_seq_content.end(), target_seq_content.begin(), ::toupper);

        // Remove existing temp file
        if (stat(tempfile.c_str(), &buffer) == 0) {
            if (std::remove(tempfile.c_str()) != 0) {
                std::cerr << "Warning: Failed to remove existing temp file: " << tempfile << std::endl;
            }
        }

        // Initialize segment pointers
        ORGC::sot = 0; // Start of target segment
        ORGC::eot = ORGC::sub_length; // End of target segment
        ORGC::sor = 0; // Start of reference segment
        ORGC::eor = ORGC::sub_length; // End of reference segment
        Position current_position; // Stores position of the last processed match

        ORGC::text = ""; // Clear global text buffer before starting segment processing

        // Main loop for processing segments
        while (true) {
            ORGC utilities; // Creates an ORGC object, which clears SequenceProcessor::hashmap
            int kmerlength_val = ORGC::kmer_length; // Use default k-mer length

            // Check boundary conditions for reference sequence
            if (ORGC::eor > reference_seq_content.length()) {
                if (ORGC::sot >= target_seq_content.size()) {
                    break;
                }
                ORGC::text = target_seq_content.substr(ORGC::sot);
                if (ORGC::text.length() <= 0) {
                    break;
                } else {
                    FileUtils::write(tempfile, ORGC::text, true);
                    break;
                }
            }
            // Check boundary conditions for target sequence
            if (ORGC::eot > target_seq_content.length()) {
                if (ORGC::sot >= target_seq_content.size()) {
                    break;
                }
                ORGC::text = target_seq_content.substr(ORGC::sot);
                if (ORGC::text.length() <= 0) {
                    break;
                } else {
                    FileUtils::write(tempfile, ORGC::text, true);
                    break;
                }
            }

            // Extract current reference and target segments
            std::string current_ref_segment = reference_seq_content.substr(ORGC::sor, ORGC::eor - ORGC::sor);
            std::string current_tar_segment = target_seq_content.substr(ORGC::sot, ORGC::eot - ORGC::sot);

            ORGC::reference = current_ref_segment;
            ORGC::target = current_tar_segment;

            std::vector<Position> list_matches = Matcher::Lmatch(current_ref_segment, current_tar_segment,
                                                                 kmerlength_val);
            // If no matches, try with a shorter k-mer length (11) as a fallback
            if (list_matches.empty()) {
                kmerlength_val = 11;
                list_matches = Matcher::Lmatch(current_ref_segment, current_tar_segment, kmerlength_val);
            }

            // If still no matches found even with shorter k-mer
            if (list_matches.empty()) {
                ORGC::mismatch++; // Increment chromosome's mismatch counter

                if (ORGC::eot >= target_seq_content.length() - 1) {
                    if (ORGC::sot >= target_seq_content.size()) {
                        break;
                    }
                    ORGC::text = target_seq_content.substr(ORGC::sot);
                    FileUtils::write(tempfile, ORGC::text, true);
                    break;
                }

                if (is_con) { // If previous segment was also unmatchable
                    controuble++;
                }
                is_con = true;

                ORGC::text += current_tar_segment + "\n"; // Append to existing ORGC::text
                FileUtils::write(tempfile, ORGC::text, true);
                ORGC::text = ""; // Reset text for next utility call or end of segment
                ORGC::sot += ORGC::sub_length;
                ORGC::eot = ORGC::sot + ORGC::sub_length;
                ORGC::eor += ORGC::sub_length;


                int difference = target_seq_content.length() - ORGC::sot;

                if (difference <= ORGC::kmer_length) {
                    if (ORGC::sot >= target_seq_content.size()) {
                        break;
                    }
                    ORGC::text = target_seq_content.substr(ORGC::sot);
                    FileUtils::write(tempfile, ORGC::text, true);
                    ORGC::text = "";
                    break;
                } else if (difference < ORGC::sub_length) {
                    ORGC::eot = target_seq_content.length();
                }

                int difference_ref = reference_seq_content.length() - ORGC::sor;

                if (difference_ref < ORGC::sub_length) {
                    ORGC::eor = reference_seq_content.length();
                }
                if (ORGC::eot >= target_seq_content.length()) { // This condition might be redundant due to earlier checks
                     // break; // Original code had a break here if eot > target_seq_content.length()
                }
                if (difference_ref <= ORGC::kmer_length) {
                    if (ORGC::sot >= target_seq_content.size()) {
                        break;
                    }
                    ORGC::text = target_seq_content.substr(ORGC::sot);
                    if (ORGC::text.length() <= 0) {
                        break;
                    } else {
                        if (ORGC::text.length() > (ORGC::sub_length * ORGC::T1)) {
                            ORGC::mismatch++;
                        }
                        if (ORGC::mismatch > ORGC::T2) {
                            ORGC::local = false; // Signal to switch to global
                            // FileUtils::write(tempfile, ORGC::text, true); // Write remaining text before failing local
                            // ORGC::text = "";
                            return; // Exit local processing
                        }
                        FileUtils::write(tempfile, ORGC::text, true);
                        ORGC::text = "";
                        break;
                    }
                }
                continue; // Proceed to next segment
            }

            // Matches found
            is_con = false; // Reset consecutive unmatchable flag
            if (controuble > 2) // If there were more than 2 consecutive unmatchable segments before this match
                ORGC::mismatch -= controuble; // Reduce mismatch count (heuristic)
            controuble = 1; // Reset counter

            // Format the found matches (this appends to static ORGC::text)
            current_position = ORGC::format_matches(list_matches);

            // If mismatch threshold T2 is exceeded after formatting current matches
            if (ORGC::mismatch > ORGC::T2) {
                ORGC::local = false; // Signal to switch to global
                // FileUtils::write(tempfile, ORGC::text, true); // Write current ORGC::text before failing
                // ORGC::text = "";
                return; // Exit local processing
            }

            ORGC::sot += current_position.getendinTar() + 1;
            ORGC::eot = ORGC::sot + ORGC::sub_length;
            ORGC::sor += ORGC::endref + 1; // ORGC::endref is updated by format_matches
            ORGC::endref = 0; 
            ORGC::eor = ORGC::sor + ORGC::sub_length;

            FileUtils::write(tempfile, ORGC::text, true); // Write the formatted matches and mismatches
            ORGC::text = ""; // Reset for the next segment or ORGC utilities call

            int difference = target_seq_content.length() - ORGC::sot;

            if (difference <= ORGC::kmer_length) {
                if (ORGC::sot < target_seq_content.size()) { // Check if there's anything left
                     ORGC::text = target_seq_content.substr(ORGC::sot);
                     if (ORGC::text.length() > 0) {
                        FileUtils::write(tempfile, ORGC::text, true);
                        ORGC::text = "";
                     }
                }
                break;
            } else if (difference < ORGC::sub_length) {
                ORGC::eot = target_seq_content.length();
            }

            int difference_ref = reference_seq_content.length() - ORGC::sor;

            if (difference_ref < ORGC::sub_length) {
                ORGC::eor = reference_seq_content.length();
            }
            if (difference_ref <= ORGC::kmer_length) {
                 if (ORGC::sot < target_seq_content.size()) {
                    ORGC::text = target_seq_content.substr(ORGC::sot);
                    if (ORGC::text.length() > 0) {
                        if (ORGC::text.length() > (ORGC::sub_length * ORGC::T1)) {
                            ORGC::mismatch++;
                        }
                        if (ORGC::mismatch > ORGC::T2) {
                            ORGC::local = false;
                            // FileUtils::write(tempfile, ORGC::text, true);
                            // ORGC::text = "";
                            return;
                        }
                        FileUtils::write(tempfile, ORGC::text, true);
                        ORGC::text = "";
                    }
                }
                break;
            }
        } // End of while(true) segment processing loop

        if (!ORGC::local) { // If loop decided local failed
            return;
        }

        // If local processing was successful up to this point
        ORGC::postprocess(tempfile, final_file_path); // Post-process the temp file into final output
        // Clean up temp file
        if (stat(tempfile.c_str(), &buffer) == 0) {
            if (std::remove(tempfile.c_str()) != 0) {
                std::cerr << "Warning: Failed to remove temp file after local processing: " << tempfile << std::endl;
            }
        }
    }

    // Processes a target genome globally against the entire reference genome
    static void process_global(const std::string& ref_genome_path, const std::string& tar_genome_path, const std::string& final_file_path, const std::string& tempfile) {
        std::cout << tar_genome_path << " is compressing (global attempt)...\n" << std::endl;
        struct stat buffer;
        // Remove existing final file, as GreadtarSeq will write headers to it
        if (stat(final_file_path.c_str(), &buffer) == 0) {
            if (std::remove(final_file_path.c_str()) != 0) {
                std::cerr << "Warning: Failed to remove existing final file before global: " << final_file_path << std::endl;
            }
        }
        // Ensure tempfile is clean (might contain data from a failed local attempt)
        if (stat(tempfile.c_str(), &buffer) == 0) {
            if (std::remove(tempfile.c_str()) != 0) {
                std::cerr << "Warning: Failed to remove existing temp file before global: " << tempfile << std::endl;
            }
        }

        // Read reference sequence for global processing (uppercase, N's removed)
        std::string reference_seq_content_global = FileUtils::readSeq(ref_genome_path, GLOBAL, REFERENCE);
        // Read target sequence: GreadtarSeq writes metadata, line length, L-list, N-list to final_file_path
        // and returns the processed sequence string (uppercase, N's removed).
        std::string target_seq_content_global = FileUtils::GreadtarSeq(tar_genome_path, final_file_path);


        if (reference_seq_content_global.empty()) {
            std::cerr << "Error: Reference sequence is empty for " << ref_genome_path << " (global)" << std::endl;
            // Write empty or minimal content to final_file_path if needed, or let it be as GreadtarSeq left it
            return;
        }
        if (target_seq_content_global.empty() && FileUtils::meta_data.empty()) { // Check if GreadtarSeq itself indicated an issue
            std::cerr << "Error: Target sequence is effectively empty for " << tar_genome_path << " (global)" << std::endl;
            // GreadtarSeq might have written headers, so final_file_path might not be empty.
            // If target_seq_content_global is empty, no matching can be done.
            return;
        }


        std::string full_reference_seq = reference_seq_content_global;
        std::string full_target_seq = target_seq_content_global;

        ORGC::reference = full_reference_seq; // For format_matches if it uses static ORGC::target
        ORGC::target = full_target_seq;       // For format_matches

        std::vector<Position> list_global_matches = Matcher::Gmatch(full_reference_seq, full_target_seq,
                                                                 ORGC::kmer_length, ORGC::limit);

        if (!full_target_seq.empty() && !list_global_matches.empty()) {
            // format_matches for global writes to tempfile
            ORGC::format_matches(list_global_matches, tempfile);
            // Post-process the tempfile and append results to final_file_path (after N/L lists)
            ORGC::postprocess(tempfile, final_file_path);
        } else if (full_target_seq.empty() && !list_global_matches.empty()){
             std::cerr << "Warning: Global matches found but target sequence content for matching was empty for " << tar_genome_path << std::endl;
        } else if (!full_target_seq.empty() && list_global_matches.empty()){
             std::cout << "Info: No global matches found for " << tar_genome_path << ". Target content might be written as is if GreadtarSeq handled it." << std::endl;
             // If GreadtarSeq wrote the sequence and there are no matches, postprocess might not be needed or might write an empty match list.
             // The current format_matches/postprocess path for global expects matches.
             // If no matches, GreadtarSeq already wrote N/L lists and metadata. The sequence itself is not written by GreadtarSeq.
             // If target is not empty but no matches, the final file will contain N/L lists and metadata, but no sequence data or match data from postprocess.
             // This seems to be the implicit behavior.
        }

        // Clean up temp file
        if (stat(tempfile.c_str(), &buffer) == 0) {
            if (std::remove(tempfile.c_str()) != 0) {
                std::cerr << "Warning: Failed to remove temp file after global processing: " << tempfile << std::endl;
            }
        }
    }

public:
    // Main function to process genomes: iterates through chromosomes and applies local or global compression.
    static void process_genome(const std::string& ref_base_path, const std::string& tar_base_path, const std::string& output_dir_path) {
        std::string final_folder = output_dir_path + "/result"; // Output directory for compressed files
        std::string misjuggements_path = final_folder + "/misjuggements.txt"; // File to log mismatch counts

        // Clear or create the misjuggements file
        std::ofstream misjuggements_file_clear(misjuggements_path, std::ios::out); // Open in overwrite mode
        if (!misjuggements_file_clear.is_open()) {
            std::cerr << "Could not open misjuggements.txt for clearing!" << std::endl;
        }
        misjuggements_file_clear.close();

        for (int i = 0; i < 24; i++) {
            ORGC::kmer_length = 21; // Reset default kmer_length for each chromosome
            ORGC::local = true;     // Assume local processing first
            ORGC::endref = ORGC::sub_length - 1; // Reset endref
            ORGC::mismatch = 0;     // Reset mismatch count for the chromosome
            // Reset T1, T2 to defaults, process_local will adjust if needed
            ORGC::T1 = 0.5;
            ORGC::T2 = 4;


            std::string current_chr_name = ORGC::chrName[i];
            // Construct full paths for reference and target chromosome files
            std::string reference_genome_path = ref_base_path + "/" + current_chr_name;
            std::string target_genome_path = tar_base_path + "/" + current_chr_name;
            // Path for the final compressed output for this chromosome
            std::string final_file_path = final_folder + "/" + current_chr_name;
            // Path for the temporary intermediate file
            std::string tempfile = final_folder + "/interim.txt";

            // Attempt local processing first
            process_local(reference_genome_path, target_genome_path, final_file_path, tempfile);

            // If local processing failed or was deemed unsuitable (ORGC::local is false)
            if (!ORGC::local) {
                // Fallback to global processing for this chromosome
                process_global(reference_genome_path, target_genome_path, final_file_path, tempfile);
            }

            // Log the number of mismatches for this chromosome
            std::ofstream misjuggements_file(misjuggements_path, std::ios::app); // Open in append mode
            if (!misjuggements_file.is_open()) {
                std::cerr << "Could not open misjuggements.txt for writing for " << current_chr_name << std::endl;
            } else {
                misjuggements_file << current_chr_name << ": " << ORGC::mismatch << std::endl;
                misjuggements_file.close();
            }
        }
    }
    
};

// Definition and initialization of static const member chrName
const std::string ORGC::chrName[] = {
        "chr1.fa", "chr2.fa", "chr3.fa", "chr4.fa", "chr5.fa", "chr6.fa", "chr7.fa", "chr8.fa", "chr9.fa",
        "chr10.fa", "chr11.fa", "chr12.fa", "chr13.fa", "chr14.fa", "chr15.fa", "chr16.fa",
        "chr17.fa", "chr18.fa", "chr19.fa", "chr20.fa",
        "chr21.fa", "chr22.fa", "chrX.fa", "chrY.fa"
};

// Initialization of static members of ORGC class
std::string ORGC::reference = "";
std::string ORGC::target = "";
std::string ORGC::text = "";
int ORGC::kmer_length = 21;      // Default k-mer length
int ORGC::sub_length = 30000;    // Default sub-sequence length for local processing
int ORGC::limit = 100;           // Default distance limit for global matching
double ORGC::T1 = 0.5;           // Default threshold for mismatch ratio in local
int ORGC::T2 = 4;                // Default threshold for mismatched segments in local
int ORGC::sot = 0;
int ORGC::eot = 0;
int ORGC::sor = 0;
int ORGC::eor = 0;
int ORGC::mismatch = 0;
int ORGC::endref = 0; // Initialized to 0. Was sub_length -1, but 0 is safer for start. Updated in process_local.
bool ORGC::local = true;         // Default to try local processing first

int FileUtils::length = 0;
std::string FileUtils::meta_data = "";

int main(int argc, char *argv[]) {

    // Check for correct number of command-line arguments
    if (argc != 4) {
        std::cout << "Usage: <program> <reference_base_path> <target_base_path> <output_dir>" << std::endl;
        std::cout << "Make sure you have inputted 3 arguments." << std::endl;
        return 1; // Return error code
    }

    // Parse command-line arguments
    std::string reference_base_path = argv[1]; // Path to directory containing reference genome files
    std::string target_base_path = argv[2];    // Path to directory containing target genome files
    std::string output_dir = argv[3];          // Path to output directory

    struct stat st = {0};
    // Create output directory if it doesn't exist
    if (stat(output_dir.c_str(), &st) == -1) { // If directory does not exist
        // Attempt to create it with read/write/execute permissions for owner, group, others
        if (mkdir(output_dir.c_str(), 0777) != 0) {
            std::cerr << "Failed to create output directory: " << output_dir << std::endl;
            return 1; // Return error code
        }
    }
    // Create "result" subdirectory within the output directory
    std::string final_folder = output_dir + "/result";
    if (stat(final_folder.c_str(), &st) == -1) { // If "result" subdirectory does not exist
        if (mkdir(final_folder.c_str(), 0777) != 0) {
            std::cerr << "Failed to create result subdirectory: " << final_folder << std::endl;
            return 1; // Return error code
        }
    }

    // Record start time and CPU time
    auto startDate = std::chrono::system_clock::now();
    long startCpuTimeNano = ORGC::getCPUTime();
    std::time_t start_time = std::chrono::system_clock::to_time_t(startDate);
    std::cout << "Start time: " << std::ctime(&start_time);

    // Pre-allocate and initialize static vectors in SequenceProcessor for global hashing
    // Maxchar seems related to total unique k-mers or hash values, maxseq to sequence positions
    SequenceProcessor::next_kmer.resize(SequenceProcessor::maxchar, -1);
    SequenceProcessor::kmer_location.resize(SequenceProcessor::maxseq, -1);

    // Start the genome processing
    ORGC::process_genome(reference_base_path, target_base_path, output_dir);

    // Compress the entire result folder using 7zip
    FileUtils::use7zip(final_folder); // final_folder is output_dir + "/result"

    std::cout << "All Done" << std::endl;
    // Calculate and print total CPU time taken
    long taskCPUTimeNano = ORGC::getCPUTime() - startCpuTimeNano;
    std::cout << "Compressed time: " << (double) taskCPUTimeNano / 1000000000.0 << " seconds." << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;

    return 0; // Successful execution
}