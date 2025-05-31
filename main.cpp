#include <iostream>
#include <string>
#include <vector>
#include <random>    // For random nucleotide generation
#include <fstream>   // For checking if files exist
#include <stdexcept> // For std::exception
#include <algorithm>
#include <iomanip> // For std::fixed and std::setprecision if printing percentages

#include "src/FastaParser.h"
#include "src/KmerIndex.h"
#include "src/MatchFinder.h"
#include "src/Encoder.h"
#include "src/Decoder.h"
#include "src/Compressor.h" // Include the new Compressor header

// Function to print usage instructions for the command-line tool.
void printUsage(const std::string& program_name) {
    std::cerr << "Usage:" << std::endl;
    std::cerr << "  " << program_name << " encode <reference.fa> <target.fa> <intermediate_output.txt> <k_mer> <k0> <L_seg> <search_range_global> <T1_mismatch_perc> <T2_mismatch_time>" << std::endl;
    std::cerr << "  " << program_name << " decode <reference.fa> <encoded_input.sccg> <decoded_output.fa>" << std::endl;
    std::cerr << "\nArguments:" << std::endl;
    std::cerr << "  mode: 'encode' or 'decode'." << std::endl;
    std::cerr << "  <reference.fa>: Path to the reference genome FASTA file." << std::endl;
    std::cerr << "  <target.fa>: Path to the target genome FASTA file (for encoding)." << std::endl;
    std::cerr << "  <encoded_output.sccg / intermediate_output.txt>: Path for output file." << std::endl;
    std::cerr << "  <encoded_input.sccg>: Path to the encoded file (for decoding)." << std::endl;
    std::cerr << "  <decoded_output.fa>: Path to save the reconstructed FASTA file (for decoding)." << std::endl;
    std::cerr << "  [kmer_size]: Optional. Integer size for k-mers (default: 12, for basic encoding only)." << std::endl;
    std::cerr << "  --- SCCG Algorithm Parameters (for full encode mode) ---" << std::endl;
    std::cerr << "  <k_mer>: Primary k-mer size." << std::endl;
    std::cerr << "  <k0>: Fallback k-mer size (k0 < k_mer)." << std::endl;
    std::cerr << "  <L_seg>: Segment length." << std::endl;
    std::cerr << "  <search_range_global>: Search range limit for global matching." << std::endl;
    std::cerr << "  <T1_mismatch_perc>: Mismatch percentage threshold for a segment." << std::endl;
    std::cerr << "  <T2_mismatch_time>: Mismatch count threshold to trigger global fallback." << std::endl;
    std::cerr << "\nExample (SCCG encode):  " << program_name << " encode ref.fa target.fa out.intermediate 12 8 1000 100 10 5" << std::endl;
    std::cerr << "Example (decode):       " << program_name << " decode ref.fa encoded.sccg decoded.fa" << std::endl;
}

// Helper function to create dummy FASTA sequences for testing purposes.
std::string generateDummySequence(size_t length) {
    std::string dummy_sequence;
    dummy_sequence.reserve(length);
    std::random_device rd;
    std::mt19937 gen(rd());
    const char nucleotides[] = {'A', 'C', 'G', 'T'};
    std::uniform_int_distribution<> dis(0, 3); // Index for nucleotides array

    for (size_t i = 0; i < length; ++i) {
        dummy_sequence += nucleotides[dis(gen)];
    }
    return dummy_sequence;
}

// void encode(const std::string& reference_path, const std::string& target_file, int kmer_size) { // This basic encode is not used by SCCG main path
//     std::string reference = FastaParser::readSequence(reference_path);
//     std::string target = FastaParser::readSequence(target_file);

//     KmerIndex ref_kmer_index(kmer_size);
//     ref_kmer_index.buildIndex(reference);

// //    MatchFinder finder(ref_kmer_index, reference, target, kmer_size); // Old MatchFinder constructor
// //    std::vector<Match> matches = finder.findSignificantMatches(); // Old Match struct and method

// //    Encoder encoder; // Old Encoder
// //    encoder.encode(reference, target, matches, kmer_size); // Old Encoder
// }

bool saveToFile(const std::string& filename, const std::string& content) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
        return false;
    }
    outFile << content;
    outFile.close();
    return true;
}

bool readFromFile(const std::string& filename, std::string& content) {
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Error: Could not open file for reading: " << filename << std::endl;
        return false;
    }
    std::stringstream buffer;
    buffer << inFile.rdbuf();
    content = buffer.str();
    inFile.close();
    return true;
}

// Helper to ensure a file exists, or create a dummy one for testing.
void ensureTestFile(const std::string& filename, const std::string& default_header, size_t default_seq_len) {
    std::ifstream f(filename.c_str());
    if (f.good()) {
        f.close();
        // std::cout << "Using existing file: " << filename << std::endl;
        return;
    }
    f.close();
    std::cout << "File " << filename << " not found. Creating dummy file for testing." << std::endl;
    std::string seq = generateDummySequence(default_seq_len);
    if (!FastaParser::writeSequence(filename, default_header, seq)) {
        std::cerr << "Warning: Could not create dummy file: " << filename << std::endl;
    }
}

// Main logic for encoding (basic version, not SCCG full algorithm)
//void run_encode(const std::string& ref_fasta_path, const std::string& target_fasta_path,
//                const std::string& encoded_output_path, int kmer_size) {
//    std::cout << "--- Starting Encoding Process ---" << std::endl;
//    ensureTestFile(ref_fasta_path, "dummy_reference_header", 500);
//    ensureTestFile(target_fasta_path, "dummy_target_header", 400);
//
//    std::cout << "Reading reference sequence from: " << ref_fasta_path << std::endl;
//    std::string reference_sequence = FastaParser::readSequence(ref_fasta_path);
//    if (reference_sequence.empty()) {
//        throw std::runtime_error("Failed to read reference sequence or file is empty: " + ref_fasta_path);
//    }
//    std::cout << "Reference sequence length: " << reference_sequence.length() << std::endl;
//
//    std::cout << "Reading target sequence from: " << target_fasta_path << std::endl;
//    std::string target_sequence = FastaParser::readSequence(target_fasta_path);
//    if (target_sequence.empty()) {
//        throw std::runtime_error("Failed to read target sequence or file is empty: " + target_fasta_path);
//    }
//    std::cout << "Target sequence length: " << target_sequence.length() << std::endl;
//
//    std::cout << "Building K-mer index for reference (k=" << kmer_size << ")..." << std::endl;
//    KmerIndex ref_kmer_idx(kmer_size);
//    ref_kmer_idx.buildIndex(reference_sequence);
//
//    std::cout << "Initializing MatchFinder..." << std::endl;
//    MatchFinder match_finder(ref_kmer_idx, reference_sequence, target_sequence, kmer_size);
//
//    std::cout << "Finding significant matches..." << std::endl;
//    std::vector<Match> significant_matches = match_finder.findSignificantMatches();
//    std::cout << "Found " << significant_matches.size() << " significant matches." << std::endl;
//
//    std::cout << "Encoding target sequence..." << std::endl;
//    Encoder encoder(target_sequence, significant_matches);
//    if (encoder.encodeToFile(encoded_output_path)) {
//        std::cout << "Encoded data successfully saved to: " << encoded_output_path << std::endl;
//    } else {
//        throw std::runtime_error("Failed to save encoded data to file: " + encoded_output_path);
//    }
//    std::cout << "--- Encoding Process Complete ---" << std::endl;
//}
//
//// Main logic for decoding
//void run_decode(const std::string& ref_fasta_path, const std::string& encoded_input_path,
//                const std::string& decoded_fasta_output_path) {
//    std::cout << "--- Starting Decoding Process ---" << std::endl;
//    ensureTestFile(ref_fasta_path, "dummy_reference_header_for_decode", 500); // Ensure ref exists
//
//    std::cout << "Reading reference sequence from: " << ref_fasta_path << std::endl;
//    std::string reference_sequence = FastaParser::readSequence(ref_fasta_path);
//    if (reference_sequence.empty()) {
//        throw std::runtime_error("Failed to read reference sequence or file is empty: " + ref_fasta_path);
//    }
//
//    std::cout << "Decoding data from: " << encoded_input_path << std::endl;
//    Decoder decoder(reference_sequence);
//    std::string reconstructed_sequence = decoder.decodeFromFile(encoded_input_path);
//    std::cout << "Reconstructed sequence length: " << reconstructed_sequence.length() << std::endl;
//
//    if (FastaParser::writeSequence(decoded_fasta_output_path, "reconstructed_target", reconstructed_sequence)) {
//        std::cout << "Reconstructed sequence successfully saved to: " << decoded_fasta_output_path << std::endl;
//    } else {
//        throw std::runtime_error("Failed to save reconstructed sequence to file: " + decoded_fasta_output_path);
//    }
//
//    // Optional: Verification if original target is available (not directly in decode mode)
//    // To verify, you'd need the original target.fa. One could add a command like:
//    // program verify <original_target.fa> <reconstructed_target.fa>
//
//    std::cout << "--- Decoding Process Complete ---" << std::endl;
//}
//
//// The 'compress' function has been moved to the Compressor class.

int main(int argc, char* argv[]) {
    if (argc < 2) { // Not enough args for even the mode
        printUsage(argv[0]);
        return 1;
    }

    std::string mode = argv[1];
    int kmer_size = 12; // Default k-mer size for basic encode

    try {
        if (mode == "encode") {
            // argc check:
            // basic: program_name encode ref target out [kmer_size] -> 5 or 6
            // sccg:  program_name encode ref target out k k0 L search_L T1 T2 -> 11
            if (argc < 5 || (argc > 6 && argc != 11)) {
                printUsage(argv[0]);
                return 1;
            }
            std::string ref_path = argv[2];
            std::string target_in_path = argv[3];
            std::string output_path = argv[4];

            if (argc == 5) { // Basic encode call
                int k_val = 12; // Default k-mer size for SCCG
                int k0_val = 8; // Default k0 for SCCG
                size_t l_val = 1000; // Default segment length
                int search_limit_val = 100; // Default search range limit
                int t1_val = 10; // Default T1 mismatch percentage
                int t2_val = 5; // Default T2 mismatch time

                Compressor sccg_compressor; // Instance created
                sccg_compressor.compress(ref_path, target_in_path, output_path, // Called on instance
                                         k_val, k0_val, l_val,
                                         search_limit_val, t1_val, t2_val);
                std::cout << "SCCG compress method finished. Intermediate output is at: " << output_path << std::endl;
            } else {
                 // This case should be caught by the initial argc check
                 printUsage(argv[0]);
                 return 1;
            }

        } else if (mode == "decode") {
            if (argc != 5) { // decode ref encoded_in decoded_out
                printUsage(argv[0]);
                return 1;
            }
            std::string ref_path = argv[2];
            std::string encoded_in_path = argv[3];
            std::string decoded_out_path = argv[4];
//            run_decode(ref_path, encoded_in_path, decoded_out_path);

        } else {
            std::cerr << "Error: Invalid mode '" << mode << "'. Use 'encode' or 'decode'." << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
