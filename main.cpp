#include <iostream>
#include <string>
#include <random>    // For random nucleotide generation
#include <stdexcept> // For std::exception
#include <algorithm>
#include <sstream>  // For std::stringstream
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
    std::cerr << "  <output_path>: Path to save the encoded output (for encoding) or decoded FASTA file (for decoding)." << std::endl;
    std::cerr << "  [kmer_size]: Optional. Integer size for k-mers (default: 12, for basic encoding only)." << std::endl;
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
                int k_val = 21; // Default k-mer size for SCCG
                int k0_val = 8; // Default k0 for SCCG
                size_t l_val = 30000; // Default segment length
                int search_limit_val = 100; // Default search range limit
                double t1_val = 0.5; // Default T1 mismatch percentage
                int t2_val = 4; // Default T2 mismatch time

                Compressor sccg_compressor; // Instance created
                sccg_compressor.compress(ref_path, target_in_path, output_path, // Called on instance
                                         k_val, k0_val, l_val,
                                         search_limit_val, t1_val, t2_val); // Pass t1_val directly as double
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
