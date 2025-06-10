#include <iostream>
#include <string>
#include <random>    // For random nucleotide generation
#include <stdexcept> // For std::exception
#include <algorithm>
#include <sstream>  // For std::stringstream
#include <fstream>  // For file operations
#include "src/Compressor.h" 
#include "src/Decompressor.h"

void printUsage(const std::string& prog) {
    std::cerr << "Usage:\n";
    std::cerr << "  " << prog << " encode <reference.fa> <target.fa> <output_dir> [k_mer k0 L_seg search_range_global T1 T2]\n";
    std::cerr << "  " << prog << " decode <reference.fa> <encoded_archive.7z> <decoded_output.fa>\n\n";
    std::cerr << "  " << prog << " check <sequence1.fa> <sequence2.fa>\n\n";

    std::cerr << "Modes:\n";
    std::cerr << "  encode   Compress <target.fa> against <reference.fa> into <output_dir>.\n";
    std::cerr << "           Optional parameters:\n";
    std::cerr << "             k_mer               = 21\n";
    std::cerr << "             k0                   = 8\n";
    std::cerr << "             L_seg                = 30000\n";
    std::cerr << "             search_range_global  = 100\n";
    std::cerr << "             T1 (mismatch ratio)  = 0.5\n";
    std::cerr << "             T2 (mismatch count)  = 4\n\n";

    std::cerr << "  decode   Reconstruct the FASTA from a single-chromosome .7z archive.\n\n";

    std::cerr << "  check    Compare two FASTA files ignoring newlines.\n";

    std::cerr << "Examples:\n";
    std::cerr << "  " << prog << " encode ref.fa target.fa ./out_dir\n";
    std::cerr << "  " << prog << " encode ref.fa target.fa ./out_dir 12 8 1000 100 0.1 2\n";
    std::cerr << "  " << prog << " decode ref.fa out_dir/chr1.7z chr1_decoded.fa\n";
    std::cerr << "  " << prog << " check seq1.fa seq2.fa\n";
}

bool areFilesEqualIgnoringNewlines(const std::string& filepath1, const std::string& filepath2) {
    std::ifstream file1(filepath1);
    std::ifstream file2(filepath2);

    if (!file1.is_open() || !file2.is_open()) {
        std::cerr << "Error opening files for comparison." << std::endl;
        return false;
    }

    std::stringstream buffer1, buffer2;
    buffer1 << file1.rdbuf();
    buffer2 << file2.rdbuf();

    std::string content1 = buffer1.str();
    std::string content2 = buffer2.str();

    content1.erase(std::remove(content1.begin(), content1.end(), '\r'), content1.end());
    content1.erase(std::remove(content1.begin(), content1.end(), '\n'), content1.end());
    content2.erase(std::remove(content2.begin(), content2.end(), '\r'), content2.end());
    content2.erase(std::remove(content2.begin(), content2.end(), '\n'), content2.end());

    if (content1 != content2) {
        std::cerr << "Files differ!" << std::endl;
        std::cerr << "File 1 (" << filepath1 << "): " << content1 << std::endl;
        std::cerr << "File 2 (" << filepath2 << "): " << content2 << std::endl;
        return false;
    }

    return true;
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
            std::string ref      = argv[2];
            std::string archive  = argv[3];
            std::string outFasta = argv[4];

            Decompressor::decompress(ref, archive, outFasta);
        } else if (mode == "check") {
            std::string seq1 = argv[2];
            std::string seq2 = argv[3];
            if (argc != 4) {
                printUsage(argv[0]);
                return 1;
            }

            if (areFilesEqualIgnoringNewlines(seq1, seq2)) {
                std::cout << "The sequences in " << seq1 << " and " << seq2 << " are equal." << std::endl;
            } else {
                std::cout << "The sequences in " << seq1 << " and " << seq2 << " are NOT equal." << std::endl;
            }
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
