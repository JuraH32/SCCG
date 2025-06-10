#include <iostream>
#include <string>
#include <random>    // For random nucleotide generation
#include <stdexcept> // For std::exception
#include <algorithm>
#include <sstream>  // For std::stringstream
#include <fstream>  // For file operations
#include "src/Compressor.h" 
#include "src/Decompressor.h"
#include <filesystem>

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
            std::string ref_path = argv[2];
            std::string target_in_path = argv[3];
            std::string output_path = argv[4];

            namespace fs = std::filesystem;
            bool ref_is_dir = fs::is_directory(ref_path);
            bool target_is_dir = fs::is_directory(target_in_path);

            auto encode_start = std::chrono::high_resolution_clock::now();

            if (ref_is_dir && target_is_dir) {
                for (const auto& entry : fs::directory_iterator(target_in_path)) {
                    if (entry.is_regular_file()) {
                        std::string filename = entry.path().filename().string();
                        fs::path ref_file = fs::path(ref_path) / filename;
                        if (fs::exists(ref_file)) {
                            auto chr_start = std::chrono::high_resolution_clock::now();
                            std::string target_file = entry.path().string();
                            std::string out_file = (fs::path(output_path) / filename).string();
                            Compressor sccg_compressor;
                            sccg_compressor.compress(ref_file.string(), target_file, out_file,
                                21, 8, 30000, 100, 0.5, 4);
                            auto chr_end = std::chrono::high_resolution_clock::now();
                            std::chrono::duration<double> chr_elapsed = chr_end - chr_start;
                            std::cout << "Compressed: " << filename << " in " << chr_elapsed.count() << " seconds." << std::endl;
                        } else {
                            std::cerr << "Reference file not found for: " << filename << std::endl;
                        }
                    }
                }
            } else if (!ref_is_dir && !target_is_dir) {
                auto chr_start = std::chrono::high_resolution_clock::now();
                Compressor sccg_compressor;
                sccg_compressor.compress(ref_path, target_in_path, output_path,
                    21, 8, 30000, 100, 0.5, 4);
                auto chr_end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> chr_elapsed = chr_end - chr_start;
                std::cout << "SCCG compress method finished. Intermediate output is at: " << output_path << std::endl;
                std::cout << "Compression time: " << chr_elapsed.count() << " seconds." << std::endl;
            } else {
                std::cerr << "Error: Reference and target must both be files or both be directories." << std::endl;
                return 1;
            }

            auto encode_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> encode_elapsed = encode_end - encode_start;
            std::cout << "Total compression time: " << encode_elapsed.count() << " seconds." << std::endl;

        } else if (mode == "decode") {
            if (argc != 5) { // decode ref encoded_in decoded_out
                printUsage(argv[0]);
                return 1;
            }
            std::string ref_path = argv[2];
            std::string archive_path = argv[3];
            std::string output_path = argv[4];

            namespace fs = std::filesystem;
            bool ref_is_dir = fs::is_directory(ref_path);
            bool archive_is_dir = fs::is_directory(archive_path);

            auto decode_start = std::chrono::high_resolution_clock::now();

            if (ref_is_dir && archive_is_dir) {
                fs::create_directories(output_path);
                for (const auto& entry : fs::directory_iterator(archive_path)) {
                    if (entry.is_regular_file()) {
                        std::string filename = entry.path().filename().string();
                        fs::path ref_file = fs::path(ref_path) / filename;
                        if (fs::exists(ref_file)) {
                            auto chr_start = std::chrono::high_resolution_clock::now();
                            std::string archive_file = entry.path().string();
                            std::string out_file = (fs::path(output_path) / filename).string();
                            Decompressor::decompress(ref_file.string(), archive_file, out_file);
                            auto chr_end = std::chrono::high_resolution_clock::now();
                            std::chrono::duration<double> chr_elapsed = chr_end - chr_start;
                            std::cout << "Decompressed: " << filename << " in " << chr_elapsed.count() << " seconds." << std::endl;
                        } else {
                            std::cerr << "Reference file not found for: " << filename << std::endl;
                        }
                    }
                }
            } else if (!ref_is_dir && !archive_is_dir) {
                auto chr_start = std::chrono::high_resolution_clock::now();
                Decompressor::decompress(ref_path, archive_path, output_path);
                auto chr_end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> chr_elapsed = chr_end - chr_start;
                std::cout << "Decompression time: " << chr_elapsed.count() << " seconds." << std::endl;
            } else {
                std::cerr << "Error: Reference and archive must both be files or both be directories." << std::endl;
                return 1;
            }

            auto decode_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> decode_elapsed = decode_end - decode_start;
            std::cout << "Total decompression time: " << decode_elapsed.count() << " seconds." << std::endl;
        }
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 1;
    }
    return 0;
}
