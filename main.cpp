#include <iostream>
#include <string>
#include <vector>
#include <random> // For random nucleotide generation


// Function to print usage instructions for the command-line tool.
void printUsage(const std::string& program_name) {
    std::cerr << "Usage: " << program_name << " <reference.fa> <target.fa> [kmer_size]" << std::endl;
    std::cerr << "  <reference.fa>: Path to the reference genome FASTA file." << std::endl;
    std::cerr << "  <target.fa>: Path to the target genome FASTA file to be compared." << std::endl;
    std::cerr << "  [kmer_size]: Optional. Integer size for k-mers (default: 12)." << std::endl;
    std::cerr << "\nExample: " << program_name << " ref.fa target.fa 10" << std::endl;
}

// Helper function to create dummy FASTA sequences for testing purposes.
std::string createDummyFastaSequence(std::int16_t sequenceLength) {
    std::string dummy_sequence;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<char> nucleotides = {'A', 'C', 'G', 'T'};
    std::uniform_int_distribution<> dis(0, 3);

    for (int i = 0; i < sequenceLength; ++i) {
        dummy_sequence += nucleotides[dis(gen)];
    }
    return dummy_sequence;
}

void encode(const std::string& reference_path, const std::string& target_file, int kmer_size) {
    //TODO: Implement the encoding logic using the Encoder class.
}

void decode(const std::string& reference_path, const std::string& target_file) {
    //TODO: Implement the decoding logic using the Decoder class.
}

int main(int argc, char* argv[]) {
    // Check for the correct number of command-line arguments.
    if (argc < 4 || argc > 6) {
        printUsage(argv[0]);
        return 1;
    }

    std::string mode = argv[1];
    std::string reference_path = argv[2];
    std::string target_file = argv[3];

    int kmer_size = 12; // Default k-mer size

    if (argc >= 5) {
        try {
            kmer_size = std::stoi(argv[4]);
            if (kmer_size <= 0) {
                throw std::invalid_argument("k-mer size must be a positive integer.");
            }
        } catch (const std::exception& e) {
            std::cerr << "Error: Invalid k-mer size. " << e.what() << std::endl;
            return 1;
        }
    }

    if (argc == 6) {
        // If a fifth argument is provided, treat it as the output file name.
        target_file = argv[5];
    }

    if (mode == "encode") {
        encode(reference_path, target_file, kmer_size);
    } else if (mode == "decode") {
        decode(reference_path, target_file);
    } else {
        std::cerr << "Error: Invalid mode. Use 'encode' or 'decode'." << std::endl;
        printUsage(argv[0]);
        return 1;
    }







    return 0;
}
