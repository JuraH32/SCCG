#include <iostream>

#define KMER_LENGTH 21

// In main parse the command line arguments and call the Matcher class
int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <reference_genome_path> <target_genome_path> <output_path>" << std::endl;
        return 1;
    }

    // Reference genome path
    std::string reference_genome_path = argv[0];
    // Target genome path
    std::string target_genome_path = argv[1];
    // Output path
    std::string output_path = argv[2];

    // Hello world message
    std::cout << "Hello, world!" << std::endl;

    return 0;
}