#include "Encoder.h"
#include <fstream>
#include <iostream>
#include <vector>
#include "MatchFinder.h"

void Encoder::encode(
    const std::string& reference,
    const std::string& target,
    const std::vector<Match>& matches,
    int kmer_size,
    const std::string& out_filename
) {
    std::ofstream out(out_filename);
    if (!out) {
        std::cerr << "Error: Cannot open output file:" << out_filename << std::endl;
        return;
    }

    out << "# Encoded file - reference_length=" << reference.size() << " target_length=" << target.size() << " k=" << kmer_size << "\n";

    size_t prev_tar_end = 0;
    for (const auto& match : matches) {
        if (match.target_pos > prev_tar_end) {
            out << "MISMATCH\t" << prev_tar_end << "\t"
                << match.target_pos - 1 << "\t"
                << target.substr(prev_tar_end, match.target_pos - prev_tar_end) << "\n";
        }
        out << "MATCH\t"
            << match.ref_pos << "\t" << match.target_pos << "\t" << match.length << "\n";
        prev_tar_end = match.target_pos + match.length;
    }
    if (prev_tar_end < target.size()) {
        out << "MISMATCH\t" << prev_tar_end << "\t" << target.size() - 1 << "\t"
            << target.substr(prev_tar_end) << "\n";
    }
    out.close();
    std::cout << "Encoding complete. Output written to " << out_filename << std::endl;
}
