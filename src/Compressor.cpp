#include "Compressor.h"
#include "MatchFinder.h" // For MatchFinder and Match struct
#include "Position.h"
#include "Utils.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm> // For std::min, std::all_of, std::toupper, std::islower
#include <vector>
#include <string>
#include <iomanip>
#include <filesystem>

Compressor::Compressor() {
}

bool Compressor::saveToFile(const std::string &filename, const std::string &content) {
    std::ofstream outFile(filename, std::ios::app);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
        return false;
    }
    outFile << content;
    outFile.close();
    return true;
}

std::string Compressor::readSequenceFromFile(const std::string &filename) {
    std::ifstream file(filename);

    //Print the absolute path of the file being read
    std::cout << "Reading sequence from file: " << std::filesystem::absolute(filename) << std::endl;

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file: " << filename << std::endl;
        return "";
    }

    std::stringstream buffer;
    std::string line;
    std::string meta_data_line;
    std::string sequence;

    // Read the first line as metadata
    if (std::getline(file, meta_data_line)) {
        meta_data = meta_data_line;
    }

    // Read the rest of the file into the buffer
    while (std::getline(file, line)) {
        buffer << line;
    }

    file.close();
    return buffer.str();
}

struct ParsedMatchInfo {
    bool is_valid_pure_match = false;
    int  ref_pos = 0;
    int  length  = 0;
};

// Helper function to parse a "MATCH R L NOMISMATCH" line
bool parsePureMatchLine(const std::string& line, ParsedMatchInfo& out_match_info) {
    auto comma = line.find(',');
    if (comma == std::string::npos) {
        out_match_info.is_valid_pure_match = false;
        return false;
    }

    // split into two pieces
    std::string left  = line.substr(0, comma);
    std::string right = line.substr(comma + 1);

    // both sides must be nonempty and all digits
    if (left.empty() || right.empty() ||
        !std::all_of(left.begin(),  left.end(),  ::isdigit) ||
        !std::all_of(right.begin(), right.end(), ::isdigit))
    {
        out_match_info.is_valid_pure_match = false;
        return false;
    }

    // parse them
    out_match_info.ref_pos = std::stoi(left);
    out_match_info.length  = std::stoi(right);
    out_match_info.is_valid_pure_match = true;
    return true;
}

static bool parseMatchLine(const std::string &ln, ParsedMatchInfo &out) {
    std::istringstream iss(ln);
    std::string kw, nom;
    int r, l;
    if (!(iss >> kw >> r >> l >> nom)) return out.is_valid_pure_match = false;
    if (kw != "MATCH" || nom != "NOMISMATCH") return out.is_valid_pure_match = false;
    out.is_valid_pure_match = true;
    out.ref_pos = r;
    out.length  = l;
    return true;
}

void Compressor::postprocess(std::string& temp_in, const std::string& final_out) {
    std::cout << "Postprocessing: merging and delta-encoding…\n";
    std::ifstream in(temp_in);
    if (!in) throw std::runtime_error("postprocess: cannot open " + temp_in);

    std::vector<std::string> merged;
    std::string line;
    bool   inMerge = false;
    int    mstart  = 0, mlen = 0;

    auto flush = [&]() {
      if (!inMerge) return;
      merged.push_back(
        "MATCH " + std::to_string(mstart) + " " + std::to_string(mlen) + " NOMISMATCH"
      );
      inMerge = false;
    };

    while (std::getline(in, line)) {
      ParsedMatchInfo pm;
      if (parseMatchLine(line, pm)) {
        if (!inMerge) {
          mstart = pm.ref_pos;
          mlen   = pm.length;
          inMerge = true;
        } else if (pm.ref_pos == mstart + mlen) {
          mlen += pm.length;
        } else {
          flush();
          mstart = pm.ref_pos;
          mlen   = pm.length;
          inMerge = true;
        }
      } else {
        flush();
        merged.push_back(line);
      }
    }
    flush();
    in.close();

    std::ofstream out(final_out);
    if (!out) throw std::runtime_error("postprocess: cannot open " + final_out);

    long long prev_end = -1;
    for (auto &l : merged) {
      ParsedMatchInfo pm;
      if (parseMatchLine(l, pm)) {
        long long delta = (long long)pm.ref_pos - prev_end;
        out << delta << "," << pm.length << "\n";
        prev_end = (long long)pm.ref_pos + pm.length - 1;
      } else {
        out << l << "\n";
      }
    }
    out.close();

    std::cout << "Postprocessing complete, result in " << final_out << "\n";
}

std::vector<Position> Compressor::getPositions(const std::string &sequence,
                                            std::function<bool(char)> is_of_interest) const {
std::vector<Position> positions;
int last_pos_begin = -1;
for (size_t i = 0; i < sequence.length(); ++i) {
    if (is_of_interest(sequence[i])) {
        if (last_pos_begin == -1) {
            last_pos_begin = i; // Start of a new segment
        }
    } else {
        if (last_pos_begin != -1) {
            // We found the end of a segment
            Position pos = Position(-1, -1, last_pos_begin, i - 1);
            positions.push_back(pos);
            last_pos_begin = -1; // Reset for the next segment
        }
    }
}
if (last_pos_begin != -1) {
    Position pos = Position(-1, -1, last_pos_begin, sequence.length() - 1);
    positions.push_back(pos);
}
return positions;
}


void Compressor::savePositionsToFile(
const std::string &filename,
const std::vector<Position> &positions,
const std::string &label_prefix
) {
std::ostringstream buffer;
buffer << label_prefix;

bool is_lowercase = (label_prefix.rfind("LOWERCASE", 0) == 0);

int prev_excl_end = 0;
bool first_segment = true;

for (const auto &pos : positions) {
    int start  = pos.startInTarget;
    int end    = pos.endInTarget;
    int length = end - start + 1;
    if (length <= 0)
        continue;

    int delta = start - prev_excl_end;

    if (is_lowercase) {
        if (first_segment) {
            length -= 1;
        } else {
            delta  += 1;
            length -= 1;
        }
    }

    std::cout << "segment: [" << start << "," << end << "]"
                << "  delta: " << delta
                << "  length: " << length
                << std::endl;

    buffer << delta << " " << length << " ";

    prev_excl_end = end + 1;
    first_segment = false;
}

buffer << "\n";

if (!saveToFile(filename, buffer.str())) {
    std::cerr << "Error: Could not save positions to file: "
                << filename << std::endl;
}
}


void Compressor::saveAlignmentSegmentsToFile(
        const std::string &filename,
        const std::vector<AlignmentSegment> &segments,
        int ref_segment_offset,
        int target_segment_offset
) {
    std::ostringstream oss;
    std::string mismatch_buf;

    for (const auto &segment : segments) {
        if (segment.is_match) {
            if (!mismatch_buf.empty()) {
                oss << mismatch_buf << "\n";
                mismatch_buf.clear();
            }

            int R0 = ref_segment_offset + segment.region.startInReference;
            int L  = segment.region.length();
            if (L > 0) {
                oss << R0 << "," << L << "\n";
            }
        } else {
            mismatch_buf += segment.mismatched_sequence;
        }
    }
    if (!mismatch_buf.empty()) {
        oss << mismatch_buf << "\n";
    }

    if (!oss.str().empty()) {
        saveToFile(filename, oss.str());
    }
}

void Compressor::compress(const std::string &ref_fasta_path,
                        const std::string &target_fasta_path,
                        const std::string &output_path,
                        int kmer_size,
                        int k0,
                        size_t segment_length,
                        int search_range_limit,
                        double mismatch_threshold_T1, // Changed from int to double
                        int mismatch_threshold_T2) {

    if (!std::filesystem::exists(output_path)) {
        std::filesystem::create_directories(output_path);
    }

    std::cout << "Loading reference sequence from: " << ref_fasta_path << std::endl;
    std::string reference_sequence = readSequenceFromFile(ref_fasta_path);

    if (reference_sequence.empty()) {
        throw std::runtime_error("Failed to read reference sequence or file is empty: " + ref_fasta_path);
    }

    std::cout << "Loading target sequence from: " << target_fasta_path << std::endl;
    std::string target_sequence = readSequenceFromFile(target_fasta_path);

    std::vector<std::string> target_path_parts = Utils::splitPath(target_fasta_path);
    std::string target_filename = target_path_parts.empty() ? "target" : target_path_parts.back();

    std::string final_filename = output_path + "/" + target_filename;
    std::string intermediate_file = output_path + "/intermediate.txt";
    Utils::removeFileIfExists(final_filename, true); 
    Utils::removeFileIfExists(intermediate_file, true); 

    if (target_sequence.empty()) {
        throw std::runtime_error("Failed to read target sequence or file is empty: " + target_fasta_path);
    }

    std::ostringstream metadata_stream;
    metadata_stream << meta_data << "\n" << target_sequence.length() << "\n";

    if (!saveToFile(intermediate_file, metadata_stream.str())) {
        throw std::runtime_error("Failed to save metadata to file: " + intermediate_file);
    }

    std::vector<Position> lowercase_positions_T;
    lowercase_positions_T = getPositions(target_sequence, [](char c) { return std::islower(c); });
    savePositionsToFile(intermediate_file, lowercase_positions_T, "LOWERCASE_POSITIONS: ");

    std::cout << "Lowercase positions in target sequence saved to file" << std::endl;

    for (char &c: reference_sequence) {
        c = std::toupper(c);
    }

    for (char &c: target_sequence) {
        c = std::toupper(c);
    }

    size_t m = (reference_sequence.length() + segment_length - 1) / segment_length; // ceiling division
    size_t n = (target_sequence.length() + segment_length - 1) / segment_length;    // ceiling division

    std::vector<std::string> ref_segments;
    std::vector<std::string> target_segments;

    int mismatch_t_ime = 0;
    bool global_fallback = false;
    int controuble = 0;
    bool is_con = false;

    for (size_t i = 0; i < m; ++i) {
        size_t start = i * segment_length;
        size_t length = std::min(segment_length, reference_sequence.length() - start);
        ref_segments.push_back(reference_sequence.substr(start, length));
    }

    for (size_t i = 0; i < n; ++i) {
        size_t start = i * segment_length;
        size_t length = std::min(segment_length, target_sequence.length() - start);
        target_segments.push_back(target_sequence.substr(start, length));
    }

    for (size_t i = 0; i < n; ++i) {
        if (i % 100 == 0) {
            std::cout << "Processing segment " << i + 1 << " of " << n << std::endl;
        }

        MatchFinder match_finder(ref_segments[i], target_segments[i], kmer_size, search_range_limit);
        std::vector<AlignmentSegment> matches = match_finder.findMatches(false);

        if (matches.empty()) {
            std::cout << "No matches found for segment " << i + 1 << ". Attempting with k0." << std::endl;
            MatchFinder matchFinder(ref_segments[i], target_segments[i], k0, search_range_limit);
            matches = matchFinder.findMatches(false);
        }

        if (matches.empty()) {
            std::cout << "No matches found for segment " << i + 1 << " even with k0." << std::endl;
            mismatch_t_ime++;
            if (is_con) {
                controuble++;
            } else {
                controuble = 1;
                is_con = true;
            }

            std::string intermediate_data = "RAW_TARGET_SEGMENT " + std::to_string(i) + " " + target_segments[i] + "\n";
            if (!saveToFile(intermediate_file, intermediate_data)) {
                throw std::runtime_error("Failed to save intermediate data stream to file: " + intermediate_file);
            }
            continue;
        }

        int current_ref_segment_offset = i * segment_length;
        int current_target_segment_offset = i * segment_length;
        saveAlignmentSegmentsToFile(intermediate_file, matches, current_ref_segment_offset, current_target_segment_offset);


        int unmatched_chars = 0;
        for (const auto &segment: matches) {
            if (!segment.is_match) {
                unmatched_chars += segment.mismatched_sequence.length();
            }
        }

        size_t total_length = target_segments[i].length();
        if (unmatched_chars / static_cast<double>(total_length) > mismatch_threshold_T1) {
            mismatch_t_ime++;
        }

                    if (is_con && controuble > 2) {
            mismatch_t_ime -= controuble;
            if (mismatch_t_ime < 0) mismatch_t_ime = 0;
        }

        controuble = 1;
        is_con = false;

        if (mismatch_t_ime > mismatch_threshold_T2) {
            std::cout << "Mismatch threshold exceeded. Triggering global fallback." << std::endl;
            global_fallback = true;
            break;
        }
    }

    if (global_fallback) {
        std::cout << "Performing global fallback due to mismatch threshold." << std::endl;

        Utils::removeFileIfExists(intermediate_file, true); // Remove intermediate file if it exists

        std::ostringstream global_metadata_stream;
        global_metadata_stream << meta_data << "\n" << target_sequence.length() << "\n";

        if (!saveToFile(intermediate_file, global_metadata_stream.str())) {
            throw std::runtime_error("Failed to save metadata to file in global fallback: " + intermediate_file);
        }

        std::vector<Position> n_fragments_T = getPositions(target_sequence, [](char c) { return c == 'N'; });

        savePositionsToFile(intermediate_file, lowercase_positions_T, "LOWERCASE_POSITIONS: ");
        savePositionsToFile(intermediate_file, n_fragments_T, "N_POSITIONS: ");

        // Remove 'N' characters from target sequence
        std::string n_free_target;
        n_free_target.reserve(target_sequence.length());
        for (char c : target_sequence) {
            if (c != 'N') {
                n_free_target += c;
            }
        }
        target_sequence = n_free_target;

        // Remove 'N' characters from reference sequence
        std::string n_free_reference;
        n_free_reference.reserve(reference_sequence.length());
        for (char c : reference_sequence) {
            if (c != 'N') {
                n_free_reference += c;
            }
        }
        reference_sequence = n_free_reference;

        std::cout << "Initializing global matching." << std::endl;

        MatchFinder global_match_finder(n_free_reference, n_free_target, kmer_size, search_range_limit);

        std::cout << "Starting global matching." << std::endl;

        std::vector<AlignmentSegment> global_matches = global_match_finder.findMatches(true);

        if (global_matches.empty()) {
            std::cout << "No global matches found. Storing raw target sequence." << std::endl;
            std::string raw_target_data = "RAW_GLOBAL_TARGET " + n_free_target + "\n";
            if (!saveToFile(intermediate_file, raw_target_data)) {
                throw std::runtime_error("Failed to save raw target sequence to file: " + intermediate_file);
            }
        } else {
            saveAlignmentSegmentsToFile(intermediate_file, global_matches, 0, 0); // Offsets are 0 for global
        }
    }

    // Postprocess the intermediate file to apply delta coding and save final output
    std::cout << "Postprocessing intermediate file: " << intermediate_file << std::endl;
    postprocess(intermediate_file, final_filename);
    std::cout << "Compression complete. Output saved to: " << final_filename << std::endl;

    Utils::compressWith7zip(final_filename);

    std::cout << "Final output compressed with 7zip." << std::endl;

    // Clean up intermediate file
//    Utils::removeFileIfExists(intermediate_file);
}
