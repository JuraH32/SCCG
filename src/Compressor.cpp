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
#include <iomanip> // For std::fixed, std::setprecision (if used for printing percentages)

Compressor::Compressor() {
    // Constructor implementation (if needed)
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

// Helper struct for postprocessing
struct ParsedMatchInfo {
    bool is_valid_pure_match = false; // True if it's a "MATCH R L NOMISMATCH" line
    int ref_pos = 0;
    int length = 0;
    // std::string original_line; // Keep the whole line if not a pure match
};

// Helper function to parse a "MATCH R L NOMISMATCH" line
bool parsePureMatchLine(const std::string& line, ParsedMatchInfo& out_match_info) {
    std::istringstream iss(line);
    std::string keyword;
    iss >> keyword;
    if (keyword == "MATCH") {
        int r_pos, l;
        std::string mismatch_token;
        if (iss >> r_pos >> l >> mismatch_token) {
            if (mismatch_token == "NOMISMATCH") {
                out_match_info.is_valid_pure_match = true;
                out_match_info.ref_pos = r_pos;
                out_match_info.length = l;
                return true;
            }
        }
    }
    out_match_info.is_valid_pure_match = false;
    return false;
}


void Compressor::postprocess(std::string& temp_input_filepath, const std::string& final_output_filepath) {
    std::cout << "Postprocessing: Merging continuous matches and applying delta encoding..." << std::endl;
    std::ifstream temp_file_stream(temp_input_filepath);
    if (!temp_file_stream.is_open()) {
        throw std::runtime_error("Postprocess: Failed to open temporary input file: " + temp_input_filepath);
    }

    std::vector<std::string> lines_after_merge_pass;
    std::string line;

    ParsedMatchInfo current_merged_segment_info;
    bool in_merged_segment_mode = false;

    // Pass 1: Merge continuous "MATCH ref_pos length NOMISMATCH" lines
    while (std::getline(temp_file_stream, line)) {
        ParsedMatchInfo parsed_line_info;
        if (parsePureMatchLine(line, parsed_line_info)) { // It's a "MATCH R L NOMISMATCH" line
            if (in_merged_segment_mode) {
                // Check for continuity: current_start == prev_start + prev_length
                if (parsed_line_info.ref_pos == (current_merged_segment_info.ref_pos + current_merged_segment_info.length)) {
                    current_merged_segment_info.length += parsed_line_info.length; // Extend current merged segment
                } else {
                    // Not continuous, finalize previous merged segment
                    lines_after_merge_pass.push_back("MATCH " + std::to_string(current_merged_segment_info.ref_pos) + " " + std::to_string(current_merged_segment_info.length) + " NOMISMATCH");
                    // Start new merged segment
                    current_merged_segment_info = parsed_line_info;
                }
            } else { // Start a new merged segment
                current_merged_segment_info = parsed_line_info;
                in_merged_segment_mode = true;
            }
        } else { // Not a pure match line (e.g., header, raw data, or MATCH with mismatch string)
            if (in_merged_segment_mode) {
                // Finalize any ongoing merged segment
                lines_after_merge_pass.push_back("MATCH " + std::to_string(current_merged_segment_info.ref_pos) + " " + std::to_string(current_merged_segment_info.length) + " NOMISMATCH");
                in_merged_segment_mode = false;
            }
            lines_after_merge_pass.push_back(line); // Pass through this line
        }
    }
    // After loop, finalize any remaining merged segment
    if (in_merged_segment_mode) {
        lines_after_merge_pass.push_back("MATCH " + std::to_string(current_merged_segment_info.ref_pos) + " " + std::to_string(current_merged_segment_info.length) + " NOMISMATCH");
    }
    temp_file_stream.close();

    // Pass 2: Delta encode "MATCH ref_pos length NOMISMATCH" lines
    std::ofstream final_output_stream(final_output_filepath);
    if (!final_output_stream.is_open()) {
        throw std::runtime_error("Postprocess: Failed to open final output file: " + final_output_filepath);
    }

    long long last_absolute_ref_end = -1; // Using long long for safety, though int might suffice. -1 indicates no previous match.
                                         // Java used prev_end_coord for delta on begin.
                                         // If prev_end_coord is inclusive end, then begin - prev_end_coord.

    for (const std::string& merged_line : lines_after_merge_pass) {
        ParsedMatchInfo parsed_delta_info;
        if (parsePureMatchLine(merged_line, parsed_delta_info)) { // It's a "MATCH R L NOMISMATCH" line
            int delta_ref_pos;
            if (last_absolute_ref_end == -1) { // First pure match encountered in this pass
                delta_ref_pos = parsed_delta_info.ref_pos;
            } else {
                // Delta is current_start - (previous_inclusive_end).
                // If adjacent, current_start = previous_inclusive_end + 1, so delta = 1.
                delta_ref_pos = parsed_delta_info.ref_pos - static_cast<int>(last_absolute_ref_end);
            }
            final_output_stream << "MATCH " << delta_ref_pos << " " << parsed_delta_info.length << " NOMISMATCH\n";
            last_absolute_ref_end = static_cast<long long>(parsed_delta_info.ref_pos) + parsed_delta_info.length -1; // Update inclusive end
        } else { // Not a pure match line, pass through
            final_output_stream << merged_line << "\n";
        }
    }
    final_output_stream.close();
    std::cout << "Postprocessing complete. Output at: " << final_output_filepath << std::endl;
}

// getPositions is passed a sequence and a function which determines whether a character is of interest.
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
    // Handle case where the last segment goes to the end of the sequence
    if (last_pos_begin != -1) {
        Position pos = Position(-1, -1, last_pos_begin, sequence.length() - 1);
        positions.push_back(pos);
    }
    return positions;
}

void Compressor::savePositionsToFile(const std::string &filename, const std::vector<Position> &positions, const std::string& label_prefix) {
    std::ostringstream buffer;
    buffer << label_prefix;
    int previous_inclusive_end = -1;

    for (const auto &pos: positions) {
        int current_start = pos.startInTarget;
        int current_end = pos.endInTarget;
        int length = current_end - current_start + 1;

        if (length <= 0) continue; // Should not happen if getPositions is correct

        int delta_start_val;
        if (previous_inclusive_end == -1) { // First position
            delta_start_val = current_start; // Absolute start for the first one
        } else {
            delta_start_val = current_start - (previous_inclusive_end + 1); // Gap from end of previous
        }
        buffer << delta_start_val << " " << length << " ";
        previous_inclusive_end = current_end;
    }

    // Add newline if the label was written, even if no positions.
    // This matches Java's behavior of writing a newline after L_list/N_list content.
    buffer << "\n";

    if (!saveToFile(filename, buffer.str())) {
        std::cerr << "Error: Could not save positions to file: " << filename << std::endl;
    }
}

void Compressor::saveAlignmentSegmentsToFile(
        const std::string &filename,
        const std::vector<AlignmentSegment> &segments,
        int ref_segment_offset,
        int target_segment_offset
) {
    std::ostringstream oss;

    for (const auto &segment: segments) {
        if (segment.is_match) {
            int absolute_ref_start = ref_segment_offset + segment.region.startInReference;
            int length = segment.region.length();
            if (length <= 0) continue;

            oss << "MATCH "
                << absolute_ref_start << " "
                << length << " NOMISMATCH\n";
        } else {
            int absolute_target_start = target_segment_offset + segment.region.startInTarget;
            if (segment.mismatched_sequence.empty()) continue;

            oss << "MISMATCH "
                << absolute_target_start << " "
                << segment.mismatched_sequence << "\n";
        }
    }

    if (!oss.str().empty()) { // Only save if there's content
        if (!saveToFile(filename, oss.str())) {
            std::cerr << "Error: Could not save alignment segments to file: " << filename << std::endl;
        }
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

    // Load the reference sequence
    std::cout << "Loading reference sequence from: " << ref_fasta_path << std::endl;
    std::string reference_sequence = readSequenceFromFile(ref_fasta_path);

    if (reference_sequence.empty()) {
        throw std::runtime_error("Failed to read reference sequence or file is empty: " + ref_fasta_path);
    }

    // Load the target sequence
    std::cout << "Loading target sequence from: " << target_fasta_path << std::endl;
    std::string target_sequence = readSequenceFromFile(target_fasta_path);

    // Split target path to get the filename
    std::vector<std::string> target_path_parts = Utils::splitPath(target_fasta_path);
    std::string target_filename = target_path_parts.empty() ? "target" : target_path_parts.back();

    std::string final_filename = output_path + "/" + target_filename;
    std::string intermediate_file = output_path + "/intermediate.txt";
    Utils::removeFileIfExists(final_filename, true); // Remove final file if it exists
    Utils::removeFileIfExists(intermediate_file, true); // Remove intermediate file if it exists

    if (target_sequence.empty()) {
        throw std::runtime_error("Failed to read target sequence or file is empty: " + target_fasta_path);
    }

    // Write the metadata and lengths to the intermediate file first
    std::ostringstream metadata_stream;
    metadata_stream << meta_data << "\n" << target_sequence.length() << "\n";

    if (!saveToFile(intermediate_file, metadata_stream.str())) {
        throw std::runtime_error("Failed to save metadata to file: " + intermediate_file);
    }

    // Line 1 (part 1): Store lowercase positions from target sequence and convert to uppercase
    std::vector<Position> lowercase_positions_T;
    lowercase_positions_T = getPositions(target_sequence, [](char c) { return std::islower(c); });
    savePositionsToFile(intermediate_file, lowercase_positions_T, "LOWERCASE_POSITIONS: ");

    std::cout << "Lowercase positions in target sequence saved to file" << std::endl;

    // Convert reference sequence to uppercase
    for (char &c: reference_sequence) {
        c = std::toupper(c);
    }

    // Calculate number of segments
    size_t m = (reference_sequence.length() + segment_length - 1) / segment_length; // ceiling division
    size_t n = (target_sequence.length() + segment_length - 1) / segment_length;    // ceiling division

    // Create segments
    std::vector<std::string> ref_segments;
    std::vector<std::string> target_segments;

    // Create mismatch_t_ime counter
    int mismatch_t_ime = 0;
    bool global_fallback = false;

    // Divide reference sequence into segments
    for (size_t i = 0; i < m; ++i) {
        size_t start = i * segment_length;
        size_t length = std::min(segment_length, reference_sequence.length() - start);
        ref_segments.push_back(reference_sequence.substr(start, length));
    }

    // Divide target sequence into segments
    for (size_t i = 0; i < n; ++i) {
        size_t start = i * segment_length;
        size_t length = std::min(segment_length, target_sequence.length() - start);
        target_segments.push_back(target_sequence.substr(start, length));
    }

    // For i to n
    for (size_t i = 0; i < n; ++i) {
//        std::cout << "Processing segment " << i + 1 << " of " << n << std::endl;
        // Print every 100th segment
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
            // Check if segments are consisted of only 'N's
            if (std::all_of(target_segments[i].begin(), target_segments[i].end(), [](char c) { return c == 'N'; })) {
                std::cout << "Segment " << i + 1 << " is all 'N's." << std::endl;
            } else if (std::all_of(ref_segments[i].begin(), ref_segments[i].end(), [](char c) { return c == 'N'; })) {
                std::cout << "Reference segment " << i + 1 << " is all 'N's." << std::endl;
            } else {
                mismatch_t_ime++;
            }

            // Store target[i] in an intermediate file using a data stream
            std::string intermediate_data = "RAW_TARGET_SEGMENT " + std::to_string(i) + " " + target_segments[i] + "\n";
            if (!saveToFile(intermediate_file, intermediate_data)) {
                throw std::runtime_error("Failed to save intermediate data stream to file: " + intermediate_file);
            }
            continue;
        }

        // If matches found, process them
//        std::cout << "Matches found for segment " << i + 1 << ". Processing matches." << std::endl;

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

        if (mismatch_t_ime > mismatch_threshold_T2) {
            std::cout << "Mismatch threshold exceeded. Triggering global fallback." << std::endl;
            global_fallback = true;
            break; // Exit the loop to perform global fallback
        }
    }

    if (global_fallback) {
        std::cout << "Performing global fallback due to mismatch threshold." << std::endl;

        Utils::removeFileIfExists(intermediate_file, true); // Remove intermediate file if it exists

        // Write the metadata and ORIGINAL target length to the intermediate file first
        // meta_data is from the original target_fasta_path. target_sequence.length() is before N-removal.
        std::ostringstream global_metadata_stream;
        global_metadata_stream << meta_data << "\n" << target_sequence.length() << "\n";

        if (!saveToFile(intermediate_file, global_metadata_stream.str())) {
            throw std::runtime_error("Failed to save metadata to file in global fallback: " + intermediate_file);
        }

        // lowercase_positions_T was already computed from the original target_sequence.
        std::vector<Position> n_fragments_T;
        n_fragments_T = getPositions(target_sequence, [](char c) { return c == 'N'; });

        // Java order: L_list then N_list.
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

        std::cout << "Performing global matching." << std::endl;

        MatchFinder global_match_finder(n_free_reference, n_free_target, kmer_size, search_range_limit);
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
