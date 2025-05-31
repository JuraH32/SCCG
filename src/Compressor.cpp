#include "Compressor.h"
#include "FastaParser.h" // For FastaParser::readSequence
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
    std::ofstream outFile(filename);
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
            // If this line breaks the sequence of pure matches, delta should restart for next pure match.
            // However, Java's `prev_end_coord` persists across unaligned text.
            // For simplicity here, if a non-pure-match line appears, it doesn't reset last_absolute_ref_end.
            // The delta is always from the *last pure match's end*.
            // If a different behavior is needed (e.g., reset delta on non-match), `last_absolute_ref_end` could be reset to -1.
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

void Compressor::savePositionsToFile(const std::string &filename, const std::vector<Position> &positions) {
    std::ostringstream buffer;
    buffer << "LOWERCASE_POSITIONS: ";
    int last_end = 0;
    for (const auto &pos: positions) {
        int start = pos.startInTarget;
        int end = pos.endInTarget;
        buffer << start - last_end << " " << end - start + 1 << " ";
    }

    if (!saveToFile(filename, buffer.str())) {
        std::cerr << "Error: Could not save positions to file: " << filename << std::endl;
    }
}

void Compressor::saveAlignmentSegmentsToFile(
        const std::string &filename,
        const std::vector<AlignmentSegment> &segments
) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
        return;
    }

    for (const auto &segment: segments) {
        if (segment.is_match) {
            outFile << "MATCH "
                    << segment.region.startInReference << " "
                    << segment.region.startInTarget << " "
                    << segment.region.length() << " NOMISMATCH\n";
        } else {
            outFile << "MISMATCH "
                    << segment.region.startInTarget << " "
                    << segment.mismatched_sequence << "\n";
        }
    }

    outFile.close();
}

//void Compressor::compress(const std::string& ref_fasta_path,
//                          const std::string& target_fasta_path,
//                          const std::string& output_path,
//                          int kmer_size,
//                          int k0,
//                          size_t segment_length,
//                          int search_range_limit,
//                          int mismatch_threshold_T1,
//                          int mismatch_threshold_T2) {
//    int controuble = 0;
//    bool is_con = false;
//    bool local = true;
//    int mismatch = 0;
//
//    std::cout << "Loading reference sequence from: " << ref_fasta_path << std::endl;
//    std::string reference_sequence_orig = FastaParser::readSequence(ref_fasta_path);
//    if (reference_sequence_orig.empty()) {
//        throw std::runtime_error("Failed to read reference sequence or file is empty: " + ref_fasta_path);
//    }
//
//    std::cout << "Loading target sequence from: " << target_fasta_path << std::endl;
//    std::string target_sequence_orig = FastaParser::readSequence(target_fasta_path);
//    if (target_sequence_orig.empty()) {
//        throw std::runtime_error("Failed to read target sequence or file is empty: " + target_fasta_path);
//    }
//
//    // Split the target path to get the filename
//    std::vector<std::string> target_path_parts = Utils::splitPath(target_fasta_path);
//    std::string target_filename = target_path_parts.empty() ? "target" : target_path_parts.back();
//    std::string final_filename = output_path + "/" + target_filename;
//
//    while(local) {
//        mismatch = 0;
//
//        Utils::removeFileIfExists(final_filename);
//        std::string reference_sequence = reference_sequence_orig;
//        std::string target_sequence = target_sequence_orig;
//        std::string temp_path = output_path + "/interim.txt";
//
//
//    }
//
//    // Line 1 (part 1): Store lowercase positions from target sequence and convert to uppercase
//    std::vector<size_t> lowercase_positions_T;
//    for (size_t i = 0; i < target_sequence.length(); ++i) {
//        if (std::islower(target_sequence[i])) {
//            lowercase_positions_T.push_back(i);
//            target_sequence[i] = std::toupper(target_sequence[i]);
//        }
//    }
//    for (char& c : reference_sequence) {
//        c = std::toupper(c);
//    }
//
//    std::string intermediate_data_stream;
//    intermediate_data_stream += "LOWERCASE_POSITIONS_T_COUNT " + std::to_string(lowercase_positions_T.size()) + "\n";
//    for (size_t pos : lowercase_positions_T) {
//        intermediate_data_stream += std::to_string(pos) + "\n";
//    }
//
//    // save the intermediate data stream to a file
//    if (!saveToFile(output_intermediate_path, intermediate_data_stream)) {
//        throw std::runtime_error("Failed to save intermediate data stream to file: " + output_intermediate_path);
//    }
//
//    // Line 1 (part 2): Divide R and T into segments
//    size_t m_ref_segments = (reference_sequence.length() + segment_length - 1) / segment_length;
//    size_t n_target_segments = (target_sequence.length() + segment_length - 1) / segment_length;
//
//    std::vector<std::string> ref_segments;
//    for (size_t i = 0; i < m_ref_segments; ++i) {
//        size_t start = i * segment_length;
//        size_t length = std::min(static_cast<size_t>(segment_length), reference_sequence.length() - start);
//        if (length > 0) ref_segments.push_back(reference_sequence.substr(start, length));
//    }
//    if (ref_segments.empty() && !reference_sequence.empty()) ref_segments.push_back(reference_sequence);
//
//
//    std::vector<std::string> target_segments;
//    for (size_t i = 0; i < n_target_segments; ++i) {
//        size_t start = i * segment_length;
//        size_t length = std::min(static_cast<size_t>(segment_length), target_sequence.length() - start);
//        if (length > 0) target_segments.push_back(target_sequence.substr(start, length));
//    }
//    if (target_segments.empty() && !target_sequence.empty()) target_segments.push_back(target_sequence);
//
//    int mismatch_t_ime = 0;
//    bool perform_global_fallback = false;
//
//    // Line 2: for i = 1 to n (loop over target_segments)
//    for (size_t i = 0; i < target_segments.size(); ++i) {
//        std::cout << "Processing target segment " << i + 1 << "/" << target_segments.size() << std::endl;
//        bool current_ti_has_match = false;
//        std::vector<Match> matches_for_ti;
//
//        // Attempt to find a match for target_segments[i] (ti) with any ref_segments[j] (rj)
//        // Line 3: Invoke Algorithm 1 with parameters (rj, ti, k, 0, false)
//        for (size_t j = 0; j < ref_segments.size(); ++j) {
//            if (ref_segments[j].length() >= static_cast<size_t>(kmer_size) && target_segments[i].length() >= static_cast<size_t>(kmer_size)) {
//                MatchFinder finder(ref_segments[j], target_segments[i], kmer_size, 0 /*search_range*/);
//                std::vector<Match> found_matches = finder.findMatches(false /*global*/);
//                for (const auto& match : found_matches) {
//                    if (match.length > 0) {
//                        current_ti_has_match = true;
//                        matches_for_ti = found_matches;
//                        goto found_match_for_ti_k; // Break from both loops for rj and go to process this match
//                    }
//                }
//            }
//        }
//        found_match_for_ti_k:;
//
//        // Line 4: if (no matched string found for (ri, ti))
//        if (!current_ti_has_match) {
//            // Line 5: Invoke Algorithm 1 with parameters (rj, ti, k0, 0, false)
//            for (size_t j = 0; j < ref_segments.size(); ++j) {
//                if (ref_segments[j].length() >= static_cast<size_t>(k0) && target_segments[i].length() >= static_cast<size_t>(k0)) {
//                    MatchFinder finder(ref_segments[j], target_segments[i], k0, 0 /*search_range*/);
//                    std::vector<Match> found_matches = finder.findMatches(false /*global*/);
//                    for (const auto& match : found_matches) {
//                        if (match.length > 0) {
//                            current_ti_has_match = true;
//                            matches_for_ti = found_matches;
//                            goto found_match_for_ti_k0; // Break from both loops for rj
//                        }
//                    }
//                }
//            }
//            found_match_for_ti_k0:;
//        }
//
//        // Line 7: if (still no matched string found for (ri, ti))
//        if (!current_ti_has_match) {
//            // Line 8: if (neither ri nor ti consists of only ‘N’ characters) - Simplified to check ti
//            bool ti_is_all_N = std::all_of(target_segments[i].begin(), target_segments[i].end(), [](char c){ return c == 'N'; });
//            if (!ti_is_all_N) {
//                // Line 9: mismatchtime ++
//                mismatch_t_ime++;
//            }
//            // Line 11: Store ti in the intermediate file
//            intermediate_data_stream += "RAW_TARGET_SEGMENT " + std::to_string(i) + " " + target_segments[i] + "\n";
//        } else { // Matched (current_ti_has_match is true)
//            // Line 13: Store (pn; ln)pairs and unmatched strings in ti
//            intermediate_data_stream += "MATCHED_TARGET_SEGMENT " + std::to_string(i) + "\n";
//            size_t total_mismatched_chars_in_ti = 0; // Will be calculated based on match coverage
//
//            size_t covered_by_matches = 0;
//            for (const auto& match : matches_for_ti) {
//                intermediate_data_stream += "MATCH " + std::to_string(match.reference_pos) + " " +
//                                            std::to_string(match.length) + " " +
//                                            (match.mismatch.empty() ? "NOMISMATCH" : match.mismatch) + "\n";
//                covered_by_matches += match.length;
//            }
//
//            // Line 14: if (the percentage of mismatching characters in ti > T1)
//            size_t total_len_ti = target_segments[i].length();
//            if (total_len_ti > 0) {
//                total_mismatched_chars_in_ti = total_len_ti - covered_by_matches;
//                double mismatch_percentage = static_cast<double>(total_mismatched_chars_in_ti) * 100.0 / total_len_ti;
//                if (mismatch_percentage > static_cast<double>(mismatch_threshold_T1)) {
//                    // Line 15: mismatchtime ++
//                    mismatch_t_ime++;
//                }
//            }
//        }
//
//        // Line 16 (second part in pseudocode, check after processing ti): if (mismatchtime > T2)
//        if (mismatch_t_ime > mismatch_threshold_T2) {
//            // Line 17: Go to line 21
//            perform_global_fallback = true;
//            break; // Break from the loop over i (target_segments)
//        }
//    } // Line 19: end for (loop over target_segments)
//
//    // After the loop (either completed or broken by T2 threshold)
//    if (perform_global_fallback) { // Came from line 17
//        std::cout << "Global fallback triggered (mismatch_t_ime=" << mismatch_t_ime << " > T2=" << mismatch_threshold_T2 << ")" << std::endl;
//        // Line 21: Clear the intermediate file (stream)
//        intermediate_data_stream.clear();
//
//        // Line 22: Convert all lowercase characters into uppercase, delete ‘N’ characters.
//        // Store info about ‘N’ fragments and lowercase characters in T into the intermediate file.
//
//        // Re-add original lowercase positions for T
//        intermediate_data_stream += "LOWERCASE_POSITIONS_T_COUNT " + std::to_string(lowercase_positions_T.size()) + "\n";
//        for (size_t pos : lowercase_positions_T) {
//            intermediate_data_stream += std::to_string(pos) + "\n";
//        }
//
//        // Use the already uppercased full sequences (reference_sequence, target_sequence)
//        std::string processed_R = reference_sequence;
//        std::string processed_T = target_sequence;
//
//        // Delete 'N' characters from T and store 'N' fragment info
//        std::string n_free_T;
//        n_free_T.reserve(processed_T.length());
//        std::vector<std::pair<size_t, size_t>> n_fragments_T; // pos, len
//        size_t current_n_start_T = std::string::npos;
//        for (size_t char_idx = 0; char_idx < processed_T.length(); ++char_idx) {
//            if (processed_T[char_idx] == 'N') {
//                if (current_n_start_T == std::string::npos) current_n_start_T = char_idx;
//            } else {
//                if (current_n_start_T != std::string::npos) {
//                    n_fragments_T.push_back({current_n_start_T, char_idx - current_n_start_T});
//                    current_n_start_T = std::string::npos;
//                }
//                n_free_T += processed_T[char_idx];
//            }
//        }
//        if (current_n_start_T != std::string::npos) {
//            n_fragments_T.push_back({current_n_start_T, processed_T.length() - current_n_start_T});
//        }
//        intermediate_data_stream += "N_FRAGMENTS_T_COUNT " + std::to_string(n_fragments_T.size()) + "\n";
//        for (const auto& frag : n_fragments_T) {
//            intermediate_data_stream += std::to_string(frag.first) + " " + std::to_string(frag.second) + "\n";
//        }
//        processed_T = n_free_T;
//
//        // Delete 'N' characters from R (no need to store N fragment info for R per pseudocode)
//        std::string n_free_R;
//        n_free_R.reserve(processed_R.length());
//        for(char c : processed_R) if(c != 'N') n_free_R += c;
//        processed_R = n_free_R;
//
//        if (processed_R.empty() || processed_T.empty() ||
//            processed_R.length() < static_cast<size_t>(kmer_size) || processed_T.length() < static_cast<size_t>(kmer_size)) {
//            std::cout << "Warning: Processed R or T is empty or too short for global matching after N removal. Storing T raw." << std::endl;
//            intermediate_data_stream += "RAW_GLOBAL_TARGET " + processed_T + "\n";
//        } else {
//            // Line 23: Invoke algorithm 1 with (processed R, processed T, k, m, true)
//            // k is kmer_size, m is search_range_limit from function args
//            std::cout << "Invoking global MatchFinder on N-free sequences (R_len=" << processed_R.length() << ", T_len=" << processed_T.length() << ")" << std::endl;
//            MatchFinder global_finder(processed_R, processed_T, kmer_size, search_range_limit);
//            std::vector<Match> global_matches = global_finder.findMatches(true /*global=true*/);
//
//            // Line 24: Store (pn; ln) pairs and unmatched sub-strings in processed T
//            intermediate_data_stream += "GLOBAL_MATCHES\n";
//            for (const auto& match : global_matches) {
//                 intermediate_data_stream += "MATCH " + std::to_string(match.reference_pos) + " " +
//                                            std::to_string(match.length) + " " +
//                                            (match.mismatch.empty() ? "NOMISMATCH" : match.mismatch) + "\n";
//            }
//        }
//    } else { // Line 20: Go to line 25 (segment-wise processing completed without fallback)
//        std::cout << "Segment-wise processing complete. No global fallback (mismatch_t_ime=" << mismatch_t_ime << ")" << std::endl;
//    }
//
//    // Line 25: Modify pn in (pn; ln) pairs in the intermediate file with delta coding;
//    // TODO: Implement merging of continuous matches before delta coding, similar to Java's postprocess().
//    std::cout << "Placeholder: Apply delta coding to pn in intermediate data." << std::endl;
//    intermediate_data_stream += "DELTA_CODING_TO_BE_APPLIED\n";
//
//    // Line 26: Compress the intermediate file using PPMd algorithm;
//    // Here we just save the intermediate data. PPMd compression is a separate step.
//    std::cout << "Saving intermediate data (intended for PPMd compression) to: " << output_intermediate_path << std::endl;
//    if (!saveToFile(output_intermediate_path, intermediate_data_stream)) { // Use non-static method
//        throw std::runtime_error("Failed to save intermediate data to " + output_intermediate_path);
//    }
//    std::cout << "SCCG algorithm processing part complete. Output: " << output_intermediate_path << std::endl;
//}

void Compressor::compress(const std::string &ref_fasta_path,
                          const std::string &target_fasta_path,
                          const std::string &output_path,
                          int kmer_size,
                          int k0,
                          size_t segment_length,
                          int search_range_limit,
                          int mismatch_threshold_T1,
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
    Utils::removeFileIfExists(final_filename);
    Utils::removeFileIfExists(intermediate_file);

    if (target_sequence.empty()) {
        throw std::runtime_error("Failed to read target sequence or file is empty: " + target_fasta_path);
    }

    // Line 1 (part 1): Store lowercase positions from target sequence and convert to uppercase
    std::vector<Position> lowercase_positions_T;
    lowercase_positions_T = getPositions(target_sequence, [](char c) { return std::islower(c); });
    savePositionsToFile(intermediate_file, lowercase_positions_T);

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
        std::cout << "Processing segment " << i + 1 << " of " << n << std::endl;
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
        std::cout << "Matches found for segment " << i + 1 << ". Processing matches." << std::endl;

        saveAlignmentSegmentsToFile(intermediate_file, matches);

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

        Utils::removeFileIfExists(intermediate_file);

        // Convert all lowercase characters into uppercase, delete 'N' characters
        std::vector<Position> n_fragments_T;
        n_fragments_T = getPositions(target_sequence, [](char c) { return c == 'N'; });
        savePositionsToFile(intermediate_file, n_fragments_T);
        savePositionsToFile(intermediate_file, lowercase_positions_T);

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

        MatchFinder global_match_finder(n_free_reference, n_free_target, kmer_size, search_range_limit);
        std::vector<AlignmentSegment> global_matches = global_match_finder.findMatches(true);

        if (global_matches.empty()) {
            std::cout << "No global matches found. Storing raw target sequence." << std::endl;
            std::string raw_target_data = "RAW_GLOBAL_TARGET " + n_free_target + "\n";
            if (!saveToFile(intermediate_file, raw_target_data)) {
                throw std::runtime_error("Failed to save raw target sequence to file: " + intermediate_file);
            }
        } else {
            saveAlignmentSegmentsToFile(intermediate_file, global_matches);
        }
    }

    // Postprocess the intermediate file to apply delta coding and save final output
    std::cout << "Postprocessing intermediate file: " << intermediate_file << std::endl;
    postprocess(intermediate_file, final_filename);
    std::cout << "Compression complete. Output saved to: " << final_filename << std::endl;

    Utils::compressWith7zip(final_filename);

    std::cout << "Final output compressed with 7zip." << std::endl;

    // Clean up intermediate file
    Utils::removeFileIfExists(intermediate_file);
}