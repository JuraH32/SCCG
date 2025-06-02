#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <list>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <limits>
#include <cctype>
#include <cmath>
#include <sys/stat.h>
#include <cstdio>
#include <cstdlib>
#include <sys/types.h>

class newkmer {
private:
    int kmerstart;
    std::string kmer;

public:
    void setkmerstart(int start) { kmerstart = start; }

    int getkmerstart() const { return kmerstart; }

    void setkmer(const std::string &k) { kmer = k; }

    std::string getkmer() const { return kmer; }
};

class Position {
private:
    int startinTar, endinTar, startinRef, endinRef;

public:
    void setstartinTar(int start) { startinTar = start; }

    int getstartinTar() const { return startinTar; }

    void setendinTar(int end) { endinTar = end; }

    int getendinTar() const { return endinTar; }

    void setstartinRef(int start) { startinRef = start; }

    int getstartinRef() const { return startinRef; }

    void setendinRef(int end) { endinRef = end; }

    int getendinRef() const { return endinRef; }
};

enum ScaleType {
    LOCAL = 0,
    GLOBAL = 1
};

enum SequenceType {
    REFERENCE = 0,
    TARGET = 1
};

class SequenceProcessor {
public:
    static std::unordered_map<int, std::vector<newkmer> > hashmap;
    static std::vector<int> next_kmer;
    static std::vector<int> kmer_location;
    static const int maxchar = 268435456; // 2^28
    static const int maxseq = 268435456 * 2; // 2^29

    static int stringHashCode(const std::string &str) {
        int hash = 0;
        for (char c: str) {
            hash = 31 * hash + static_cast<int>(c);
        }
        return hash;
    }

    static void buildLhashtable(const std::string &read,
                                int kmer_length_param) {
        int current_length = read.length(), i = kmer_length_param;
        std::string Nkmer = "";
        while (i > 0) {
            Nkmer += "N";
            i--;
        }

        int Nkey = stringHashCode(Nkmer);
        i = 0;

        int count = 0;

        while (i < current_length - kmer_length_param + 1) {

            std::string kmer_str = read.substr(i, kmer_length_param);
            newkmer newKMer;
            newKMer.setkmerstart(i);
            newKMer.setkmer(kmer_str);
            int key = stringHashCode(kmer_str);
            if (hashmap.find(key) != hashmap.end()) {
                std::vector<newkmer> &list = hashmap[key];
                list.push_back(newKMer);
            } else {
                std::vector<newkmer> list;
                list.push_back(newKMer);
                hashmap[key] = list;
            }
            i++;
            if (key == Nkey) {
                while (i < current_length - kmer_length_param + 1 &&
                       (read.substr(i, 1) == "n" || read.substr(i, 1) == "N")) {
                    i++;
                }
            }
        }
    }

    static void buildGhashtable(const std::string &read, int kmer_length_param) {
        int iteration = read.length() - kmer_length_param + 1;

        for (int i = 0; i < maxseq; i++) {
            kmer_location[i] = -1;
        }
        for (int i = 0; i < iteration; i++) {

            std::string kmer_str = read.substr(i, kmer_length_param);
            long key = std::abs(static_cast<long>(stringHashCode(kmer_str)));

            if (key == -2147483648LL) {
                key = 0;
            }

            while (key > maxseq - 1) {
                key = key / 2;
            }
            next_kmer[i] = kmer_location[static_cast<int>(key)];
            kmer_location[static_cast<int>(key)] = i;
        }
    }

    static std::vector<Position> lowercase_position(const std::string &sequence) {
        std::vector<Position> list;
        bool successive = false;
        int start = 0, current_end = 0;

        for (int i = 0; i < sequence.length(); i++) {
            if (std::islower(sequence[i])) {
                if (successive) {
                    current_end += 1;
                } else {
                    start = i;
                    current_end += 1;
                    successive = true;
                }
            } else {
                if (successive) {
                    Position position;
                    position.setstartinTar(start);
                    position.setendinTar(current_end - 1);
                    list.push_back(position);
                }
                successive = false;
                start = 0;
                current_end = i + 1;
            }
        }

        if (successive) {
            Position position;
            position.setstartinTar(start);
            position.setendinTar(current_end);
            list.push_back(position);
        }
        return list;
    }
};

std::unordered_map<int, std::vector<newkmer> > SequenceProcessor::hashmap;
std::vector<int> SequenceProcessor::next_kmer;
std::vector<int> SequenceProcessor::kmer_location;

class Matcher {
public:
    static std::vector<Position>
    Lmatch(const std::string &ref_seq, const std::string &tar_seq, int kmer_len) {
        std::vector<Position> list;
        int index = 0, startIndex;
        int increment, most_incre, key_val;
        int kmerstart_val, endinRef_val, endinTar_val, Refendindex, Tarendindex;
        std::string kmer_str;

        SequenceProcessor::buildLhashtable(ref_seq, kmer_len);

        while (true) {
            increment = 0;
            most_incre = 0;
            if (index + kmer_len > tar_seq.length()) {
                break;
            }
            kmer_str = tar_seq.substr(index, kmer_len);
            key_val = SequenceProcessor::stringHashCode(kmer_str);

            if (SequenceProcessor::hashmap.find(key_val) == SequenceProcessor::hashmap.end()) {
                startIndex = std::numeric_limits<int>::max();
                index = index + 1;
                if (index + kmer_len > tar_seq.length()) {
                    break;
                }
                continue;
            }

            std::vector<newkmer> &klist = SequenceProcessor::hashmap[key_val];
            startIndex = std::numeric_limits<int>::max();
            most_incre = 0;
            for (const newkmer &newKMer: klist) {
                if (newKMer.getkmer() == kmer_str) {
                    kmerstart_val = newKMer.getkmerstart();
                    endinRef_val = kmerstart_val + kmer_len - 1;
                    endinTar_val = index + kmer_len - 1;
                    Refendindex = ref_seq.length() - 1;
                    Tarendindex = tar_seq.length() - 1;
                    increment = get_incre(ref_seq, tar_seq, endinRef_val, endinTar_val, Refendindex,
                                          Tarendindex);

                    if (klist.size() > 1) {
                        if (increment == most_incre) {
                            if (list.size() >
                                1) { // Should be list.size() > 0 or !list.empty()
                                int lastEIR = list.back().getendinRef();
                                if (std::abs(kmerstart_val - lastEIR) <
                                    std::abs(startIndex - lastEIR))
                                    startIndex = kmerstart_val;
                            } else if (list.empty()) {
                                startIndex = kmerstart_val;
                            }
                        } else if (increment > most_incre) {
                            most_incre = increment;
                            startIndex = kmerstart_val;
                        }
                    } else {
                        most_incre = increment;
                        startIndex = kmerstart_val;
                        break;
                    }
                }
            }

            if (startIndex == std::numeric_limits<int>::max()) {
                index = index + 1;
                if (index + kmer_len > tar_seq.length()) {
                    break;
                }
                continue;
            }

            Position position_obj;
            position_obj.setstartinTar(index);
            position_obj.setendinTar(index + kmer_len + most_incre - 1);
            position_obj.setstartinRef(startIndex);
            position_obj.setendinRef(startIndex + kmer_len + most_incre - 1);
            list.push_back(position_obj);
            index = index + kmer_len + most_incre +
                    1;
            if (index + kmer_len > tar_seq.length()) {
                break;
            }
        }
        return list;
    }

    static std::vector<Position>
    Gmatch(const std::string &ref_seq, const std::string &tar_seq, int kmer_len, int limit_param) {
        std::vector<Position> list;
        int index = 0, startIndex, lastEIR = 0;
        std::string kmer_str;

        SequenceProcessor::buildGhashtable(ref_seq, kmer_len);

        while (true) {
            int increment = 0, most_incre = 0;
            if (index + kmer_len > tar_seq.length()) {
                break;
            }
            kmer_str = tar_seq.substr(index, kmer_len);
            int key_val = std::abs(SequenceProcessor::stringHashCode(kmer_str));

            if (key_val == -2147483648) { // Integer.MIN_VALUE
                key_val = 0;
            }
            while (key_val > SequenceProcessor::maxseq - 1) {
                key_val = key_val / 2;
            }

            if (SequenceProcessor::kmer_location[key_val] == -1) {
                startIndex = std::numeric_limits<int>::max();
                index = index + 1;
                if (index + kmer_len > tar_seq.length()) {
                    break;
                }
                continue;
            }

            startIndex = std::numeric_limits<int>::max();
            most_incre = 0;
            bool match_found = false;

            for (int k = SequenceProcessor::kmer_location[key_val]; k != -1; k = SequenceProcessor::next_kmer[k]) {
                increment = 0;
                if (k + kmer_len > ref_seq.length()) continue;
                std::string Rkmer = ref_seq.substr(k, kmer_len);
                if (kmer_str != Rkmer) {
                    continue;
                }
                try {
                    if (!list.empty()) {
                        lastEIR = list.back().getendinRef();
                    } else {
                        lastEIR = 0;
                    }
                } catch (...) {
                    lastEIR = 0;
                }
                if (std::abs(k - lastEIR) > limit_param) { // Use limit_param
                    continue;
                }

                match_found = true;
                int ref_idx = k + kmer_len;
                int tar_idx = index + kmer_len;
                while (ref_idx < ref_seq.length() && tar_idx < tar_seq.length()
                       && (ref_seq.substr(ref_idx, 1) == tar_seq.substr(tar_idx, 1))) {
                    ref_idx++;
                    tar_idx++;
                    increment++;
                }
                if (increment == most_incre) {
                    if (!list.empty()) {
                        if (std::abs(k - lastEIR) < std::abs(startIndex - lastEIR))
                            startIndex = k;
                    } else if (list.empty() && startIndex == std::numeric_limits<int>::max()) {
                        startIndex = k;
                    }
                } else if (increment > most_incre) {
                    most_incre = increment;
                    startIndex = k;
                }
            }
            if (!match_found) {
                for (int k = SequenceProcessor::kmer_location[key_val]; k != -1; k = SequenceProcessor::next_kmer[k]) {
                    increment = 0;
                    if (k + kmer_len > ref_seq.length()) continue;
                    std::string Rkmer = ref_seq.substr(k, kmer_len);
                    if (kmer_str != Rkmer) {
                        continue;
                    }

                    int ref_idx = k + kmer_len;
                    int tar_idx = index + kmer_len;

                    while (ref_idx < ref_seq.length() && tar_idx < tar_seq.length()
                           && (ref_seq.substr(ref_idx, 1) == tar_seq.substr(tar_idx, 1))) {
                        ref_idx++;
                        tar_idx++;
                        increment++;
                    }
                    if (increment == most_incre) {
                        if (!list.empty()) {
                            if (std::abs(k - lastEIR) < std::abs(startIndex - lastEIR))
                                startIndex = k;
                        } else if (list.empty() && startIndex == std::numeric_limits<int>::max()) {
                            startIndex = k;
                        }
                    } else if (increment > most_incre) {
                        most_incre = increment;
                        startIndex = k;
                    }
                }
            }

            if (startIndex == std::numeric_limits<int>::max()) {
                index = index + 1;
                if (index + kmer_len > tar_seq.length()) {
                    break;
                }
                continue;
            }

            Position position_obj;
            position_obj.setstartinTar(index);
            position_obj.setendinTar(index + kmer_len + most_incre - 1);
            position_obj.setstartinRef(startIndex);
            position_obj.setendinRef(startIndex + kmer_len + most_incre - 1);
            list.push_back(position_obj);
            index = index + kmer_len + most_incre + 1;
            if (index + kmer_len > tar_seq.length()) {
                break;
            }
        }
        return list;
    }

private:
    static int get_incre(const std::string& reference_seq, const std::string& target_seq,
                         int endinRef_param, int endinTar_param, int Refendindex, int Tarendindex) {
        int position = 0;
        int endIndex;

        if (Refendindex - endinRef_param <= Tarendindex - endinTar_param) {
            endIndex = Refendindex - endinRef_param + 1;
        } else {
            endIndex = Tarendindex - endinTar_param + 1;
        }

        for (int i = 1; i < endIndex; i++) {
            if ((endinTar_param + i) >= target_seq.length() || (endinRef_param + i) >= reference_seq.length())
                std::cerr << "OUT OF BOUNDS: endinTar_param + i = " << (endinTar_param + i)
                          << " target_seq.length() = " << target_seq.length()
                          << ", endinRef_param + i = " << (endinRef_param + i)
                          << " reference_seq.length() = " << reference_seq.length() << std::endl;
            break;
        }
        return position;
    }
};

class FileUtils {
public:
    static std::string meta_data;
    static int length;
    static std::string readSeq(const std::string &sequenceFileName, ScaleType scaleType = ScaleType::LOCAL, SequenceType sequenceType = SequenceType::TARGET) {
        std::ifstream file(sequenceFileName);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + sequenceFileName);
        }

        std::stringstream stringbuilder;
        std::string line;
        int line_length;
        char temp_ch;
        bool base = false;

        std::getline(file, line);

        // Meta data is important only from target sequence when scale is local
        if (sequenceType == SequenceType::TARGET && scaleType == ScaleType::LOCAL) {
            meta_data = line; // Store the first line as metadata
        }

        while (std::getline(file, line)) {
            line_length = line.length();
            if (sequenceType == SequenceType::REFERENCE && scaleType == ScaleType::GLOBAL) {
                for (int i = 0; i < line_length; i++) {
                    temp_ch = std::toupper(line[i]);
                    if (temp_ch != 'N') {
                        stringbuilder << temp_ch;
                    }
                }
            } else {
                stringbuilder << line;
            }
            if (!base && scaleType == ScaleType::LOCAL) {
                length = line_length;
                base = true;
            }
        }
        file.close();
        return stringbuilder.str();
    }

    static std::string GreadtarSeq(const std::string &sequenceFileName, const std::string &fileName) {
        std::ifstream file(sequenceFileName);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + sequenceFileName);
        }

        std::ofstream out(fileName, std::ios::app);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot open output file: " + fileName);
        }

        std::stringstream Llist;
        std::stringstream Nlist;

        std::string meta_data_local, line;

        int line_length, totallength = 0, Llen = 0, Nlen = 0, last_L = 0, last_N = 0, start_L = 0,
                start_N = 0;
        char temp_ch;
        bool is_low = false;
        bool is_N = false;
        bool linelen = true;

        std::stringstream stringbuilder;

        std::getline(file, meta_data_local);
        out << meta_data_local << '\n';

        while (std::getline(file, line)) {
            line_length = line.length();

            if (linelen) {
                out << std::to_string(line_length) << '\n';
                linelen = false;
            }

            for (int i = 0; i < line_length; i++) {
                temp_ch = line[i];

                if (std::islower(temp_ch)) {
                    if (is_low) {
                        Llen++;
                    } else {
                        is_low = true;
                        start_L = i + totallength;
                        Llist << (start_L - last_L) << ' ';
                        Llen++;
                    }
                    temp_ch = std::toupper(temp_ch);
                } else {
                    if (is_low) {
                        Llen--;
                        Llist << Llen << ' ';
                        last_L = start_L + Llen;
                    }
                    is_low = false;
                    Llen = 0;
                }

                if (temp_ch == 'N') {
                    if (is_N) {
                        Nlen++;
                    } else {
                        is_N = true;
                        start_N = i + totallength;
                        Nlist << (start_N - last_N) << ' ';
                        Nlen++;
                    }
                } else {
                    if (is_N) {
                        Nlist << Nlen << ' ';
                        last_N = start_N + Nlen;
                    }
                    is_N = false;
                    Nlen = 0;

                    stringbuilder << temp_ch;
                }
            }
            totallength += line_length;
        }

        std::string LlistStr = Llist.str();
        std::string NlistStr = Nlist.str();

        if (LlistStr.length() > 0) {
            if (is_low) {
                LlistStr += std::to_string(Llen);
            } else {
                LlistStr.pop_back();
            }
        }

        if (NlistStr.length() > 0) {
            if (is_N) {
                NlistStr += std::to_string(Nlen);
            } else {
                NlistStr.pop_back();
            }
        }

        file.close();

        out << LlistStr << '\n';
        out << NlistStr << '\n';
        out.flush();
        out.close();
        return stringbuilder.str();
    }

    static void
    write(const std::string &filename, const std::string &text_content, bool append) {
        std::ofstream output;
        if (append) {
            output.open(filename, std::ios::app);
        } else {
            output.open(filename);
        }
        if (!output.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        output << text_content;
        output.flush();
        output.close();
    }

    static void
    write(std::string filename, std::vector<Position> list_param, bool append, std::string auxiliary) {
        std::ofstream output;
        if (append) {
            output.open(filename, std::ios::app);
        } else {
            output.open(filename);
        }

        if (!output.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        std::ostringstream text_stream;
        int temp_end = 0;
        output << auxiliary;

        for (Position position_obj: list_param) {
            int start_val = position_obj.getstartinTar();
            int end_val = position_obj.getendinTar();
            text_stream << (start_val - temp_end) << " " << (end_val - start_val) << " ";
            temp_end = end_val;
        }
        output << text_stream.str();
        output << "\n";
        output.flush();
        output.close();
    }

    static void use7zip(std::string filename) {
        struct stat buffer;
        if (stat(filename.c_str(), &buffer) != 0) {
            std::cerr << "Warning: File " << filename << " not found for 7zip." << std::endl;
            return;
        }

        std::string exec = "./7za a " + filename + ".7z " + filename + " -m0=PPMD";
        try {
            int result = std::system(exec.c_str());
            if (result != 0) {
                std::cerr << "Warning: 7zip command failed with result " << result << std::endl;
            }
        } catch (const std::exception &e) {
            std::cerr << "Exception in use7zip: " << e.what() << std::endl;
        }
    }
};

class ORGC {
public:
    static std::string reference, target;
    static std::string text;
    static int kmer_length;
    static int sub_length, limit;
    static double T1;
    static int T2;
    static int sot, eot, sor, eor;
    static int mismatch, endref;
    static bool local;

    // Define chrName array here
    static const std::string chrName[];


    ORGC() {
        text = "";
        SequenceProcessor::hashmap = std::unordered_map<int, std::vector<newkmer> >();
    }

    static long getCPUTime() {
        auto now = std::chrono::high_resolution_clock::now();
        auto duration = now.time_since_epoch();
        auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);
        return nanoseconds.count();
    }

    static Position format_matches(const std::vector<Position> &list_param) {
        int startinTar_val, startinRef_val, endinRef_val;
        int trouble = 0;
        if (list_param.empty()) {
            return Position();
        }

        for (int i = 0; i < list_param.size(); i++) {

            if (i == 0) {
                startinTar_val = list_param[i].getstartinTar();
                startinRef_val = list_param[i].getstartinRef();
                endinRef_val = list_param[i].getendinRef();
                if (endinRef_val >= endref) {
                    endref = endinRef_val;
                }

                if (startinTar_val > 0) {
                    std::string preamble = target.substr(0, startinTar_val);
                    text += preamble + "\n";
                    trouble += preamble.length();
                }
                text += "" + std::to_string(startinRef_val + sor) + "," + std::to_string(endinRef_val + sor) +
                        "\n";
                continue;
            }

            startinTar_val = list_param[i].getstartinTar();
            startinRef_val = list_param[i].getstartinRef();
            endinRef_val = list_param[i].getendinRef();
            if (endinRef_val >= endref) {
                endref = endinRef_val;
            }

            int endinTarPrev = list_param[i - 1].getendinTar();
            // Bounds check for substr
            if (endinTarPrev + 1 < startinTar_val) {
                std::string mismatch_str = target.substr(endinTarPrev + 1, startinTar_val - (endinTarPrev +
                                                                                             1));
                if (mismatch_str.length() > 0) {
                    text += mismatch_str + "\n";
                    trouble += mismatch_str.length();
                }
            }

            text += std::to_string(startinRef_val + sor) + "," + std::to_string(endinRef_val + sor) +
                    "\n";
        }
        if (trouble > (sub_length * T1)) // ORGC::sub_length, ORGC::T1
        {
            mismatch++;
        }

        return list_param.back();
    }

    static void format_matches(const std::vector<Position> &list_param, const std::string &fileName) {
        std::stringstream stringbuilder;

        int startinTar_val, startinRef_val, endinRef_val, endinTar_val = 0;
        for (int i = 0; i < list_param.size(); i++) {
            if (i == 0) {
                startinTar_val = list_param[i].getstartinTar();
                endinTar_val = list_param[i].getendinTar();
                startinRef_val = list_param[i].getstartinRef();
                endinRef_val = list_param[i].getendinRef();
                if (startinTar_val > 0) {
                    std::string preamble = target.substr(0, startinTar_val);
                    stringbuilder << preamble << "\n";
                }
                stringbuilder << startinRef_val << "," << endinRef_val << "\n";
                continue;
            }

            startinTar_val = list_param[i].getstartinTar();
            startinRef_val = list_param[i].getstartinRef();
            endinRef_val = list_param[i].getendinRef();
            endinTar_val = list_param[i].getendinTar();
            int endinTarPrev = list_param[i - 1].getendinTar();
            if (endinTarPrev + 1 < startinTar_val) {
                std::string mismatch_str = target.substr(endinTarPrev + 1, startinTar_val - (endinTarPrev +
                                                                                             1));
                if (mismatch_str.length() > 0) {
                    stringbuilder << mismatch_str << "\n";
                }
            } else if (endinTarPrev + 1 == startinTar_val) {
                // No mismatch
            } else {
                // Overlapping
            }
            stringbuilder << startinRef_val << "," << endinRef_val << "\n";
        }
        if (endinTar_val < (target.length() - 1)) {
            if (endinTar_val + 1 < target.length()) {
                stringbuilder
                        << target.substr(endinTar_val + 1, (target.length()) - (endinTar_val + 1));
            }
        }
        FileUtils::write(fileName, stringbuilder.str(), true);
    }

    static void postprocess(std::string filename, std::string final_file) {
        std::ifstream input(filename);
        if (!input.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        std::string line;
        std::ostringstream stringbuilder;
        std::vector<int> num_list;

        while (std::getline(input, line)) {
            if (line.find(",") != std::string::npos) {
                size_t comma_pos = line.find(",");
                int begin = std::stoi(line.substr(0, comma_pos));
                int end_val = std::stoi(line.substr(comma_pos + 1));

                if (num_list.size() > 0) {
                    int prevEnd = num_list.back();
                    if (begin != prevEnd + 1) {
                        stringbuilder << num_list[0] << "," << num_list.back() << "\n";
                        num_list.clear();
                    }
                }
                num_list.push_back(begin);
                num_list.push_back(end_val);
            } else if (line.length() > 0) {

                if (!num_list.empty()) {
                    stringbuilder << num_list[0] << "," << num_list.back() << "\n";
                }

                if (line.find("^") == std::string::npos) {
                    stringbuilder << line << "\n";
                }
                num_list.clear();
            }
        }
        if (num_list.size() > 0) {
            stringbuilder << num_list[0] << "," << num_list.back() << "\n";
        }
        input.close();

        std::istringstream inputStream(stringbuilder.str());
        stringbuilder.str("");
        stringbuilder.clear();
        num_list.clear();
        std::vector<std::string> stringList;

        while (std::getline(inputStream, line)) {
            stringList.push_back(line);
        }

        int prev = 0;
        bool successive = false;
        for (size_t i = 0; i < stringList.size(); i++) {
            std::string str = stringList[i];

            if (str.find(",") != std::string::npos) {
                size_t comma_pos = str.find(",");
                int begin = std::stoi(str.substr(0, comma_pos));
                int end_val = std::stoi(str.substr(comma_pos + 1));
                if (!successive) {
                    num_list.push_back(begin);
                    num_list.push_back(end_val - begin);
                    prev = end_val;

                    if (!num_list.empty()) stringbuilder << num_list[0] << "," << num_list.back() << "\n";

                    successive = true;
                } else {
                    num_list.push_back(begin - prev);
                    num_list.push_back(end_val - begin);
                    prev = end_val;
                    if (!num_list.empty()) stringbuilder << num_list[0] << "," << num_list.back() << "\n";
                }
                num_list.clear();
            } else if (str.length() > 0) {
                stringbuilder << str << "\n";
            }
        }

        if (!num_list.empty()) {
            stringbuilder << num_list[0] << "," << num_list.back() << "\n";
        }

        FileUtils::write(final_file, stringbuilder.str(), true);
    }

private:
    static void process_local(const std::string& ref_genome_path, const std::string& tar_genome_path, const std::string& final_file_path, const std::string& tempfile) {
        int controuble = 0;
        bool is_con = false;

        std::cout << tar_genome_path << " is compressing (local attempt)...\n" << std::endl;

        struct stat buffer;
        if (stat(final_file_path.c_str(), &buffer) == 0) {
            if (std::remove(final_file_path.c_str()) != 0) {
                std::cerr << "Warning: Failed to remove existing final file: " << final_file_path << std::endl;
            }
        }

        std::string reference_seq_content = FileUtils::readSeq(ref_genome_path, LOCAL, REFERENCE);
        std::string target_seq_content = FileUtils::readSeq(tar_genome_path, LOCAL, TARGET);
        if (reference_seq_content.empty()) {
            std::cerr << "Error: Reference sequence is empty for " << ref_genome_path << std::endl;
            ORGC::local = false; // Mark local as failed
            return;
        }
        if (target_seq_content.empty()) {
            std::cerr << "Error: Target sequence is empty for " << tar_genome_path << std::endl;
            ORGC::local = false; // Mark local as failed
            return;
        }

        int target_length_val = target_seq_content.length();

        if (target_length_val < ORGC::sub_length * 5) {
            ORGC::local = false;
            return; // Exit if local processing is not suitable
        } else if (target_length_val < ORGC::sub_length * 1333) {
            ORGC::T1 = 0.1;
            ORGC::T2 = 0;
        }

        std::string auxiliary = FileUtils::meta_data + "\n" + std::to_string(FileUtils::length) + "\n";
        std::vector<Position> L_list = SequenceProcessor::lowercase_position(target_seq_content);
        FileUtils::write(final_file_path, L_list, false, auxiliary);
        FileUtils::write(final_file_path, "\n", true);

        std::transform(reference_seq_content.begin(), reference_seq_content.end(), reference_seq_content.begin(),
                       ::toupper);
        std::transform(target_seq_content.begin(), target_seq_content.end(), target_seq_content.begin(), ::toupper);

        if (stat(tempfile.c_str(), &buffer) == 0) {
            if (std::remove(tempfile.c_str()) != 0) {
                 std::cerr << "Warning: Failed to remove existing temp file: " << tempfile << std::endl;
            }
        }

        ORGC::sot = 0;
        ORGC::eot = ORGC::sub_length;
        ORGC::sor = 0;
        ORGC::eor = ORGC::sub_length;
        Position current_position;

        ORGC::text = ""; // Reset text for the first segment

        while (true) {
            ORGC utilities; // Calls constructor, resets ORGC::text and SequenceProcessor::hashmap
            int kmerlength_val = ORGC::kmer_length;

            if (ORGC::eor > reference_seq_content.length()) {
                if (ORGC::sot >= target_seq_content.size()) {
                    break;
                }
                ORGC::text = target_seq_content.substr(ORGC::sot);
                if (ORGC::text.length() <= 0) {
                    break;
                } else {
                    FileUtils::write(tempfile, ORGC::text, true);
                    break;
                }
            }
            if (ORGC::eot > target_seq_content.length()) {
                if (ORGC::sot >= target_seq_content.size()) {
                    break;
                }
                ORGC::text = target_seq_content.substr(ORGC::sot);
                if (ORGC::text.length() <= 0) {
                    break;
                } else {
                    FileUtils::write(tempfile, ORGC::text, true);
                    break;
                }
            }
            std::string current_ref_segment = reference_seq_content.substr(ORGC::sor, ORGC::eor - ORGC::sor);
            std::string current_tar_segment = target_seq_content.substr(ORGC::sot, ORGC::eot - ORGC::sot);

            ORGC::reference = current_ref_segment;
            ORGC::target = current_tar_segment;

            std::vector<Position> list_matches = Matcher::Lmatch(current_ref_segment, current_tar_segment,
                                                              kmerlength_val);

            if (list_matches.empty()) {
                kmerlength_val = 11;
                list_matches = Matcher::Lmatch(current_ref_segment, current_tar_segment, kmerlength_val);
            }

            if (list_matches.empty()) {
                ORGC::mismatch++;

                if (ORGC::eot >= target_seq_content.length() - 1) {
                    if (ORGC::sot >= target_seq_content.size()) {
                        break;
                    }
                    ORGC::text = target_seq_content.substr(ORGC::sot);
                    FileUtils::write(tempfile, ORGC::text, true);
                    break;
                }

                if (is_con) {
                    controuble++;
                }
                is_con = true;

                ORGC::text += current_tar_segment + "\n"; // Append to existing ORGC::text
                FileUtils::write(tempfile, ORGC::text, true);
                ORGC::text = ""; // Reset text for next utility call or end of segment
                ORGC::sot += ORGC::sub_length;
                ORGC::eot = ORGC::sot + ORGC::sub_length;
                ORGC::eor += ORGC::sub_length;


                int difference = target_seq_content.length() - ORGC::sot;

                if (difference <= ORGC::kmer_length) {
                    if (ORGC::sot >= target_seq_content.size()) {
                        break;
                    }
                    ORGC::text = target_seq_content.substr(ORGC::sot);
                    FileUtils::write(tempfile, ORGC::text, true);
                    ORGC::text = "";
                    break;
                } else if (difference < ORGC::sub_length) {
                    ORGC::eot = target_seq_content.length();
                }

                int difference_ref = reference_seq_content.length() - ORGC::sor;

                if (difference_ref < ORGC::sub_length) {
                    ORGC::eor = reference_seq_content.length();
                }
                if (ORGC::eot >= target_seq_content.length()) { // This condition might be redundant due to earlier checks
                     // break; // Original code had a break here if eot > target_seq_content.length()
                }
                if (difference_ref <= ORGC::kmer_length) {
                    if (ORGC::sot >= target_seq_content.size()) {
                        break;
                    }
                    ORGC::text = target_seq_content.substr(ORGC::sot);
                    if (ORGC::text.length() <= 0) {
                        break;
                    } else {
                        if (ORGC::text.length() > (ORGC::sub_length * ORGC::T1)) {
                            ORGC::mismatch++;
                        }
                        if (ORGC::mismatch > ORGC::T2) {
                            ORGC::local = false; // Signal to switch to global
                            // FileUtils::write(tempfile, ORGC::text, true); // Write remaining text before failing local
                            // ORGC::text = "";
                            return; // Exit local processing
                        }
                        FileUtils::write(tempfile, ORGC::text, true);
                        ORGC::text = "";
                        break;
                    }
                }
                continue;
            }

            is_con = false;
            if (controuble > 2)
                ORGC::mismatch -= controuble;
            controuble = 1;

            current_position = ORGC::format_matches(list_matches); // Appends to ORGC::text

            if (ORGC::mismatch > ORGC::T2) {
                ORGC::local = false; // Signal to switch to global
                // FileUtils::write(tempfile, ORGC::text, true); // Write current ORGC::text before failing
                // ORGC::text = "";
                return; // Exit local processing
            }

            ORGC::sot += current_position.getendinTar() + 1;
            ORGC::eot = ORGC::sot + ORGC::sub_length;
            ORGC::sor += ORGC::endref + 1; // ORGC::endref is updated by format_matches
            ORGC::endref = 0; 
            ORGC::eor = ORGC::sor + ORGC::sub_length;

            FileUtils::write(tempfile, ORGC::text, true); // Write the formatted matches and mismatches
            ORGC::text = ""; // Reset for the next segment or ORGC utilities call

            int difference = target_seq_content.length() - ORGC::sot;

            if (difference <= ORGC::kmer_length) {
                if (ORGC::sot < target_seq_content.size()) { // Check if there's anything left
                     ORGC::text = target_seq_content.substr(ORGC::sot);
                     if (ORGC::text.length() > 0) {
                        FileUtils::write(tempfile, ORGC::text, true);
                        ORGC::text = "";
                     }
                }
                break;
            } else if (difference < ORGC::sub_length) {
                ORGC::eot = target_seq_content.length();
            }

            int difference_ref = reference_seq_content.length() - ORGC::sor;

            if (difference_ref < ORGC::sub_length) {
                ORGC::eor = reference_seq_content.length();
            }
            if (difference_ref <= ORGC::kmer_length) {
                 if (ORGC::sot < target_seq_content.size()) {
                    ORGC::text = target_seq_content.substr(ORGC::sot);
                    if (ORGC::text.length() > 0) {
                        if (ORGC::text.length() > (ORGC::sub_length * ORGC::T1)) {
                            ORGC::mismatch++;
                        }
                        if (ORGC::mismatch > ORGC::T2) {
                            ORGC::local = false;
                            // FileUtils::write(tempfile, ORGC::text, true);
                            // ORGC::text = "";
                            return;
                        }
                        FileUtils::write(tempfile, ORGC::text, true);
                        ORGC::text = "";
                    }
                }
                break;
            }
        } // End of while(true) segment processing loop

        if (!ORGC::local) { // If loop decided local failed
            return;
        }

        // If local processing was successful up to this point
        ORGC::postprocess(tempfile, final_file_path);
        if (stat(tempfile.c_str(), &buffer) == 0) {
            if (std::remove(tempfile.c_str()) != 0) {
                std::cerr << "Warning: Failed to remove temp file after local processing: " << tempfile << std::endl;
            }
        }
    }

    static void process_global(const std::string& ref_genome_path, const std::string& tar_genome_path, const std::string& final_file_path, const std::string& tempfile) {
        std::cout << tar_genome_path << " is compressing (global attempt)...\n" << std::endl;
        struct stat buffer;
        if (stat(final_file_path.c_str(), &buffer) == 0) {
            if (std::remove(final_file_path.c_str()) != 0) {
                 std::cerr << "Warning: Failed to remove existing final file before global: " << final_file_path << std::endl;
            }
        }
        // Ensure tempfile is clean if a failed local attempt left it.
        if (stat(tempfile.c_str(), &buffer) == 0) {
            if (std::remove(tempfile.c_str()) != 0) {
                std::cerr << "Warning: Failed to remove existing temp file before global: " << tempfile << std::endl;
            }
        }


        std::string reference_seq_content_global = FileUtils::readSeq(ref_genome_path, GLOBAL, REFERENCE);
        // GreadtarSeq writes metadata and N/L lists directly to final_file_path, then returns sequence content
        std::string target_seq_content_global = FileUtils::GreadtarSeq(tar_genome_path, final_file_path);


        if (reference_seq_content_global.empty()) {
            std::cerr << "Error: Reference sequence is empty for " << ref_genome_path << " (global)" << std::endl;
            // Write empty or minimal content to final_file_path if needed, or let it be as GreadtarSeq left it
            return;
        }
        if (target_seq_content_global.empty() && FileUtils::meta_data.empty()) { // Check if GreadtarSeq itself indicated an issue
            std::cerr << "Error: Target sequence is effectively empty for " << tar_genome_path << " (global)" << std::endl;
            // GreadtarSeq might have written headers, so final_file_path might not be empty.
            // If target_seq_content_global is empty, no matching can be done.
            return;
        }


        std::string full_reference_seq = reference_seq_content_global;
        std::string full_target_seq = target_seq_content_global;

        ORGC::reference = full_reference_seq; // For format_matches if it uses static ORGC::target
        ORGC::target = full_target_seq;       // For format_matches

        std::vector<Position> list_global_matches = Matcher::Gmatch(full_reference_seq, full_target_seq,
                                                                 ORGC::kmer_length, ORGC::limit);

        if (!full_target_seq.empty() && !list_global_matches.empty()) {
            // format_matches for global writes to tempfile
            ORGC::format_matches(list_global_matches, tempfile);
            // postprocess reads tempfile and appends to final_file_path
            ORGC::postprocess(tempfile, final_file_path);
        } else if (full_target_seq.empty() && !list_global_matches.empty()){
             std::cerr << "Warning: Global matches found but target sequence content for matching was empty for " << tar_genome_path << std::endl;
        } else if (!full_target_seq.empty() && list_global_matches.empty()){
             std::cout << "Info: No global matches found for " << tar_genome_path << ". Target content might be written as is if GreadtarSeq handled it." << std::endl;
             // If GreadtarSeq wrote the sequence and there are no matches, postprocess might not be needed or might write an empty match list.
             // The current format_matches/postprocess path for global expects matches.
             // If no matches, GreadtarSeq already wrote N/L lists and metadata. The sequence itself is not written by GreadtarSeq.
             // If target is not empty but no matches, the final file will contain N/L lists and metadata, but no sequence data or match data from postprocess.
             // This seems to be the implicit behavior.
        }


        if (stat(tempfile.c_str(), &buffer) == 0) {
             if (std::remove(tempfile.c_str()) != 0) {
                 std::cerr << "Warning: Failed to remove temp file after global processing: " << tempfile << std::endl;
             }
        }
    }

public:
    static void process_genome(const std::string& ref_base_path, const std::string& tar_base_path, const std::string& output_dir_path) {
        std::string final_folder = output_dir_path + "/result";
        std::string misjuggements_path = final_folder + "/misjuggements.txt";

        std::ofstream misjuggements_file_clear(misjuggements_path, std::ios::out);
        if (!misjuggements_file_clear.is_open()) {
            std::cerr << "Could not open misjuggements.txt for clearing!" << std::endl;
        }
        misjuggements_file_clear.close();

        for (int i = 0; i < 24; i++) {
            ORGC::kmer_length = 21; // Reset default kmer_length for each chromosome
            ORGC::local = true;     // Assume local processing first
            ORGC::endref = ORGC::sub_length - 1; // Reset endref
            ORGC::mismatch = 0;     // Reset mismatch count for the chromosome
            // Reset T1, T2 to defaults, process_local will adjust if needed
            ORGC::T1 = 0.5;
            ORGC::T2 = 4;


            std::string current_chr_name = ORGC::chrName[i];
            std::string reference_genome_path = ref_base_path + "/" + current_chr_name;
            std::string target_genome_path = tar_base_path + "/" + current_chr_name;
            std::string final_file_path = final_folder + "/" + current_chr_name;
            std::string tempfile = final_folder + "/interim.txt"; // Common temp file name

            process_local(reference_genome_path, target_genome_path, final_file_path, tempfile);

            if (!ORGC::local) { // If local processing failed or decided it's not suitable
                process_global(reference_genome_path, target_genome_path, final_file_path, tempfile);
            }

            std::ofstream misjuggements_file(misjuggements_path, std::ios::app);
            if (!misjuggements_file.is_open()) {
                std::cerr << "Could not open misjuggements.txt for writing for " << current_chr_name << std::endl;
            } else {
                misjuggements_file << current_chr_name << ": " << ORGC::mismatch << std::endl;
                misjuggements_file.close();
            }
        }
    }
    
};

const std::string ORGC::chrName[] = {
        "chr1.fa", "chr2.fa", "chr3.fa", "chr4.fa", "chr5.fa", "chr6.fa", "chr7.fa", "chr8.fa", "chr9.fa",
        "chr10.fa", "chr11.fa", "chr12.fa", "chr13.fa", "chr14.fa", "chr15.fa", "chr16.fa",
        "chr17.fa", "chr18.fa", "chr19.fa", "chr20.fa",
        "chr21.fa", "chr22.fa", "chrX.fa", "chrY.fa"
};

std::string ORGC::reference = "";
std::string ORGC::target = "";
std::string ORGC::text = "";
int ORGC::kmer_length = 21;
int ORGC::sub_length = 30000;
int ORGC::limit = 100;
double ORGC::T1 = 0.5;
int ORGC::T2 = 4;
int ORGC::sot = 0;
int ORGC::eot = 0;
int ORGC::sor = 0;
int ORGC::eor = 0;
int ORGC::mismatch = 0;
int ORGC::endref = ORGC::sub_length - 1;
bool ORGC::local = true;

int FileUtils::length = 0;
std::string FileUtils::meta_data = "";

int main(int argc, char *argv[]) {

    if (argc != 4) {
        std::cout << "Make sure you have inputted 3 arguments." << std::endl;
        return 0;
    }

    std::string reference_base_path = argv[1];
    std::string target_base_path = argv[2];
    std::string output_dir = argv[3];

    struct stat st = {0};
    if (stat(output_dir.c_str(), &st) == -1) {
        if (mkdir(output_dir.c_str(), 0777) != 0) {
            std::cerr << "Failed to create output directory: " << output_dir << std::endl;
            return 1;
        }
    }
    std::string final_folder = output_dir + "/result";
    if (stat(final_folder.c_str(), &st) == -1) {
        if (mkdir(final_folder.c_str(), 0777) != 0) { // Changed from final_folder.c_str() to output_dir for the first mkdir
            std::cerr << "Failed to create result subdirectory: " << final_folder << std::endl;
            return 1;
        }
    }


    auto startDate = std::chrono::system_clock::now();
    long startCpuTimeNano = ORGC::getCPUTime();
    std::time_t start_time = std::chrono::system_clock::to_time_t(startDate);
    std::cout << "Start time: " << std::ctime(&start_time);

    SequenceProcessor::next_kmer.resize(SequenceProcessor::maxchar, -1);
    SequenceProcessor::kmer_location.resize(SequenceProcessor::maxseq, -1);

    ORGC::process_genome(reference_base_path, target_base_path, output_dir);

    FileUtils::use7zip(final_folder); // final_folder is output_dir + "/result"
    std::cout << "All Done" << std::endl;
    long taskCPUTimeNano = ORGC::getCPUTime() - startCpuTimeNano;
    std::cout << "Compressed time: " << (double) taskCPUTimeNano / 1000000000.0 << " seconds." << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;

    return 0;
}
