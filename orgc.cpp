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

class newkmer {
private:
    int kmerstart;
    std::string kmer;
    
public:
    void setkmerstart(int start) { kmerstart = start; }
    int getkmerstart() const { return kmerstart; }
    void setkmer(const std::string& k) { kmer = k; }
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

class ORGC {
public:
    static std::string reference, target;
    static std::string meta_data;
    static std::string text;
    static int kmer_length;
    static int sub_length, limit;
    static double T1;
    static int T2;
    static int sot, eot, sor, eor; // start of target,end of target,start of reference,end of reference
    static int length, mismatch, endref;
    static const int maxchar = 268435456; // 2^28
    static const int maxseq = 268435456*2; // 2^29
    static bool local;
    static std::unordered_map<int, std::vector<newkmer>> hashmap;
    static std::vector<int> next_kmer;
    static std::vector<int> kmer_location;

    // Constructor
    ORGC() {
        text = "";
        hashmap = std::unordered_map<int, std::vector<newkmer>>();
    }

    // Java String.hashCode() equivalent
    static int stringHashCode(const std::string& str) {
        int hash = 0;
        for (char c : str) {
            hash = 31 * hash + static_cast<int>(c);
        }
        return hash;
    }

    static long getCPUTime() {
        auto now = std::chrono::high_resolution_clock::now();
        auto duration = now.time_since_epoch();
        auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);
        return nanoseconds.count();
    }

    // read sequence for local matching
    static std::string LreadSeq(const std::string& sequenceFileName) {
        std::ifstream file(sequenceFileName);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + sequenceFileName);
        }
        
        std::stringstream stringbuilder;
        std::string line;
        bool base = false;
        
        std::getline(file, meta_data); // Chromosome label
        while (std::getline(file, line)) {
            stringbuilder << line;
            if (!base) {
                length = line.length();
                base = true;
            }
        }
        file.close();
        return stringbuilder.str();
    }

    // read reference sequence for global matching
    static std::string GreadrefSeq(const std::string& sequenceFileName) {
        std::ifstream file(sequenceFileName);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + sequenceFileName);
        }
        
        std::string line;
        int line_length;
        char temp_ch;

        std::stringstream stringbuilder;
        std::getline(file, line);
        while (std::getline(file, line)) {
            line_length = line.length();
            for (int i = 0; i < line_length; i++) {
                temp_ch = line[i];
                if (!std::isupper(temp_ch)) {
                    temp_ch = std::toupper(temp_ch);
                }
                if (temp_ch != 'N') {
                    stringbuilder << temp_ch;
                }
            }
        }
        file.close();
        return stringbuilder.str();
    }

    // read target sequence for global matching
    static std::string GreadtarSeq(const std::string& sequenceFileName, const std::string& fileName) {
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

        std::string meta_data_local, line; // Renamed to avoid conflict with static meta_data

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

    // build local hash table
    static void buildLhashtable(const std::string& read, int kmer_length_param) { // Renamed kmer_length to avoid shadowing static member
        int current_length = read.length(), i = kmer_length_param; // Use param
        std::string Nkmer = "";
        while (i > 0) {
            Nkmer += "N";
            i--;
        }

        int Nkey = stringHashCode(Nkmer);
        i = 0;

        std::cout << "buildLhashtable: current_length = " << current_length << ", kmer_length_param = " << kmer_length_param << std::endl;
         int count = 0;

        while (i < current_length - kmer_length_param + 1) {
            if (++count % 5000 == 0) std::cout << "Adding kmer #" << count << std::endl;

            std::string kmer_str = read.substr(i, kmer_length_param); // Use param, renamed kmer to kmer_str
            newkmer newKMer;
            newKMer.setkmerstart(i);
            newKMer.setkmer(kmer_str);
            int key = stringHashCode(kmer_str);
            if (hashmap.find(key) != hashmap.end()) {
                std::vector<newkmer>& list = hashmap[key];
                list.push_back(newKMer);
            } else {
                std::vector<newkmer> list;
                list.push_back(newKMer);
                hashmap[key] = list;
            }
            i++;
            if (key == Nkey) {
                while (i < current_length - kmer_length_param + 1 && (read.substr(i, 1) == "n" || read.substr(i, 1) == "N")) {
                    i++;
                }
            }
        }
        std::cout << "Exiting buildLhashtable normally" << std::endl;
    }

    // build global hash table
    static void buildGhashtable(const std::string& read, int kmer_length_param) { // Renamed kmer_length
        int iteration = read.length() - kmer_length_param + 1; // Use param

        for (int i = 0; i < maxseq; i++) {
            kmer_location[i] = -1;
        }
        for (int i = 0; i < iteration; i++) {
            if (i % 1000000 == 0) {
               std::cout << "buildGhashtable: Processing kmer at index " << i << std::endl;
            }
            std::string kmer_str = read.substr(i, kmer_length_param); // Use param, renamed kmer
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

    static std::vector<Position> lowercase_position(const std::string& sequence) {
        std::vector<Position> list;
        bool successive = false;
        int start = 0, current_end = 0; // Renamed end to current_end

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
            position.setendinTar(current_end); // Corrected: was end, should be current_end
            list.push_back(position);
        }
        return list;
    }

    static int get_incre(int endinRef_param, int endinTar_param, int Refendindex, int Tarendindex) { // Renamed params
      std::cout << "get_incre called with endinRef_param: " << endinRef_param
          << " endinTar_param: " << endinTar_param
          << " Refendindex: " << Refendindex
          << " Tarendindex: " << Tarendindex << std::endl;

        int position = 0;
        int endIndex;

        if (Refendindex - endinRef_param <= Tarendindex - endinTar_param) {
            endIndex = Refendindex - endinRef_param + 1;
        } else {
            endIndex = Tarendindex - endinTar_param + 1;
        }

        for (int i = 1; i < endIndex; i++) {
         // BOUNDS CHECK
         if ((endinTar_param + i) >= target.size() || (endinRef_param + i) >= reference.size())
            std::cerr << "OUT OF BOUNDS: endinTar_param + i = " << (endinTar_param + i)
               << " target.size() = " << target.size()
               << ", endinRef_param + i = " << (endinRef_param + i)
               << " reference.size() = " << reference.size() << std::endl;
            break; // Prevent out-of-bounds
         if (target[endinTar_param + i] == reference[endinRef_param + i]) {
            position++;
         } else {
            break;
         }
      }

        return position;
    }

    // local match
    static std::vector<Position> Lmatch(const std::string& ref_seq, const std::string& tar_seq, int kmer_len) { // Renamed params
        std::vector<Position> list;
        int index = 0, startIndex;
        int increment, most_incre, key_val; // Renamed key to key_val
        int kmerstart_val, endinRef_val, endinTar_val, Refendindex, Tarendindex; // Renamed params
        std::string kmer_str; // Renamed kmer

        std::cout << "In Lmatch, ref_seq.size(): " << ref_seq.size() << " tar_seq.size(): " << tar_seq.size() << " kmer_len: " << kmer_len << std::endl;
        buildLhashtable(ref_seq, kmer_len);
        std::cout << "buildLhashtable OK" << std::endl;
        
        while (true) {
            increment = 0; most_incre = 0;
            if (index + kmer_len > tar_seq.length()) { // Check before substr
                break;
            }
            kmer_str = tar_seq.substr(index, kmer_len);
            key_val = stringHashCode(kmer_str);

            if (hashmap.find(key_val) == hashmap.end()) {
                startIndex = std::numeric_limits<int>::max();
                index = index + 1;
                if (index + kmer_len > tar_seq.length()) {
                    break;
                }
                continue;
            }                       
           
            std::vector<newkmer>& klist = hashmap[key_val];
            startIndex = std::numeric_limits<int>::max();
            most_incre = 0;
            for (const newkmer& newKMer : klist) {
                if (newKMer.getkmer() == kmer_str) {
                    kmerstart_val = newKMer.getkmerstart();
                    endinRef_val = kmerstart_val + kmer_len - 1;
                    endinTar_val = index + kmer_len - 1;
                    Refendindex = ref_seq.length() - 1; // Use ref_seq
                    Tarendindex = tar_seq.length() - 1; // Use tar_seq                   	
                    increment = get_incre(endinRef_val, endinTar_val, Refendindex, Tarendindex); // Pass local endinRef/Tar
                    
                    if (klist.size() > 1) {                  		                   	
                        if (increment == most_incre) {
                            if (list.size() > 1) { // Should be list.size() > 0 for safety, or check if list is not empty
                                int lastEIR = list.back().getendinRef(); // Use back() for last element
                                if (std::abs(kmerstart_val - lastEIR) < std::abs(startIndex - lastEIR)) // Compare absolute differences
                                    startIndex = kmerstart_val;
                            } else if (list.empty()){ // If list is empty, this is the first potential match
                                startIndex = kmerstart_val;
                            }
                        }
                        else if (increment > most_incre) {
                            most_incre = increment;
                            startIndex = kmerstart_val;
                        }
                    }
                    else { // klist.size() == 1
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
            
            Position position_obj; // Renamed position
            position_obj.setstartinTar(index);
            position_obj.setendinTar(index + kmer_len + most_incre - 1);
            position_obj.setstartinRef(startIndex);
            position_obj.setendinRef(startIndex + kmer_len + most_incre - 1);
            list.push_back(position_obj);
            index = index + kmer_len + most_incre + 1; // Should be most_incre, not most_incre + 1, if get_incre returns count of matching chars
            if (index + kmer_len > tar_seq.length()) {
                break;
            }      
        }
        return list;
    }

    static std::vector<Position> Gmatch(const std::string& ref_seq, const std::string& tar_seq, int kmer_len) { // Renamed params
        std::vector<Position> list;
        int index = 0, startIndex, lastEIR = 0;
        std::string kmer_str; // Renamed kmer

        buildGhashtable(ref_seq, kmer_len); // Pass ref_seq and kmer_len

        while (true) {
            int increment = 0, most_incre = 0;
            if (index + kmer_len > tar_seq.length()) { // Check before substr
                break;
            }
            kmer_str = tar_seq.substr(index, kmer_len);
            int key_val = std::abs(stringHashCode(kmer_str)); // Renamed key
            
            if (key_val == -2147483648) { // This condition will never be true due to std::abs
                key_val = 0;
            } 			
            while (key_val > maxseq - 1) {
                key_val = key_val / 2;
            }

            if (kmer_location[key_val] == -1) {
                startIndex = std::numeric_limits<int>::max();
                index = index + 1;
                if (index + kmer_len > tar_seq.length()) {
                    break;
                }
                continue;
            }

            startIndex = std::numeric_limits<int>::max();
            most_incre = 0;
            bool match_found = false; // Renamed match to match_found
            
            for (int k = kmer_location[key_val]; k != -1; k = next_kmer[k]) {
                increment = 0;
                if (k + kmer_len > ref_seq.length()) continue; // Bounds check for Rkmer
                std::string Rkmer = ref_seq.substr(k, kmer_len);
                if (kmer_str != Rkmer) {
                    continue;
                }
                try {
                    if (!list.empty()) { // Check if list is not empty before accessing
                        lastEIR = list.back().getendinRef();
                    } else {
                        lastEIR = 0; // Default if list is empty
                    }
                } catch (...) { // Catch specific exceptions if possible
                    lastEIR = 0;
                }
                if (std::abs(k - lastEIR) > limit) { // Check absolute difference against limit
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
                     if (!list.empty()) { // Check if list is not empty
                        if (std::abs(k - lastEIR) < std::abs(startIndex - lastEIR)) // Compare absolute differences
                            startIndex = k;
                    } else if (list.empty() && startIndex == std::numeric_limits<int>::max()){
                         startIndex = k; // If list is empty, this is the first potential match
                    }
                } else if (increment > most_incre) {
                    most_incre = increment;
                    startIndex = k;
                }
            }
            if (!match_found) { // Redundant loop, logic should be combined or re-evaluated
                for (int k = kmer_location[key_val]; k != -1; k = next_kmer[k]) {
                    increment = 0;
                    if (k + kmer_len > ref_seq.length()) continue; // Bounds check
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
                        if (!list.empty()) { // Check if list is not empty
                             if (std::abs(k - lastEIR) < std::abs(startIndex - lastEIR)) // Compare absolute differences
                                startIndex = k;
                        } else if (list.empty() && startIndex == std::numeric_limits<int>::max()){
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

            Position position_obj; // Renamed position
            position_obj.setstartinTar(index);
            position_obj.setendinTar(index + kmer_len + most_incre - 1);
            position_obj.setstartinRef(startIndex);
            position_obj.setendinRef(startIndex + kmer_len + most_incre - 1);
            list.push_back(position_obj);
            index = index + kmer_len + most_incre + 1; // As in LMatch, consider if +1 is needed
            if (index + kmer_len > tar_seq.length()) {
                break;
            }
        }
        return list;
    }
    
    static Position format_matches(const std::vector<Position>& list_param) { // Renamed list
        // This function uses static `text`, `target`, `sor`, `sub_length`, `T1`, `mismatch`, `endref`
        int startinTar_val, startinRef_val, endinRef_val; // Renamed
        int trouble = 0;
        if (list_param.empty()) { // Handle empty list
             // Decide what to return or throw an error. Returning an empty/default Position might be an option.
             // For now, let's assume it implies no match and doesn't modify static `text`.
             // However, the original code would crash here: list[list.size() - 1]
             // Returning a default constructed Position if list is empty.
             return Position(); 
        }

        for (int i = 0; i < list_param.size(); i++) {

            if (i == 0) {
                startinTar_val = list_param[i].getstartinTar();
                startinRef_val = list_param[i].getstartinRef();
                endinRef_val = list_param[i].getendinRef();
                if (endinRef_val >= endref) { // ORGC::endref
                    endref = endinRef_val;    // ORGC::endref
                }
                
                if (startinTar_val > 0) {
                    std::string preamble = target.substr(0, startinTar_val); // ORGC::target
                    text += preamble + "\n"; // ORGC::text
                    trouble += preamble.length();
                }
                text += "" + std::to_string(startinRef_val + sor) + "," + std::to_string(endinRef_val + sor) + "\n"; // ORGC::text, ORGC::sor
                continue;
            }

            startinTar_val = list_param[i].getstartinTar();
            startinRef_val = list_param[i].getstartinRef();
            endinRef_val = list_param[i].getendinRef();
            if (endinRef_val >= endref) { // ORGC::endref
                endref = endinRef_val;    // ORGC::endref
            }                    

            int endinTarPrev = list_param[i - 1].getendinTar();
            // Bounds check for substr
            if (endinTarPrev + 1 < startinTar_val) {
                 std::string mismatch_str = target.substr(endinTarPrev + 1, startinTar_val - (endinTarPrev + 1)); // ORGC::target, renamed mismatch
                 if (mismatch_str.length() > 0) {
                    text += mismatch_str + "\n"; // ORGC::text
                    trouble += mismatch_str.length();
                }
            } else if (endinTarPrev + 1 == startinTar_val) {
                // No mismatch, string would be empty
            } else {
                // Overlapping or invalid sequence of matches, handle error or log
            }


            text += std::to_string(startinRef_val + sor) + "," + std::to_string(endinRef_val + sor) + "\n"; // ORGC::text, ORGC::sor
        }
        if (trouble > (sub_length * T1)) // ORGC::sub_length, ORGC::T1
        {
            mismatch++; // ORGC::mismatch
        }
        
        if (!list_param.empty()) {
            int endinTar_last = list_param.back().getendinTar();
            if (endinTar_last < target.length() - 1) {
                std::string tail = target.substr(endinTar_last + 1, target.length() - (endinTar_last + 1));
                if (!tail.empty()) {
                    // For text accumulation version:
                    text += tail + "\n";
                    // For stringstream/file version:
                    // stringbuilder << tail << "\n";
                }
            }
        }

        return list_param.back(); // Returns last element
    }
    
    static void format_matches(const std::vector<Position>& list_param, const std::string& fileName) { // Renamed list
        // This function uses static `target`
        std::stringstream stringbuilder;

        int startinTar_val, startinRef_val, endinRef_val, endinTar_val = 0; // Renamed
        for (int i = 0; i < list_param.size(); i++) {
            if (i == 0) {
                startinTar_val = list_param[i].getstartinTar();
                endinTar_val = list_param[i].getendinTar();
                startinRef_val = list_param[i].getstartinRef();
                endinRef_val = list_param[i].getendinRef();
                if (startinTar_val > 0) {
                    std::string preamble = target.substr(0, startinTar_val); // ORGC::target
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
            // Bounds check for substr
            if (endinTarPrev + 1 < startinTar_val) {
                std::string mismatch_str = target.substr(endinTarPrev + 1, startinTar_val - (endinTarPrev + 1)); // ORGC::target, renamed mismatch
                if (mismatch_str.length() > 0) {
                    stringbuilder << mismatch_str << "\n";
                }
            } else if (endinTarPrev + 1 == startinTar_val) {
                // No mismatch
            } else {
                // Overlapping, handle error
            }
            stringbuilder << startinRef_val << "," << endinRef_val << "\n";
        }
        if (endinTar_val < (target.length() - 1)) { // ORGC::target
            // Bounds check for substr
            if (endinTar_val + 1 < target.length()) {
                stringbuilder << target.substr(endinTar_val + 1, (target.length()) - (endinTar_val + 1)) << "\n"; // ORGC::target
            }
        }
        write(fileName, stringbuilder.str(), true); // ORGC::write
        if (!list_param.empty()) {
            int endinTar_last = list_param.back().getendinTar();
            if (endinTar_last < target.length() - 1) {
                std::string tail = target.substr(endinTar_last + 1, target.length() - (endinTar_last + 1));
                if (!tail.empty()) {
                    // For text accumulation version:
                    text += tail + "\n";
                    // For stringstream/file version:
                    // stringbuilder << tail << "\n";
                }
            }
        }
    }
    
    static std::vector<Position> format_matches_N(const std::string& sequence) { // Renamed to avoid overload ambiguity if types were closer
        std::vector<Position> list;
        bool successive = false;
        int start = 0, current_end = 0; // Renamed end

        for (int i = 0; i < sequence.length(); i++) {
            if (sequence[i] == 'N' || sequence[i] == 'n') {
                if (successive) {
                    current_end += 1;
                } else {
                    start = i;
                    current_end += 1; // Should be start + 1 or i + 1 if current_end tracks length
                    successive = true;
                }
            } else {
                if (successive) {
                    Position position_obj; // Renamed
                    position_obj.setstartinTar(start);
                    position_obj.setendinTar(current_end -1); // If current_end is next position, then -1 is correct end index
                    list.push_back(position_obj);
                }
                successive = false;
                start = 0; // Reset start
                current_end = i + 1; // current_end should track the end of the segment or next start
            }
        }

        if (successive) {
            Position position_obj; // Renamed
            position_obj.setstartinTar(start);
            position_obj.setendinTar(current_end -1); // Assuming current_end was length or next pos
            list.push_back(position_obj);
        }
        return list;
    }
    
    static void write(const std::string& filename, const std::string& text_content, bool append) { // Renamed text to text_content
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

// REMOVE THE FOLLOWING DUPLICATE CLASS DEFINITIONS (approx. lines 575-633 in original)
// class newkmer { ... };
// static std::unordered_map<int, std::vector<newkmer>> hashmap; // This was line 635
// class Position { ... };

    static void write(std::string filename, std::vector<Position> list_param, bool append, std::string auxiliary) { // Renamed list
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
        
        for (Position position_obj : list_param) { // Renamed position
            int start_val = position_obj.getstartinTar(); // Renamed start
            int end_val = position_obj.getendinTar();     // Renamed end
            text_stream << (start_val - temp_end) << " " << (end_val - start_val) << " ";
            temp_end = end_val;
        }
        output << text_stream.str();
        output << "\n";
        output.flush();
        output.close();
    }

// REMOVE THE FOLLOWING DUPLICATE write METHOD (approx. lines 660-673 in original)
// static void write(std::string filename, std::string content, bool append) { ... }

    static void postprocess(std::string filename, std::string final_file) {
        std::ifstream input(filename);
        if (!input.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        
        std::string line;
        std::ostringstream stringbuilder;
        std::vector<int> num_list; // Renamed list to num_list
        
        while (std::getline(input, line)) {
            if (line.find(",") != std::string::npos) {
                size_t comma_pos = line.find(",");
                int begin = std::stoi(line.substr(0, comma_pos));
                int end_val = std::stoi(line.substr(comma_pos + 1)); // Renamed end
                
                if (num_list.size() > 0) {
                    int prevEnd = num_list.back(); // Use back() for last element
                    if (begin != prevEnd + 1) {
                        stringbuilder << num_list[0] << "," << num_list.back() << "\n";
                        num_list.clear();
                    }
                }
                num_list.push_back(begin);
                num_list.push_back(end_val);
            } else if (line.length() > 0) {
                // Original code had a try-catch for an out_of_range exception from list[0] or list[list.size()-1]
                // This implies num_list could be empty here.
                if (!num_list.empty()) { // Check if num_list is not empty
                    stringbuilder << num_list[0] << "," << num_list.back() << "\n";
                }
                // The original code had a condition `if (line.find("^") == std::string::npos)`
                // This check is maintained.
                if (line.find("^") == std::string::npos) {
                     stringbuilder << line << "\n";
                }
                num_list.clear();
            }
        }
        if (num_list.size() > 0) { // Check if num_list is not empty
            stringbuilder << num_list[0] << "," << num_list.back() << "\n";
        }
        input.close();

        std::istringstream inputStream(stringbuilder.str());
        stringbuilder.str(""); // Clear stringbuilder
        stringbuilder.clear(); // Clear error flags
        num_list.clear();      // Clear num_list again
        std::vector<std::string> stringList;

        while (std::getline(inputStream, line)) {
            stringList.push_back(line);
        }

        int prev = 0;
        bool successive = false;
        for (size_t i = 0; i < stringList.size(); i++) { // Use size_t for loop counter
            std::string str = stringList[i];

            if (str.find(",") != std::string::npos) {
                size_t comma_pos = str.find(",");
                int begin = std::stoi(str.substr(0, comma_pos));
                int end_val = std::stoi(str.substr(comma_pos + 1)); // Renamed end
                if (!successive) {
                    num_list.push_back(begin);
                    num_list.push_back(end_val - begin);
                    prev = end_val;
                    // stringbuilder << num_list[0] << "," << num_list.back() << "\n"; // Original logic was complex here
                    // The original code printed list[0],list[list.size()-1] which would be begin, (end-begin)
                    // This seems to be transforming coordinates. Let's stick to the original output format if it was intended.
                    // The original code here was: stringbuilder << list[0] << "," << list[list.size() - 1] << "\n";
                    // If list contains [begin, end-begin], then list.size()-1 is index 1.
                    // So it prints begin, (end-begin).
                    if (!num_list.empty()) stringbuilder << num_list[0] << "," << num_list.back() << "\n";

                    successive = true;
                } else {
                    num_list.push_back(begin - prev);
                    num_list.push_back(end_val - begin);
                    prev = end_val;
                    // Same logic for printing as above
                    if (!num_list.empty()) stringbuilder << num_list[0] << "," << num_list.back() << "\n";
                }
                num_list.clear(); // Cleared after each processed coordinate line
            } else if (str.length() > 0) {
                stringbuilder << str << "\n";
                successive = false; // A non-coordinate line might break succession
            }
        }

        if (!num_list.empty()) { // Check if num_list is not empty
            stringbuilder << num_list[0] << "," << num_list.back() << "\n";
        }

        write(final_file, stringbuilder.str(), true); // ORGC::write
    }

    // 7zip
    static void use7zip(std::string filename) {
        struct stat buffer;
        if (stat(filename.c_str(), &buffer) != 0) {
            // File does not exist or error, perhaps log or handle
            std::cerr << "Warning: File " << filename << " not found for 7zip." << std::endl;
            return;
        }
        
        std::string exec = "./7za a " + filename + ".7z " + filename + " -m0=PPMD";
        try {
            int result = std::system(exec.c_str());
            if (result != 0) {
                std::cerr << "Warning: 7zip command failed with result " << result << std::endl;
            }
            // (void)result; // suppress unused variable warning if not checking result
        } catch (const std::exception& e) { // Catch standard exceptions
            std::cerr << "Exception in use7zip: " << e.what() << std::endl;
        }
    }

// REMOVE THE ENTIRE BLOCK OF DUPLICATE/STUB HELPER METHODS HERE
// (approx. lines 827-915 in original file)
// This includes LreadSeq, lowercase_position, Lmatch, format_matches (all overloads),
// GreadrefSeq, GreadtarSeq, Gmatch, getCPUTime.

}; // End of ORGC class

// Static member definitions must be outside the class
std::vector<int> ORGC::next_kmer;
std::vector<int> ORGC::kmer_location;
std::string ORGC::reference = "";
std::string ORGC::target = "";
std::string ORGC::meta_data = "";
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
int ORGC::length = 0;
int ORGC::mismatch = 0;
int ORGC::endref = ORGC::sub_length - 1; // ORGC::sub_length should be initialized before this
// int ORGC::next_kmer[ORGC::maxchar]; // Definition for arrays is tricky, often done in .cpp if header/source split
// int ORGC::kmer_location[ORGC::maxseq];
// For single file, these are fine, but ensure maxchar/maxseq are usable at this point.
// They are const int, so it's okay.
bool ORGC::local = true;
std::unordered_map<int, std::vector<newkmer>> ORGC::hashmap;


// main function should be global, not a static member of ORGC
int main(int argc, char* argv[]) {
    
    if (argc != 4) {
        std::cout << "Make sure you have inputted 3 arguments." << std::endl;
        return 0;
    }
    
    std::string reference_genome_path; // Renamed to avoid conflict with ORGC::reference
    std::string target_genome_path;    // Renamed to avoid conflict with ORGC::target
    std::string final_folder = std::string(argv[3]) + "/result"; // output file folder
    std::string mkdir_cmd = "mkdir \"" + final_folder + "\"";  // Enclose path in quotes for spaces
    std::system(mkdir_cmd.c_str()); // Use c_str()
    std::string final_file_path; // Renamed
    int controuble;
    bool is_con;
    
    std::string chrName[] = { "chr1.fa", "chr2.fa", "chr3.fa", "chr4.fa", "chr5.fa", "chr6.fa", "chr7.fa", "chr8.fa", "chr9.fa",
            "chr10.fa", "chr11.fa", "chr12.fa", "chr13.fa", "chr14.fa", "chr15.fa", "chr16.fa", "chr17.fa", "chr18.fa", "chr19.fa", "chr20.fa",
            "chr21.fa", "chr22.fa", "chrX.fa", "chrY.fa" };

    // Get system memory info (simplified)
    std::cout << "MAX RAM: Available" << std::endl;

    auto startDate = std::chrono::system_clock::now();
    long startCpuTimeNano = ORGC::getCPUTime(); // ORGC scope
    std::time_t start_time = std::chrono::system_clock::to_time_t(startDate);
    std::cout << "Start time: " << std::ctime(&start_time);
   
   std::cout << "Resizing next_kmer..." << std::endl;
   ORGC::next_kmer.resize(ORGC::maxchar, -1);
   std::cout << "next_kmer done. Resizing kmer_location..." << std::endl;
   ORGC::kmer_location.resize(ORGC::maxseq, -1);
   std::cout << "Both vectors resized." << std::endl;

    
    for (int i = 0; i < 24; i++) {
        ORGC::kmer_length = 21; // ORGC scope

        controuble = 0;
        is_con = false;
        ORGC::local = true;      // ORGC scope
        ORGC::endref = ORGC::sub_length - 1; // ORGC scope
        
        reference_genome_path = std::string(argv[1]) + "/" + chrName[i]; 
        target_genome_path = std::string(argv[2]) + "/" + chrName[i]; 
        final_file_path = final_folder + "/" + chrName[i];

        while (ORGC::local) { // ORGC scope
            ORGC::mismatch = 0; // ORGC scope
            std::cout << target_genome_path << " is compressing...\n" << std::endl;

            struct stat buffer;
            if (stat(final_file_path.c_str(), &buffer) == 0) {
                std::remove(final_file_path.c_str());
            }

            std::string greference_path = reference_genome_path; // Renamed
            std::string gtarget_path = target_genome_path;       // Renamed
            std::string tempfile = final_folder + "/interim.txt";

            std::string reference_seq_content = ORGC::LreadSeq(greference_path); // ORGC scope, renamed var
            std::string target_seq_content = ORGC::LreadSeq(gtarget_path);       // ORGC scope, renamed var
            if (reference_seq_content.empty()) {
               std::cerr << "Error: Reference sequence is empty for " << greference_path << std::endl;
               break;
            }
            if (target_seq_content.empty()) {
               std::cerr << "Error: Target sequence is empty for " << gtarget_path << std::endl;
               break;
            }
            
            std::cout << "Reference: " << greference_path << " length: " << reference_seq_content.size() << std::endl;
            std::cout << "Target: " << gtarget_path << " length: " << target_seq_content.size() << std::endl;

            int target_length_val = target_seq_content.length(); // Renamed var
            
            if (target_length_val < ORGC::sub_length * 5) { // ORGC scope
                ORGC::local = false; // ORGC scope
                break;
            } else if (target_length_val < ORGC::sub_length * 1333) { // ORGC scope
                ORGC::T1 = 0.1;  // ORGC scope
                ORGC::T2 = 0;    // ORGC scope
            }
            
            std::string auxiliary = ORGC::meta_data + "\n" + std::to_string(ORGC::length) + "\n"; // ORGC scope

            std::cout << "Read sequences OK" << std::endl;
            std::vector<Position> L_list = ORGC::lowercase_position(target_seq_content); // ORGC scope
            std::cout << "lowercase_position OK" << std::endl;
            ORGC::write(final_file_path, L_list, false, auxiliary); // ORGC scope
            std::cout << "write 1 OK" << std::endl;
            ORGC::write(final_file_path, "\n", true);               // ORGC scope
            std::cout << "write 2 OK" << std::endl;

            std::cout << "Before std::transform upper" << std::endl;
            std::transform(reference_seq_content.begin(), reference_seq_content.end(), reference_seq_content.begin(), ::toupper);
            std::cout << "After upper ref" << std::endl;
            std::transform(target_seq_content.begin(), target_seq_content.end(), target_seq_content.begin(), ::toupper);
            std::cout << "After upper target" << std::endl;

            if (stat(tempfile.c_str(), &buffer) == 0) {
                std::remove(tempfile.c_str());
            }
            
            ORGC::sot = 0;                         // ORGC scope
            ORGC::eot = ORGC::sub_length;          // ORGC scope
            ORGC::sor = 0;                         // ORGC scope
            ORGC::eor = ORGC::sub_length;          // ORGC scope
            Position current_position;             // Renamed var

            ORGC::text = ""; // Explicitly reset text for this local attempt if constructor doesn't run here.
                             // The ORGC utilities object below will reset it.

            while (true) {
                std::cout << "New segment: sot=" << ORGC::sot << " eot=" << ORGC::eot << " sor=" << ORGC::sor << " eor=" << ORGC::eor << std::endl;
                ORGC utilities; // Calls constructor, resets ORGC::text and ORGC::hashmap
                int kmerlength_val = ORGC::kmer_length; // ORGC scope, Renamed var

                if (ORGC::eor > reference_seq_content.length()) { // ORGC scope
                  if (ORGC::sot >= target_seq_content.size()) {
    std::cerr << "sot (" << ORGC::sot << ") >= target_seq_content.size() (" << target_seq_content.size() << ")" << std::endl;
    break;
}
                    ORGC::text = target_seq_content.substr(ORGC::sot); // ORGC scope
                    if (ORGC::text.length() <= 0) { // ORGC scope
                        break;
                    } else {
                        ORGC::write(tempfile, ORGC::text, true); // ORGC scope
                        break;
                    }
                }
                if (ORGC::eot > target_seq_content.length()) { // ORGC scope
                  if (ORGC::sot >= target_seq_content.size()) {
    std::cerr << "sot (" << ORGC::sot << ") >= target_seq_content.size() (" << target_seq_content.size() << ")" << std::endl;
    break;
}
                    ORGC::text = target_seq_content.substr(ORGC::sot); // ORGC scope
                    if (ORGC::text.length() <= 0) { // ORGC scope
                        break;
                    } else {
                        ORGC::write(tempfile, ORGC::text, true); // ORGC scope
                        break;
                    }
                }
                // Assign to static members ORGC::reference and ORGC::target for Lmatch/Gmatch if they use these static members
                // The Lmatch/Gmatch functions were changed to take string& parameters, so this is better:
                std::cout << "Before substrings" << std::endl;
                std::string current_ref_segment = reference_seq_content.substr(ORGC::sor, ORGC::eor - ORGC::sor);
                std::string current_tar_segment = target_seq_content.substr(ORGC::sot, ORGC::eot - ORGC::sot);
                std::cout << "Created substrings, ref_segment.size(): " << current_ref_segment.size() << " tar_segment.size(): " << current_tar_segment.size() << std::endl;

               ORGC::reference = current_ref_segment;
               ORGC::target = current_tar_segment;


                std::vector<Position> list_matches = ORGC::Lmatch(current_ref_segment, current_tar_segment, kmerlength_val); // ORGC scope
                std::cout << "Lmatch OK, list_matches.size(): " << list_matches.size() << std::endl;

                if (list_matches.empty()) { // Check .empty() instead of size() <=0
                    kmerlength_val = 11; 
                    list_matches = ORGC::Lmatch(current_ref_segment, current_tar_segment, kmerlength_val); // ORGC scope
                }

                if (list_matches.empty()) {
                    ORGC::mismatch++; // ORGC scope
                    
                    if (ORGC::eot >= target_seq_content.length() - 1) { // ORGC scope
                     if (ORGC::sot >= target_seq_content.size()) {
    std::cerr << "sot (" << ORGC::sot << ") >= target_seq_content.size() (" << target_seq_content.size() << ")" << std::endl;
    break;
}
                        ORGC::text = target_seq_content.substr(ORGC::sot); // ORGC scope
                        ORGC::write(tempfile, ORGC::text, true); // ORGC scope
                        break;
                    }
                    
                    if (is_con) {
                        controuble++;
                    }
                    is_con = true;
                    
                    // ORGC::text += ORGC::target + "\n"; // ORGC::target was the segment, use current_tar_segment
                    ORGC::text += current_tar_segment + "\n"; 
                    ORGC::write(tempfile, ORGC::text, true); // ORGC scope
                    ORGC::sot += ORGC::sub_length; // ORGC scope
                    ORGC::eot = ORGC::sot + ORGC::sub_length; // ORGC scope
                    ORGC::eor += ORGC::sub_length; // ORGC scope


                    int difference = target_seq_content.length() - ORGC::sot; // ORGC scope

                    if (difference <= ORGC::kmer_length) { // ORGC scope
                     if (ORGC::sot >= target_seq_content.size()) {
    std::cerr << "sot (" << ORGC::sot << ") >= target_seq_content.size() (" << target_seq_content.size() << ")" << std::endl;
    break;
}
                        ORGC::text = target_seq_content.substr(ORGC::sot); // ORGC scope
                        ORGC::write(tempfile, ORGC::text, true); // ORGC scope
                        break;
                    } else if (difference < ORGC::sub_length) { // ORGC scope
                        ORGC::eot = target_seq_content.length() -1; 
                         ORGC::eot = target_seq_content.length(); 
                    }

                    int difference_ref = reference_seq_content.length() - ORGC::sor; // ORGC scope

                    if (difference_ref < ORGC::sub_length) { // ORGC scope
                        ORGC::eor = reference_seq_content.length(); // ORGC scope, assuming exclusive
                    }
                    if (ORGC::eot >= target_seq_content.length()) { // ORGC scope
                        break;
                    }
                    if (difference_ref <= ORGC::kmer_length) { // ORGC scope
                     if (ORGC::sot >= target_seq_content.size()) {
    std::cerr << "sot (" << ORGC::sot << ") >= target_seq_content.size() (" << target_seq_content.size() << ")" << std::endl;
    break;
}
                        ORGC::text = target_seq_content.substr(ORGC::sot); // ORGC scope
                        if (ORGC::text.length() <= 0) { // ORGC scope
                            break;
                        } else {
                            if (ORGC::text.length() > (ORGC::sub_length * ORGC::T1)) { // ORGC scope
                                ORGC::mismatch++; // ORGC scope
                            }
                            if (ORGC::mismatch > ORGC::T2) { // ORGC scope
                                ORGC::local = false; // ORGC scope
                                break;
                            }
                            ORGC::write(tempfile, ORGC::text, true); // ORGC scope
                            break;
                        }
                    }
                    continue;
                }

                is_con = false;
                if (controuble > 2)
                    ORGC::mismatch -= controuble; // ORGC scope
                controuble = 1;
                
                // format_matches modifies static ORGC::text and ORGC::endref
                current_position = ORGC::format_matches(list_matches); // ORGC scope

                if (ORGC::mismatch > ORGC::T2) { // ORGC scope
                    ORGC::local = false; // ORGC scope
                    break;
                }

                ORGC::sot += current_position.getendinTar() + 1; 
                ORGC::eot = ORGC::sot + ORGC::sub_length;
                ORGC::sor += ORGC::endref + 1;                  
                ORGC::endref = 0; 
                ORGC::eor = ORGC::sor + ORGC::sub_length;  
                ORGC::write(tempfile, ORGC::text, true); // ORGC scope
                ORGC::text = ""; 


                int difference = target_seq_content.length() - ORGC::sot; // ORGC scope

                if (difference <= ORGC::kmer_length) { // ORGC scope
                  if (ORGC::sot >= target_seq_content.size()) {
    std::cerr << "sot (" << ORGC::sot << ") >= target_seq_content.size() (" << target_seq_content.size() << ")" << std::endl;
    break;
}
                    ORGC::text = target_seq_content.substr(ORGC::sot); // ORGC scope
                    if (ORGC::text.length() <= 0) { // ORGC scope
                        break;
                    } else {
                        ORGC::write(tempfile, ORGC::text, true); // ORGC scope
                        break;
                    }

                } else if (difference < ORGC::sub_length) { // ORGC scope
                    ORGC::eot = target_seq_content.length(); 

                int difference_ref = reference_seq_content.length() - ORGC::sor; // ORGC scope

                if (difference_ref < ORGC::sub_length) { 
                    ORGC::eor = reference_seq_content.length(); 
                }
                if (difference_ref <= ORGC::kmer_length) { 
                  if (ORGC::sot >= target_seq_content.size()) {
    std::cerr << "sot (" << ORGC::sot << ") >= target_seq_content.size() (" << target_seq_content.size() << ")" << std::endl;
    break;
}
                    ORGC::text = target_seq_content.substr(ORGC::sot); 
                    if (ORGC::text.length() <= 0) { 
                        break;
                    } else {
                        if (ORGC::text.length() > (ORGC::sub_length * ORGC::T1)) { // ORGC scope
                            ORGC::mismatch++; // ORGC scope
                        }
                        if (ORGC::mismatch > ORGC::T2) { // ORGC scope
                            ORGC::local = false; // ORGC scope
                            break;
                        }
                        ORGC::write(tempfile, ORGC::text, true); // ORGC scope
                        break;
                    }
                }
            } 

            if (!ORGC::local) { // ORGC scope
                break;
            }
            
            ORGC::postprocess(tempfile, final_file_path); // ORGC scope
            std::remove(tempfile.c_str());
            break; 
        } 
        if (!ORGC::local) { // ORGC scope
        struct stat buffer;
        if (stat(final_file_path.c_str(), &buffer) == 0) {
            std::remove(final_file_path.c_str());
        }
        
        std::string greference_path = reference_genome_path; // Renamed
        std::string gtarget_path = target_genome_path;       // Renamed
        std::string tempfile = final_folder + "/interim.txt";

        std::string reference_seq_content_global = ORGC::GreadrefSeq(greference_path); // ORGC scope
        std::string target_seq_content_global = ORGC::GreadtarSeq(gtarget_path, final_file_path); // ORGC scope
        if (reference_seq_content_global.empty()) {
            std::cerr << "Error: Reference sequence is empty for " << greference_path << std::endl;
            continue;
         }
         if (target_seq_content_global.empty()) {
            std::cerr << "Error: Target sequence is empty for " << gtarget_path << std::endl;
            continue;
         }

        if (stat(tempfile.c_str(), &buffer) == 0) {
            std::remove(tempfile.c_str());
        }
        
        std::string full_reference_seq = reference_seq_content_global; 
        std::string full_target_seq = target_seq_content_global;

        ORGC::reference = full_reference_seq;
        ORGC::target = full_target_seq;
        ORGC::target = full_target_seq; // Set ORGC::target for the context of Gmatch results

        std::vector<Position> list_global_matches = ORGC::Gmatch(full_reference_seq, full_target_seq, ORGC::kmer_length); // ORGC scope

        if (!ORGC::target.empty() && !list_global_matches.empty()) {
            ORGC::format_matches(list_global_matches, tempfile); // ORGC scope
            ORGC::postprocess(tempfile, final_file_path);      // ORGC scope
        } else {
            std::cerr << "Warning: target sequence or matches are empty for " << chrName[i] << std::endl;
        }
        
        std::remove(tempfile.c_str());
    }
    }

    ORGC::use7zip(final_folder); // ORGC scope
    std::cout << "All Done" << std::endl;
    long taskCPUTimeNano = ORGC::getCPUTime() - startCpuTimeNano; // ORGC scope
    std::cout << "Compressed time: " << (double) taskCPUTimeNano / 1000000000.0 << " seconds." << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;
    
    return 0;
}