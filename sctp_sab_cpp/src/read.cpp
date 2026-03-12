#include "read.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <regex>

// Trim whitespace from both ends
static std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

// Split string by delimiter
static std::vector<std::string> split(const std::string& s, char delim) {
    std::vector<std::string> tokens;
    std::istringstream iss(s);
    std::string token;
    while (std::getline(iss, token, delim)) {
        std::string t = trim(token);
        if (!t.empty()) tokens.push_back(t);
    }
    return tokens;
}

void read_net(const std::string& network, const std::string& base_dir,
              LinkMapD& Capacity, LinkMapD& Length, LinkMapD& Fftime,
              AdjList& Innode, AdjList& Outnode, std::vector<Link>& Link_list,
              int& node_count, int& link_count) {

    std::string filepath = base_dir + "/network/" + network + "/" + network + "_net.txt";
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(file, line)) {
        if (!trim(line).empty()) {
            lines.push_back(line);
        }
    }
    file.close();

    // Parse header: extract node_count and link_count
    node_count = 0;
    link_count = 0;
    int start_line = -1;
    for (size_t i = 0; i < lines.size(); ++i) {
        if (lines[i].find("NUMBER OF NODES") != std::string::npos) {
            // Extract last number from the line
            auto parts = split(lines[i], ' ');
            node_count = std::stoi(parts.back());
        }
        if (lines[i].find("NUMBER OF LINKS") != std::string::npos) {
            auto parts = split(lines[i], ' ');
            link_count = std::stoi(parts.back());
        }
        if (lines[i].find('~') != std::string::npos) {
            start_line = static_cast<int>(i) + 1;
            break;
        }
    }

    if (start_line < 0) {
        throw std::runtime_error("Could not find data start marker '~' in " + filepath);
    }

    // Initialize adjacency lists
    for (int i = 1; i <= node_count; ++i) {
        Innode[i] = {};
        Outnode[i] = {};
    }

    // Parse link data
    for (int i = start_line; i < start_line + link_count && i < (int)lines.size(); ++i) {
        std::string l = lines[i];
        // Remove trailing semicolon and whitespace
        while (!l.empty() && (l.back() == ';' || l.back() == '\t' || l.back() == ' ' || l.back() == '\r')) {
            l.pop_back();
        }
        // Remove leading whitespace
        l = trim(l);

        auto tokens = split(l, '\t');
        if (tokens.size() < 5) {
            // Try space-separated
            tokens = split(l, ' ');
        }
        if (tokens.size() < 5) continue;

        int ii = std::stoi(tokens[0]);
        int jj = std::stoi(tokens[1]);
        double cap = std::stod(tokens[2]);
        double len = std::stod(tokens[3]);
        double fft = std::stod(tokens[4]) / 60.0;  // Convert to hours (same as Python)

        Link lk = {ii, jj};
        Capacity[lk] = cap;
        Length[lk] = len;
        Fftime[lk] = fft;
        Innode[jj].push_back(ii);
        Outnode[ii].push_back(jj);
        Link_list.push_back(lk);
    }
}

void read_trp(const std::string& network, const std::string& base_dir,
              OdTree& Odtree, DemandMap& Demand, double& flow_count) {

    std::string filepath = base_dir + "/network/" + network + "/" + network + "_trp.txt";
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(file, line)) {
        if (!trim(line).empty()) {
            lines.push_back(line);
        }
    }
    file.close();

    // Find "Origin" lines
    std::vector<int> origin_lines;
    for (size_t i = 0; i < lines.size(); ++i) {
        if (lines[i].find("Origin") != std::string::npos) {
            origin_lines.push_back(static_cast<int>(i));
        }
    }

    flow_count = 0.0;

    for (size_t oi = 0; oi < origin_lines.size(); ++oi) {
        int origin_idx = origin_lines[oi];
        int next_idx = (oi + 1 < origin_lines.size()) ? origin_lines[oi + 1] : (int)lines.size();

        // Parse origin ID
        std::string origin_line = lines[origin_idx];
        // Remove "Origin" and extract the number
        std::regex num_re("[0-9]+");
        std::smatch match;
        int o_id = 0;
        // Find the last number on the Origin line
        std::string temp = origin_line;
        while (std::regex_search(temp, match, num_re)) {
            o_id = std::stoi(match[0]);
            temp = match.suffix().str();
        }

        std::vector<int> d_list;

        // Parse demand entries from lines after Origin until next Origin
        for (int li = origin_idx + 1; li < next_idx; ++li) {
            std::string dline = lines[li];
            // Remove spaces, split by semicolon to get "d:demand" pairs
            dline.erase(std::remove(dline.begin(), dline.end(), ' '), dline.end());

            // Remove trailing semicolon
            while (!dline.empty() && dline.back() == ';') dline.pop_back();

            auto pairs = split(dline, ';');
            for (auto& pair : pairs) {
                auto kv = split(pair, ':');
                if (kv.size() >= 2) {
                    int d_id = std::stoi(kv[0]);
                    double demand = std::stod(kv[1]);
                    if (d_id != o_id && std::abs(demand) > 1e-8) {
                        d_list.push_back(d_id);
                        Demand[{o_id, d_id}] = demand;
                        flow_count += demand;
                    }
                }
            }
        }

        Odtree[o_id] = d_list;
    }
}

NetworkData read_notoll(const std::string& network, const std::string& base_dir) {
    NetworkData data;

    read_net(network, base_dir,
             data.Capacity, data.Length, data.Fftime,
             data.Innode, data.Outnode, data.Link_list,
             data.node_count, data.link_count);

    read_trp(network, base_dir,
             data.Odtree, data.Demand, data.flow_count);

    return data;
}
