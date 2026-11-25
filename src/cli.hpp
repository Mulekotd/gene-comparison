#pragma once

#include <cstddef>
#include <filesystem>
#include <iosfwd>
#include <string>

struct ProgramOptions {
    std::string gene = "recA";
    std::filesystem::path input_directory;
    std::filesystem::path output_directory;
    std::size_t requested_threads = 1;
    bool protein_mode = true;
    bool benchmark = false;
    bool help = false;
};

ProgramOptions parse_arguments(int argc, char* argv[], std::size_t available_threads);
void print_usage(std::ostream& output);
