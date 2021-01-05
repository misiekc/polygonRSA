#include "PackingGenerator.h"
#include "shape/ShapeFactory.h"
#include "Utils.h"
#include "ProgramArguments.h"


#include <iostream>
#include <unistd.h>
#include <sys/wait.h>
#include <cstring>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <memory>


//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0
void process_mem_usage(double &vm_usage, double &resident_set) {
    using std::ios_base;
    using std::ifstream;
    using std::string;

    vm_usage = 0.0;
    resident_set = 0.0;

    // 'file' stat seems to give the most reliable results
    //
    ifstream stat_stream("/proc/self/stat", ios_base::in);

    // dummy vars for leading entries in stat that we don't care about
    //
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // the two fields we want
    //
    unsigned long vsize;
    long rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
                >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
                >> utime >> stime >> cutime >> cstime >> priority >> nice
                >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

    stat_stream.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE); // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0 / 1024.0;
    resident_set = rss * page_size_kb / 1024.0 / 1024.0;
}

void makeDatFileForPackingsInDirectory(const Parameters *params, const std::string &dirName) {
    std::string sFile = params->getPackingSignature() + ".dat";
    std::ofstream dataFile(sFile);
    if (!dataFile)
        die("Cannot open file " + sFile + " to store packing info");
    dataFile.precision(std::numeric_limits<double>::digits10 + 1);

    auto filenames = PackingGenerator::findPackingsInDir(dirName);
    for (const auto &filename : filenames) {
        int no1 = lastIndexOf(filename, '_');
        int no2 = lastIndexOf(filename, '.');
        std::string seed = filename.substr(no1 + 1, no2 - no1 - 1);

        Packing packing;
        packing.restore(filename);

        dataFile << seed << "\t" << packing.back()->no << "\t" << packing.back()->time << std::endl;
        std::cout << ".";
        std::cout.flush();
    }
    std::cout << std::endl;
}

void runSingleSimulation(int seed, Parameters *params, std::ofstream &dataFile) {
    double vm, rss;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    PackingGenerator pg(seed, params);
    pg.run();
    const Packing &packing = pg.getPacking();

    if (params->storePackings) {
        std::string sPackingFile = params->getPackingSignature() + "_" + std::to_string(seed) + ".bin";
        pg.getPacking().store(sPackingFile);
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    process_mem_usage(vm, rss);
    dataFile << seed << "\t" << packing.size() << "\t" << packing.back()->time << std::endl;
    dataFile.flush();
    std::cout << "T:" << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "; VM: " << vm
              << "; RSS: " << rss << std::endl;
}

int simulate(const ProgramArguments &arguments) {
    Parameters params = arguments.getParameters();

    std::string sFile = params.getPackingSignature() + ".dat";
    std::ofstream file(sFile);
    if (!file)
        die("Cannot open file " + sFile + " to store packing info");
    file.precision(std::numeric_limits<double>::digits10 + 1);

    for (std::size_t i = 0; i < params.collectors; i++) {
        runSingleSimulation(static_cast<int>(params.from + i), &params, file);
    }
    file.close();
    return 1;
}


void test(const ProgramArguments &arguments) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() < 2)
        die(arguments.formatUsage("<file in> <max time>"));

    Packing packing;
    packing.restore(positionalArguments[0]);
    PackingGenerator pg(1, &arguments.getParameters());
    pg.testPacking(packing, std::stod(positionalArguments[1]));
}

void debug(const ProgramArguments &arguments) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.empty())
        die(arguments.formatUsage("<packing generator file>"));

    Parameters params = arguments.getParameters();
    const char *cfile = positionalArguments[0].c_str();
    std::string filename(cfile);
    char buf[20];
    sprintf(buf, "%.0f", pow(params.surfaceSize, 2));
    std::string size(buf);

    std::ifstream file(filename, std::ios::binary);
    if (!file)
        die("Cannot open file " + filename + " to restore packing generator");

    PackingGenerator pg(0, &params);
    pg.restore(file);
    pg.run();
}

void dat(const ProgramArguments &arguments) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.empty())
        die(arguments.formatUsage("<directory>"));

    std::string directory = positionalArguments[0];
    makeDatFileForPackingsInDirectory(&arguments.getParameters(), directory);
}

void povray(const ProgramArguments &arguments) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.empty())
        die(arguments.formatUsage("<file in>"));

    std::string file(positionalArguments[0]);
    Packing packing;
    packing.restore(file);
    PackingGenerator::toPovray(packing, arguments.getParameters().surfaceSize, nullptr, file + ".pov");
}

void wolfram(const ProgramArguments &arguments) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.empty())
        die(arguments.formatUsage("<file in>"));

    std::string file(positionalArguments[0]);
    Packing packing;
    packing.restore(file);
    PackingGenerator::toWolfram(packing, nullptr, file + ".nb");
}

void bc_expand(const ProgramArguments &arguments) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.empty())
        die(arguments.formatUsage("<file in> (file out = file in)"));

    std::string fileIn(positionalArguments[0]);
    std::string fileOut = (positionalArguments.size() == 1 ? fileIn : positionalArguments[1]);
    Packing packing;
    packing.restore(fileIn);
    packing.expandOnPBC(arguments.getParameters().surfaceSize, 0.1);
    packing.store(fileOut);
}

int main(int argc, char **argv) {
    std::unique_ptr<ProgramArguments> arguments;
    try {
        arguments = std::make_unique<ProgramArguments>(argc, argv);
    } catch (InvalidArgumentsException &e) {
        die(e.what());
    }

    Parameters params = arguments->getParameters();
    ShapeFactory::initShapeClass(params.particleType, params.particleAttributes);
    
    std::string mode = arguments->getMode();
    if (mode == "simulate")
        simulate(*arguments);
    else if (mode == "test")
        test(*arguments);
    else if (mode == "debug")
        debug(*arguments);
    else if (mode == "dat")
        dat(*arguments);
    else if (mode == "povray")
        povray(*arguments);
    else if (mode == "wolfram")
        wolfram(*arguments);
    else if (mode == "bc_expand")
        bc_expand(*arguments);
    else
        die("Unknown mode: " + mode);

    return EXIT_SUCCESS;
}
