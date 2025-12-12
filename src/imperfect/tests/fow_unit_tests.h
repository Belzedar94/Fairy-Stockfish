#ifndef FOW_UNIT_TESTS_H_INCLUDED
#define FOW_UNIT_TESTS_H_INCLUDED

#include <ostream>

namespace Stockfish {
namespace FogOfWar {

// Runs the C++ unit-style FoW suites and prints a compact report.
// Returns 0 on success, non-zero on any failure.
int run_fow_unit_suites(std::ostream& out);

// Runs the FoW micro-benchmark and comparison probes; returns non-zero on failures.
int run_fow_benchmark_suites(std::ostream& out);

// Runs sanitizer-friendly quick probes intended for ASan/TSan builds.
int run_fow_sanitizer_suites(std::ostream& out);

}  // namespace FogOfWar
}  // namespace Stockfish

#endif
