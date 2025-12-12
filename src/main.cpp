/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2022 The Stockfish developers (see AUTHORS file)

  Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Stockfish is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>

#include "bitboard.h"
#include "endgame.h"
#include "position.h"
#include "psqt.h"
#include "search.h"
#include "syzygy/tbprobe.h"
#include "thread.h"
#include "tt.h"
#include "uci.h"
#include "imperfect/tests/fow_unit_tests.h"

#include "piece.h"
#include "variant.h"
#include "xboard.h"


using namespace Stockfish;

int main(int argc, char* argv[]) {

  std::cout << engine_info() << std::endl;

  pieceMap.init();
  variants.init();
  CommandLine::init(argc, argv);
  UCI::init(Options);
  Tune::init();
  PSQT::init(variants.find(Options["UCI_Variant"])->second);
  Bitboards::init();
  Position::init();
  Bitbases::init();
  Endgames::init();
  Threads.set(size_t(Options["Threads"]));
  Search::clear(); // After threads are up
  Eval::NNUE::init();

  bool runFoWUnits = false;
  bool runFoWSanitizers = false;
  bool runFoWBenchmarks = false;
  for (int i = 1; i < argc; ++i)
  {
      std::string arg(argv[i]);
      if (arg == "--fow-unittests")
          runFoWUnits = true;
      else if (arg == "--fow-sanitizers")
          runFoWSanitizers = true;
      else if (arg == "--fow-benchmarks")
          runFoWBenchmarks = true;
  }

  if (runFoWUnits || runFoWSanitizers || runFoWBenchmarks)
  {
      int status = 0;
      if (runFoWUnits)
          status |= FogOfWar::run_fow_unit_suites(std::cout);
      if (runFoWSanitizers)
          status |= FogOfWar::run_fow_sanitizer_suites(std::cout);
      if (runFoWBenchmarks)
          status |= FogOfWar::run_fow_benchmark_suites(std::cout);
      Threads.set(0);
      variants.clear_all();
      pieceMap.clear_all();
      delete XBoard::stateMachine;
      return status;
  }

  UCI::loop(argc, argv);

  Threads.set(0);
  variants.clear_all();
  pieceMap.clear_all();
  delete XBoard::stateMachine;
  return 0;
}
