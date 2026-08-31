/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2026 The Stockfish developers (see AUTHORS file)

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

#include "engine.h"

#include <algorithm>
#include <cassert>
#include <charconv>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <deque>
#include <iosfwd>
#include <memory>
#include <ostream>
#include <sstream>
#include <string_view>
#include <utility>
#include <vector>

#include "evaluate.h"
#include "koth.h"
#include "misc.h"
#include "nnue/network.h"
#include "nnue/nnue_common.h"
#include "numa.h"
#include "perft.h"
#include "position.h"
#include "search.h"
#include "shm.h"
#include "types.h"
#include "uci.h"
#include "ucioption.h"

namespace Stockfish {

namespace NN = Eval::NNUE;

int MaxThreads = std::max(1024, 4 * int(get_hardware_concurrency()));

// The default configuration will attempt to group L3 domains up to 32 threads.
// This size was found to be a good balance between the Elo gain of increased
// history sharing and the speed loss from more cross-cache accesses (see
// PR#6526). The user can always explicitly override this behavior.
constexpr NumaAutoPolicy DefaultNumaPolicy = BundledL3Policy{32};

namespace {

bool parse_decimal(std::string_view token, int minimum, int maximum) {
    if (token.empty())
        return false;

    int value               = 0;
    const auto [end, error] = std::from_chars(token.data(), token.data() + token.size(), value);
    return error == std::errc{} && end == token.data() + token.size() && value >= minimum
        && value <= maximum;
}

std::optional<PositionSetError> canonical_koth_fen(const std::string& input,
                                                   std::string&       canonical) {
    std::istringstream       stream(input);
    std::vector<std::string> fields;
    std::string              field;

    while (stream >> field)
        fields.push_back(field);

    if (fields.size() != 6)
        return PositionSetError("code=INVALID_FEN_FIELD_COUNT expected=6 actual="
                                + std::to_string(fields.size()));

    if (fields[1] != "w" && fields[1] != "b")
        return PositionSetError("code=INVALID_FEN_SIDE");

    if (fields[2] != "-")
    {
        constexpr std::string_view order    = "KQkq";
        std::size_t                previous = 0;
        bool                       first    = true;
        for (char right : fields[2])
        {
            const auto index = order.find(right);
            if (index == std::string_view::npos || (!first && index <= previous))
                return PositionSetError("code=INVALID_FEN_CASTLING_RIGHTS");
            first    = false;
            previous = index;
        }
    }

    if (fields[3] != "-"
        && (fields[3].size() != 2 || fields[3][0] < 'a' || fields[3][0] > 'h'
            || (fields[3][1] != '3' && fields[3][1] != '6')))
        return PositionSetError("code=INVALID_FEN_EP_FIELD");

    if (!parse_decimal(fields[4], 0, 32767))
        return PositionSetError("code=INVALID_FEN_HALFMOVE");

    if (!parse_decimal(fields[5], 1, 100000))
        return PositionSetError("code=INVALID_FEN_FULLMOVE");

    canonical = fields[0] + " " + fields[1] + " " + fields[2] + " " + fields[3] + " " + fields[4]
              + " " + fields[5];
    return std::nullopt;
}

}  // namespace

Engine::Engine(std::optional<std::filesystem::path> path) :
    binaryDirectory(path ? CommandLine::get_binary_directory(*path) : std::filesystem::path{}),
    numaContext(NumaConfig::from_system(DefaultNumaPolicy)),
    states(new std::deque<StateInfo>(1)),
    threads(),
    networkFile{std::nullopt, ""},
    network(numaContext, get_default_network()) {

    pos.set(StartFEN, false, &states->back());

    options.add(  //
      "Debug Log File", Option("", [](const Option& o) {
          start_logger(path_from_utf8(std::string(o)));
          return std::nullopt;
      }));

    options.add(  //
      "NumaPolicy", Option("auto", [this](const Option& o) {
          if (!set_numa_config_from_option(o))
              return "NumaPolicy: invalid value '" + std::string(o) + "', keeping previous config.";
          return numa_config_information_as_string() + "\n"
               + thread_allocation_information_as_string();
      }));

    options.add(  //
      "Threads", Option(1, 1, MaxThreads, [this](const Option&) {
          resize_threads();
          return thread_allocation_information_as_string();
      }));

    options.add(  //
      "Hash", Option(16, 1, MaxHashMB, [this](const Option& o) {
          set_tt_size(o);
          return std::nullopt;
      }));

    options.add(  //
      "Clear Hash", Option([this](const Option&) {
          search_clear();
          return std::nullopt;
      }));

    options.add(  //
      "Ponder", Option(false));

    options.add(  //
      "MultiPV", Option(1, 1, MAX_MOVES));

    options.add("Skill Level", Option(20, 0, 20));

    options.add("Move Overhead", Option(10, 0, 5000));

    options.add("nodestime", Option(0, 0, 10000));

    options.add("EvalFile", Option(EvalFileDefaultName, [this](const Option& o) {
                    load_network(path_from_utf8(std::string(o)));
                    return std::nullopt;
                }));

    options.add("UCI_Variant", Option("var kingofthehill", "kingofthehill"));

    options.add("UCI_LimitStrength", Option(false));

    options.add("UCI_Elo",
                Option(Stockfish::Search::Skill::LowestElo, Stockfish::Search::Skill::LowestElo,
                       Stockfish::Search::Skill::HighestElo));

    threads.clear();
    threads.ensure_network_replicated();
    resize_threads();
}

std::variant<u64, PositionSetError> Engine::perft(Depth depth, bool gameDomain) {
    StateInfo state;
    Position  copy;
    copy.clone_from(pos, state);
    return Benchmark::perft<true>(copy, depth, gameDomain);
}

void Engine::go(Search::LimitsType& limits) {
    assert(limits.perft == 0);

    const auto status = Koth::classify(pos);
    if (status.terminal())
    {
        if (onVerifyNetwork)
            onVerifyNetwork("koth " + Koth::serialize(status));

        if (updateContext.onBestmove)
            updateContext.onBestmove(UCIEngine::move(Move::none()), "");
        return;
    }

    if (limits.mate)
    {
        if (onVerifyNetwork)
            onVerifyNetwork("error code=UNSUPPORTED_KOTH_GO_LIMIT limit=mate");

        if (updateContext.onBestmove)
            updateContext.onBestmove(UCIEngine::move(Move::none()), "");
        return;
    }

    verify_network();
    threads.start_thinking(options, pos, states, limits);
}
void Engine::stop() { threads.stop = true; }

void Engine::search_clear() {
    wait_for_search_finished();

    tt.clear(threads);
    threads.clear();
}

void Engine::set_on_update_no_moves(std::function<void(const Engine::InfoShort&)>&& f) {
    updateContext.onUpdateNoMoves = std::move(f);
}

void Engine::set_on_update_full(std::function<void(const Engine::InfoFull&)>&& f) {
    updateContext.onUpdateFull = std::move(f);
}

void Engine::set_on_iter(std::function<void(const Engine::InfoIter&)>&& f) {
    updateContext.onIter = std::move(f);
}

void Engine::set_on_bestmove(std::function<void(std::string_view, std::string_view)>&& f) {
    updateContext.onBestmove = std::move(f);
}

void Engine::set_on_start(std::function<void()>&& f) { updateContext.onStart = std::move(f); }

void Engine::set_on_verify_network(std::function<void(std::string_view)>&& f) {
    onVerifyNetwork = std::move(f);
}

void Engine::wait_for_search_finished() { threads.main_thread()->wait_for_search_finished(); }

std::optional<PositionSetError> Engine::set_position(const std::string&              fen,
                                                     const std::vector<std::string>& moves) {
    std::string canonicalFen;
    if (auto err = canonical_koth_fen(fen, canonicalFen))
        return err;

    auto     candidateStates = StateListPtr(new std::deque<StateInfo>(1));
    Position candidate;
    auto     err = candidate.set(canonicalFen, false, &candidateStates->back());
    if (err.has_value())
        return err;

    if (candidate.fen() != canonicalFen)
        return PositionSetError("code=NONCANONICAL_OR_INCONSISTENT_FEN");

    if (Koth::any_king_on_hill(candidate))
        return PositionSetError("code=AMBIGUOUS_GOAL_ROOT");

    if (const auto rootStatus =
          Koth::classify(candidate, Koth::AdjudicationContext::RAW_ROOT_VALIDATION);
        rootStatus.terminal())
        return PositionSetError("code=TERMINAL_RAW_ROOT primary="
                                + Koth::primary_reason_name(rootStatus.primary));

    for (std::size_t index = 0; index < moves.size(); ++index)
    {
        const auto& move = moves[index];
        if (!Koth::is_canonical_uci_move(move))
            return PositionSetError(
              "code=NONCANONICAL_TRAJECTORY_MOVE index=" + std::to_string(index) + " move=" + move);

        if (Koth::classify(candidate).terminal())
            return PositionSetError(
              "code=POSTTERMINAL_TRAJECTORY_MOVE index=" + std::to_string(index) + " move=" + move);

        auto m = UCIEngine::to_move(candidate, move);

        if (m == Move::none())
            return PositionSetError("code=ILLEGAL_TRAJECTORY_MOVE index=" + std::to_string(index)
                                    + " move=" + move);

        candidateStates->emplace_back();
        candidate.do_move(m, candidateStates->back());

        if (Koth::classify(candidate).terminal() && index + 1 != moves.size())
            return PositionSetError("code=POSTTERMINAL_TRAJECTORY_TAIL terminal_index="
                                    + std::to_string(index));
    }

    states = std::move(candidateStates);
    pos.clone_from(candidate, states->back());

    return std::nullopt;
}

// modifiers

bool Engine::set_numa_config_from_option(const std::string& o) {
    if (o == "auto" || o == "system")
    {
        numaContext.set_numa_config(NumaConfig::from_system(DefaultNumaPolicy));
    }
    else if (o == "hardware")
    {
        // Don't respect affinity set in the system.
        numaContext.set_numa_config(NumaConfig::from_system(DefaultNumaPolicy, false));
    }
    else if (o == "none")
    {
        numaContext.set_numa_config(NumaConfig{});
    }
    else
    {
        auto parsed = NumaConfig::from_string(o);
        if (!parsed.has_value())
            return false;
        numaContext.set_numa_config(std::move(*parsed));
    }

    // Force reallocation of threads in case affinities need to change.
    resize_threads();
    threads.ensure_network_replicated();
    return true;
}

void Engine::resize_threads() {
    threads.wait_for_search_finished();
    threads.set(numaContext.get_numa_config(), {options, threads, tt, sharedHists, network},
                updateContext);

    // Reallocate the hash with the new threadpool size
    set_tt_size(options["Hash"]);
    threads.ensure_network_replicated();
}

void Engine::set_tt_size(usize mb) {
    wait_for_search_finished();
    tt.resize(mb, threads);
}

void Engine::set_ponderhit(bool b) { threads.main_manager()->ponder = b; }

// network related

void Engine::verify_network() {
    const auto requested = path_from_utf8(std::string(options["EvalFile"]));
    if (!network->is_initialized())
        load_network(requested);
    network->verify(onVerifyNetwork, networkFile, requested);
}

std::unique_ptr<Eval::NNUE::Network> Engine::get_default_network() {
    return std::make_unique<NN::Network>();
}

void Engine::load_network(const std::filesystem::path& file) {
    auto candidate     = std::make_unique<NN::Network>();
    auto candidateFile = networkFile;
    if (auto error = candidate->load(binaryDirectory, file, candidateFile))
    {
        const std::string message = "error " + *error;
        if (onVerifyNetwork)
            onVerifyNetwork(message);
        else
            sync_cout << "info string " << message << sync_endl;
        std::exit(EXIT_FAILURE);
    }

    network     = std::move(candidate);
    networkFile = std::move(candidateFile);
    threads.clear();
    threads.ensure_network_replicated();
}

void Engine::save_network(const std::optional<std::filesystem::path>& file) {
    network.modify_and_replicate(
      [&file, this](NN::Network& network_) { network_.save(networkFile, file); });
}

// utility functions

void Engine::trace_eval() {
    verify_network();
    if (const auto status = Koth::classify(pos); status.terminal())
    {
        if (onVerifyNetwork)
            onVerifyNetwork("koth " + Koth::serialize(status));
        return;
    }
    sync_cout << Eval::trace(pos, *network) << sync_endl;
}

const OptionsMap& Engine::get_options() const { return options; }
OptionsMap&       Engine::get_options() { return options; }

std::string Engine::fen() const { return pos.fen(); }

std::string Engine::koth_status() const {
    std::ostringstream out;
    const auto&        transition = pos.transition();

    out << Koth::serialize(Koth::classify(pos)) << " repetition_count=" << pos.repetition_count()
        << " transition=";

    switch (transition.kind)
    {
    case Koth::TransitionKind::ROOT :
        out << "ROOT";
        break;
    case Koth::TransitionKind::ACCEPTED_MOVE :
        out << "ACCEPTED_MOVE";
        break;
    case Koth::TransitionKind::SEARCH_NULL :
        out << "SEARCH_NULL";
        break;
    }

    out << " king_entered_hill=" << (transition.kingEnteredHill ? "true" : "false") << " fen_ep="
        << (pos.fen_ep_square() == SQ_NONE ? "-" : UCIEngine::square(pos.fen_ep_square()))
        << " legal_ep=" << (pos.ep_square() == SQ_NONE ? "-" : UCIEngine::square(pos.ep_square()));
    return out.str();
}

std::string Engine::koth_moves() const {
    const auto         physical = Koth::physical_moves(pos);
    const auto         game     = Koth::game_moves(pos);
    std::ostringstream out;

    out << "physical_count=" << physical.size() << " game_count=" << game.size() << " physical=";
    for (Move move : physical)
        out << UCIEngine::move(move) << ',';

    out << " game=";
    for (Move move : game)
        out << UCIEngine::move(move) << ',';

    return out.str();
}

std::string Engine::koth_selftest() const { return Koth::state_selftest(); }

std::string Engine::koth_network_status() const { return network->status(networkFile); }

std::string Engine::koth_network_eval() {
    verify_network();
    if (const auto status = Koth::classify(pos); status.terminal())
        return "terminal=true " + Koth::serialize(status);

    const auto raw = network->evaluate_raw(pos);
    const bool lazy =
      std::abs(raw.psqt) > NN::Network::LegacyLazyThreshold * NN::Network::LegacyOutputScale;
    const i32 psqt = raw.psqt / NN::Network::LegacyOutputScale;
    const i32 total =
      lazy ? psqt : i32((i64(raw.psqt) + raw.positional) / NN::Network::LegacyOutputScale);
    const i32 positional = total - psqt;
    return "terminal=false bucket=" + std::to_string(raw.bucket)
         + " lazy=" + (lazy ? "true" : "false") + " psqt_raw=" + std::to_string(raw.psqt)
         + " positional_raw=" + std::to_string(raw.positional) + " psqt=" + std::to_string(psqt)
         + " positional=" + std::to_string(positional) + " total=" + std::to_string(total);
}

std::optional<PositionSetError> Engine::flip() {
    return PositionSetError("code=UNSUPPORTED_KOTH_COMMAND command=flip");
}

std::string Engine::visualize() const {
    std::stringstream ss;
    ss << pos;
    return ss.str();
}

int Engine::get_hashfull(int maxAge) const { return tt.hashfull(maxAge); }

std::vector<std::pair<usize, usize>> Engine::get_bound_thread_count_by_numa_node() const {
    auto                                 counts = threads.get_bound_thread_count_by_numa_node();
    const NumaConfig&                    cfg    = numaContext.get_numa_config();
    std::vector<std::pair<usize, usize>> ratios;
    NumaIndex                            n = 0;
    for (; n < counts.size(); ++n)
        ratios.emplace_back(counts[n], cfg.num_cpus_in_numa_node(n));
    if (!counts.empty())
        for (; n < cfg.num_numa_nodes(); ++n)
            ratios.emplace_back(0, cfg.num_cpus_in_numa_node(n));
    return ratios;
}

std::string Engine::get_numa_config_as_string() const {
    return numaContext.get_numa_config().to_string();
}

std::string Engine::numa_config_information_as_string() const {
    auto cfgStr = get_numa_config_as_string();
    return "Available processors: " + cfgStr;
}

std::string Engine::thread_binding_information_as_string() const {
    auto              boundThreadsByNode = get_bound_thread_count_by_numa_node();
    std::stringstream ss;
    if (boundThreadsByNode.empty())
        return ss.str();

    bool isFirst = true;

    for (auto&& [current, total] : boundThreadsByNode)
    {
        if (!isFirst)
            ss << ":";
        ss << current << "/" << total;
        isFirst = false;
    }

    return ss.str();
}

std::string Engine::thread_allocation_information_as_string() const {
    std::stringstream ss;

    usize threadsSize = threads.size();
    ss << "Using " << threadsSize << (threadsSize > 1 ? " threads" : " thread");

    auto boundThreadsByNodeStr = thread_binding_information_as_string();
    if (boundThreadsByNodeStr.empty())
        return ss.str();

    ss << " with NUMA node thread binding: ";
    ss << boundThreadsByNodeStr;

    return ss.str();
}
}
