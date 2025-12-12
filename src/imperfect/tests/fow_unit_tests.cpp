#include "fow_unit_tests.h"

#include <algorithm>
#include <array>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "../Belief.h"
#include "../CFR.h"
#include "../Evaluator.h"
#include "../Selection.h"
#include "../Subgame.h"
#include "../Visibility.h"
#include "../Planner.h"
#include "../../uci.h"
#include "../../movegen.h"
#include "../../position.h"
#include "../../thread.h"
#include "../../variant.h"

namespace Stockfish {
namespace FogOfWar {

namespace {

struct TestCase {
    std::string name;
    bool passed = false;
    std::string detail;
};

struct TestSuiteResult {
    std::string name;
    std::vector<TestCase> cases;
};

Variant const* chess_variant() {
    return variants.find("chess")->second;
}

Position make_position(const std::string& fen) {
    Position pos;
    StateInfo st;
    pos.set(chess_variant(), fen, false, &st, nullptr);
    return pos;
}

TestCase visibility_castling_and_ep() {
    Position pos = make_position("r3k2r/8/8/8/8/8/8/R3K2R w KQkq d6 0 1");
    VisibilityInfo vi = compute_visibility(pos);

    const bool epVisible = is_visible(pos, SQ_D6, vi);
    const bool kingHomeVisible = is_visible(pos, SQ_E1, vi) && is_visible(pos, SQ_A1, vi);

    return {"visibility: castling + en-passant", epVisible && kingHomeVisible,
            epVisible ? "ep target surfaced" : "ep target missing"};
}

TestCase visibility_pawn_masking() {
    Position pos = make_position("8/8/8/8/8/3pP3/8/8 w - - 0 10");
    VisibilityInfo vi = compute_visibility(pos);

    const bool blockerHidden = !(vi.visible & SQ_D6) && (vi.visible & SQ_E6);
    return {"visibility: pawn blocker hidden", blockerHidden,
            blockerHidden ? "blocked pawn respected" : "blocker leaked"};
}

TestSuiteResult run_visibility_suite() {
    TestSuiteResult suite{"visibility"};
    suite.cases.push_back(visibility_castling_and_ep());
    suite.cases.push_back(visibility_pawn_masking());
    return suite;
}

ObservationHistory seed_observation_history(const Position& pos) {
    ObservationHistory hist;
    hist.add_observation(create_observation(pos));
    return hist;
}

TestCase belief_enumeration_cap() {
    Position pos = make_position("rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR b KQkq - 0 1");
    BeliefState belief;
    ObservationHistory history = seed_observation_history(pos);
    belief.rebuild_from_observations(history, pos);
    belief.compress(128);

    const bool nonEmpty = belief.size() > 0 && belief.size() <= 128;
    std::ostringstream detail;
    detail << "belief count=" << belief.size();
    return {"belief: enumeration + cap", nonEmpty, detail.str()};
}

TestCase belief_incremental_filter() {
    Position pos = make_position("8/2k5/8/8/8/8/2K5/8 w - - 0 1");
    BeliefState belief;
    ObservationHistory history = seed_observation_history(pos);
    belief.rebuild_from_observations(history, pos);
    const auto initial = belief.size();

    // Add a constrained observation (king visibility only)
    Observation tighter = create_observation(pos);
    tighter.visible = pos.pieces(WHITE, KING);
    history.add_observation(tighter);
    belief.update_incrementally(history, pos);

    const bool filtered = belief.size() <= initial;
    std::ostringstream detail;
    detail << "before=" << initial << " after=" << belief.size();
    return {"belief: incremental filter", filtered, detail.str()};
}

TestSuiteResult run_belief_suite() {
    TestSuiteResult suite{"belief"};
    suite.cases.push_back(belief_enumeration_cap());
    suite.cases.push_back(belief_incremental_filter());
    return suite;
}

TestCase selection_purification_support() {
    InfosetNode infoset;
    infoset.actions = { MOVE_NONE, MOVE_NONE, MOVE_NONE, MOVE_NONE };
    infoset.regrets = { 1.0f, 0.5f, -0.25f, 0.0f };
    infoset.qValues = { 0.25f, 0.1f, -0.5f, 0.0f };
    infoset.variances = { 0.01f, 0.01f, 0.01f, 0.01f };
    infoset.strategy = { 0.4f, 0.3f, 0.2f, 0.1f };

    ActionSelection selector;
    auto margins = selector.compute_margins(&infoset);
    auto purified = selector.purify_strategy(infoset.strategy, margins, false);

    size_t support = 0;
    for (float p : purified)
        if (p > 0.0f)
            ++support;

    return {"selection: purified support<=3", support <= 3, "support=" + std::to_string(support)};
}

TestCase selection_resolve_determinism() {
    InfosetNode infoset;
    infoset.actions = { MOVE_NONE, MOVE_NONE };
    infoset.regrets = { 0.0f, 1.0f };
    infoset.qValues = { 0.01f, 0.1f };
    infoset.variances = { 0.0f, 0.0f };
    infoset.strategy = { 0.5f, 0.5f };

    ActionSelection selector;
    selector.set_max_support(1);
    auto margins = selector.compute_margins(&infoset);
    auto purified = selector.purify_strategy(infoset.strategy, margins, true);

    const bool deterministic = purified[0] == 0.0f || purified[1] == 0.0f;
    return {"selection: resolve deterministic", deterministic, deterministic ? "collapsed" : "mixed"};
}

TestSuiteResult run_selection_suite() {
    TestSuiteResult suite{"selection"};
    suite.cases.push_back(selection_purification_support());
    suite.cases.push_back(selection_resolve_determinism());
    return suite;
}

TestCase kluss_sequence_id_stability() {
    std::vector<Move> sequence = { MOVE_NONE, make_move(SQ_E2, SQ_E4), make_move(SQ_E7, SQ_E5) };
    auto idA = compute_sequence_id_from_moves(sequence);
    auto idB = compute_sequence_id_from_moves(sequence);

    return {"kluss: sequence id stable", idA == idB, idA == idB ? "stable" : "drift"};
}

TestSuiteResult run_kluss_suite() {
    TestSuiteResult suite{"kluss"};
    suite.cases.push_back(kluss_sequence_id_stability());
    return suite;
}

TestCase cfr_reset_and_stop() {
    CFRSolver solver;
    solver.reset();
    const bool zeroed = solver.get_iterations() == 0;
    solver.stop();
    const bool stopped = !solver.is_running();
    return {"cfr: reset+stop", zeroed && stopped, stopped ? "stopped" : "running"};
}

TestSuiteResult run_cfr_suite() {
    TestSuiteResult suite{"cfr"};
    suite.cases.push_back(cfr_reset_and_stop());
    return suite;
}

TestSuiteResult run_benchmark_suite() {
    TestSuiteResult suite{"benchmarks"};

    Position pos = make_position("rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1");
    ObservationHistory history = seed_observation_history(pos);
    BeliefState belief;
    belief.rebuild_from_observations(history, pos);

    // quick one-ply evaluations for FoW vs baseline consistency
    float baseline = evaluate(pos);
    float fowLeaf = evaluate_belief_state(belief, pos.variant());

    suite.cases.push_back({"benchmark: evaluator parity", std::abs(baseline - fowLeaf) < 1.5f, "baseline=" + std::to_string(baseline)});
    return suite;
}

TestSuiteResult run_sanitizer_suite() {
    TestSuiteResult suite{"sanitizers"};
    // Minimal construction/destruction loops to flag leaks/races in ASan/TSan runs
    for (int i = 0; i < 3; ++i) {
        Position pos = make_position("8/8/8/8/8/8/8/8 w - - 0 1");
        ObservationHistory hist = seed_observation_history(pos);
        BeliefState belief;
        belief.rebuild_from_observations(hist, pos);
    }
    suite.cases.push_back({"asan/tsan: init/teardown", true, ""});
    return suite;
}

int emit_report(const std::vector<TestSuiteResult>& suites, std::ostream& out) {
    bool ok = true;
    for (const auto& suite : suites) {
        out << "[FoW] " << suite.name << "\n";
        for (const auto& tc : suite.cases) {
            out << "  - " << std::left << std::setw(32) << tc.name << (tc.passed ? "OK" : "FAIL");
            if (!tc.detail.empty())
                out << " :: " << tc.detail;
            out << "\n";
            ok &= tc.passed;
        }
    }
    return ok ? 0 : 1;
}

}  // namespace

int run_fow_unit_suites(std::ostream& out) {
    std::vector<TestSuiteResult> suites;
    suites.push_back(run_visibility_suite());
    suites.push_back(run_belief_suite());
    suites.push_back(run_cfr_suite());
    suites.push_back(run_selection_suite());
    suites.push_back(run_kluss_suite());
    return emit_report(suites, out);
}

int run_fow_benchmark_suites(std::ostream& out) {
    std::vector<TestSuiteResult> suites;
    suites.push_back(run_benchmark_suite());
    suites.push_back(run_selection_suite());
    return emit_report(suites, out);
}

int run_fow_sanitizer_suites(std::ostream& out) {
    std::vector<TestSuiteResult> suites;
    suites.push_back(run_sanitizer_suite());
    suites.push_back(run_belief_suite());
    return emit_report(suites, out);
}

}  // namespace FogOfWar
}  // namespace Stockfish
