/*
  KOTH-Stockfish, a UCI chess playing engine derived from Stockfish
  Copyright (C) 2004-2026 The Stockfish developers (see AUTHORS file)

  KOTH-Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  KOTH-Stockfish is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef NETWORK_H_INCLUDED
#define NETWORK_H_INCLUDED

#include <array>
#include <filesystem>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <vector>

#include "../misc.h"
#include "../types.h"
#include "nnue_common.h"
#include "nnue_misc.h"

namespace Stockfish {
class Position;
}

namespace Stockfish::Eval::NNUE {

class AccumulatorStack;
struct AccumulatorCaches;

using NetworkOutput = std::tuple<Value, Value>;

struct LegacyRawOutput {
    i32   psqt;
    i32   positional;
    usize bucket;
};

// Compatibility evaluator for the exact owner-authored legacy network. Its
// physical feature geometry and quantized container are the official Stockfish
// HalfKAv2 architecture introduced by upstream commit e8d64af1. The project
// deliberately evaluates it through a scalar reference lane first; search is a
// separate gate.
//
// The object remains in-line and trivially copyable because Stockfish places
// Network instances in content-addressed shared memory.
class Network {
   public:
    static constexpr usize LegacyFileBytes       = 47721371;
    static constexpr u32   LegacyVersion         = 0x7AF32F20u;
    static constexpr u32   LegacyArchitecture    = 0x3C103E72u;
    static constexpr u32   LegacyFeatureHash     = 0x5F2348B8u;
    static constexpr u32   LegacyLayerHash       = 0x633376CAu;
    static constexpr usize LegacyInputDimensions = 45056;
    static constexpr usize LegacyL1              = 512;
    static constexpr usize LegacyLayerStacks     = 8;
    static constexpr usize LegacyPsqtBuckets     = 8;
    static constexpr i32   LegacyOutputScale     = 16;
    static constexpr i32   LegacyLazyThreshold   = 1400;

    static constexpr std::string_view LegacyCanonicalName = "kingofthehill-978b86d0e6a4.nnue";
    static constexpr std::string_view LegacyAliasName     = "KOTH_v1.nnue";
    static constexpr std::string_view LegacySha256 =
      "978B86D0E6A45E05F9F1375DCED129CEA0ACEA13041EA65960691632EDC47AF7";

    Network();

    Network(const Network&) = default;
    Network(Network&&)      = default;

    Network& operator=(const Network&) = default;
    Network& operator=(Network&&)      = default;

    // A failed load leaves both this object and evalFile unchanged.
    std::optional<std::string> load(const std::filesystem::path& rootDirectory,
                                    std::filesystem::path        evalfilePath,
                                    EvalFile&                    evalFile);

    bool save(const EvalFile& evalFile, const std::optional<std::filesystem::path>& filename) const;

    usize get_content_hash() const;
    bool  is_initialized() const { return initialized; }

    NetworkOutput evaluate(const Position&    pos,
                           AccumulatorStack&  accumulatorStack,
                           AccumulatorCaches& cache) const;

    LegacyRawOutput evaluate_raw(const Position& pos) const;

    void verify(const std::function<void(std::string_view)>& f,
                const EvalFile&                              evalFile,
                std::filesystem::path                        evalfilePath) const;

    std::string status(const EvalFile& evalFile) const;

    NnueEvalTrace trace_evaluate(const Position&    pos,
                                 AccumulatorStack&  accumulatorStack,
                                 AccumulatorCaches& cache) const;

   private:
    static constexpr usize LegacyNetworkInput = LegacyL1 * 2;
    static constexpr usize LegacyFc0Outputs   = 16;
    static constexpr usize LegacyFc1Outputs   = 32;
    static constexpr usize LegacyFc1PaddedIn  = 32;

    struct alignas(CacheLineSize) LegacyFeatureTransformer {
        std::array<i16, LegacyL1>                                  biases;
        std::array<i16, LegacyInputDimensions * LegacyL1>          weights;
        std::array<i32, LegacyInputDimensions * LegacyPsqtBuckets> psqtWeights;
    };

    struct alignas(CacheLineSize) LegacyLayerStack {
        std::array<i32, LegacyFc0Outputs>                     fc0Biases;
        std::array<i8, LegacyFc0Outputs * LegacyNetworkInput> fc0Weights;
        std::array<i32, LegacyFc1Outputs>                     fc1Biases;
        std::array<i8, LegacyFc1Outputs * LegacyFc1PaddedIn>  fc1Weights;
        std::array<i32, 1>                                    outputBias;
        std::array<i8, LegacyFc1Outputs>                      outputWeights;
    };

    std::optional<std::string> load_single_file(const std::filesystem::path& file,
                                                std::string&                 description);
    std::optional<std::string> parse_exact_bytes(const std::vector<u8>& bytes,
                                                 std::string&           description);
    LegacyRawOutput            evaluate_raw_bucket(const Position& pos, usize bucket) const;

    LegacyFeatureTransformer                        featureTransformer;
    std::array<LegacyLayerStack, LegacyLayerStacks> layerStacks;
    bool                                            initialized;
};

static_assert(std::is_trivially_destructible_v<Network>);
static_assert(std::is_trivially_move_constructible_v<Network>);
static_assert(std::is_trivially_copy_constructible_v<Network>);

}  // namespace Stockfish::Eval::NNUE

template<>
struct std::hash<Stockfish::Eval::NNUE::Network> {
    Stockfish::usize operator()(const Stockfish::Eval::NNUE::Network& network) const noexcept {
        return network.get_content_hash();
    }
};

#endif
