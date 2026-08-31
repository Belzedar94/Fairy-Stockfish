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

#include "network.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "../bitboard.h"
#include "../misc.h"
#include "../position.h"
#include "../types.h"
#include "nnue_accumulator.h"

namespace Stockfish::Eval::NNUE {

namespace fs = std::filesystem;

namespace {

constexpr std::string_view LegacyDescription =
  "Network trained with the https://github.com/glinscott/nnue-pytorch trainer.";
constexpr usize LegacyDescriptionBytes = 75;
constexpr usize LegacyFirstLayerOffset = 47580251;
constexpr usize LegacyLayerBytes       = 17640;

constexpr std::array<std::array<usize, PIECE_NB>, COLOR_NB> PieceSquareIndex = {{
  {{0, 0, 128, 256, 384, 512, 640, 0, 0, 64, 192, 320, 448, 576, 640, 0}},
  {{0, 64, 192, 320, 448, 576, 640, 0, 0, 0, 128, 256, 384, 512, 640, 0}},
}};

constexpr u32 rotate_right(u32 value, unsigned count) {
    return (value >> count) | (value << (32 - count));
}

std::string sha256_hex(const std::vector<u8>& bytes) {
    constexpr std::array<u32, 64> constants = {
      0x428A2F98u, 0x71374491u, 0xB5C0FBCFu, 0xE9B5DBA5u, 0x3956C25Bu, 0x59F111F1u, 0x923F82A4u,
      0xAB1C5ED5u, 0xD807AA98u, 0x12835B01u, 0x243185BEu, 0x550C7DC3u, 0x72BE5D74u, 0x80DEB1FEu,
      0x9BDC06A7u, 0xC19BF174u, 0xE49B69C1u, 0xEFBE4786u, 0x0FC19DC6u, 0x240CA1CCu, 0x2DE92C6Fu,
      0x4A7484AAu, 0x5CB0A9DCu, 0x76F988DAu, 0x983E5152u, 0xA831C66Du, 0xB00327C8u, 0xBF597FC7u,
      0xC6E00BF3u, 0xD5A79147u, 0x06CA6351u, 0x14292967u, 0x27B70A85u, 0x2E1B2138u, 0x4D2C6DFCu,
      0x53380D13u, 0x650A7354u, 0x766A0ABBu, 0x81C2C92Eu, 0x92722C85u, 0xA2BFE8A1u, 0xA81A664Bu,
      0xC24B8B70u, 0xC76C51A3u, 0xD192E819u, 0xD6990624u, 0xF40E3585u, 0x106AA070u, 0x19A4C116u,
      0x1E376C08u, 0x2748774Cu, 0x34B0BCB5u, 0x391C0CB3u, 0x4ED8AA4Au, 0x5B9CCA4Fu, 0x682E6FF3u,
      0x748F82EEu, 0x78A5636Fu, 0x84C87814u, 0x8CC70208u, 0x90BEFFFAu, 0xA4506CEBu, 0xBEF9A3F7u,
      0xC67178F2u};

    std::array<u32, 8> state = {0x6A09E667u, 0xBB67AE85u, 0x3C6EF372u, 0xA54FF53Au,
                                0x510E527Fu, 0x9B05688Cu, 0x1F83D9ABu, 0x5BE0CD19u};

    auto transform = [&](const u8* block) {
        std::array<u32, 64> schedule{};
        for (usize i = 0; i < 16; ++i)
            schedule[i] = (u32(block[i * 4]) << 24) | (u32(block[i * 4 + 1]) << 16)
                        | (u32(block[i * 4 + 2]) << 8) | u32(block[i * 4 + 3]);

        for (usize i = 16; i < schedule.size(); ++i)
        {
            const u32 s0 = rotate_right(schedule[i - 15], 7) ^ rotate_right(schedule[i - 15], 18)
                         ^ (schedule[i - 15] >> 3);
            const u32 s1 = rotate_right(schedule[i - 2], 17) ^ rotate_right(schedule[i - 2], 19)
                         ^ (schedule[i - 2] >> 10);
            schedule[i] = schedule[i - 16] + s0 + schedule[i - 7] + s1;
        }

        u32 a = state[0];
        u32 b = state[1];
        u32 c = state[2];
        u32 d = state[3];
        u32 e = state[4];
        u32 f = state[5];
        u32 g = state[6];
        u32 h = state[7];

        for (usize i = 0; i < schedule.size(); ++i)
        {
            const u32 sum1     = rotate_right(e, 6) ^ rotate_right(e, 11) ^ rotate_right(e, 25);
            const u32 choose   = (e & f) ^ (~e & g);
            const u32 temp1    = h + sum1 + choose + constants[i] + schedule[i];
            const u32 sum0     = rotate_right(a, 2) ^ rotate_right(a, 13) ^ rotate_right(a, 22);
            const u32 majority = (a & b) ^ (a & c) ^ (b & c);
            const u32 temp2    = sum0 + majority;

            h = g;
            g = f;
            f = e;
            e = d + temp1;
            d = c;
            c = b;
            b = a;
            a = temp1 + temp2;
        }

        state[0] += a;
        state[1] += b;
        state[2] += c;
        state[3] += d;
        state[4] += e;
        state[5] += f;
        state[6] += g;
        state[7] += h;
    };

    const usize complete = bytes.size() / 64;
    for (usize i = 0; i < complete; ++i)
        transform(bytes.data() + i * 64);

    std::array<u8, 128> tail{};
    const usize         remaining = bytes.size() % 64;
    std::memcpy(tail.data(), bytes.data() + complete * 64, remaining);
    tail[remaining] = 0x80;

    const usize tailBytes = remaining < 56 ? 64 : 128;
    const u64   bitLength = u64(bytes.size()) * 8;
    for (usize i = 0; i < 8; ++i)
        tail[tailBytes - 1 - i] = u8(bitLength >> (i * 8));

    transform(tail.data());
    if (tailBytes == 128)
        transform(tail.data() + 64);

    std::ostringstream out;
    out << std::uppercase << std::hex << std::setfill('0');
    for (u32 word : state)
        out << std::setw(8) << word;
    return out.str();
}

std::string hex_u32(u32 value) {
    std::ostringstream out;
    out << "0x" << std::uppercase << std::hex << std::setw(8) << std::setfill('0') << value;
    return out.str();
}

class ByteCursor {
   public:
    explicit ByteCursor(const std::vector<u8>& source) :
        bytes(source) {}

    usize offset() const { return cursor; }

    bool read_u32(u32& value) {
        if (remaining() < 4)
            return false;
        value = u32(bytes[cursor]) | (u32(bytes[cursor + 1]) << 8) | (u32(bytes[cursor + 2]) << 16)
              | (u32(bytes[cursor + 3]) << 24);
        cursor += 4;
        return true;
    }

    bool read_string(usize count, std::string& value) {
        if (remaining() < count)
            return false;
        value.assign(reinterpret_cast<const char*>(bytes.data() + cursor), count);
        cursor += count;
        return true;
    }

    template<typename T>
    bool read_array(T* output, usize count) {
        static_assert(std::is_integral_v<T>);
        if (count > remaining() / sizeof(T))
            return false;

        if constexpr (sizeof(T) == 1)
            std::memcpy(output, bytes.data() + cursor, count * sizeof(T));
        else if (IsLittleEndian)
            std::memcpy(output, bytes.data() + cursor, count * sizeof(T));
        else
        {
            using U = std::make_unsigned_t<T>;
            for (usize i = 0; i < count; ++i)
            {
                U value = 0;
                for (usize j = 0; j < sizeof(T); ++j)
                    value |= U(bytes[cursor + i * sizeof(T) + j]) << (j * 8);
                std::memcpy(output + i, &value, sizeof(T));
            }
        }

        cursor += count * sizeof(T);
        return true;
    }

   private:
    usize remaining() const { return bytes.size() - cursor; }

    const std::vector<u8>& bytes;
    usize                  cursor = 0;
};

u8 clipped_relu(i64 value) {
    if (value <= 0)
        return 0;
    return u8(std::min<i64>(127, value / 64));
}

i32 checked_i32(i64 value) {
    assert(value >= std::numeric_limits<i32>::min());
    assert(value <= std::numeric_limits<i32>::max());
    return i32(value);
}

template<typename T>
T wrapping_add(T left, T right) {
    static_assert(std::is_signed_v<T> && std::is_integral_v<T>);
    using U = std::make_unsigned_t<T>;
    U leftBits;
    U rightBits;
    std::memcpy(&leftBits, &left, sizeof(T));
    std::memcpy(&rightBits, &right, sizeof(T));
    const U sumBits = U(leftBits + rightBits);
    T       result;
    std::memcpy(&result, &sumBits, sizeof(T));
    return result;
}

}  // namespace

Network::Network() {
    // Shared-memory copies include padding, so initialize every byte rather
    // than leaving non-semantic padding indeterminate.
    std::memset(static_cast<void*>(this), 0, sizeof(*this));
}

std::optional<std::string>
Network::load(const fs::path& rootDirectory, fs::path evalfilePath, EvalFile& evalFile) {
    if (evalfilePath.empty())
        evalfilePath = fs::path(EvalFile::defaultName);

    const std::string basename = evalfilePath.filename().string();
    if (basename != LegacyCanonicalName && basename != LegacyAliasName)
        return "code=KOTH_NET_BASENAME expected=kingofthehill-978b86d0e6a4.nnue_or_KOTH_v1.nnue actual="
             + (basename.empty() ? std::string("<empty>") : basename);

    std::vector<fs::path> candidates;
    if (evalfilePath.is_absolute())
        candidates.push_back(evalfilePath);
    else
    {
        candidates.push_back(evalfilePath);
        const fs::path rooted = rootDirectory / evalfilePath;
        if (rooted.lexically_normal() != evalfilePath.lexically_normal())
            candidates.push_back(rooted);
    }

    std::error_code ec;
    for (const fs::path& candidatePath : candidates)
    {
        if (!fs::is_regular_file(candidatePath, ec))
        {
            ec.clear();
            continue;
        }

        auto        candidate = std::make_unique<Network>();
        std::string description;
        if (auto error = candidate->load_single_file(candidatePath, description))
            return error;

        *this                   = *candidate;
        evalFile.current        = evalfilePath;
        evalFile.netDescription = std::move(description);
        return std::nullopt;
    }

    return "code=KOTH_NET_MISSING basename=" + basename;
}

std::optional<std::string> Network::load_single_file(const fs::path& file,
                                                     std::string&    description) {
    std::ifstream stream(file, std::ios::binary | std::ios::ate);
    if (!stream)
        return "code=KOTH_NET_OPEN_FAILED basename=" + file.filename().string();

    const auto end = stream.tellg();
    if (end < 0)
        return "code=KOTH_NET_SIZE_QUERY_FAILED basename=" + file.filename().string();

    const auto size = u64(end);
    if (size != LegacyFileBytes)
        return "code=KOTH_NET_SIZE expected=" + std::to_string(LegacyFileBytes)
             + " actual=" + std::to_string(size);

    std::vector<u8> bytes(LegacyFileBytes);
    stream.seekg(0);
    stream.read(reinterpret_cast<char*>(bytes.data()), std::streamsize(bytes.size()));
    if (!stream || usize(stream.gcount()) != bytes.size())
        return "code=KOTH_NET_SHORT_READ expected=" + std::to_string(LegacyFileBytes)
             + " actual=" + std::to_string(stream.gcount());

    return parse_exact_bytes(bytes, description);
}

std::optional<std::string> Network::parse_exact_bytes(const std::vector<u8>& bytes,
                                                      std::string&           description) {
    if (bytes.size() != LegacyFileBytes)
        return "code=KOTH_NET_SIZE expected=" + std::to_string(LegacyFileBytes)
             + " actual=" + std::to_string(bytes.size());

    ByteCursor cursor(bytes);
    u32        version          = 0;
    u32        architecture     = 0;
    u32        descriptionBytes = 0;
    if (!cursor.read_u32(version) || !cursor.read_u32(architecture)
        || !cursor.read_u32(descriptionBytes))
        return "code=KOTH_NET_TRUNCATED_HEADER";

    if (version != LegacyVersion)
        return "code=KOTH_NET_VERSION expected=" + hex_u32(LegacyVersion)
             + " actual=" + hex_u32(version);
    if (architecture != LegacyArchitecture)
        return "code=KOTH_NET_ARCHITECTURE expected=" + hex_u32(LegacyArchitecture)
             + " actual=" + hex_u32(architecture);
    if (descriptionBytes != LegacyDescriptionBytes)
        return "code=KOTH_NET_DESCRIPTION_SIZE expected=" + std::to_string(LegacyDescriptionBytes)
             + " actual=" + std::to_string(descriptionBytes);
    if (!cursor.read_string(descriptionBytes, description))
        return "code=KOTH_NET_TRUNCATED_DESCRIPTION";
    if (description != LegacyDescription)
        return "code=KOTH_NET_DESCRIPTION_MISMATCH";

    u32 featureHash = 0;
    if (!cursor.read_u32(featureHash))
        return "code=KOTH_NET_TRUNCATED_FEATURE_HEADER";
    if (featureHash != LegacyFeatureHash)
        return "code=KOTH_NET_FEATURE_HASH expected=" + hex_u32(LegacyFeatureHash)
             + " actual=" + hex_u32(featureHash);

    if (!cursor.read_array(featureTransformer.biases.data(), featureTransformer.biases.size())
        || !cursor.read_array(featureTransformer.weights.data(), featureTransformer.weights.size())
        || !cursor.read_array(featureTransformer.psqtWeights.data(),
                              featureTransformer.psqtWeights.size()))
        return "code=KOTH_NET_TRUNCATED_FEATURE_PARAMETERS offset="
             + std::to_string(cursor.offset());

    if (cursor.offset() != LegacyFirstLayerOffset)
        return "code=KOTH_NET_FEATURE_SECTION_SIZE expected_offset="
             + std::to_string(LegacyFirstLayerOffset)
             + " actual_offset=" + std::to_string(cursor.offset());

    for (usize bucket = 0; bucket < LegacyLayerStacks; ++bucket)
    {
        const usize expectedOffset = LegacyFirstLayerOffset + bucket * LegacyLayerBytes;
        if (cursor.offset() != expectedOffset)
            return "code=KOTH_NET_LAYER_OFFSET bucket=" + std::to_string(bucket) + " expected="
                 + std::to_string(expectedOffset) + " actual=" + std::to_string(cursor.offset());

        u32 layerHash = 0;
        if (!cursor.read_u32(layerHash))
            return "code=KOTH_NET_TRUNCATED_LAYER_HEADER bucket=" + std::to_string(bucket);
        if (layerHash != LegacyLayerHash)
            return "code=KOTH_NET_LAYER_HASH bucket=" + std::to_string(bucket)
                 + " expected=" + hex_u32(LegacyLayerHash) + " actual=" + hex_u32(layerHash);

        auto& layer = layerStacks[bucket];
        if (!cursor.read_array(layer.fc0Biases.data(), layer.fc0Biases.size())
            || !cursor.read_array(layer.fc0Weights.data(), layer.fc0Weights.size())
            || !cursor.read_array(layer.fc1Biases.data(), layer.fc1Biases.size())
            || !cursor.read_array(layer.fc1Weights.data(), layer.fc1Weights.size())
            || !cursor.read_array(layer.outputBias.data(), layer.outputBias.size())
            || !cursor.read_array(layer.outputWeights.data(), layer.outputWeights.size()))
            return "code=KOTH_NET_TRUNCATED_LAYER_PARAMETERS bucket=" + std::to_string(bucket)
                 + " offset=" + std::to_string(cursor.offset());
    }

    if (cursor.offset() != bytes.size())
        return "code=KOTH_NET_TRAILING_BYTES expected_eof=" + std::to_string(cursor.offset())
             + " actual_size=" + std::to_string(bytes.size());

    const std::string digest = sha256_hex(bytes);
    if (digest != LegacySha256)
        return "code=KOTH_NET_SHA256 expected=" + std::string(LegacySha256) + " actual=" + digest;

    initialized = true;
    return std::nullopt;
}

bool Network::save(const EvalFile&, const std::optional<fs::path>&) const {
    sync_cout << "Failed to export a net: code=KOTH_NET_EXPORT_LICENSE_UNRESOLVED" << sync_endl;
    return false;
}

usize Network::get_content_hash() const { return initialized ? usize(0x978B86D0E6A45E05ULL) : 0; }

LegacyRawOutput Network::evaluate_raw_bucket(const Position& pos, usize bucket) const {
    if (!initialized || bucket >= LegacyLayerStacks)
        std::abort();

    std::array<std::array<i16, LegacyL1>, COLOR_NB>          accumulation{};
    std::array<std::array<i32, LegacyPsqtBuckets>, COLOR_NB> psqt{};

    for (Color perspective : {WHITE, BLACK})
    {
        for (usize j = 0; j < LegacyL1; ++j)
            accumulation[perspective][j] = featureTransformer.biases[j];

        const Square kingSquare =
          Square(int(pos.square<KING>(perspective)) ^ (perspective == BLACK ? 56 : 0));
        Bitboard occupied = pos.pieces();
        while (occupied)
        {
            const Square square   = pop_lsb(occupied);
            const Piece  piece    = pos.piece_on(square);
            const usize  oriented = usize(int(square) ^ (perspective == BLACK ? 56 : 0));
            const usize  index =
              oriented + PieceSquareIndex[perspective][piece] + 704 * usize(kingSquare);
            assert(index < LegacyInputDimensions);

            const usize weightOffset = index * LegacyL1;
            for (usize j = 0; j < LegacyL1; ++j)
                accumulation[perspective][j] = wrapping_add(
                  accumulation[perspective][j], featureTransformer.weights[weightOffset + j]);

            const usize psqtOffset = index * LegacyPsqtBuckets;
            for (usize j = 0; j < LegacyPsqtBuckets; ++j)
                psqt[perspective][j] = wrapping_add(psqt[perspective][j],
                                                    featureTransformer.psqtWeights[psqtOffset + j]);
        }
    }

    std::array<u8, LegacyNetworkInput> transformed{};
    const std::array<Color, COLOR_NB>  perspectives = {pos.side_to_move(), ~pos.side_to_move()};
    for (usize p = 0; p < COLOR_NB; ++p)
        for (usize j = 0; j < LegacyL1; ++j)
            transformed[p * LegacyL1 + j] =
              u8(std::clamp<i64>(accumulation[perspectives[p]][j], 0, 127));

    const i64 psqtRaw = (psqt[perspectives[0]][bucket] - psqt[perspectives[1]][bucket]) / 2;

    const auto&                       layer = layerStacks[bucket];
    std::array<i64, LegacyFc0Outputs> fc0{};
    for (usize output = 0; output < LegacyFc0Outputs; ++output)
    {
        i64         sum = layer.fc0Biases[output];
        const usize row = output * LegacyNetworkInput;
        for (usize input = 0; input < LegacyNetworkInput; ++input)
            sum += i64(layer.fc0Weights[row + input]) * transformed[input];
        fc0[output] = sum;
    }

    std::array<u8, LegacyFc0Outputs> ac0{};
    for (usize i = 0; i < ac0.size(); ++i)
        ac0[i] = clipped_relu(fc0[i]);

    std::array<i64, LegacyFc1Outputs> fc1{};
    for (usize output = 0; output < LegacyFc1Outputs; ++output)
    {
        i64         sum = layer.fc1Biases[output];
        const usize row = output * LegacyFc1PaddedIn;
        for (usize input = 0; input < LegacyFc0Outputs; ++input)
            sum += i64(layer.fc1Weights[row + input]) * ac0[input];
        fc1[output] = sum;
    }

    std::array<u8, LegacyFc1Outputs> ac1{};
    for (usize i = 0; i < ac1.size(); ++i)
        ac1[i] = clipped_relu(fc1[i]);

    i64 positionalRaw = layer.outputBias[0];
    for (usize input = 0; input < LegacyFc1Outputs; ++input)
        positionalRaw += i64(layer.outputWeights[input]) * ac1[input];

    return {checked_i32(psqtRaw), checked_i32(positionalRaw), bucket};
}

LegacyRawOutput Network::evaluate_raw(const Position& pos) const {
    const int pieces = pos.count<ALL_PIECES>();
    if (pieces < 1 || pieces > 32)
        std::abort();
    return evaluate_raw_bucket(pos, usize((pieces - 1) / 4));
}

NetworkOutput Network::evaluate(const Position&    pos,
                                AccumulatorStack&  accumulatorStack,
                                AccumulatorCaches& cache) const {
    (void) accumulatorStack;
    (void) cache;
    const auto  raw  = evaluate_raw(pos);
    const Value psqt = Value(raw.psqt / LegacyOutputScale);
    if (std::abs(raw.psqt) > LegacyLazyThreshold * LegacyOutputScale)
        return {psqt, VALUE_ZERO};
    const Value total = Value((i64(raw.psqt) + raw.positional) / LegacyOutputScale);
    return {psqt, total - psqt};
}

void Network::verify(const std::function<void(std::string_view)>& f,
                     const EvalFile&                              evalFile,
                     fs::path                                     evalfilePath) const {
    if (evalfilePath.empty())
        evalfilePath = fs::path(EvalFile::defaultName);

    if (!initialized || !evalFile.current.has_value() || evalFile.current != evalfilePath)
    {
        if (f)
            f("error code=KOTH_NET_NOT_LOADED expected=" + evalfilePath.filename().string());
        std::exit(EXIT_FAILURE);
    }

    if (f)
        f("koth network loaded=true basename=" + evalfilePath.filename().string()
          + " bytes=" + std::to_string(LegacyFileBytes) + " sha256=" + std::string(LegacySha256)
          + " container=" + hex_u32(LegacyVersion) + " architecture=" + hex_u32(LegacyArchitecture)
          + " feature_hash=" + hex_u32(LegacyFeatureHash)
          + " topology=45056x512-16-32-1 buckets=8 evaluator=scalar_reference");
}

std::string Network::status(const EvalFile& evalFile) const {
    if (!initialized || !evalFile.current.has_value())
        return "loaded=false default=kingofthehill-978b86d0e6a4.nnue";
    return "loaded=true basename=" + evalFile.current->filename().string()
         + " bytes=" + std::to_string(LegacyFileBytes) + " sha256=" + std::string(LegacySha256)
         + " architecture=" + hex_u32(LegacyArchitecture) + " evaluator=scalar_reference";
}

NnueEvalTrace Network::trace_evaluate(const Position&    pos,
                                      AccumulatorStack&  accumulatorStack,
                                      AccumulatorCaches& cache) const {
    (void) accumulatorStack;
    (void) cache;
    NnueEvalTrace trace{};
    const int     pieces = pos.count<ALL_PIECES>();
    if (pieces < 1 || pieces > 32)
        std::abort();
    trace.correctBucket = usize((pieces - 1) / 4);

    for (usize bucket = 0; bucket < LegacyLayerStacks; ++bucket)
    {
        const auto  raw          = evaluate_raw_bucket(pos, bucket);
        const Value psqt         = Value(raw.psqt / LegacyOutputScale);
        const Value total        = std::abs(raw.psqt) > LegacyLazyThreshold * LegacyOutputScale
                                   ? psqt
                                   : Value((i64(raw.psqt) + raw.positional) / LegacyOutputScale);
        trace.psqt[bucket]       = psqt;
        trace.positional[bucket] = total - psqt;
    }
    return trace;
}

}  // namespace Stockfish::Eval::NNUE
