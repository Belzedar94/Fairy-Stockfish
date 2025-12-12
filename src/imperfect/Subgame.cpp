/*
  Fairy-Stockfish, a UCI chess variant playing engine derived from Stockfish
  Copyright (C) 2018-2024 Fabian Fichter

  Fairy-Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Fairy-Stockfish is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Subgame.h"
#include "../movegen.h"
#include <algorithm>
#include <numeric>
#include <stack>
#include <tuple>

namespace Stockfish {
namespace FogOfWar {

/// compute_sequence_id_from_moves() generates a hash for a move sequence
SequenceId compute_sequence_id_from_moves(const std::vector<Move>& moves) {
    SequenceId hash = 0xcbf29ce484222325ULL; // FNV-1a offset basis
    constexpr SequenceId prime = 0x100000001b3ULL;

    for (Move m : moves) {
        hash ^= static_cast<SequenceId>(m);
        hash *= prime;
    }

    return hash;
}

SequenceId Subgame::compute_sequence_id(const std::vector<Move>& moves) {
    return compute_sequence_id_from_moves(moves);
}

SequenceId Subgame::extend_sequence_id(SequenceId base, Move move) const {
    SequenceId hash = base;
    constexpr SequenceId prime = 0x100000001b3ULL;
    hash ^= static_cast<SequenceId>(move);
    hash *= prime;
    return hash;
}

std::unique_ptr<GameTreeNode> Subgame::acquire_node() {
    std::unique_ptr<GameTreeNode> node;
    if (!nodePool.empty()) {
        node = std::move(nodePool.back());
        nodePool.pop_back();
        *node = GameTreeNode();
    } else {
        node = std::make_unique<GameTreeNode>();
    }

    ++liveNodeCount;
    return node;
}

void Subgame::release_subtree(std::unique_ptr<GameTreeNode>& node) {
    if (!node)
        return;

    for (auto& child : node->children)
        release_subtree(child);

    node->children.clear();
    nodePool.push_back(std::move(node));
    node = nullptr;
    if (liveNodeCount > 0)
        --liveNodeCount;
}

std::shared_ptr<InfosetNode> Subgame::get_infoset(SequenceId seqId, Color player) {
    auto it = infosets.find(seqId);
    if (it != infosets.end())
        return it->second;

    auto iset = std::make_shared<InfosetNode>();
    iset->sequenceId = seqId;
    iset->player = player;
    infosets[seqId] = iset;
    return iset;
}

void Subgame::construct(const std::vector<std::string>& sampledStateFens,
                        int minInfosetSize) {
    (void)minInfosetSize; // Parameter currently unused

    // Clear existing tree
    if (rootNode)
        release_subtree(rootNode);

    rootNode = acquire_node();
    infosets.clear();
    nodeIdCounter = 0;
    resolveEntered = false;
    liveNodeCount = rootNode ? 1 : 0;

    // Build tree from sampled states
    build_tree_from_samples(sampledStateFens);

    // Compute KLUSS region (2-KLUSS: order-2 neighborhood, unfrozen at distance 1)
    compute_kluss_region(sampledStateFens);
}

void Subgame::build_tree_from_samples(const std::vector<std::string>& sampledStateFens) {
    std::unique_lock<std::shared_mutex> lock(treeMutex);

    if (sampledStateFens.empty())
        return;

    // Initialize root node with first sampled state
    rootNode->nodeId = nodeIdCounter++;
    rootNode->stateFen = sampledStateFens.front();
    // Store FEN instead of Position
    // TODO: Parse FEN when needed for move generation
    rootNode->ourSequence = 0;
    rootNode->theirSequence = 0;
    rootNode->depth = 0;
    rootNode->inKLUSS = true;

    // For now, create a simple root infoset without parsing FENs
    // TODO: Parse FENs to determine correct player and actions
    auto rootInfoset = get_infoset(0, WHITE); // Default to WHITE

    // Initialize with empty actions (will be filled during expansion)
    if (rootInfoset) {
        std::lock_guard<std::mutex> guard(rootInfoset->infosetMutex);
        rootInfoset->regrets.clear();
        rootInfoset->strategy.clear();
        rootInfoset->cumulativeStrategy.clear();
        rootInfoset->visitCounts.clear();
        rootInfoset->qValues.clear();
        rootInfoset->variances.clear();
    }

    prune_outside_kluss();
}

void Subgame::compute_kluss_region(const std::vector<std::string>& sampledStateFens) {
    // For 2-KLUSS: nodes are in the knowledge region if they are reachable
    // within 2 moves from any sampled state
    // Simplified implementation: mark root and immediate children as in KLUSS
    (void)sampledStateFens; // Suppress unused parameter warning

    std::unique_lock<std::shared_mutex> lock(treeMutex);

    if (!rootNode)
        return;

    // Root is always in KLUSS
    std::stack<GameTreeNode*> stack;
    stack.push(rootNode.get());

    while (!stack.empty()) {
        GameTreeNode* node = stack.top();
        stack.pop();

        if (!node)
            continue;

        node->inKLUSS = node->depth <= 2;

        // Track frozen/unfrozen infosets based on distance
        mark_frozen_state(node);

        for (auto& child : node->children) {
            child->depth = node->depth + 1;
            if (child->depth <= 2)
                stack.push(child.get());
        }
    }
}

bool Subgame::is_in_kluss(const GameTreeNode* node) const {
    return node && node->inKLUSS;
}

GameTreeNode* Subgame::expand_node(GameTreeNode* leaf, Position& pos) {
    if (!leaf || leaf->expanded)
        return nullptr;

    // Generate children for this leaf
    StateInfo st;
    std::vector<Move> legalMoves;

    for (const auto& m : MoveList<LEGAL>(pos))
        legalMoves.push_back(m);

    if (legalMoves.empty()) {
        // Terminal node (checkmate or stalemate)
        leaf->terminal = true;
        leaf->terminalValue = pos.checkers() ? -1.0f : 0.0f; // Checkmate or stalemate
        return leaf;
    }

    // Create child nodes
    for (Move m : legalMoves) {
        auto child = acquire_node();
        child->nodeId = nodeIdCounter++;
        child->parent = leaf;
        child->depth = leaf->depth + 1;

        // Initialize sequences from parent so both perspectives remain valid
        child->ourSequence = leaf->ourSequence;
        child->theirSequence = leaf->theirSequence;

        // Update sequences for the mover
        Color mover = pos.side_to_move();
        if (mover == WHITE)
            child->ourSequence = extend_sequence_id(child->ourSequence, m);
        else
            child->theirSequence = extend_sequence_id(child->theirSequence, m);

        // Make move to get child state FEN
        pos.do_move(m, st);
        child->stateFen = pos.fen();
        pos.undo_move(m);

        child->inKLUSS = child->depth <= 2;
        mark_frozen_state(child.get());

        leaf->children.push_back(std::move(child));
    }

    leaf->expanded = true;
    return leaf;
}

size_t Subgame::count_nodes() const {
    std::shared_lock<std::shared_mutex> lock(treeMutex);

    if (!rootNode)
        return 0;

    size_t count = 1;
    std::vector<const GameTreeNode*> stack = {rootNode.get()};

    while (!stack.empty()) {
        const GameTreeNode* node = stack.back();
        stack.pop_back();

        for (const auto& child : node->children) {
            count++;
            stack.push_back(child.get());
        }
    }

    return count;
}

int Subgame::average_depth() const {
    std::shared_lock<std::shared_mutex> lock(treeMutex);

    if (!rootNode)
        return 0;

    int totalDepth = 0;
    int nodeCount = 0;
    std::vector<const GameTreeNode*> stack = {rootNode.get()};

    while (!stack.empty()) {
        const GameTreeNode* node = stack.back();
        stack.pop_back();

        totalDepth += node->depth;
        nodeCount++;

        for (const auto& child : node->children)
            stack.push_back(child.get());
    }

    return nodeCount > 0 ? totalDepth / nodeCount : 0;
}

std::vector<std::shared_ptr<InfosetNode>> Subgame::snapshot_infosets() const {
    std::shared_lock<std::shared_mutex> lock(treeMutex);
    std::vector<std::shared_ptr<InfosetNode>> result;
    result.reserve(infosets.size());
    for (const auto& kv : infosets)
        result.push_back(kv.second);
    return result;
}

void Subgame::mark_frozen_state(GameTreeNode* node) {
    if (!node)
        return;

    // Determine which sequence to use based on depth parity
    Color nodePlayer = (node->depth % 2 == 0) ? WHITE : BLACK;
    SequenceId seqId = nodePlayer == WHITE ? node->ourSequence : node->theirSequence;

    auto infoset = get_infoset(seqId, nodePlayer);
    if (!infoset)
        return;

    std::lock_guard<std::mutex> guard(infoset->infosetMutex);
    infoset->unfrozen = node->depth <= 1;

    if (!infoset->unfrozen && infoset->trunkStrategy.empty()) {
        if (!infoset->strategy.empty())
            infoset->trunkStrategy = infoset->strategy;
        else if (!infoset->actions.empty())
            infoset->trunkStrategy.assign(infoset->actions.size(), 1.0f / infoset->actions.size());
    }
}

void Subgame::prune_outside_kluss() {
    if (!rootNode)
        return;

    std::unique_lock<std::shared_mutex> lock(treeMutex);
    std::vector<GameTreeNode*> stack = {rootNode.get()};

    while (!stack.empty()) {
        GameTreeNode* node = stack.back();
        stack.pop_back();

        for (size_t i = 0; i < node->children.size();) {
            if (!node->children[i]->inKLUSS) {
                auto pruned = std::move(node->children[i]);
                node->children.erase(node->children.begin() + i);
                release_subtree(pruned);
                continue;
            }

            stack.push_back(node->children[i].get());
            ++i;
        }
    }
}

void Subgame::enforce_node_limit() {
    if (!rootNode)
        return;

    std::unique_lock<std::shared_mutex> lock(treeMutex);

    auto select_prunable_leaf = [this](GameTreeNode* root) {
        GameTreeNode* target = nullptr;
        GameTreeNode* targetParent = nullptr;
        size_t targetIndex = 0;
        int bestDepth = -1;
        bool preferOutside = false;

        std::vector<GameTreeNode*> stack = {root};

        while (!stack.empty()) {
            GameTreeNode* node = stack.back();
            stack.pop_back();

            for (size_t idx = 0; idx < node->children.size(); ++idx) {
                GameTreeNode* child = node->children[idx].get();
                if (!child)
                    continue;

                bool isLeaf = child->children.empty();
                bool outside = !child->inKLUSS || child->depth > 2;

                if (isLeaf && (!target || outside > preferOutside ||
                               (outside == preferOutside && child->depth > bestDepth))) {
                    target = child;
                    targetParent = node;
                    targetIndex = idx;
                    bestDepth = child->depth;
                    preferOutside = outside;
                }

                if (!isLeaf)
                    stack.push_back(child);
            }
        }

        return std::tuple<GameTreeNode*, GameTreeNode*, size_t, bool>(target, targetParent, targetIndex, preferOutside);
    };

    while (liveNodeCount > nodeLimit) {
        auto [leaf, parent, index, found] = select_prunable_leaf(rootNode.get());
        if (!leaf || !parent)
            break;

        (void)found;
        auto removed = std::move(parent->children[index]);
        parent->children.erase(parent->children.begin() + index);
        release_subtree(removed);
    }
}

/// compute_alternative_value() for Resolve gadget (Appendix B.3.1)
/// Uses current (x,y) instead of best-response values for stability
float compute_alternative_value(const InfosetNode* infoset,
                                 const std::vector<float>& currentX,
                                 const std::vector<float>& currentY) {
    (void)currentX; // Suppress unused parameter warning
    if (!infoset)
        return 0.0f;

    // Alternative value uses Resolve prior
    std::vector<float> prior = compute_resolve_prior(infoset, currentY);
    float altValue = 0.0f;

    const size_t actions = infoset->actions.size();
    for (size_t i = 0; i < actions; ++i) {
        float childVal = (i < infoset->qValues.size()) ? infoset->qValues[i] : infoset->value;
        float weight = i < prior.size() ? prior[i] : 0.0f;
        altValue += weight * childVal;
    }

    // Blend with current estimate for stability
    float currentValue = infoset->value;
    return 0.5f * altValue + 0.5f * currentValue;
}

/// compute_gift() for Resolve gadget (Appendix B.3.1)
float compute_gift(const InfosetNode* infoset,
                   const std::vector<float>& currentX,
                   const std::vector<float>& currentY) {
    // Gift is the value opponent forfeits by playing into the subgame
    float altValue = compute_alternative_value(infoset, currentX, currentY);
    float currentValue = infoset ? infoset->value : 0.0f;
    return altValue - currentValue;
}

std::vector<float> compute_resolve_prior(const InfosetNode* infoset,
                                         const std::vector<float>& opponentStrategy) {
    if (!infoset)
        return {};

    size_t n = infoset->actions.size();
    std::vector<float> prior(n, 0.0f);

    float uniform = n ? 1.0f / static_cast<float>(n) : 0.0f;

    for (size_t i = 0; i < n; ++i) {
        float opp = (i < opponentStrategy.size()) ? opponentStrategy[i] : 0.0f;
        prior[i] = 0.5f * uniform + 0.5f * opp;
    }

    // Renormalize to guard against zero-sum issues
    float sum = std::accumulate(prior.begin(), prior.end(), 0.0f);
    if (sum > 0.0f) {
        for (float& p : prior)
            p /= sum;
    }

    return prior;
}

} // namespace FogOfWar
} // namespace Stockfish
