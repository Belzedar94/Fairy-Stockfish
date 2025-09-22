/*
  Fairy-Stockfish, a UCI chess variant playing engine derived from Stockfish
  Copyright (C) 2018-2022 Fabian Fichter

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

#include <algorithm>
#include <array>
#include <cctype>
#include <deque>
#include <iostream>
#include <sstream>
#include <string>

#include "partner.h"
#include "thread.h"
#include "uci.h"

namespace Stockfish {

PartnerHandler Partner; // Global object

namespace {

struct PieceMaskParseResult {
    uint64_t mask = 0;
    bool clear = false;
};

PieceType piece_type_from_token(const Position& pos, char symbol) {
    const std::string& mapping = pos.piece_to_char();
    const std::string& synonyms = pos.piece_to_char_synonyms();
    char lowered = static_cast<char>(std::tolower(static_cast<unsigned char>(symbol)));

    for (int pt = PAWN; pt < KING; ++pt)
    {
        PieceType pieceType = PieceType(pt);
        size_t whiteIndex = size_t(make_piece(WHITE, pieceType));
        size_t blackIndex = size_t(make_piece(BLACK, pieceType));
        auto matches = [&](char candidate) {
            return candidate && candidate != ' ' && std::tolower(static_cast<unsigned char>(candidate)) == lowered;
        };

        if ((whiteIndex < mapping.size() && matches(mapping[whiteIndex]))
            || (blackIndex < mapping.size() && matches(mapping[blackIndex]))
            || (whiteIndex < synonyms.size() && matches(synonyms[whiteIndex]))
            || (blackIndex < synonyms.size() && matches(synonyms[blackIndex])))
            return pieceType;
    }

    return NO_PIECE_TYPE;
}

PieceMaskParseResult parse_piece_list(const Position& pos, const std::string& text) {
    PieceMaskParseResult result;
    std::string lower = text;
    std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return char(std::tolower(c)); });

    if (lower.find("none") != std::string::npos || lower.find("clear") != std::string::npos)
        result.clear = true;

    for (char raw : text)
    {
        unsigned char c = static_cast<unsigned char>(raw);
        if (std::isspace(c) || raw == ',' || raw == ';' || raw == '.')
            continue;
        if (raw == '-' || raw == '0')
        {
            result.clear = true;
            continue;
        }

        PieceType pt = piece_type_from_token(pos, raw);
        if (pt != NO_PIECE_TYPE && pt != KING)
        {
            result.mask |= 1ULL << pt;
        }
    }

    if (!result.mask && !result.clear)
    {
        std::string trimmed = text;
        trimmed.erase(std::remove_if(trimmed.begin(), trimmed.end(), [](unsigned char c) { return std::isspace(c); }), trimmed.end());
        if (trimmed.empty())
            result.clear = true;
    }

    if (result.clear && !result.mask)
        result.mask = 0;

    return result;
}

std::string mask_to_string(const Position& pos, uint64_t mask, Color color) {
    if (!mask)
        return std::string();

    std::string result;
    const std::string& mapping = pos.piece_to_char();
    for (int pt = PAWN; pt < KING; ++pt)
    {
        if (!(mask & (1ULL << pt)))
            continue;

        PieceType pieceType = PieceType(pt);
        size_t index = size_t(make_piece(color, pieceType));
        if (index < mapping.size() && mapping[index] && mapping[index] != ' ')
            result += mapping[index];
    }
    return result;
}

std::string describe_mask(const Position& pos, uint64_t mask) {
    std::string description = mask_to_string(pos, mask, WHITE);
    if (description.empty())
        description = mask_to_string(pos, mask, BLACK);
    return description;
}

struct PieceFlowInfo {
    std::array<int, PIECE_TYPE_NB> toPartner{};
    std::array<int, PIECE_TYPE_NB> consume{};
    std::array<int, PIECE_TYPE_NB> request{};
};

uint64_t counts_to_mask(const std::array<int, PIECE_TYPE_NB>& counts) {
    uint64_t mask = 0;
    for (int pt = PAWN; pt < KING; ++pt)
        if (counts[pt] > 0)
            mask |= 1ULL << pt;
    return mask;
}

std::string format_piece_counts(const Position& pos, const std::array<int, PIECE_TYPE_NB>& counts,
                                Color color, uint64_t allowedMask = ~0ULL) {
    std::string result;
    const std::string& mapping = pos.piece_to_char();

    for (int pt = PAWN; pt < KING; ++pt)
    {
        if (!(allowedMask & (1ULL << pt)))
            continue;

        int count = counts[pt];
        if (count <= 0)
            continue;

        PieceType pieceType = PieceType(pt);
        size_t index = size_t(make_piece(color, pieceType));
        if (index >= mapping.size())
            continue;

        char symbol = mapping[index];
        if (!symbol || symbol == ' ')
            continue;

        result.append(count, symbol);
    }

    return result;
}

PieceFlowInfo analyze_piece_flow(const Position& rootPos, const std::vector<Move>& pv, size_t maxPlies) {
    PieceFlowInfo info;

    if (!rootPos.two_boards() || pv.empty() || pv[0] == MOVE_NONE)
        return info;

    StateListPtr states(new std::deque<StateInfo>(1));
    states->front() = StateInfo();

    Position pvPos;
    pvPos.set(rootPos.variant(), rootPos.fen(), rootPos.is_chess960(), &states->back(), rootPos.this_thread());

    Color us = rootPos.side_to_move();
    std::array<int, PIECE_TYPE_NB> available{};
    for (int pt = PAWN; pt < KING; ++pt)
        available[pt] = rootPos.count_in_hand(us, PieceType(pt));

    size_t limit = std::min(pv.size(), maxPlies);
    for (size_t ply = 0; ply < limit; ++ply)
    {
        Move m = pv[ply];
        if (m == MOVE_NONE)
            break;

        states->emplace_back();
        pvPos.do_move(m, states->back());
        const StateInfo& st = states->back();
        bool ourMove = (ply & 1) == 0;

        if (!ourMove)
            continue;

        if (st.capturedPiece)
        {
            Piece captured = st.capturedPiece;
            if (st.capturedpromoted)
                captured = st.unpromotedCapturedPiece ? st.unpromotedCapturedPiece
                                                      : make_piece(color_of(captured), rootPos.main_promotion_pawn_type(color_of(captured)));

            PieceType pt = type_of(captured);
            if (pt != KING)
                info.toPartner[pt]++;
        }

        if (type_of(m) == DROP)
        {
            PieceType dropType = in_hand_piece_type(m);
            if (dropType != NO_PIECE_TYPE && dropType != KING)
            {
                info.consume[dropType]++;
                if (available[dropType] > 0)
                    --available[dropType];
                else
                    info.request[dropType]++;
            }
        }
    }

    return info;
}

} // namespace

void PartnerHandler::reset() {
    fast = sitRequested = partnerDead = weDead = weWin = weVirtualWin = weVirtualLoss = false;
    time = opptime = 0;
    moveRequested = MOVE_NONE;
    partnerNeedsMask = 0;
    partnerBansMask = 0;
    ourNeedsMask = 0;
    lastFlowSummary.clear();
    lastNeedSummary.clear();
    lastDeliveredMask = 0;
    lastBlockedMask = 0;
    lastSitAnnounce = 0;
}

template <PartnerType p>
void PartnerHandler::ptell(const std::string& message) {
    if (p == ALL_PARTNERS || (p == FAIRY && isFairy) || (p == HUMAN && !isFairy))
        sync_cout << "tellics ptell " << message << sync_endl;
}

void PartnerHandler::parse_partner(std::istringstream& is) {
    std::string token;
    if (is >> token)
        // handshake to identify Fairy-Stockfish
        ptell("partner Fairy-Stockfish is an engine. Ask it 'help' for supported commands.");
    else
        isFairy = false;
}

void PartnerHandler::parse_ptell(std::istringstream& is, const Position& pos) {
    std::string token;
    is >> token;
    if (token == "partner")
    {
        // handshake to identify Fairy-Stockfish
        if (is >> token && token == "Fairy-Stockfish")
            isFairy = true;
    }
    else if (token == "help")
    {
        if (!(is >> token))
        {
            ptell<HUMAN>("I listen to the commands help, sit, go, move, fast, slow, dead, x, time, otim, need, ban, and allow.");
            ptell<HUMAN>("Tell 'help sit', etc. for details.");
        }
        else if (token == "sit")
            ptell<HUMAN>("After receiving 'sit', I stop moving. Also see 'go'.");
        else if (token == "go")
            ptell<HUMAN>("After receiving 'go', I will no longer sit.");
        else if (token == "move")
        {
            ptell<HUMAN>("After receiving 'move', I will move immediately." );
            ptell<HUMAN>("If you specify a valid move, e.g., 'move e2e4', I will play it.");
        }
        else if (token == "fast")
            ptell<HUMAN>("After receiving 'fast', I will play fast.");
        else if (token == "slow")
            ptell<HUMAN>("After receiving 'slow', I will play at normal speed.");
        else if (token == "dead")
            ptell<HUMAN>("After receiving 'dead', I assume you are dead and I play fast.");
        else if (token == "x")
            ptell<HUMAN>("After receiving 'x', I assume I can play normally again.");
        else if (token == "time")
        {
            ptell<HUMAN>("'time' together with your time in centiseconds allows me to consider your time.");
            ptell<HUMAN>("E.g., 'time 1000' for 10 seconds.");
        }
        else if (token == "otim")
            ptell<HUMAN>("'otim' together with your opponent's time in centiseconds allows me to consider his time.");
        else if (token == "need")
        {
            ptell<HUMAN>("Use 'need <pieces>' to tell me what material you require. '-' clears the request.");
            ptell<HUMAN>("Examples: 'need np' or 'need -'.");
        }
        else if (token == "ban")
        {
            ptell<HUMAN>("Use 'ban <pieces>' to tell me which drops I should not rely on. '-' clears bans.");
        }
        else if (token == "allow")
        {
            ptell<HUMAN>("Use 'allow <pieces>' to remove bans for individual pieces or '-' to clear all bans.");
        }
    }
    else if (!pos.two_boards())
        return;
    else if (token == "sit")
    {
        // Avoid deadlocking sit
        if (!isFairy || !weWin)
            sitRequested = true;
        ptell<HUMAN>("I sit, tell me 'go' to continue");
    }
    else if (token == "go")
    {
        sitRequested = false;
        Threads.stop = true;
    }
    else if (token == "move")
    {
        if (is >> token)
        {
            // if the given move is valid and we can still abort the search, play it
            Move move = UCI::to_move(pos, token);
            if (move && !Threads.abort.exchange(true))
                moveRequested = move;
            else
                ptell<HUMAN>("sorry, not possible");
        }
        else
            // Move immediately on request
            Threads.stop = true;
    }
    else if (token == "fast")
    {
        fast = true;
        ptell<HUMAN>("I play fast, tell me 'slow' to play normally again");
    }
    else if (token == "slow")
    {
        fast = false;
        ptell<HUMAN>("I play at normal speed again.");
    }
    else if (token == "dead")
    {
        partnerDead = true;
        ptell<HUMAN>("I play fast, tell me 'x' if you are no longer dead.");
    }
    else if (token == "x")
    {
        partnerDead = false;
        sitRequested = false;
        ptell<HUMAN>("I play normally again");
    }
    else if (token == "time")
    {
        int value;
        time = (is >> value) ? value * 10 : 0;
    }
    else if (token == "otim")
    {
        int value;
        opptime = (is >> value) ? value * 10 : 0;
    }
    else if (token == "need")
    {
        std::string rest;
        std::getline(is, rest);
        PieceMaskParseResult parsed = parse_piece_list(pos, rest);
        uint64_t mask = parsed.clear ? 0 : parsed.mask;
        partnerNeedsMask = mask;

        if (mask)
        {
            std::string pieces = describe_mask(pos, mask);
            if (!pieces.empty())
                ptell<HUMAN>("Acknowledged request for " + pieces + ".");
            else
                ptell<HUMAN>("Acknowledged piece request.");
        }
        else
            ptell<HUMAN>("Cleared piece requests.");
    }
    else if (token == "ban")
    {
        std::string rest;
        std::getline(is, rest);
        PieceMaskParseResult parsed = parse_piece_list(pos, rest);
        if (parsed.clear)
        {
            partnerBansMask = 0;
            ptell<HUMAN>("Cleared piece bans.");
        }
        else if (parsed.mask)
        {
            uint64_t mask = partnerBansMask.load();
            mask |= parsed.mask;
            partnerBansMask = mask;
            std::string pieces = describe_mask(pos, parsed.mask);
            if (!pieces.empty())
                ptell<HUMAN>("Will avoid relying on " + pieces + ".");
            else
                ptell<HUMAN>("Recorded piece bans.");
        }
    }
    else if (token == "allow")
    {
        std::string rest;
        std::getline(is, rest);
        PieceMaskParseResult parsed = parse_piece_list(pos, rest);
        if (parsed.clear)
        {
            partnerBansMask = 0;
            ptell<HUMAN>("Cleared piece bans.");
        }
        else if (parsed.mask)
        {
            uint64_t mask = partnerBansMask.load();
            mask &= ~parsed.mask;
            partnerBansMask = mask;
            std::string pieces = describe_mask(pos, parsed.mask);
            if (!pieces.empty())
                ptell<HUMAN>("Removed bans for " + pieces + ".");
            else
                ptell<HUMAN>("Updated piece bans.");
        }
    }
}

void PartnerHandler::update_piece_flow(const Position& rootPos, const std::vector<Move>& pv) {
    if (!rootPos.two_boards() || pv.empty() || pv[0] == MOVE_NONE)
    {
        if (!lastFlowSummary.empty())
        {
            ptell("flow -");
            lastFlowSummary.clear();
        }
        if (!lastNeedSummary.empty())
        {
            ptell("need -");
            lastNeedSummary.clear();
        }
        if (lastDeliveredMask)
        {
            ptell("feed -");
            lastDeliveredMask = 0;
        }
        lastBlockedMask = 0;
        return;
    }

    constexpr size_t MaxPlies = 12;
    PieceFlowInfo info = analyze_piece_flow(rootPos, pv, MaxPlies);

    Color us = rootPos.side_to_move();
    std::string incoming = format_piece_counts(rootPos, info.toPartner, ~us);
    std::string consume = format_piece_counts(rootPos, info.consume, us);

    std::string flowMessage;
    if (!incoming.empty() || !consume.empty())
    {
        flowMessage = "flow";
        if (!incoming.empty())
            flowMessage += " +" + incoming;
        if (!consume.empty())
            flowMessage += " -" + consume;

        if (flowMessage != lastFlowSummary)
        {
            ptell(flowMessage);
            lastFlowSummary = flowMessage;
        }
    }
    else if (!lastFlowSummary.empty())
    {
        ptell("flow -");
        lastFlowSummary.clear();
    }

    uint64_t requestMask = counts_to_mask(info.request);
    ourNeedsMask = requestMask;
    uint64_t bannedMask = partnerBansMask.load();
    uint64_t allowedMask = requestMask & ~bannedMask;

    std::array<int, PIECE_TYPE_NB> allowedCounts = info.request;
    if (bannedMask)
        for (int pt = PAWN; pt < KING; ++pt)
            if (bannedMask & (1ULL << pt))
                allowedCounts[pt] = 0;

    std::string needMessage;
    if (allowedMask)
    {
        std::string neededPieces = format_piece_counts(rootPos, allowedCounts, us, allowedMask);
        if (!neededPieces.empty())
            needMessage = "need " + neededPieces;
    }

    if (!needMessage.empty())
    {
        if (needMessage != lastNeedSummary)
        {
            ptell(needMessage);
            lastNeedSummary = needMessage;
        }
    }
    else if (!lastNeedSummary.empty())
    {
        ptell("need -");
        lastNeedSummary.clear();
    }

    uint64_t blockedMask = requestMask & bannedMask;
    if (blockedMask && blockedMask != lastBlockedMask)
    {
        std::string blockedPieces = format_piece_counts(rootPos, info.request, us, blockedMask);
        if (!blockedPieces.empty())
            ptell<HUMAN>("Cannot request " + blockedPieces + " (banned).");
        lastBlockedMask = blockedMask;
    }
    else if (!blockedMask)
        lastBlockedMask = 0;

    uint64_t incomingMask = counts_to_mask(info.toPartner);
    uint64_t deliverMask = incomingMask & partnerNeedsMask.load();
    if (deliverMask != lastDeliveredMask)
    {
        if (deliverMask)
        {
            std::string deliverPieces = format_piece_counts(rootPos, info.toPartner, ~us, deliverMask);
            if (!deliverPieces.empty())
                ptell("feed " + deliverPieces);
        }
        else if (lastDeliveredMask)
            ptell("feed -");

        lastDeliveredMask = deliverMask;
    }
}

template void PartnerHandler::ptell<HUMAN>(const std::string&);
template void PartnerHandler::ptell<FAIRY>(const std::string&);
template void PartnerHandler::ptell<ALL_PARTNERS>(const std::string&);

} // namespace Stockfish
