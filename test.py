# -*- coding: utf-8 -*-

import faulthandler
import unittest

import pyffish as sf

faulthandler.enable()


class TestPyffish(unittest.TestCase):

    @staticmethod
    def _gating_info(move: str):
        if "," in move:
            prefix, suffix = move.split(",", 1)
            if "@" in prefix:
                tag, square = prefix.split("@", 1)
                if len(tag) == 1 and tag.isalpha() and square and square.isalnum():
                    return suffix, tag.lower(), square
        if len(move) < 5:
            return None
        rank_last = move[-1]
        if rank_last not in "0123456789":
            return None
        if len(move) >= 6 and move[-2] in "0123456789":
            square = move[-3:]
            gate_index = -4
        else:
            square = move[-2:]
            gate_index = -3
        if square[0] not in "abcdefghijklmnopqrstuvwxyz":
            return None
        gate_char = move[gate_index]
        if not gate_char.isalpha():
            return None
        return move[:gate_index], gate_char.lower(), square

    @staticmethod
    def _first_normal_move(moves):
        for move in moves:
            if TestPyffish._gating_info(move) is None and "@" not in move:
                return move
        raise AssertionError("No normal move available")

    @classmethod
    def _filter_potion_moves(cls, moves, kind):
        result = []
        for move in moves:
            info = cls._gating_info(move)
            if info and info[1] == kind:
                result.append(move)
        return result

    @classmethod
    def _has_potion_move(cls, moves, kind):
        for move in moves:
            info = cls._gating_info(move)
            if info and info[1] == kind:
                return True
        return False

    @staticmethod
    def _board_state(fen: str):
        board = {}
        board_part = fen.split()[0]
        rows = board_part.split("/")
        height = len(rows)
        width = 0
        for ch in rows[0]:
            if ch.isdigit():
                width += int(ch)
            elif ch.isalpha():
                width += 1
        for rank_index, row in enumerate(rows):
            file_index = 0
            for ch in row:
                if ch.isdigit():
                    file_index += int(ch)
                elif ch.isalpha():
                    square = chr(ord("a") + file_index) + str(height - rank_index)
                    board[square] = ch
                    file_index += 1
        return board, width, height

    @staticmethod
    def _square_to_coords(square: str):
        return ord(square[0]) - ord("a"), int(square[1:]) - 1

    @staticmethod
    def _coords_to_square(file_index: int, rank_index: int):
        return chr(ord("a") + file_index) + str(rank_index + 1)

    @classmethod
    def _freeze_zone_has_enemy(cls, board, width, height, square, enemy_color):
        file_index, rank_index = cls._square_to_coords(square)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                nx = file_index + dx
                ny = rank_index + dy
                if not (0 <= nx < width and 0 <= ny < height):
                    continue
                target = cls._coords_to_square(nx, ny)
                piece = board.get(target)
                if not piece:
                    continue
                if enemy_color == "w" and piece.isupper():
                    return True
                if enemy_color == "b" and piece.islower():
                    return True
        return False

    @classmethod
    def _squares_between(cls, origin: str, target: str):
        ox, oy = cls._square_to_coords(origin)
        tx, ty = cls._square_to_coords(target)
        dx = tx - ox
        dy = ty - oy
        step_x = (dx > 0) - (dx < 0)
        step_y = (dy > 0) - (dy < 0)
        if step_x and step_y and abs(dx) != abs(dy):
            return []
        steps = max(abs(dx), abs(dy))
        squares = []
        for step in range(1, steps):
            square = cls._coords_to_square(ox + step_x * step, oy + step_y * step)
            squares.append(square)
        return squares

    @classmethod
    def _jump_uses_gate(cls, move: str, board, width):
        info = cls._gating_info(move)
        if not info:
            return False
        base = info[0]
        gate = info[2]
        if len(base) < 4:
            return False
        origin = base[:2]
        target = base[2:4]
        if gate not in board:
            return False
        return gate in cls._squares_between(origin, target)

    def test_version(self):
        result = sf.version()
        self.assertEqual(len(result), 3)

    def test_info(self):
        result = sf.info()
        self.assertTrue(result.startswith("Fairy-Stockfish"))

    def test_variants_loaded(self):
        variants = sf.variants()
        self.assertIn("spell-chess", variants)

    def test_set_option(self):
        result = sf.set_option("UCI_Variant", "spell-chess")
        self.assertIsNone(result)

    def test_start_fen(self):
        fen = sf.start_fen("spell-chess")
        self.assertTrue(fen)

    def test_validate_fen(self):
        fen = sf.start_fen("spell-chess")
        self.assertEqual(sf.validate_fen(fen, "spell-chess"), sf.FEN_OK)

    def test_spell_chess_freeze_blocks_origin(self):
        start = sf.start_fen("spell-chess")
        moves = sf.legal_moves("spell-chess", start, [])
        freeze_moves = self._filter_potion_moves(moves, "f")
        self.assertTrue(freeze_moves)
        for move in freeze_moves:
            self.assertIn(",", move)
            self.assertEqual(move[0].lower(), "f")
            self.assertEqual(move[1], "@")

        freeze_on_e2 = [m for m in freeze_moves if self._gating_info(m)[2] == "e2"]
        self.assertTrue(freeze_on_e2)

        for move in freeze_on_e2:
            base, _, _ = self._gating_info(move)
            self.assertNotEqual(base[:2], "e2")

    def test_spell_chess_jump_allows_leap(self):
        fen = "4k3/8/8/8/8/8/P7/R3K3[JJFFFFFjjfffff] w - - 0 1"
        moves = sf.legal_moves("spell-chess", fen, [])
        jump_moves = self._filter_potion_moves(moves, "j")
        self.assertTrue(jump_moves)
        for move in jump_moves:
            self.assertIn(",", move)
            self.assertEqual(move[0].lower(), "j")
            self.assertEqual(move[1], "@")

        leap_moves = [m for m in jump_moves if self._gating_info(m)[0] == "a1a3" and self._gating_info(m)[2] == "a2"]
        self.assertTrue(leap_moves)
        self.assertNotIn("a1a3", moves)

    def test_spell_chess_has_useless_freeze_targets(self):
        start = sf.start_fen("spell-chess")
        moves = sf.legal_moves("spell-chess", start, [])
        freeze_moves = self._filter_potion_moves(moves, "f")
        self.assertTrue(freeze_moves)
        board, width, height = self._board_state(start)
        enemy = "b" if start.split()[1] == "w" else "w"
        useless = [
            move
            for move in freeze_moves
            if not self._freeze_zone_has_enemy(board, width, height, self._gating_info(move)[2], enemy)
        ]
        self.assertTrue(useless)

    def test_spell_chess_has_useless_jump_targets(self):
        start = sf.start_fen("spell-chess")
        history = [
            "e2e4",
            "e7e5",
            "f1c4",
            "f8c5",
            "f@d7,c4f7",
            "f@e6,c5f2",
            "e1f2",
        ]
        moves = sf.legal_moves("spell-chess", start, history)
        jump_moves = self._filter_potion_moves(moves, "j")
        self.assertTrue(jump_moves)
        fen_after = sf.get_fen("spell-chess", start, history)
        board, width, height = self._board_state(fen_after)
        useless = [m for m in jump_moves if not self._jump_uses_gate(m, board, width)]
        self.assertIn("j@a1,e8f7", useless)

    def test_spell_chess_freeze_cooldown(self):
        start = "4k3/8/8/8/8/8/8/4K1N1[JJFFFFFjjfffff] w - - 0 1"
        history = []

        moves = sf.legal_moves("spell-chess", start, history)
        freeze_moves = self._filter_potion_moves(moves, "f")
        self.assertTrue(freeze_moves)
        history.append(freeze_moves[0])

        cooldown_full_moves = 3
        for _ in range(max(cooldown_full_moves - 1, 0)):
            moves = sf.legal_moves("spell-chess", start, history)
            history.append(self._first_normal_move(moves))

            moves = sf.legal_moves("spell-chess", start, history)
            self.assertFalse(self._has_potion_move(moves, "f"))
            history.append(self._first_normal_move(moves))

        moves = sf.legal_moves("spell-chess", start, history)
        history.append(self._first_normal_move(moves))

        moves = sf.legal_moves("spell-chess", start, history)
        self.assertTrue(self._has_potion_move(moves, "f"))

    def test_spell_chess_freeze_cooldown_regression(self):
        history = [
            "d2d4",
            "e7e6",
            "c2c3",
            "c7c6",
            "e2e4",
            "g8f6",
            "f1d3",
            "f8e7",
            "f@e7,e4e5",
            "b8a6",
            "e5f6",
            "e7f6",
            "j@b2,c1a3",
            "d8a5",
            "e1f1",
            "d7d6",
            "d1e2",
            "f6e7",
            "a3d6",
            "a5g5",
            "g1f3",
            "g5f6",
            "d6e5",
            "f@f1,f6f3",
            "f@f4,b1d2",
            "j@h7,h8h2",
            "e2f3",
            "h2h1",
            "j@f1,a1h1",
            "c8d7",
            "h1h7",
            "e8d8",
            "f3f7",
            "d7e8",
            "f@d8,f7e7",
            "f@f6,a8b8",
            "h7g7",
            "d8c8",
            "e7e6",
            "e8d7",
            "e6d7",
        ]

        moves = sf.legal_moves("spell-chess", "startpos", history)
        self.assertIn("f@f6,c8d7", moves)

    def test_spell_chess_jump_capture_wins_immediately(self):
        fen = "5rk1/1p2ppb1/2p1q1p1/3p2Np/3P2n1/3BP3/PPP2PPP/R1B1K2R[JFFFFjffff] b KQ - 5 12"
        result = sf.game_result("spell-chess", fen, ["j@e3,e6e1"])
        self.assertEqual(result, -sf.VALUE_MATE)

    def test_spell_chess_freeze_check_does_not_win(self):
        fen = "rnbqk1nr/pppp1ppp/8/1N2p3/1b6/8/PPPPPPPP/R1BQKBNR[JJFFFFFjjfffff] w KQkq - 2 3"
        result = sf.game_result("spell-chess", fen, ["f@d7,b5c7"])
        self.assertNotEqual(result, sf.VALUE_MATE)
        self.assertNotEqual(result, -sf.VALUE_MATE)

    def test_spell_chess_freeze_check_has_evasion(self):
        fen = "rnbqkbnr/pp1p1ppp/2p1p3/8/4P3/5Q2/PPPP1PPP/RNB1KBNR[JJFFFFFjjfffff] {F@-:0,J@-:0,f@-:0,j@-:0} w KQkq - 0 3"
        moves = sf.legal_moves("spell-chess", fen, ["f@d7,f3f7"])
        self.assertNotIn("f@f7,g8h6", moves)
        self.assertIn("f@f6,g8h6", moves)

    def test_spell_chess_freeze_check_blocks_frozen_knight(self):
        fen = "rnbqkbnr/pp1p1Qpp/2p1p3/8/4P3/8/PPPP1PPP/RNB1KBNR[JJFFFFjjfffff] {F@f8:2,J@-:0,f@-:0,j@-:0} b KQkq - 0 3"
        moves = sf.legal_moves("spell-chess", fen, [])
        self.assertNotIn("g8h6", moves)

    def test_spell_chess_freeze_zone_does_not_force_mate(self):
        fen = "rnbqkbnr/pp1p1Qpp/2p1p3/8/4P3/8/PPPP1PPP/RNB1KBNR[JJFFFFjjfffff] {F@e8:2,J@-:0,f@-:0,j@-:0} b KQkq - 0 3"
        moves = sf.legal_moves("spell-chess", fen, [])
        self.assertIn("g8h6", moves)
        self.assertNotIn("f@f7,g8h6", moves)
        result = sf.game_result("spell-chess", fen, [])
        self.assertNotEqual(result, sf.VALUE_MATE)
        self.assertNotEqual(result, -sf.VALUE_MATE)

    def test_spell_chess_freeze_zone_defense_prevents_mate(self):
        fen = "rnbqkbnr/pp1p1Qpp/2p5/4p3/4P3/8/PPPP1PPP/RNB1KBNR[JJFFFFjjfffff] {F@f7:2,J@-:0,f@-:0,j@-:0} b KQkq - 0 3"
        moves = sf.legal_moves("spell-chess", fen, [])
        self.assertNotIn("g8h6", moves)
        self.assertIn("f@f6,d8f6", moves)
        result = sf.game_result("spell-chess", fen, [])
        self.assertNotEqual(result, sf.VALUE_MATE)
        self.assertNotEqual(result, -sf.VALUE_MATE)

    def test_spell_chess_freeze_mate_in_two_threat(self):
        fen = "rnbqkb1r/pppppppp/8/1n2P3/8/8/PPPP1PPP/R1BQK1NR[JJFFFFFjjffff] {F@-:0,J@-:0,f@-:2,j@-:0} w KQkq - 0 5"
        moves = sf.legal_moves("spell-chess", fen, [])
        self.assertIn("f@g7,d1h5", moves)
        result = sf.game_result("spell-chess", fen, ["f@g7,d1h5", "a7a6", "h5f7"])
        self.assertGreaterEqual(result, sf.VALUE_MATE)

    def test_spell_chess_freeze_zone_blocks_starting_inside(self):
        fen = "r3kbnr/pp1n1ppp/2p1p3/3pP1B1/3P4/1Q3N2/PqP2PPP/RN3RK1[JJFFFFjjffff] {F@-:0,J@-:0,f@-:0,j@-:0} b kq - 1 9"
        moves = sf.legal_moves("spell-chess", fen, [])
        self.assertNotIn("f@b1,b2b3", moves)

    def test_spell_chess_freeze_zone_blocks_diagonal_start(self):
        history = [
            "f@d7,e2e4", "a7a6", "e4e5", "e7e6", "f1e2", "b8c6",
            "b1c3", "c6e5", "d2d4", "f@d2,f8b4", "f@a4,d4e5", "d8h4",
            "e1f1", "j@h2,h4h1", "e2h5", "g7g6", "h5f3", "h1h2",
            "g1h3", "j@h7,h8h3", "g2h3", "b4c5", "c3e4", "c5b6",
            "d1e2", "h2h3", "f3g2", "h3h4", "f@d8,e4d6", "f@d5,g8e7",
            "e2g4", "c7d6", "g4h4", "d6e5", "h4f6", "d7d5",
        ]
        moves = sf.legal_moves("spell-chess", "startpos", history)
        self.assertNotIn("f@e7,f6h8", moves)

    def test_spell_chess_jump_gate_capture_allowed(self):
        history = [
            "e2e4", "b8c6", "f1b5", "e7e5", "f@c7,b5c6", "f@c5,f8d6", "d1h5", "e8e7",
            "j@h7,h5h8", "j@g8,d8h8", "c6d5", "g8f6", "c2c3", "f6d5", "e4d5", "h8e8",
            "d2d3", "f7f6", "g1f3", "f@d2,e7f8", "f@e7,f3d2", "e5e4", "e1f1",
            "e4d3", "d2c4", "e8e2", "j@f1,h1e1"
        ]
        moves = sf.legal_moves("spell-chess", "startpos", history)
        self.assertIn("e2f1", moves)

    def test_spell_chess_castling_illegal_while_in_check(self):
        moves = [
            "e2e4", "e7e5", "f1c4", "f8c5", "d1e2",
            "g8f6", "c2c3", "d7d5", "f@d7,c4b5",
            "f@a5,a7a6", "j@a2,a1a6"
        ]
        fen = sf.get_fen("spell-chess", "startpos", moves)
        legal = sf.legal_moves("spell-chess", fen, [])
        self.assertNotIn("e8g8", legal)

    def test_spell_chess_potion_consumes_hand(self):
        start = sf.start_fen("spell-chess")
        moves = sf.legal_moves("spell-chess", start, [])
        freeze_moves = self._filter_potion_moves(moves, "f")
        self.assertTrue(freeze_moves)

        start_pocket = start.split()[0]
        start_pocket = start_pocket[start_pocket.index('[') + 1:start_pocket.index(']')]
        start_f = start_pocket.count('F')
        start_j = start_pocket.count('J')
        start_f_black = start_pocket.count('f')
        start_j_black = start_pocket.count('j')

        fen_after = sf.get_fen("spell-chess", start, [freeze_moves[0]])
        after_pocket = fen_after.split()[0]
        after_pocket = after_pocket[after_pocket.index('[') + 1:after_pocket.index(']')]
        self.assertEqual(after_pocket.count('F'), max(start_f - 1, 0))
        self.assertEqual(after_pocket.count('J'), start_j)
        self.assertEqual(after_pocket.count('f'), start_f_black)
        self.assertEqual(after_pocket.count('j'), start_j_black)

    def test_spell_chess_jump_pawn_double_step(self):
        fen = "r1b5/pp1p1k2/1qpPp1pB/8/2Bb4/2P5/PPQ2PPP/R3K2R[JJFFjjffff] {F@-:0,J@-:0,f@-:0,j@-:0} b KQ - 3 16"
        moves = sf.legal_moves("spell-chess", fen, [])
        self.assertIn("j@d6,d7d5", moves)

    def test_spell_chess_castling_blocked_by_unfreeze(self):
        fen = "r7/pp1b1k2/2pPp1pB/3p4/1qBP4/P7/1PQ2PPP/R3K2R[JJFjfff] {F@-:2,J@-:0,f@b2:2,j@-:0} w KQ - 1 18"
        moves = sf.legal_moves("spell-chess", fen, [])
        self.assertNotIn("e1g1", moves)

    def test_spell_chess_jump_zone_persists_two_plies(self):
        fen = "r4kr1/pp3p1p/3P2p1/2np4/8/5N1b/PPP2PPP/R2Q1BK1[JFFjjfff] {F@-:0,J@-:0,f@-:0,j@-:0} w - - 1 19"
        moves = sf.legal_moves("spell-chess", fen, ["j@g2,f1h3"])
        self.assertIn("j@g6,g8g1", moves)

    def test_spell_chess_jump_zone_expires_after_two_plies(self):
        fen = ("2rqk2r/pp2nppp/8/1p1nQ3/1P1nP3/3PBP2/PP3P1P/"
               "R3K2R[JFFFjjfff] {F@-:1,J@d4:2,f@-:0,j@-:0} b KQk - 1 13")
        moves = sf.legal_moves("spell-chess", fen, ["g7g6"])
        self.assertIn("e3d4", moves)

    def test_spell_chess_jump_zone_expired_allows_capture(self):
        fen = "2rqk2r/pp2nppp/8/1p2Q3/1P1nP3/3PBP2/PP3P1P/R3K2R[JFFFjjfff] {F@-:1,J@d4:2,f@-:0,j@-:0} w KQk - 1 14"
        moves = sf.legal_moves("spell-chess", fen, [])
        self.assertIn("e3d4", moves)

    def test_spell_chess_freeze_zone_expires_after_two_plies(self):
        fen = ("rnbqk2r/1pp2ppp/R4n2/1Bbpp3/4P3/2P5/PP1PQPPP/"
               "1NB1K1NR[JFFFFjjffff] {F@-:1,J@a2:2,f@a4:3,j@-:0} w Kkq - 0 6")
        fen_after = sf.get_fen("spell-chess", fen, ["h2h3"])
        state = fen_after[fen_after.index('{') + 1:fen_after.index('}')]
        self.assertIn("f@-:2", state)
        moves = sf.legal_moves("spell-chess", fen_after, [])
        self.assertNotIn("e8g8", moves)

    def test_spell_chess_freeze_zone_history_allows_mate(self):
        fen = "rnbqkbnr/pp1p1Qpp/2p5/4p3/4P3/8/PPPP1PPP/RNB1KBNR[JJFFFFjjfffff] {F@f8:2,J@-:0,f@-:0,j@-:0} b KQkq - 0 3"
        moves = sf.legal_moves("spell-chess", fen, ["d7d5"])
        self.assertIn("f7e8", moves)
        result = sf.game_result("spell-chess", fen, ["d7d5", "f7e8"])
        self.assertEqual(result, -sf.VALUE_MATE)

    def test_spell_chess_capture_commoner_in_check(self):
        fen = "4k3/4b3/8/8/8/8/4R3/4K3[JJFFFFjjffff] b - - 0 1"
        moves = sf.legal_moves("spell-chess", fen, [])
        self.assertIn("e7b4", moves)
        moves = sf.legal_moves("spell-chess", fen, ["e7b4"])
        self.assertIn("e2e8", moves)
        result = sf.game_result("spell-chess", fen, ["e7b4", "e2e8"])
        self.assertEqual(result, -sf.VALUE_MATE)


if __name__ == '__main__':
    unittest.main()
