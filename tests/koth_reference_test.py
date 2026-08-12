import importlib.util
import json
import pathlib
import subprocess
import sys
import unittest

import chess


ROOT = pathlib.Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "koth_reference", ROOT / "tools" / "koth_reference.py"
)
assert SPEC and SPEC.loader
REFERENCE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = REFERENCE
SPEC.loader.exec_module(REFERENCE)


def request(root_fen, *moves):
    return {
        "schema": "koth-reference-input-v1",
        "rule_profile": "KOTH_LICHESS_V1",
        "root_fen": root_fen,
        "moves": list(moves),
    }


class KothReferenceTest(unittest.TestCase):
    def test_every_goal_and_color(self):
        cases = (
            ("7k/8/8/8/8/2K5/8/8 w - - 0 1", "c3d4", "white"),
            ("7k/8/8/8/8/3K4/8/8 w - - 0 1", "d3e4", "white"),
            ("7k/8/8/8/2K5/8/8/8 w - - 0 1", "c4d5", "white"),
            ("7k/8/8/8/5K2/8/8/8 w - - 0 1", "f4e5", "white"),
            ("8/8/2k5/8/8/8/8/K7 b - - 0 1", "c6d5", "black"),
            ("8/8/3k4/8/8/8/8/K7 b - - 0 1", "d6e5", "black"),
            ("8/8/8/2k5/8/8/8/K7 b - - 0 1", "c5d4", "black"),
            ("8/8/8/5k2/8/8/8/K7 b - - 0 1", "f5e4", "black"),
        )
        for fen, move, winner in cases:
            with self.subTest(move=move):
                result = REFERENCE.replay(request(fen, move))
                self.assertTrue(result["accepted"])
                self.assertEqual(result["terminal"]["primary_reason"], "HILL")
                self.assertEqual(result["terminal"]["winner"], winner)
                self.assertEqual(result["plies"][0]["san"][-1], "#")
                self.assertEqual(result["plies"][0]["game_legal_move_count"], 0)

    def test_every_goal_entry_attacked_for_both_colors_is_illegal(self):
        cases = (
            ("3r3k/8/8/8/8/2K5/8/8 w - - 0 1", "c3d4"),
            ("4r2k/8/8/8/8/3K4/8/8 w - - 0 1", "d3e4"),
            ("3r3k/8/8/8/2K5/8/8/8 w - - 0 1", "c4d5"),
            ("4r2k/8/8/8/5K2/8/8/8 w - - 0 1", "f4e5"),
            ("8/8/2k5/8/8/8/8/K2R4 b - - 0 1", "c6d5"),
            ("8/8/3k4/8/8/8/8/K3R3 b - - 0 1", "d6e5"),
            ("8/8/8/2k5/8/8/8/K2R4 b - - 0 1", "c5d4"),
            ("8/8/8/5k2/8/8/8/K3R3 b - - 0 1", "f5e4"),
        )
        for fen, move in cases:
            with self.subTest(move=move), self.assertRaises(
                REFERENCE.ReferenceError
            ) as caught:
                REFERENCE.replay(request(fen, move))
            self.assertEqual(caught.exception.code, "KREF_E_MOVE_ILLEGAL")

    def test_non_king_on_goal_is_not_terminal(self):
        result = REFERENCE.replay(
            request("7k/8/8/8/8/8/2N5/K7 w - - 0 1", "c2d4")
        )
        self.assertFalse(result["terminal"]["terminal"])
        self.assertEqual(result["plies"][0]["san"], "Nd4")

    def test_loaded_goal_is_rejected_for_both_turns(self):
        for turn in ("w", "b"):
            with self.subTest(turn=turn), self.assertRaises(
                REFERENCE.ReferenceError
            ) as caught:
                REFERENCE.replay(
                    request(f"7k/8/8/8/3K4/8/8/8 {turn} - - 1 1")
                )
            self.assertEqual(caught.exception.code, "KREF_E_AMBIGUOUS_GOAL_ROOT")

    def test_move_after_terminal_is_rejected(self):
        with self.assertRaises(REFERENCE.ReferenceError) as caught:
            REFERENCE.replay(
                request(
                    "7k/8/8/8/8/2K5/8/8 w - - 0 1",
                    "c3d4",
                    "h8h7",
                )
            )
        self.assertEqual(caught.exception.code, "KREF_E_MOVE_AFTER_TERMINAL")

    def test_fifty_move_root_is_rejected(self):
        with self.assertRaises(REFERENCE.ReferenceError) as caught:
            REFERENCE.replay(request("7k/8/8/8/8/2K5/8/8 w - - 100 1"))
        self.assertEqual(caught.exception.code, "KREF_E_TERMINAL_ROOT")

    def test_checkmate_and_hill_predicates_are_both_retained(self):
        result = REFERENCE.replay(
            request("2k2N2/5N2/8/8/8/2K5/8/1RR5 w - - 0 1", "c3d4")
        )
        terminal = result["terminal"]
        self.assertEqual(terminal["predicates"], ("CHECKMATE", "HILL"))
        self.assertEqual(terminal["primary_reason"], "CHECKMATE")
        self.assertEqual(result["plies"][0]["san"], "Kd4#")

    def test_hill_precedes_stalemate_and_fifty_move(self):
        stalemate = REFERENCE.replay(
            request("k1B5/8/2N5/8/8/2K5/8/8 w - - 0 1", "c3d4")
        )
        self.assertEqual(stalemate["terminal"]["predicates"], ("HILL", "STALEMATE"))
        self.assertEqual(stalemate["terminal"]["primary_reason"], "HILL")

        fifty = REFERENCE.replay(
            request("7k/8/8/8/8/2K5/8/8 w - - 99 1", "c3d4")
        )
        self.assertEqual(fifty["terminal"]["predicates"], ("HILL", "FIFTY_MOVE"))
        self.assertEqual(fifty["terminal"]["primary_reason"], "HILL")

    def test_threefold_ends_an_honest_replay(self):
        result = REFERENCE.replay(
            request(
                chess.STARTING_FEN,
                "g1f3",
                "g8f6",
                "f3g1",
                "f6g8",
                "g1f3",
                "g8f6",
                "f3g1",
                "f6g8",
            )
        )
        self.assertEqual(result["terminal"]["predicates"], ("THREEFOLD",))
        self.assertEqual(result["terminal"]["primary_reason"], "AUTOMATIC_DRAW")

    def test_synthetic_fivefold_sets_both_repetition_predicates(self):
        board = chess.Board()
        for _ in range(4):
            for token in ("g1f3", "g8f6", "f3g1", "f6g8"):
                board.push_uci(token)
        terminal = REFERENCE.terminal_status(
            board, hill_transition=False, mover=chess.BLACK
        )
        self.assertEqual(terminal.predicates, ("THREEFOLD", "FIVEFOLD"))
        self.assertEqual(terminal.primary_reason, "AUTOMATIC_DRAW")

    def test_cli_consumes_utf8_json_and_emits_canonical_json(self):
        payload = json.dumps(
            request("7k/8/8/8/8/2K5/8/8 w - - 0 1", "c3d4")
        ).encode("utf-8")
        completed = subprocess.run(
            [sys.executable, str(ROOT / "tools" / "koth_reference.py"), "replay"],
            input=payload,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr.decode("utf-8"))
        decoded = json.loads(completed.stdout.decode("utf-8"))
        self.assertTrue(decoded["accepted"])
        canonical = json.dumps(decoded, sort_keys=True, separators=(",", ":")) + "\n"
        self.assertEqual(completed.stdout.decode("utf-8"), canonical)


if __name__ == "__main__":
    unittest.main()
