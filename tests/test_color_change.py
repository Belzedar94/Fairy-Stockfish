import subprocess
import unittest
from pathlib import Path


ENGINE_PATH = Path(__file__).resolve().parent.parent / "src" / "stockfish"


def run_engine(commands):
    """Run the engine with the provided list of UCI commands."""

    assert ENGINE_PATH.exists(), f"Expected engine binary at {ENGINE_PATH}"
    script = "\n".join(commands) + "\n"
    result = subprocess.run(
        [str(ENGINE_PATH)],
        input=script,
        text=True,
        capture_output=True,
        check=True,
    )
    return result.stdout


class ColorChangeRegressionTests(unittest.TestCase):

    def test_benedict_conversion_no_false_mate(self):
        output = run_engine(
            [
                "uci",
                "setoption name UCI_Variant value benedict",
                "isready",
                "position startpos moves e2e3 e7e6 d1e2 h7h5 e2b5 f7f5 b5e5",
                "go depth 2",
                "quit",
            ]
        )

        self.assertIn("bestmove", output)
        self.assertNotIn("score mate", output)

    def test_benedict_mass_conversion_stability(self):
        output = run_engine(
            [
                "uci",
                "setoption name UCI_Variant value benedict",
                "isready",
                "position fen rnbqkbnr/pppppppp/pppppppp/pppppppp/3Q4/PPPPPPPP/PPPPPPPP/RNBQKBNR w KQkq - 0 1 moves d4e5",
                "go depth 1",
                "quit",
            ]
        )

        self.assertIn("bestmove", output)
        self.assertNotIn("bestmove (none)", output)
        self.assertNotIn("Unknown command", output)

    def test_antiandernach_short_search(self):
        output = run_engine(
            [
                "uci",
                "setoption name UCI_Variant value antiandernach",
                "isready",
                "position startpos moves e2e4 e7e5 g1f3 b8c6 f1c4 g8f6",
                "go depth 2",
                "quit",
            ]
        )

        self.assertIn("bestmove", output)
        self.assertNotIn("bestmove (none)", output)


if __name__ == "__main__":
    unittest.main()
