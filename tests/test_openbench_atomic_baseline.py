#!/usr/bin/env python3

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
MAKEFILE = (ROOT / "src" / "Makefile").read_text(encoding="utf-8")
BENCHMARK = (ROOT / "src" / "benchmark.cpp").read_text(encoding="utf-8")
EVALUATE = (ROOT / "src" / "evaluate.h").read_text(encoding="utf-8")
HARNESS = MAKEFILE.split("### Section 6. Frozen Atomic OpenBench harness", 1)[1]


class OpenBenchAtomicBaselineContractTests(unittest.TestCase):

    def test_evalfile_only_replaces_the_default_goal_for_openbench(self):
        self.assertRegex(
            HARNESS,
            r"ifneq \(\$\(strip \$\(EVALFILE\)\),\)\s*"
            r"\.DEFAULT_GOAL := openbench-atomic-baseline\s*endif\s*$",
        )

    def test_windows_and_linux_use_the_normative_compilers(self):
        mapping = re.search(
            r"ifeq \(\$\(OS\),Windows_NT\)\s*"
            r"OPENBENCH_ATOMIC_COMP = (?P<windows>\S+).*?\s*else\s*"
            r"OPENBENCH_ATOMIC_COMP = (?P<linux>\S+)",
            HARNESS,
            re.DOTALL,
        )
        self.assertIsNotNone(mapping)
        self.assertEqual(
            mapping.groupdict(), {"windows": "mingw", "linux": "gcc"}
        )

    def test_native_windows_replays_the_frozen_lto_link_contract(self):
        windows = re.search(
            r"ifeq \(\$\(OS\),Windows_NT\)(?P<body>.*?)\s*else\s*"
            r"OPENBENCH_ATOMIC_COMP = gcc",
            HARNESS,
            re.DOTALL,
        )
        self.assertIsNotNone(windows)
        body = windows.group("body")
        self.assertIn(
            "OPENBENCH_ATOMIC_CXXFLAGS = -flto -flto-partition=one", body
        )
        self.assertIn(
            "OPENBENCH_ATOMIC_LDFLAGS = -flto -flto-partition=one "
            "-save-temps -Wl,--no-insert-timestamp",
            body,
        )
        self.assertIn(
            "EXTRACXXFLAGS='$(EXTRACXXFLAGS) "
            "$(OPENBENCH_ATOMIC_CXXFLAGS) -DOPENBENCH_ATOMIC_BASELINE'",
            HARNESS,
        )
        self.assertIn(
            "EXTRALDFLAGS='$(EXTRALDFLAGS) $(OPENBENCH_ATOMIC_LDFLAGS)'",
            HARNESS,
        )

    def test_public_build_is_the_frozen_small_bmi2_configuration(self):
        self.assertRegex(
            HARNESS,
            r"\+\$\(MAKE\) all EXE=\"\$\(EXE\)\" CXX=\"\$\(CXX\)\"\s*\\\s*"
            r"ARCH=x86-64-bmi2 COMP=\$\(OPENBENCH_ATOMIC_COMP\)\s*\\\s*"
            r"all=no largeboards=no nnue=yes",
        )
        self.assertIn("-DOPENBENCH_ATOMIC_BASELINE", HARNESS)

    def test_authenticated_network_uses_fairys_canonical_embedded_name(self):
        self.assertRegex(
            EVALUATE,
            r"#ifdef OPENBENCH_ATOMIC_BASELINE\s*"
            r'#define EvalFileDefaultName\s+"atomic_run3b_e202_l05\.nnue"\s*'
            r"#else\s*"
            r'#define EvalFileDefaultName\s+"nn-3475407dc199\.nnue"\s*'
            r"#endif",
        )
        canonical = "atomic_run3b_e202_l05.nnue"
        self.assertIn(f"OPENBENCH_ATOMIC_NET = {canonical}", HARNESS)
        self.assertIn('cp "$(EVALFILE)" "$(OPENBENCH_ATOMIC_NET).tmp"', HARNESS)
        self.assertIn(
            'mv "$(OPENBENCH_ATOMIC_NET).tmp" "$(OPENBENCH_ATOMIC_NET)"',
            HARNESS,
        )

    def test_harness_only_changes_missing_bench_defaults(self):
        self.assertIn(
            """#ifdef OPENBENCH_ATOMIC_BASELINE
      // The frozen OpenBench harness must validate the same Atomic NNUE path
      // used by its games. Explicit benchmark arguments remain unchanged.
      varname = "atomic";
#else""",
            BENCHMARK,
        )
        self.assertIn(
            """#ifdef OPENBENCH_ATOMIC_BASELINE
                                   : "NNUE";
#else""",
            BENCHMARK,
        )


if __name__ == "__main__":
    unittest.main()
