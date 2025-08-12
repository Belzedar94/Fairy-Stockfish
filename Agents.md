1. For most variants, you will only have to add a new entry to variants.ini
2. Check to see if the options to create your variant already exist.
   * For rule variants, see variants.ini -> "Rule definition options"
   * For custom pieces, see variants.ini -> "Custom pieces"
3. Most variants are added to variants.ini, but some are added to variant.cpp (generally those that are very popular or the root of a large family of games)
4. If your game rules aren't supported:
    1. Break it down into what changes will be required. For instance, "Connect 4" can be broken down into settings such as:
            enclosingDrop = top
            connectN = 4

        rather than:
            playConnect4 = true

    2. In variant.h -> "struct Variant", add a variable for each new setting
    3. In parser.cpp -> "parse(Variant* v)", add "parse_attribute" to read from variants.ini into the variable.
    4. In position.h, most settings have a getter with a snake_case name which reads the camelCase variable. Remember to declare the getter under "class Position".
    5. If your rule changes how moves occur, see:
        1. movegen.cpp->“make_move_and_gating(const Position& pos, ExtMove* moveList, Color us, Square from, Square to, PieceType pt = NO_PIECE_TYPE)”
        2. position.cpp -> "do_move(Move m, StateInfo& newSt, bool givesCheck)"
        3. position.cpp->”legal(Move m)”
        4. position.cpp->”pseudo_legal(const Move m)”
        5. position.cpp->”undo_move(Move m)”
    6. If your rule includes a new Betza modifier, piece.cpp->”from_betza(const std::string& betza, const std::string& name)”
    7. If your rule adjudicates win, loss, or draw, see position.cpp -> "is_immediate_game_end(Value& result, int ply)" or "is_optional_game_end(Value& result, int ply, int countStarted)"
    8. Document your new setting in the comment section at the beginning on variants.ini
    9. See https://github.com/fairy-stockfish/Fairy-Stockfish/wiki/Understanding-the-code for more descriptions.
    10. Do not break backwards compatibility with previous .ini files that people may have. If you are making a more flexible version of a previous configuration option, read the previous option into the new option with defaults that reflect the old behaviour.
    11. Performance and portability are top priorities. Try to avoid excessive dependencies, or unnecessary intensive computation inside main loops. It is ok in the main loop to have code that is a bit harder to read if there is better performance.
    12. C++ Code is indented 2 spaces for the first function level, then 4 spaces each level after that.
    13. Don't change variable names unless specifically required. We know that in some cases local variables overshadow global ones; we are ok with that.
    14. You must compile and test that your feature is read all the way from variants.ini (or your test .ini) and used by the engine.
    15. Comments should be appropriate for experienced developers. Don't change copyright year in comments.
    16. Bitwise operators are overloaded between Squares and Bitboards in bitboard.h; you don't have to explicitly convert in most cases.
5. Testing your variant
    1. Compile using “make”. Type “make help” to see important options. Remember that you need to run "make" from the "src" folder, the executable will be created there, and the default location for "variants.ini" is in there.
        1. A simple compile looks like: "make -j build ARCH=x86-64-modern"
        2. “largeboards=yes” if your board is greater than 8x8. Boards greater than 10x12 are not supported.
        3. “all=yes” if your variant has a large branching factor such as “Duck Chess” or “Game of the Amazons”
        4. “debug=yes” and “optimize=no” if you’re troubleshooting
    2. Once you launch the executable, here are useful options:
        1. “setoption name VariantPath value variants.ini”
        2. “setoption name UCI_Variant value [your_variant]”
    3. To play moves, start with “position”, either “startpos” or “fen” then use coordinates to play your moves. ie.
        1. “position startpos moves e2e4 e7e5”
        2. “position fen 4k3/8/8/8/8/8/p7/4K2R w K - 0 1 moves e1g1 a2a1q”
        3. Castling is encoded as “king moves two spaces”
        4. Promotion piece is indicated with just letter (no “=”)
        5. Drops are encoded with “@”. For instance, tictactoe looks like:
            position startpos moves P@b2 P@a1 P@c1
        6. “d” (prints current board using text characters)
        7. “go movetime [milliseconds]”
        8. “go depth [ply]”
        9. “quit”
    4. If you want to automate your testing, you can use redirection from the command line. There are other acceptable ways to test as well: 
        Create test.txt:

            position startpos moves e2e4 d7d5
            go movetime 100
            d
            quit

        Run stockfish:

            stockfish < test.txt > output.txt
    5. Run "stockfish check variants.ini" if you edited variants.ini
    6. Common mistakes in testing:
        1. If your variant uses standard kings, make sure your FEN includes both kings, and that they are not already in checkmate or stalemate. Also, if a player is in check, the move sequence may be forcing, rather than allow your intended move.
