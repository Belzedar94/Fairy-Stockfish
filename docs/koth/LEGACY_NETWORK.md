# Official Legacy Compatibility Network

The project owner trained and authorized the following exact network for local
KOTH-Stockfish engineering and compatibility testing:

| Field | Value |
| --- | --- |
| Historic filename | `kingofthehill-978b86d0e6a4.nnue` |
| Reserved release alias | `KOTH_v1.nnue` |
| Size | `47,721,371` bytes |
| SHA-256 | `978B86D0E6A45E05F9F1375DCED129CEA0ACEA13041EA65960691632EDC47AF7` |
| Default network | No |
| Bytes in this repository | No |

The alias is valid only when its bytes are identical to the historic network.
Renaming another file to `KOTH_v1.nnue` must fail authentication. A matching
container header or file length is not sufficient.

## Distribution status

No standard public redistribution or modification license has been assigned to
these weights. Public availability, a historic download, and the engine's GPL
license do not substitute for a network license. Consequently:

- the network is not committed, vendored, cached, or bundled here;
- CI and release jobs must not download it from an unofficial mirror;
- contributors may test only with a copy they are independently authorized to
  use;
- a release cannot bundle the bytes until the license gate is resolved.

## Loader status

The pinned upstream loader does not yet satisfy the project contract. In the
discovery binary, the reserved alias and unrelated basenames could silently
fall back to classical evaluation, while a same-size payload mutation was not
authenticated. The engineering gate therefore remains open.

Certification requires positive loading of both authorized names and
fail-closed behavior for missing, wrong-name, wrong-hash, corrupt, truncated,
incompatible, and extra-byte inputs. It also requires full-byte authentication
before activation and proof that the alias does not mutate the default.
