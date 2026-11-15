#!/usr/bin/env python3
"""FeatherCrypt differential-characteristic CVC generator for STP.

This script models reduced-round FeatherCrypt under the same-key XOR
 differential setting.  In this setting the initial key addition and every
round-key XOR cancel in the difference, so the differential state transition is

    X_r -- 16 parallel 4-bit S-boxes --> S_r -- bit permutation --> X_{r+1}

The FeatherCrypt paper gives the bit permutation family

    pi(i) = (r * i) mod 64,  i in [0, 63],

and gives r=17 as its concrete example.  Because the paper uses wording such as
"e.g., r=17" rather than a separate fixed permutation table, this generator uses
17 by default but exposes --perm-multiplier so the value can be changed if an
official implementation or author clarification specifies another value.  The
paper also does not explicitly state whether bit 0 is the LSB or MSB, so this
generator supports both interpretations through --bit-numbering.

The S-box DDT is generated automatically from the published S-box.  Impossible
DDT transitions are explicitly forbidden; possible transitions receive exact
integer weights w = -log2(p).  For the published 4-bit S-box the only DDT
counts are 0, 2, 4 and 16, hence the weights are 3, 2 and 0.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable, List, Optional, Sequence


# -----------------------------------------------------------------------------
# FeatherCrypt specification (IJETT, Vol. 74, Issue 7, July 2026)
# -----------------------------------------------------------------------------
BLOCK_BITS = 64
SBOX_BITS = 4
NUM_SBOXES = BLOCK_BITS // SBOX_BITS
FULL_ROUNDS = 16
DEFAULT_PERM_MULTIPLIER = 17
WEIGHT_BITS = 16

# Published mapping:
# 0->C, 1->5, 2->6, 3->B, 4->9, 5->0, 6->A, 7->D,
# 8->3, 9->E, A->F, B->8, C->4, D->7, E->1, F->2
SBOX: List[int] = [
    0xC, 0x5, 0x6, 0xB,
    0x9, 0x0, 0xA, 0xD,
    0x3, 0xE, 0xF, 0x8,
    0x4, 0x7, 0x1, 0x2,
]


# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------
def compute_ddt(sbox: Sequence[int]) -> List[List[int]]:
    """Return the XOR difference-distribution table of an n-bit S-box."""
    size = len(sbox)
    n = int(math.log2(size))
    if 1 << n != size:
        raise ValueError("S-box size must be a power of two")

    ddt = [[0 for _ in range(size)] for _ in range(size)]
    for din in range(size):
        for x in range(size):
            dout = sbox[x] ^ sbox[x ^ din]
            ddt[din][dout] += 1
    return ddt


def validate_spec(ddt: Sequence[Sequence[int]], perm_multiplier: int) -> None:
    """Sanity-check the published S-box and the permutation."""
    if sorted(SBOX) != list(range(16)):
        raise ValueError("Published FeatherCrypt S-box is not bijective")

    if any(sum(row) != 16 for row in ddt):
        raise ValueError("Invalid DDT: each row must sum to 16")

    differential_uniformity = max(max(row) for row in ddt[1:])
    if differential_uniformity != 4:
        raise ValueError(
            f"Unexpected differential uniformity: {differential_uniformity} (expected 4)"
        )

    counts = sorted({v for row in ddt for v in row})
    if counts != [0, 2, 4, 16]:
        raise ValueError(
            "This generator currently expects DDT counts {0,2,4,16}; "
            f"observed {counts}"
        )

    perm = paper_permutation(perm_multiplier)
    if sorted(perm) != list(range(BLOCK_BITS)):
        raise ValueError(f"pi(i)=({perm_multiplier}*i) mod 64 is not a permutation")


def paper_permutation(perm_multiplier: int) -> List[int]:
    """Return scatter permutation pi: input paper-bit i -> output paper-bit pi(i)."""
    if math.gcd(perm_multiplier, BLOCK_BITS) != 1:
        raise ValueError(
            f"perm-multiplier={perm_multiplier} is not coprime to {BLOCK_BITS}"
        )
    return [(perm_multiplier * i) % BLOCK_BITS for i in range(BLOCK_BITS)]


def inverse_permutation(perm: Sequence[int]) -> List[int]:
    """For scatter perm[out=perm[in]], return src_for_output[out] = in."""
    inv = [-1] * len(perm)
    for src, dst in enumerate(perm):
        if inv[dst] != -1:
            raise ValueError("Permutation contains duplicate destinations")
        inv[dst] = src
    if any(v < 0 for v in inv):
        raise ValueError("Permutation is incomplete")
    return inv


def stp_source_for_output(bit_numbering: str, perm_multiplier: int) -> List[int]:
    """Map each STP output bit index to the STP source bit index.

    FeatherCrypt says input paper-bit i is moved to paper-bit pi(i).
    STP numbers bits from the LSB: bit 0 is the rightmost bit.

    bit_numbering='lsb0': paper bit i == STP bit i.
    bit_numbering='msb0': paper bit i == STP bit 63-i.
    """
    inv = inverse_permutation(paper_permutation(perm_multiplier))

    if bit_numbering == "lsb0":
        return inv

    if bit_numbering == "msb0":
        src_for_out = [-1] * BLOCK_BITS
        for stp_out in range(BLOCK_BITS):
            paper_out = (BLOCK_BITS - 1) - stp_out
            paper_in = inv[paper_out]
            stp_in = (BLOCK_BITS - 1) - paper_in
            src_for_out[stp_out] = stp_in
        return src_for_out

    raise ValueError(f"Unknown bit-numbering convention: {bit_numbering}")


def cvc_hex(value: int, bits: int) -> str:
    """Format a non-negative integer as an STP CVC hexadecimal bitvector."""
    if value < 0 or value >= (1 << bits):
        raise ValueError(f"Value {value} does not fit in {bits} bits")
    digits = (bits + 3) // 4
    return f"0hex{value:0{digits}x}"


def parse_diff(text: Optional[str]) -> Optional[int]:
    """Parse an optional 64-bit hexadecimal/decimal difference."""
    if text is None:
        return None
    value = int(text, 0)
    if not 0 <= value < (1 << BLOCK_BITS):
        raise argparse.ArgumentTypeError("difference must fit in 64 bits")
    return value


def concat_nibbles(names: Iterable[str]) -> str:
    """Concatenate nibble variables supplied from MS nibble to LS nibble."""
    return "@".join(names)


def permutation_expression(input_var: str, bit_numbering: str, perm_multiplier: int) -> str:
    """Build an STP expression implementing the FeatherCrypt P-layer.

    @ places the left operand in the more-significant side, so output bits are
    emitted from bit 63 down to bit 0.
    """
    src_for_out = stp_source_for_output(bit_numbering, perm_multiplier)
    return "@".join(
        f"{input_var}[{src_for_out[out_bit]}:{src_for_out[out_bit]}]"
        for out_bit in range(BLOCK_BITS - 1, -1, -1)
    )


def ddt_weight(count: int) -> int:
    """Convert a non-zero 4-bit S-box DDT count to -log2(count/16)."""
    if count <= 0:
        raise ValueError("DDT=0 is an impossible transition and has no finite weight")
    if count & (count - 1):
        raise ValueError(
            "Non-power-of-two DDT count encountered; integer-weight model is insufficient"
        )
    return SBOX_BITS - int(math.log2(count))


def print_ddt(ddt: Sequence[Sequence[int]]) -> None:
    print("FeatherCrypt S-box DDT (rows=Delta_in, columns=Delta_out):")
    print("     " + " ".join(f"{i:X}" for i in range(16)))
    for i, row in enumerate(ddt):
        print(f"{i:X}: " + " ".join(f"{v:2d}" for v in row))


def print_permutation(bit_numbering: str, perm_multiplier: int) -> None:
    perm = paper_permutation(perm_multiplier)
    src = stp_source_for_output(bit_numbering, perm_multiplier)
    print("Paper scatter mapping i -> pi(i):")
    for i in range(0, BLOCK_BITS, 8):
        print("  " + ", ".join(f"{j:02d}->{perm[j]:02d}" for j in range(i, i + 8)))
    print(f"STP interpretation: {bit_numbering}")
    print("STP output bit -> STP input bit:")
    for i in range(0, BLOCK_BITS, 8):
        print("  " + ", ".join(f"{j:02d}<-{src[j]:02d}" for j in range(i, i + 8)))


# -----------------------------------------------------------------------------
# CVC generation
# -----------------------------------------------------------------------------
def generate_cvc(
    output_path: Path,
    rounds: int,
    weight_bound: int,
    bound_mode: str,
    bit_numbering: str,
    perm_multiplier: int,
    input_diff: Optional[int],
    output_diff: Optional[int],
    ddt: Sequence[Sequence[int]],
) -> None:
    if not 1 <= rounds <= FULL_ROUNDS:
        raise ValueError(f"rounds must be between 1 and {FULL_ROUNDS}")
    if not 0 <= weight_bound < (1 << WEIGHT_BITS):
        raise ValueError(f"weight-bound must fit in {WEIGHT_BITS} bits")
    if bound_mode not in {"le", "eq"}:
        raise ValueError("bound_mode must be 'le' or 'eq'")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    possible_counts = sorted({v for row in ddt for v in row if v != 0})

    with output_path.open("w", encoding="ascii", newline="\n") as fp:
        # DDT count lookup: key = Delta_in @ Delta_out (8 bits).
        fp.write("diff_dist: ARRAY BITVECTOR(8) OF BITVECTOR(8) ;\n\n")

        # Round differential states. X_0 is the plaintext difference after the
        # initial key XOR (same as plaintext difference under same-key DC).
        for r in range(rounds + 1):
            fp.write(f"X_{r}: BITVECTOR(64) ;\n")
        fp.write("\n")

        # Per-round S-box difference variables and weights.
        for r in range(rounds):
            fp.write(f"SOUT_{r}: BITVECTOR(64) ;\n")
            for i in range(NUM_SBOXES):
                fp.write(f"DIN_{r}_{i}: BITVECTOR(4) ;\n")
                fp.write(f"DOUT_{r}_{i}: BITVECTOR(4) ;\n")
                fp.write(f"W_{r}_{i}: BITVECTOR({WEIGHT_BITS}) ;\n")
            fp.write("\n")

        fp.write(f"total_weight: BITVECTOR({WEIGHT_BITS}) ;\n\n")

        # Store the exact DDT counts published S-box -> automatically computed.
        for din in range(16):
            for dout in range(16):
                key = (din << 4) | dout
                fp.write(
                    f"ASSERT( diff_dist[{cvc_hex(key, 8)}] = "
                    f"{cvc_hex(ddt[din][dout], 8)} ) ;\n"
                )
        fp.write("\n")

        # Non-trivial characteristic: plaintext difference must be non-zero,
        # unless the user fixes it explicitly (the fixed value is still checked).
        if input_diff is None:
            fp.write(
                f"ASSERT NOT( X_0 = {cvc_hex(0, BLOCK_BITS)} ) ;\n"
            )
        else:
            if input_diff == 0:
                raise ValueError("input-diff=0 gives the trivial characteristic")
            fp.write(f"ASSERT( X_0 = {cvc_hex(input_diff, BLOCK_BITS)} ) ;\n")
        fp.write("\n")

        # Round model: S-box layer -> P-layer. Round-key XOR is omitted because
        # it has zero difference in a same-key differential characteristic.
        for r in range(rounds):
            for i in range(NUM_SBOXES):
                lo = SBOX_BITS * i
                hi = lo + SBOX_BITS - 1
                key_expr = f"DIN_{r}_{i}@DOUT_{r}_{i}"

                fp.write(
                    f"ASSERT( DIN_{r}_{i} = X_{r}[{hi}:{lo}] ) ;\n"
                )

                # DDT=0 is not merely penalized: it is forbidden.
                fp.write(
                    f"ASSERT NOT( diff_dist[{key_expr}] = {cvc_hex(0, 8)} ) ;\n"
                )

                # Every possible DDT count is mapped to its exact integer weight.
                for count in possible_counts:
                    weight = ddt_weight(count)
                    fp.write(
                        f"ASSERT( ( diff_dist[{key_expr}] = {cvc_hex(count, 8)} ) "
                        f"=> ( W_{r}_{i} = {cvc_hex(weight, WEIGHT_BITS)} ) ) ;\n"
                    )

            # Reassemble 16 output nibbles into a 64-bit difference state.
            fp.write(
                f"ASSERT( SOUT_{r} = "
                + concat_nibbles(f"DOUT_{r}_{i}" for i in range(NUM_SBOXES - 1, -1, -1))
                + " ) ;\n"
            )

            # FeatherCrypt P-layer: input bit i -> output bit (r*i) mod 64.
            fp.write(
                f"ASSERT( X_{r + 1} = "
                f"{permutation_expression(f'SOUT_{r}', bit_numbering, perm_multiplier)} ) ;\n\n"
            )

        if output_diff is not None:
            fp.write(
                f"ASSERT( X_{rounds} = {cvc_hex(output_diff, BLOCK_BITS)} ) ;\n\n"
            )

        weights = [
            f"W_{r}_{i}" for r in range(rounds) for i in range(NUM_SBOXES)
        ]
        fp.write(
            f"ASSERT( total_weight = BVPLUS({WEIGHT_BITS}, "
            + ", ".join(weights)
            + " ) ) ;\n"
        )

        bound_const = cvc_hex(weight_bound, WEIGHT_BITS)
        if bound_mode == "le":
            fp.write(f"ASSERT( BVLE(total_weight, {bound_const}) ) ;\n")
        else:
            fp.write(f"ASSERT( total_weight = {bound_const} ) ;\n")

        # With X_0 != 0 and a bijective S-box, total_weight cannot be 0, but this
        # additional guard makes accidental trivial models obvious.
        fp.write(
            f"ASSERT( BVGT(total_weight, {cvc_hex(0, WEIGHT_BITS)}) ) ;\n\n"
        )

        fp.write("QUERY(FALSE) ;\n")
        fp.write("COUNTEREXAMPLE ;\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate an STP/CVC model for FeatherCrypt differential characteristics."
    )
    parser.add_argument(
        "-r", "--rounds", type=int, default=1,
        help=f"number of reduced rounds to model (1..{FULL_ROUNDS}, default: 1)",
    )
    parser.add_argument(
        "-w", "--weight-bound", type=int, default=64,
        help="target total differential weight (default: 64)",
    )
    parser.add_argument(
        "--bound-mode", choices=("le", "eq"), default="le",
        help="'le': total_weight <= bound; 'eq': total_weight == bound (default: le)",
    )
    parser.add_argument(
        "--bit-numbering", choices=("lsb0", "msb0"), default="lsb0",
        help=(
            "interpretation of paper bit index 0. The paper does not state this "
            "explicitly; default is lsb0 to match STP's native indexing."
        ),
    )
    parser.add_argument(
        "--perm-multiplier", type=int, default=DEFAULT_PERM_MULTIPLIER,
        help=(
            "permutation multiplier r in pi(i)=(r*i) mod 64. The paper gives "
            "r=17 as an example; default: 17"
        ),
    )
    parser.add_argument(
        "--input-diff", type=parse_diff, default=None,
        help="optional fixed 64-bit input difference, e.g. 0x0000000000000001",
    )
    parser.add_argument(
        "--output-diff", type=parse_diff, default=None,
        help="optional fixed 64-bit output difference",
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=None,
        help="output CVC path; default: FeatherCrypt_<r>r_w<bound>_<numbering>.cvc",
    )
    parser.add_argument(
        "--show-ddt", action="store_true", help="print the generated S-box DDT"
    )
    parser.add_argument(
        "--show-permutation", action="store_true",
        help="print paper and STP permutation mappings",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    ddt = compute_ddt(SBOX)
    validate_spec(ddt, args.perm_multiplier)

    if args.show_ddt:
        print_ddt(ddt)
        print()
    if args.show_permutation:
        print_permutation(args.bit_numbering, args.perm_multiplier)
        print()

    differential_uniformity = max(max(row) for row in ddt[1:])
    max_sbox_dp = differential_uniformity / 16.0
    max_sbox_weight = -math.log2(max_sbox_dp)

    output = args.output
    if output is None:
        output = Path(
            f"FeatherCrypt_{args.rounds}r_w{args.weight_bound}_{args.bit_numbering}.cvc"
        )

    generate_cvc(
        output_path=output,
        rounds=args.rounds,
        weight_bound=args.weight_bound,
        bound_mode=args.bound_mode,
        bit_numbering=args.bit_numbering,
        perm_multiplier=args.perm_multiplier,
        input_diff=args.input_diff,
        output_diff=args.output_diff,
        ddt=ddt,
    )

    print("FeatherCrypt STP/CVC model generated successfully.")
    print(f"  output            : {output}")
    print(f"  rounds            : {args.rounds}")
    print(f"  bound             : {args.bound_mode} {args.weight_bound}")
    print(f"  bit numbering     : {args.bit_numbering}")
    print(f"  perm multiplier   : {args.perm_multiplier}")
    print(f"  S-box DU          : {differential_uniformity}")
    print(f"  max single S-box DP: {max_sbox_dp:g} = 2^-{max_sbox_weight:g}")
    print("  DDT counts        : 0, 2, 4, 16")
    if args.bit_numbering == "lsb0":
        print(
            "  note              : the paper does not explicitly define whether bit 0 "
            "is LSB or MSB; use --bit-numbering msb0 to test the opposite convention."
        )
    if args.perm_multiplier == DEFAULT_PERM_MULTIPLIER:
        print(
            "  note              : r=17 is the only concrete permutation multiplier "
            "given in the paper, but it is written as an example (e.g., r=17)."
        )


if __name__ == "__main__":
    main()
