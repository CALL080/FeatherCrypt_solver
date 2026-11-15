# Code for "Four 16-Bit Ciphers in One: Full-Round Structural Attacks on the Lightweight Block Cipher FeatherCrypt"

This repository reproduces every number in the paper.
FeatherCrypt is the 64-bit block cipher of S. Krishnapriya (IJETT 74(7), 2026) with the bit permutation
`pi(i) = r * i mod 64`; the paper analyses `r = 17`, the only value given in the specification.

## Requirements

- Python 3.10 or later. `attack_demo.py` needs NumPy; the other scripts use only the standard library.
- For the directory `stp/`: the solver STP (we used version 2.3.3) in the `PATH`.

The round keys in all experiments are independent uniformly random 64-bit words.
This covers every key schedule, because the attacks do not use the key schedule.

## Claims of the paper and how to check them

| Paper | Claim | Command | Expected output (file in `outputs/`) |
|---|---|---|---|
| Theorem 1, Sect. 3.1 | The cipher is four independent 16-bit ciphers (2000 random plaintexts and keys) | `python verify_claims.py` | `== 1`: `ok: E_K = F0 x F1 x F2 x F3 holds in all trials` |
| Sect. 3.2, Table 4 | One flipped plaintext bit changes 8.00 ciphertext bits on average (12.5%), at most 16 | `python attack_demo.py` | `== D` |
| Sect. 4.4 | Distinguisher: 400000 of 400000 pairs | `python attack_demo.py` | `== A` |
| Sect. 4.4 | Tables from 2^16 chosen plaintexts; 10^6 of 10^6 ciphertexts decrypted | `python attack_demo.py` | `== B` (the times depend on the machine) |
| Sect. 4.4 | Known plaintexts: mean 2^19.70, max 2^20.21, 99.5% within 2^20.1 (200 trials) | `python attack_demo.py` | `== C` |
| Sect. 4.4 | The split also holds if bit 0 is the most significant bit | `python attack_demo.py` | `== E` |
| Sect. 5.1 | One-round linear approximation with correlation 1/2 | `python explore_multipliers.py` | first line (`E4`) |
| Sect. 5.2, Table 5 | Minimum weights 2, 4, 8, ..., 53 for 1 to 15 rounds, and the numbers of optimal characteristics 384, ..., 352 | `python verify_claims.py` | `== 5` (exhaustive dynamic program, no solver) |
| Sect. 5.2, Table 5 | The same minimum weights with STP | `cd stp && ./run_stp.sh 2 15` | `stp/results_r17.txt` (weight W: characteristic found; weight W-1: none) |
| Sect. 5.3 | The rotation by four bits commutes with the permutation | `python verify_claims.py` | `== 2` |
| Sect. 5.3 | No right pair for the characteristic of Fig. 4 in 60 keys; 0.47 pairs per key with its output difference | `python verify_claims.py` | `== 4` |
| Sect. 5.3 | Only 8 of the 56 "zero" bits filter: pass rate 2^-8.01 | `python verify_claims.py` | `== 3` |
| Sect. 6, Table 6 | Components and avalanche effect for all 32 multipliers | `python explore_multipliers.py` | table `E1-E3` (columns `comps`, `aval16%`) |
| Sect. 6 | 16-round characteristics of weight 39 (r=43), 44 (r=3), 46 (r=31, 63) | `python trail_other_multipliers.py` | the characteristics, re-checked with the DDT |
| Sect. 6 | STP confirms that characteristics of these weights exist | see `stp/results_other_multipliers.txt` | |

Run times on one PC: `verify_claims.py` about 15 s, `attack_demo.py` about 2 min,
`explore_multipliers.py` about 20 s, `trail_other_multipliers.py` about 1 min.
`stp/run_stp.sh 2 9` takes about 7 min; the higher rounds take longer (results of all rounds: `stp/results_r17.txt`).
The random seeds are fixed, so the outputs are the same on every run (except the times).

## Files

- `verify_claims.py` — the split, the rotation symmetry, the pre-filter, right pairs, minimum weights by dynamic programming.
- `attack_demo.py` — the three attacks of Table 1 on the full 16 rounds, the avalanche effect, the bit numbering.
- `explore_multipliers.py` — all 32 admissible multipliers.
- `trail_other_multipliers.py` — 16-round characteristics for multipliers that do not split the cipher.
- `stp/FeatherCrypt_stp.py` — generator of the STP model (CVC format) for differential characteristics.
- `stp/run_stp.sh` — generates and solves the models for weight W and W-1 for the given rounds.
- `outputs/` — the outputs of the commands above on our machine.
