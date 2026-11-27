// goldbach_v6_4_gap.cpp -- the twin/cousin-starter engine of goldbach-starters/
// (goldbach_v6_4_tc.cpp) plus review fixes (32)-(36) and the single-gap build
// (37)-(38). Default build = TC anchors; -DGAP=g = g-starter anchors (6 | g).
// NOT the campaign source of the arXiv paper: that file lives, unchanged, in
// goldbach-starters/ and must stay byte-identical there.
//
// goldbach_v6_4_tc.cpp -- TC-ANCHOR VARIANT of goldbach_v6_4 (minp build).
// [SX] Build with -DSEXY_ANCHORS for the SISTER conjecture: anchors are the
// sexy starters only (p prime with p + 6 prime). Because 3 | 6, a sexy
// starter can sit in either nonzero class mod 3, so there is no mod-3
// gearbox: one gear should serve every class of N. Under -DSEXY_ANCHORS the
// version string, artifact names, starter test and MINP columns all switch
// (minp type 'S'; column p_twin holds the minimal sexy partner, p_cousin 0).
// ONE SEMANTIC CHANGE: the anchor table holds ONLY twin/cousin starter primes --
// odd primes p with p + 2 prime (twin starter) or p + 4 prime (cousin starter).
// The claim verified is therefore STRONGER than Goldbach:
//   every even N in [start, end] admits N = p + q with q prime and p a
//   twin/cousin starter, p <= p_anchor_limit.
// Everything downstream inherits the restriction from base_primes: the
// anchor_lookup table, the odd-indexed anchor bitset (abits) used by the
// batched QHot ring, the cold-path ascending scan, and the MINP log (whose p
// is now the minimal twin/cousin anchor partner of its N). For consistency,
// escalated_search chases misses over twin/cousin starters only, and
// validate_misses re-checks the starter property on restored
// ANCHOR_EXHAUSTED rows. Artifact file names and PROG_VERSION are distinct
// from the unrestricted build so campaigns can never mix checkpoints.
// Mechanism note (why both +2 and +4 are needed): twin starters have
// p == 2 (mod 3), so for N == 2 (mod 3) the partner q = N - p would be
// divisible by 3; cousin starters have p == 1 (mod 3) and serve exactly that
// class. {+2, +4} is a complete system mod 3 (p = 3 serves as the edge case
// in both). Expect ~zero misses at modest p_anchor_limit if the conjecture
// holds. Empirically the minimal starter is small but NOT bounded by 10^4:
// over sampled windows its maximum was 6,197 up to 10^13 and 6,299 near
// 4*10^18, while exhaustive records give at least 20,747 below 10^13 (at
// N = 6,527,240,154,856) and at least 21,599 below 10^14 (manuscript, Sec. 5).
// Choose p_anchor_limit well above these (campaigns used 10^7 or 2*10^7).
//
// Review fixes (tc-r1):
//  (32) [SX] N = 6 has no sexy-starter partition (6 - 5 = 1; 3 is not a sexy
//       starter because 9 is composite), so every SX run starting at 6 ended
//       with MISS(UNRESOLVED) N=6 and exit code 4. SX now starts at N = 8 with
//       an explicit NOTE, the same way N = 4 is excluded for all builds.
//  (33) RESUME CHECKS THE PRODUCING BUILD. --resume previously accepted any
//       checkpoint that merely carried a cfg_version key. It now requires
//       cfg_version == PROG_VERSION, so a checkpoint from the unrestricted,
//       TC or SX build can never be resumed by a different build (e.g. after
//       a copy or rename). Restricted builds never wrote version-less
//       checkpoints, so a checkpoint without cfg_version is rejected too.
//  (34) Artifact wording follows the anchor mode: the scan-order provenance
//       line lists the actual first anchors (3, 5, 7, 11, ... for TC; 5, 7,
//       11, 13, ... for SX), and miss texts say "starter partition".
//  (35) -Wunused-but-set-variable on `ty` in the SX build silenced.
//  (36) -DN3FIRST (TC only, experimental, off by default): for N == 2 or 4
//       (mod 6), test N = 3 + (N - 3) before the ring or any other probe
//       (Lemmas 1-2: 3 is the only partner the mismatched gear can use in
//       those classes). Hits are certified directly and excluded from the
//       QHot counters; their Miller-Rabin calls are included in mr_calls.
//       Distinct PROG_VERSION, so (33) keeps its checkpoints separate.
//       The claim verified and the miss records are unchanged.
//  (37) -DGAP=g (g a positive multiple of 6): single-gap anchors, the
//       g-starters (p prime with p + g prime). -DSEXY_ANCHORS is the historical
//       name of GAP=6 and keeps its version string and artifact names. Each
//       gap gets its own version string and artifact files, so (33) never lets
//       one gap resume another's checkpoint. N = 6 is skipped (3 is never a
//       g-starter when 6 | g); MINP type column is G.
//  (38) Minimal-partner window statistics: count, max minimal partner (and its
//       smallest N), mean, and N without a partner, merged serially in block
//       order. Printed in the summary and, for a run that finished its range,
//       appended as one row to <prefix>_minpstats.csv. With --minp-limit=0 no
//       per-N rows are built, so a window of 10^8 N fits in memory.
//
// ---- Original v6.4 header follows ----
// goldbach_v6_3_test.cpp -- EXPERIMENT for v6.3: -DCARRYIDX (carried index and shifts in the
// AVX-512 and AVX2 prefixes) and -DCOLDLEAN (no ring lookup, no repeated trial division in the cold path).
// Based on goldbach_v6_2_idx.cpp: -DCARRYIDX also carries the gather index and shifts
// of the AVX-512 prefix between batches (no q load or index arithmetic on the carried path).
// goldbach_v6_2.cpp -- v6.2: v6.1 plus words carried between consecutive batches in the
// prefix. On by default for AVX-512 and AVX2 (-DNO_CARRY disables); available for NEON and
// the scalar prefix with -DCARRY. Tagged [CARRY].
// Carrying halves the gathers of the prefix: when a batch advances by 128, each prefix
// entry's low word is the previous batch's high word, provided the ring is unchanged.
// Based on goldbach_v6_1.cpp -- v6.1: v6 plus the fully-covered-batch shortcut (always on, except
// in -DDUMP_ALL builds). KFIX defaults to 56 with -DSIMD_AVX512 or -DSIMD_NEON (v6.2), else 32.
// Vector prefix: -DSIMD_AVX512 | -DSIMD_AVX2 | -DSIMD_NEON | none (scalar).
//
// goldbach_v6_simd.cpp -- EXPERIMENT: KFIX prefix (default 32) with optional
// -DKACC4 (four independent accumulators, portable) or -DSIMD_NEON (ARM NEON, two
// entries per vector). KFIX must be a multiple of 4 for those variants.
// Based on goldbach_v5_kfix.cpp: sorted ring + fixed-length unrolled prefix of KFIX
// probes (compile with -DKFIX=n; 0 = unchanged walk). Tagged [KFIX].
// goldbach_v5_sorted.cpp -- EXPERIMENT: v5 with the QHot ring kept sorted by q
// (descending), evicting the smallest q, and a ring walk that stops at the
// first stale entry. Changes tagged [SORT]. Test build, not for campaigns.
//
// goldbach_v5_batch.cpp -- PROTOTYPE of batched QHot. Fork of goldbach_v4.cpp.
//
// Identical to v4 except the per-N loop inside each block (tagged [B]):
//   * the QHot ring is evaluated 64 consecutive even N at a time, using one
//     64-bit slice of an odd-indexed anchor-prime bitset per ring entry;
//   * N not covered by the ring go through v4's cold path unchanged (ordered
//     anchor scan, trial-division prefilter, deterministic Miller-Rabin), and a
//     newly found q is also applied to the later N of the same batch;
//   * misses and escalation are unchanged.
// A hot certificate is (p, q) with p an anchor prime (bitset built from the
// same sieve as v4's anchor table) and q a ring entry, which by the v4 ring
// invariant was confirmed prime. Checkpoint and miss files use distinct names.
// Status: prototype. Not reviewed. Do not use for campaigns.
//
// ---- Original v4 header follows ----
// goldbach_v4.cpp — Goldbach verification engine
//
// Verifies that every even N in [start, end] admits a representation
// N = p + q with p, q odd primes and p <= p_anchor_limit.
//
// Build:  g++ -O3 -march=native -fopenmp -std=c++17 goldbach_v4.cpp -o goldbach_v4
//
// Changes from goldbach_deadclean.cpp:
//   (1) all dead/experimental code removed
//   (2) no exception can escape the OpenMP structured block
//   (3) trial-division prefilter renamed to reflect what it actually tests
//   (4) separate cost counters: QHot probes, cold candidates, prefilter
//       survivors, Miller-Rabin calls
//   (5) --resume accepted in any argv position
//   (6) misses are logged with N and escalated past the anchor limit, so
//       "anchor exhausted" is distinguished from "no partition found"
//   (+) portable includes; single contains_q call on the accept path;
//       QHot-path witnesses are now sampled alongside cold-path witnesses
//   (+) witness source tag SIEVE renamed to COLD -- no build of this
//       program has used a segmented sieve; the cold path is Miller-Rabin.
//       SIEVE is still accepted on resume as a legacy alias.
//   (7) REPRODUCIBLE SAMPLING. The witness sample is now a pure function of
//       (start, end, block_bits, p_anchor_limit, qhot_size, sample_limit).
//       It does not depend on thread count, on OpenMP scheduling, or on
//       where a run was interrupted and resumed. Four changes make this hold:
//         a. the QHot ring is cleared at the start of every block, so a
//            block's accepted witnesses depend only on that block;
//         b. sample candidates are stored per block index, not per thread,
//            so every block contributes regardless of who ran it;
//         c. bucket ownership and the sample_limit cutoff are applied in a
//            serial pass over blocks in ascending N;
//         d. misses are reported and written in the same serial order.
//       Note: elapsed_sec and the throughput figure are wall-clock and will
//       still differ between runs. The N,p,q,source rows will not.
//
// Changes in v4 (this file):
//   (8) DURABLE CHECKPOINTING. save_checkpoint verifies the stream after
//       every write and after close, fsyncs the tmp file and the containing
//       directory, and only then renames. A failed write (e.g. disk full)
//       leaves the previous checkpoint untouched instead of clobbering it
//       with a truncated file. rename() itself is checked. write_miss_file
//       gets the same treatment. Both return success; callers warn loudly
//       on failure but keep computing (the old checkpoint stays valid).
//   (9) RESUME VALIDATES CONFIG. The checkpoint records block_bits,
//       p_anchor_limit, sample_limit and qhot_size; --resume now refuses to
//       run if the command line disagrees with any of them, because the (7)
//       reproducibility guarantee is a function of those parameters and a
//       mismatched resume would silently produce a mixed artifact.
//       (threads may differ: it does not affect the sample by design.)
//       Checkpoints from older builds without cfg keys resume with a NOTE.
//  (10) The "(2) no exception can escape the OpenMP structured block" claim
//       is now enforced rather than asserted: every allocation site inside
//       the parallel region (QHot::init, per-block miss push_back) is wrapped
//       in try/catch that sets a fatal flag and requests shutdown. All
//       threads still reach the omp-for barrier, so no deadlock. A fatal
//       chunk is discarded exactly like an interrupted one; exit code 5.
//  (11) Argument parsing rejects negative numbers ("-4" no longer wraps to
//       ~2^64 via stoull). Ranges starting below 6 print an explicit NOTE
//       that N < 6 is excluded (p = 2 is excluded by design, so 4 = 2 + 2
//       is out of scope). Ranges within one block-span of 2^64 are rejected
//       up front instead of silently skipping the wrapped final block.
//  (12) MISS RECORDS SURVIVE RESUME. --resume reloads goldbach_misses.csv
//       into memory and cross-checks the row count against the checkpoint's
//       total_misses. Previously the first post-resume checkpoint rewrote
//       the miss file from an empty vector, destroying earlier miss details
//       while the counters still claimed them.
//  (13) A fresh (non-resume) run removes a stale miss file from a previous
//       run, so an empty miss set can never be shadowed by old rows.
//  (14) Checkpoint loading validates required fields (original_start/end,
//       last_verified) and invariants (start <= last_verified <= end,
//       verified + misses == total, unresolved <= misses) and states the
//       reason for rejection, instead of accepting any file where
//       last_verified > 0.
//  (15) Hardening: u64 sieve loop counter (p_anchor_limit = 2^32-1 no longer
//       wraps into an infinite loop); at most 7 positional args; per-chunk
//       block count capped at 2^24 (block_bits too small for the chunk size
//       is a clean error, not a multi-GB bad_alloc); out-of-memory in setup
//       or per-chunk bookkeeping is caught and reported, never uncaught.
//  (16) Diagnostics track "a valid checkpoint exists" (loaded or saved)
//       rather than "this process saved one", so a resumed run interrupted
//       before its first new save no longer claims no checkpoint exists.
//  (17) MISSES LIVE IN THE UNIFIED CHECKPOINT. Miss records are written as
//       # miss=N,p,q,status rows in the checkpoint file, committed by the
//       same atomic rename as the counters -- the crash window between the
//       checkpoint write and a separate miss-file write is gone. --resume
//       restores misses from the checkpoint; the counter-vs-records match
//       is enforced, and a shortfall blocks the resume unless
//       --allow-miss-loss is given. Transitional checkpoints (cfg_version
//       present, no miss rows) fall back to the separate report.
//  (18) Current-format checkpoints (cfg_version present) are validated for
//       completeness: all 18 required fields must be present exactly once
//       and parse cleanly, so a damaged file cannot slip through on
//       zero-defaulted counters. The counter invariant uses the
//       overflow-safe subtractive form (near-2^64 corrupted values cannot
//       wrap an addition and pass).
//  (19) goldbach_misses.csv is now a derived report, rewritten atomically at
//       every checkpoint (header-only when there are no misses) and stamped
//       with version, range and expected count for reconciliation. The
//       fresh-run early deletion is gone: a previous run's report survives
//       argument typos and setup failures, and is replaced only when this
//       run first checkpoints.
//  (20) Sampling-bucket allocation and resume artifact loading are under
//       bad_alloc handling like the other large allocations.
//  (21) Resuming an already-complete checkpoint returns the checkpointed
//       outcome (4 for unresolved, 2 for anchor-exhausted misses, 0 clean)
//       instead of unconditional 0, so automation sees the same status
//       from a fresh run and from a resume of its checkpoint.
//  (22) Restored miss records are validated semantically, not just counted:
//       N even, in range, strictly increasing; UNRESOLVED rows carry zero
//       p/q and tally with the unresolved counter; ANCHOR_EXHAUSTED rows
//       satisfy p + q = N (overflow-safe), p > p_anchor_limit, and both
//       witnesses re-pass is_prime64. Corruption is never overridable.
//  (23) Field completeness uses a per-field presence bitmask: missing and
//       duplicated keys are detected individually and named, where a bare
//       count of 18 could not tell one-removed-plus-one-duplicated from
//       intact.
//  (24) The derived miss report is written only after a successful
//       checkpoint save and stamped with that checkpoint's last_verified
//       and miss count; the final rewrite from in-memory state is gone.
//       The report can therefore never describe progress the authoritative
//       checkpoint does not hold.
//  (25) EXPLICIT FORMAT MARKER. Checkpoints declare cfg_format=2 (embedded
//       misses). An embedded-format checkpoint with counters claiming
//       misses but no # miss= rows, or a row/counter mismatch, is rejected
//       as corruption -- it can no longer masquerade as a transitional
//       checkpoint and silently fall back to the separate file. Genuine
//       transitional files (cfg_version, no cfg_format) keep the fallback.
//       Formats newer than this build are rejected.
//  (26) Artifact fields parse with full consumption ("49998junk" is a
//       parse failure, not 49998) and elapsed_sec must be finite and
//       non-negative.
//  (27) The miss report records checkpoint_total_misses (authoritative
//       counter) and available_miss_records (rows on hand) separately, so
//       loss accepted via --allow-miss-loss is visible in the artifact.
//  (28) The final summary names the miss report only if this process wrote
//       it successfully; otherwise it points at the checkpoint as the
//       authoritative record.
//  (29) strict_u64 is digits-only: stoull accepts a minus sign with full
//       consumption, so "-1" parsed as 2^64-1 and passed the trailing-
//       character check. strict_int had the same hole via stoi. Both now
//       reject any non-digit, the CLI parsers are wrappers over the same
//       implementation, and a PRESENT cfg_format key must equal exactly the
//       supported format (absent = transitional; explicit 1 or -1 is
//       rejected, not tolerated).
//  (30) TRUTHFUL SAVE DIAGNOSTICS. save_checkpoint and write_miss_file
//       return a three-state SaveResult: FAILED_BEFORE_RENAME (previous
//       artifact untouched and authoritative), RENAMED_NOT_DURABLE (new
//       artifact already visible, directory sync failed, crash durability
//       unconfirmed), SUCCESS. The caller reports each case accurately --
//       "previous checkpoint kept" is only claimed when it is true. Resume
//       checkpoint loading is also guarded against bad_alloc/length_error.
//  (31) MINP WINDOWS (v6.4). --minp-per-window=K logs, for the first K even
//       N of every SIEVE_SAMPLE_INTERVAL stretch of the range, the MINIMAL
//       anchor partner: the smallest anchor prime p (ascending scan over
//       base_primes; p = 2 excluded by design) with N - p prime per
//       is_prime64. The scan is independent of the QHot ring and of build
//       flags, so MINP rows are a pure function of (original_start, N,
//       p_anchor_limit): unconditioned minimal partners for EVERY N in the
//       window, not just for ring misses. Rows go to MINP_FILE, a derived
//       artifact rewritten atomically after every successful checkpoint
//       (same tmp+fsync+rename pattern as the miss report), capped by
//       --minp-limit (default 1000000). The checkpoint records
//       cfg_minp_per_window / cfg_minp_limit; --resume validates them when
//       present (older checkpoints lack them: NOTE, not error). A crash
//       between the checkpoint rename and the minp write loses at most the
//       last chunk's minp rows; checkpoint_last_verified in the minp header
//       exposes any such gap. The witness sample of (7) is unchanged.

#include <algorithm>
#if defined(SIMD_NEON)
#include <arm_neon.h>
#endif
#if defined(SIMD_AVX512) || defined(SIMD_AVX2)
#include <immintrin.h>
#endif
#ifdef SIMD_CHECK
#include <atomic>
static std::atomic<unsigned long long> simd_checks{0};
#endif
#if !defined(DUMP_ALL) && !defined(KACC4) && !defined(NO_CARRY) && !defined(CARRY) && \
    (defined(SIMD_AVX512) || defined(SIMD_AVX2))
#define CARRY 1   // [CARRY] on by default for AVX-512 and AVX2 (-DNO_CARRY disables).
                  // NEON and scalar: off by default (-DCARRY enables); measured neutral
                  // on Apple M4 at K=56 and slower for the scalar prefix on x86.
#endif
#ifndef KFIX
#if defined(SIMD_AVX512) || defined(SIMD_NEON) || defined(V61_PORTABLE)
#define KFIX 56   // measured best: AVX-512 with carried words (EPYC 9J14), NEON (Apple M4)
#else
#define KFIX 32   // scalar and AVX2 (AVX2 not yet re-swept after carried words)
#endif
#endif
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif
#include <csignal>

// (8) POSIX fsync for durable checkpointing. On non-POSIX builds the sync
// degrades to a no-op; atomicity of rename is then platform-dependent.
#if defined(__unix__) || defined(__APPLE__)
#include <fcntl.h>
#include <unistd.h>
static bool fsync_path(const char* path) {
    int fd = open(path, O_RDONLY);
    if (fd < 0) return false;
    bool ok = (fsync(fd) == 0);
    close(fd);
    return ok;
}
#else
static bool fsync_path(const char*) { return true; }
#endif

using namespace std;

// (37) [GAP] Single-gap anchor mode. -DGAP=g (g a positive multiple of 6)
// makes the anchors the g-starters: primes p with p + g prime. -DSEXY_ANCHORS
// is kept as the historical name of GAP=6 (same version string and artifact
// names as before, so existing SX checkpoints still resume). SINGLE_GAP marks
// every code path shared by the two.
#if defined(SEXY_ANCHORS) && defined(GAP) && (GAP != 6)
#error "SEXY_ANCHORS means GAP=6; do not combine it with another -DGAP"
#endif
#if defined(SEXY_ANCHORS) && !defined(GAP)
#define GAP 6
#endif
#ifdef GAP
#define SINGLE_GAP 1
static_assert(GAP >= 6 && GAP % 6 == 0,
              "-DGAP must be a positive multiple of 6: for 6 \xe2\x88\xa4 g a whole residue class "
              "mod 3 has no g-starter partition, so every N in it would be escalated as a miss");
static_assert((unsigned long long)GAP < 1000000000ULL, "-DGAP too large for the anchor sieve");
#define GB_STR2(x) #x
#define GB_STR(x) GB_STR2(x)
#endif
#if defined(N3FIRST) && defined(SINGLE_GAP)
#error "N3FIRST needs p = 3 as an anchor; 3 is never a g-starter when 6 | g (3 + g is divisible by 3)"
#endif
#if defined(N3FIRST) && defined(DUMP_ALL)
#error "N3FIRST is not wired into the DUMP_ALL certificate dump"
#endif
#ifdef SEXY_ANCHORS
static const char* PROG_VERSION = "goldbach_v6_4_sx_anchor";
#elif defined(SINGLE_GAP)
static const char* PROG_VERSION = "goldbach_v6_4_gap" GB_STR(GAP) "_anchor";   // (37)
#elif defined(N3FIRST)
// (36) Distinct version: N3FIRST changes which N reach the ring and the cold
// path, hence the witness sample, so (33) must keep its checkpoints apart.
static const char* PROG_VERSION = "goldbach_v6_4_tc_anchor_n3first";
#else
static const char* PROG_VERSION = "goldbach_v6_4_tc_anchor";
#endif
#ifdef N3FIRST
// (36) Run-level counters for the lemma-first check (printed in the summary;
// not checkpointed, so they cover this process only).
static std::atomic<unsigned long long> g_n3_tests{0}, g_n3_hits{0};
#endif
// (25) Checkpoint format: 2 = miss records embedded in the checkpoint.
// Format 1 (no cfg_format key) is the transitional layout with misses in a
// separate file.
static constexpr int CKPT_FORMAT = 2;


using u64 = uint64_t;
using u32 = uint32_t;
using u8  = uint8_t;

// (26)/(29) Full-consumption, digits-only parsing for unsigned artifact
// fields. stoull("12junk") happily returns 12, and stoull("-1") accepts the
// minus sign, consumes the whole string, and wraps to 2^64-1 -- so both the
// trailing check AND the digits-only check are needed. Every numeric field
// this program writes is a non-negative decimal.
static u64 strict_u64(const string& s) {
    if (s.empty() || s.find_first_not_of("0123456789") != string::npos)
        throw invalid_argument("expected unsigned decimal integer: " + s);
    size_t used = 0;
    u64 v = stoull(s, &used);
    if (used != s.size()) throw invalid_argument("trailing characters: " + s);
    return v;
}
static int strict_int(const string& s) {
    u64 v = strict_u64(s);   // (29) same hole existed here: stoi("-1") == -1
    if (v > (u64)numeric_limits<int>::max())
        throw out_of_range("value out of int range: " + s);
    return (int)v;
}
static double strict_nonneg_double(const string& s) {
    size_t used = 0;
    double v = stod(s, &used);
    if (used != s.size() || !isfinite(v) || v < 0.0)
        throw invalid_argument("not a finite non-negative number: " + s);
    return v;
}

// ------------------ Config ------------------

#ifndef CKPT_EVENS
#define CKPT_EVENS 100000000000ULL   // even numbers per checkpoint chunk (-DCKPT_EVENS=... to change)
#endif
static constexpr u64 CKPT_INTERVAL         = CKPT_EVENS;
static constexpr u64 SIEVE_SAMPLE_INTERVAL = 10000000000ULL;  // one witness per this much N-distance
static constexpr u64 ESCALATION_LIMIT      = 100000000ULL;    // (6) how far past the anchor limit a miss is chased

// (31) Window membership for the minimal-partner log: the first
// minp_per_window even N of every SIEVE_SAMPLE_INTERVAL stretch, measured
// from the range's original_start. A pure function of N: no dependence on
// threads, scheduling, block partition, build flags or the QHot ring.
static inline bool minp_window_member(u64 N, u64 orig_start, u64 per_window) {
    if (per_window == 0 || N < orig_start) return false;
    return (N - orig_start) % SIEVE_SAMPLE_INTERVAL < 2ULL * per_window;
}

// (30) Three-state save outcome. A plain bool cannot distinguish "the old
// artifact is still authoritative" from "the new artifact is already in
// place but its survival across a power loss is unconfirmed" -- and the
// caller's diagnostics must not claim the former when the latter happened.
enum class SaveResult {
    SUCCESS,              // written, renamed, directory sync confirmed
    FAILED_BEFORE_RENAME, // previous artifact untouched and authoritative
    RENAMED_NOT_DURABLE   // new artifact visible; crash durability unconfirmed
};

#ifdef SEXY_ANCHORS
static const char* CKPT_FILE = "goldbach_sx_checkpoint.csv";   // [SX] never mixes with TC or unrestricted
static const char* MISS_FILE = "goldbach_sx_misses.csv";        // [SX]
static const char* MINP_FILE = "goldbach_sx_minp.csv";          // [SX]
#elif defined(SINGLE_GAP)
// (37) One artifact set per gap: a gap-30 run can never resume a gap-210 checkpoint.
static const char* CKPT_FILE = "goldbach_gap" GB_STR(GAP) "_checkpoint.csv";
static const char* MISS_FILE = "goldbach_gap" GB_STR(GAP) "_misses.csv";
static const char* MINP_FILE = "goldbach_gap" GB_STR(GAP) "_minp.csv";
#else
static const char* CKPT_FILE = "goldbach_tc_checkpoint.csv";   // [TC] never the unrestricted campaign file
static const char* MISS_FILE = "goldbach_tc_misses.csv";        // [TC]
static const char* MINP_FILE = "goldbach_tc_minp.csv";          // (31) derived minimal-partner log
#endif
// (38) Append-only summary of the minimal-partner window(s) of each completed
// run: one row per run, independent of --minp-limit (so a window of 10^8 N
// can be summarised without holding 10^8 rows in memory).
#ifdef SEXY_ANCHORS
static const char* MINPSTATS_FILE = "goldbach_sx_minpstats.csv";
#elif defined(SINGLE_GAP)
static const char* MINPSTATS_FILE = "goldbach_gap" GB_STR(GAP) "_minpstats.csv";
#else
static const char* MINPSTATS_FILE = "goldbach_tc_minpstats.csv";
#endif

static atomic<int> g_shutdown{0};
// (10) Set when an exception is caught inside the parallel region; the run
// aborts, discarding the in-flight chunk, instead of letting the exception
// escape the OpenMP structured block (undefined behaviour).
static atomic<int> g_fatal{0};

static void signal_handler(int) {
    g_shutdown.store(1, memory_order_relaxed);
}

// ------------------ Deterministic 64-bit Miller-Rabin ------------------

static inline u64 mod_mul64(u64 a, u64 b, u64 mod) {
    return (u64)((__uint128_t)a * b % mod);
}

static inline u64 mod_pow64(u64 a, u64 d, u64 mod) {
    u64 r = 1;
    while (d) {
        if (d & 1ULL) r = mod_mul64(r, a, mod);
        a = mod_mul64(a, a, mod);
        d >>= 1ULL;
    }
    return r;
}

// (3) Renamed. This is NOT a primality test: it returns true for any n that
// survives trial division by the primes <= 17 (e.g. 361 = 19^2 returns true).
// It is only ever used as a cheap rejection filter ahead of Miller-Rabin.
// Survival rate for random odd n is phi(510510)/510510 ~= 0.171.
static inline bool survives_small_trial_division(u64 n) {
    if (n < 2) return false;
    static const u64 small[] = {2, 3, 5, 7, 11, 13, 17};
    for (u64 p : small) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }
    return true;
}

// Deterministic for all n < 2^64 (Sinclair 7-base set, verified against the
// Feitsma base-2 pseudoprime database). Convention: bases a >= n are skipped
// rather than reduced mod n. This is safe here because every n at or below
// the largest base (1795265022) either falls to the trial division above or
// is decided by the smaller bases, a region that is exhaustively checkable;
// for all larger n every base is < n and the published determinism result
// applies directly. State this convention explicitly when citing the set.
static bool is_prime64(u64 n) {
    if (n < 2) return false;
    static const u64 small[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (u64 p : small) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }

    u64 d = n - 1;
    int s = 0;
    while ((d & 1ULL) == 0) { d >>= 1ULL; ++s; }

    static const u64 bases[] = {2ULL, 325ULL, 9375ULL, 28178ULL,
                                450775ULL, 9780504ULL, 1795265022ULL};
    for (u64 a : bases) {
        if (a >= n) continue;
        u64 x = mod_pow64(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool witness = true;
        for (int r = 1; r < s; ++r) {
            x = mod_mul64(x, x, n);
            if (x == n - 1) { witness = false; break; }
        }
        if (witness) return false;
    }
    return true;
}

#ifdef COLDLEAN
// [COLDLEAN] is_prime64 for numbers that already passed survives_small_trial_division
// (no factor 2..17): skip that part of the trial division, keep 19..37, then the same
// deterministic Miller-Rabin. Same result as is_prime64 for such n.
static bool is_prime64_after_td(u64 n) {
    static const u64 small[] = {19, 23, 29, 31, 37};
    for (u64 p : small) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }
    u64 d = n - 1;
    int s = 0;
    while ((d & 1ULL) == 0) { d >>= 1ULL; ++s; }
    static const u64 bases[] = {2ULL, 325ULL, 9375ULL, 28178ULL,
                                450775ULL, 9780504ULL, 1795265022ULL};
    for (u64 a : bases) {
        if (a >= n) continue;
        u64 x = mod_pow64(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool witness = true;
        for (int r = 1; r < s; ++r) {
            x = mod_mul64(x, x, n);
            if (x == n - 1) { witness = false; break; }
        }
        if (witness) return false;
    }
    return true;
}
#endif

// ------------------ Base primes (anchor table) ------------------
//
// [TC] The anchor table holds ONLY twin/cousin starters: odd primes p with
// p + 2 prime or p + 4 prime. anchor_lookup and abits are filled from this
// table in main(), so the QHot ring, the batch slices, the cold path and the
// MINP scan all inherit the restriction with no further change.

static vector<u32> base_primes;
// [TC-INST] starter_types[i] classifies base_primes[i]:
// bit 0 = twin starter (p + 2 prime), bit 1 = cousin starter (p + 4 prime).
// Filled by build_base_primes from the same sieve, so it is exact (no MR).
static vector<u8> starter_types;

// [TC] Starter test by Miller-Rabin, for the paths that step past the sieved
// table (escalated_search, validate_misses). p + 4 cannot overflow for any p
// this program passes here (p <= p_anchor_limit + ESCALATION_LIMIT << 2^64).
static inline bool is_tc_starter64(u64 p) {
#ifdef SINGLE_GAP
    return is_prime64(p) && is_prime64(p + (u64)GAP);   // [SX]/(37) g-starter
#else
    return is_prime64(p) && (is_prime64(p + 2) || is_prime64(p + 4));
#endif
}

static void build_base_primes(u32 limit) {
    // [TC] Sieve 4 past the limit so the starter test for p = limit itself
    // reads sieved bits, not special cases. u64 size guards limit = 2^32-1.
#ifdef SINGLE_GAP
    const u64 top = (u64)limit + (u64)GAP;   // [SX]/(37) sieve GAP past the limit
#else
    const u64 top = (u64)limit + 4ULL;
#endif
    vector<bool> is_prime((size_t)top + 1, true);
    is_prime[0] = false;
    is_prime[1] = false;

    for (u64 i = 2; i * i <= top; ++i) {
        if (is_prime[(size_t)i]) {
            for (u64 j = i * i; j <= top; j += i)
                is_prime[(size_t)j] = false;
        }
    }

    base_primes.clear();
    starter_types.clear();
    // p = 2 is deliberately excluded: N is even, so q = N - 2 is even and > 2.
    // (15) u64 counter: with a u32 counter and limit = 2^32-1 (which the
    // validation range [3, 2^32) permits), i += 2 wraps and loops forever.
    // [TC] Keep only twin/cousin starters: p + 2 or p + 4 prime.
    for (u64 i = 3; i <= limit; i += 2) {
        if (!is_prime[(size_t)i]) continue;
#ifdef SINGLE_GAP
        // [SX]/(37) bit 0 = g-starter; the MINP scan reports it in column p_twin.
        const u8 t = (u8)(is_prime[(size_t)(i + (u64)GAP)] ? 1 : 0);
#else
        const u8 t = (u8)((is_prime[(size_t)(i + 2)] ? 1 : 0) |
                          (is_prime[(size_t)(i + 4)] ? 2 : 0));
#endif
        if (t) {
            base_primes.push_back((u32)i);
            starter_types.push_back(t);   // [TC-INST]
        }
    }
}

// ------------------ QHot: per-thread ring of recently confirmed primes q ------------------
//
// Invariant: every slot in buf[0 .. count-1] holds an odd prime that was
// confirmed either by is_prime64 or by already being present in the ring.
// Therefore a hit (p from the sieved anchor table, q from the ring,
// p + q == N by construction) is a complete certificate for N.

struct QHot {
    vector<u64> buf;
    int pos   = 0;
    int count = 0;
    int size  = 1024;
    u64 ver   = 0;      // [CARRY] bumped whenever the ring's contents or order change

    // (2) Precondition: sz is a positive power of two. Validated in main().
    // No throw here: this runs inside an OpenMP structured block, and an
    // exception escaping one is undefined behaviour.
    void init(int sz) {
        size = sz;
        pos = 0;
        count = 0;
        buf.assign((size_t)size, 0);
    }

    // (7) Cheap logical clear. Entries beyond count-1 are never read
    // (contains_q scans [0, count); hit() walks pos-1 down to pos-count),
    // so no memset is needed. Called at the top of every block so that a
    // block's ring state depends only on that block, not on which blocks
    // this thread happened to process earlier.
    inline void reset() {
        pos   = 0;
        count = 0;
        ++ver;          // [CARRY]
    }

    // [SORT] buf[0..count) is kept in descending order of q (ascending d).
    // A full ring evicts its smallest q, i.e. the entry with the largest
    // partner d = N - q, which is the first to become useless.
    inline void record(u64 q) {
        if (q < 3 || ((q & 1ULL) == 0)) return;
        int n = count;
        if (n == size) {
            if (q <= buf[(size_t)(n - 1)]) return;   // stalest of all: not worth keeping
            --n;                                      // drop the current smallest
        }
        int i = n;
        while (i > 0 && buf[(size_t)(i - 1)] < q) { buf[(size_t)i] = buf[(size_t)(i - 1)]; --i; }
        buf[(size_t)i] = q;
        count = n + 1;
        ++ver;          // [CARRY]
    }

    inline bool contains_q(u64 q) const {
        for (int i = 0; i < count; ++i)
            if (buf[(size_t)i] == q) return true;
        return false;
    }

    // (4) Returns the probe count as a signed ordinal: > 0 means a hit was
    // found after that many probes, <= 0 means no hit after -r probes. Returning
    // it keeps the caller's counter in a register -- writing through a u64&
    // reference forces it to the stack and costs ~7% in the hot loop.
    inline int hit(u64 N, u64 p_limit, const vector<u8>& anchor_lookup,
                   u64& p_out, u64& q_out) const {
        if (count == 0) return 0;

        const int mask = size - 1;

#define QHOT_CHECK(IDX, ORD)                                     \
        {                                                        \
            u64 q_ = buf[(size_t)(IDX)];                         \
            if (q_ < N) {                                        \
                u64 p_ = N - q_;                                 \
                if (p_ <= p_limit && anchor_lookup[(size_t)p_]) {\
                    p_out = p_; q_out = q_;                      \
                    return (ORD);                                \
                }                                                \
            }                                                    \
        }

        (void)mask;
        for (int i = 0; i < count; ++i)                 // [SORT] largest q first
            QHOT_CHECK(i, i + 1);

#undef QHOT_CHECK
        return -count;
    }
};

// [B] 64-bit slice of the odd-indexed anchor bitset starting at bit b
// (bit k <-> p = 2k + 1). The bitset carries two zero words of padding, so
// w + 1 is always in range for b <= (p_anchor_limit - 1)/2.
static inline u64 anchor_slice(const vector<u64>& a, u64 b) {
    const size_t w = (size_t)(b >> 6);
    const unsigned s = (unsigned)(b & 63ULL);
    return s ? ((a[w] >> s) | (a[w + 1] << (64U - s))) : a[w];
}

// ------------------ (6) Miss handling ------------------

enum MissStatus {
    MISS_ANCHOR_EXHAUSTED = 0, // a partition exists, but needs p > p_anchor_limit
    MISS_UNRESOLVED       = 1  // no partition found up to the escalation limit
};

struct MissRecord {
    u64 N = 0;
    u64 p = 0;   // 0 when unresolved
    u64 q = 0;
    int status = MISS_UNRESOLVED;
};

// (7) One sample candidate slot per block, per path. Written by whichever
// thread owns the block; read back in block-index order after the parallel
// region. Block index -- not thread id, not completion order -- decides which
// witness reaches the artifact.
struct BlockWitness {
    u8  have_qhot = 0, have_cold = 0;
    u64 qN = 0, qp = 0, qq = 0;
    u64 cN = 0, cp = 0, cq = 0;
};

// Chase a miss past the anchor limit with an unbounded (well, escalation-limited)
// ordered scan. Cost is irrelevant because this path is expected never to run.
// [TC] The escalation keeps the variant's claim coherent: it searches
// twin/cousin starters only. ANCHOR_EXHAUSTED then means "a starter partition
// exists but needs p > p_anchor_limit"; UNRESOLVED means "no starter
// partition up to the escalation limit" -- the conjecture-breaking case.
static bool escalated_search(u64 N, u64 p_from, u64 p_to, u64& p_out, u64& q_out) {
    u64 p = p_from | 1ULL;
    if (p < 3) p = 3;
    for (; p <= p_to && p < N; p += 2) {
        if (!survives_small_trial_division(p)) continue;
        if (!is_tc_starter64(p)) continue;   // [TC] prime AND twin/cousin starter
        u64 q = N - p;
        if (is_prime64(q)) { p_out = p; q_out = q; return true; }
    }
    return false;
}

// (8) Returns false if the file could not be written completely. Writes to a
// tmp file and renames, so a failed write never truncates an existing miss
// file from an earlier checkpoint.
// (12) Load previously recorded misses so a resumed run carries them forward.
// Without this, the first post-resume checkpoint rewrites the miss file from
// an empty vector, destroying every miss recorded before the interruption
// while the checkpoint counters still claim them.
static vector<MissRecord> load_miss_file(u64& malformed_out) {
    vector<MissRecord> misses;
    malformed_out = 0;
    ifstream f(MISS_FILE);
    if (!f.is_open()) return misses;

    string line;
    while (getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (line.rfind("N,", 0) == 0) continue;  // header
        size_t c1 = line.find(',');
        size_t c2 = (c1 == string::npos) ? string::npos : line.find(',', c1 + 1);
        size_t c3 = (c2 == string::npos) ? string::npos : line.find(',', c2 + 1);
        if (c3 == string::npos) { ++malformed_out; continue; }
        try {
            MissRecord m;
            m.N = strict_u64(line.substr(0, c1));
            m.p = strict_u64(line.substr(c1 + 1, c2 - c1 - 1));
            m.q = strict_u64(line.substr(c2 + 1, c3 - c2 - 1));
            const string status = line.substr(c3 + 1);
            if      (status == "ANCHOR_EXHAUSTED") m.status = MISS_ANCHOR_EXHAUSTED;
            else if (status == "UNRESOLVED")       m.status = MISS_UNRESOLVED;
            else { ++malformed_out; continue; }
            misses.push_back(m);
        } catch (...) { ++malformed_out; }
    }
    return misses;
}

// (22) Semantic validation of restored miss records. A corrupted checkpoint
// can carry the RIGHT NUMBER of miss rows with wrong content; row count
// alone must not admit it. Returns an empty string when valid, else the
// reason. skip_unresolved_count is set when a shortfall was accepted via
// --allow-miss-loss (the per-status tally is then meaningless).
static string validate_misses(const vector<MissRecord>& misses,
                              u64 original_start, u64 last_verified,
                              u64 expected_unresolved, u64 p_anchor_limit,
                              bool skip_unresolved_count) {
    u64 unresolved_rows = 0;
    u64 prev_N = 0;
    for (const auto& m : misses) {
        ostringstream o;
        if ((m.N & 1ULL) || m.N < original_start || m.N > last_verified) {
            o << "miss N=" << m.N << " is odd or outside"
                 " [original_start, last_verified]";
            return o.str();
        }
        if (prev_N != 0 && m.N <= prev_N) {
            o << "miss rows unordered or duplicated at N=" << m.N;
            return o.str();
        }
        prev_N = m.N;
        if (m.status == MISS_UNRESOLVED) {
            ++unresolved_rows;
            if (m.p != 0 || m.q != 0) {
                o << "UNRESOLVED miss N=" << m.N << " carries nonzero p/q";
                return o.str();
            }
        } else {
            // Overflow-safe p + q == N, then the escalation contract
            // (p > p_anchor_limit), then primality of both witnesses.
            if (m.q > m.N || m.N - m.q != m.p) {
                o << "miss N=" << m.N << ": p + q != N";
                return o.str();
            }
            if (p_anchor_limit != 0 && m.p <= p_anchor_limit) {
                o << "ANCHOR_EXHAUSTED miss N=" << m.N
                  << " has p <= p_anchor_limit";
                return o.str();
            }
            if (!is_prime64(m.p) || !is_prime64(m.q)) {
                o << "miss N=" << m.N << ": p or q is not prime";
                return o.str();
            }
            // [TC] The escalation contract of this variant: the witness p
            // must itself be a twin/cousin starter.
            if (!is_tc_starter64(m.p)) {
                o << "miss N=" << m.N << ": p is not a twin/cousin starter";
                return o.str();
            }
        }
    }
    if (!skip_unresolved_count && unresolved_rows != expected_unresolved) {
        ostringstream o;
        o << "unresolved miss rows (" << unresolved_rows
          << ") do not match checkpoint unresolved counter ("
          << expected_unresolved << ")";
        return o.str();
    }
    return "";
}

// (8)/(19)/(24) Derived, human-readable report. Written ONLY after a
// successful checkpoint save and stamped with that checkpoint's
// last_verified and miss count, so it can never describe state the
// authoritative artifact does not hold. Header-only when there are zero
// misses; tmp+rename, so a stale file from an earlier run is replaced by
// this run's first successful checkpoint rather than deleted up front (13).
static SaveResult write_miss_file(const vector<MissRecord>& misses,
                                  u64 original_start, u64 original_end,
                                  u64 last_verified, u64 total_misses_counter) {
    string tmp = string(MISS_FILE) + ".tmp";
    ofstream f(tmp, ios::out | ios::trunc);
    if (!f.is_open()) return SaveResult::FAILED_BEFORE_RENAME;
    f << "# Goldbach verification misses (derived report; authoritative copy"
         " = # miss= rows in " << CKPT_FILE << ")\n";
    f << "# generated_by=" << PROG_VERSION << "\n";
    f << "# range_start=" << original_start << "\n";
    f << "# range_end="   << original_end   << "\n";
    f << "# checkpoint_last_verified=" << last_verified << "\n";
    // (27) Both figures on purpose: after --allow-miss-loss the available
    // records can be fewer than the authoritative counter, and the report
    // must show that shortfall rather than present itself as complete.
    f << "# checkpoint_total_misses="  << total_misses_counter << "\n";
    f << "# available_miss_records="   << misses.size() << "\n";
    // (34) Wording follows the anchor mode: p must be a starter.
    f << "# status ANCHOR_EXHAUSTED: a starter partition exists but requires p > p_anchor_limit\n";
    f << "# status UNRESOLVED: no starter partition found up to the escalation limit "
      << ESCALATION_LIMIT << " -- INVESTIGATE\n";
    f << "N,p,q,status\n";
    for (const auto& m : misses) {
        f << m.N << ","
          << (m.status == MISS_ANCHOR_EXHAUSTED ? m.p : 0) << ","
          << (m.status == MISS_ANCHOR_EXHAUSTED ? m.q : 0) << ","
          << (m.status == MISS_ANCHOR_EXHAUSTED ? "ANCHOR_EXHAUSTED" : "UNRESOLVED")
          << "\n";
    }
    f.flush();
    if (!f.good()) return SaveResult::FAILED_BEFORE_RENAME;
    f.close();
    if (!f.good()) return SaveResult::FAILED_BEFORE_RENAME;
    if (!fsync_path(tmp.c_str())) return SaveResult::FAILED_BEFORE_RENAME;
    if (rename(tmp.c_str(), MISS_FILE) != 0) return SaveResult::FAILED_BEFORE_RENAME;
    // (30) Past this point the new file is already visible; a directory-sync
    // failure downgrades the outcome, it does not undo it.
    if (!fsync_path(".")) return SaveResult::RENAMED_NOT_DURABLE;
    return SaveResult::SUCCESS;
}

// (31) Derived minimal-partner log. Same safety pattern as write_miss_file:
// tmp + fsync + rename, written only after a successful checkpoint save and
// stamped with that checkpoint's last_verified. A crash between the
// checkpoint rename and this write loses at most the last chunk's rows;
// the checkpoint_last_verified stamp exposes the gap on inspection.
static SaveResult write_minp_file(const vector<string>& rows,
                                  u64 original_start, u64 original_end,
                                  u64 last_verified,
                                  u64 minp_per_window, u64 minp_limit,
                                  u64 p_anchor_limit) {
    string tmp = string(MINP_FILE) + ".tmp";
    ofstream f(tmp, ios::out | ios::trunc);
    if (!f.is_open()) return SaveResult::FAILED_BEFORE_RENAME;
    f << "# Goldbach minimal-partner log (derived artifact; see (31))\n";
    f << "# generated_by=" << PROG_VERSION << "\n";
    f << "# range_start=" << original_start << "\n";
    f << "# range_end="   << original_end   << "\n";
    f << "# checkpoint_last_verified=" << last_verified << "\n";
    f << "# minp_per_window=" << minp_per_window << "\n";
    f << "# minp_limit="      << minp_limit      << "\n";
    f << "# p_anchor_limit="  << p_anchor_limit  << "\n";
    f << "# window: the first minp_per_window even N of every "
      << SIEVE_SAMPLE_INTERVAL << " stretch, measured from range_start\n";
    f << "# p is the MINIMAL anchor partner of its N: the smallest p among\n";
#ifdef SEXY_ANCHORS
    f << "# base_primes -- [SX] sexy starters only (p prime with p+6 prime),\n";
    f << "# ascending (5, 7, 11, 13, 17, ...; 3 is excluded since 9 is composite) with\n";
#elif defined(SINGLE_GAP)
    f << "# base_primes -- (37) g-starters only, g = " << GAP << " (p prime with p+g prime),\n";
    f << "# ascending (3 is never one, since 3 | 3+g) with\n";
#else
    f << "# base_primes -- [TC] twin/cousin starters only (p prime with p+2\n";
    f << "# or p+4 prime), ascending (3, 5, 7, 11, 13, ...; p > 2 by design) with\n";
#endif
    f << "# N - p prime per is_prime64 (deterministic 7-base Miller-Rabin,\n";
    f << "# Sinclair base set, bases a >> n skipped as documented there).\n";
    f << "# Independent of the QHot ring: rows are a pure function of\n";
    f << "# (range_start, N, p_anchor_limit). N with no anchor partner\n";
    f << "# produce no row here and appear in the miss records instead.\n";
    f << "# minp_format=2 [TC-INST]: three extra columns --\n";
    f << "#   type: T = p is a twin starter only (p+2 prime), C = cousin only\n";
    f << "#         (p+4 prime), B = both (p, p+2, p+4 all prime)\n";
    f << "#   p_twin:   minimal TWIN-starter partner of N (0 = none <= limit)\n";
    f << "#   p_cousin: minimal COUSIN-starter partner of N (0 = none <= limit)\n";
    f << "# p_twin == 0 marks an N that NEEDS the +4 gear (mod-3 blocked or\n";
    f << "# unlucky); min(p_twin!=0, p_cousin!=0) == p by construction.\n";
#ifdef SEXY_ANCHORS
    f << "# [SX] SEXY_ANCHORS build: base_primes are sexy starters (p, p+6 prime);\n";
    f << "#   type is always S, p_twin = minimal sexy partner (== p), p_cousin = 0.\n";
#elif defined(SINGLE_GAP)
    f << "# (37) GAP build, g = " << GAP << ": base_primes are g-starters;\n";
    f << "#   type is always G, p_twin = minimal g-starter partner (== p), p_cousin = 0.\n";
#endif
    f << "N,p,q,type,p_twin,p_cousin\n";
    for (const auto& r : rows) f << r << "\n";
    f.flush();
    if (!f.good()) return SaveResult::FAILED_BEFORE_RENAME;
    f.close();
    if (!f.good()) return SaveResult::FAILED_BEFORE_RENAME;
    if (!fsync_path(tmp.c_str())) return SaveResult::FAILED_BEFORE_RENAME;
    if (rename(tmp.c_str(), MINP_FILE) != 0) return SaveResult::FAILED_BEFORE_RENAME;
    if (!fsync_path(".")) return SaveResult::RENAMED_NOT_DURABLE;
    return SaveResult::SUCCESS;
}

// ------------------ Checkpoint (unified CSV) ------------------

struct CkptState {
    u64 original_start = 0;
    u64 original_end   = 0;
    u64 last_verified  = 0;
    u64 total_verified = 0;
    u64 total_total    = 0;
    u64 qhot_probes    = 0;
    u64 cold_cand      = 0;
    u64 cold_surv      = 0;
    u64 mr_calls       = 0;
    u64 total_qhot     = 0;
    u64 total_misses   = 0;
    u64 unresolved     = 0;
    double elapsed_sec = 0.0;
    bool valid = false;

    // (9) Run configuration recorded in the checkpoint. 0 = key absent
    // (checkpoint written by an older build); validation is then skipped.
    // These four parameters, plus (start, end), fully determine the witness
    // sample under guarantee (7), so --resume must not change them.
    int cfg_block_bits     = 0;
    u64 cfg_p_anchor_limit = 0;
    u64 cfg_sample_limit   = 0;   // 0 also means "absent": a sample_limit=0
                                  // run writes no samples worth protecting.
    int cfg_qhot_size      = 0;
    // (31) Minimal-partner logging config. has_minp_cfg distinguishes "keys
    // absent (older build)" from "feature off (recorded 0)".
    u64 cfg_minp_per_window = 0;
    u64 cfg_minp_limit      = 0;
    bool has_minp_cfg       = false;

    // (17) Miss records read back from # miss= rows in this file.
    vector<MissRecord> misses;
    // (18) True when the checkpoint declares cfg_version (current format);
    // such checkpoints are validated strictly for field completeness.
    bool has_version = false;
    string cfg_version_str;   // (33) value of cfg_version, checked on resume
    // (25) Checkpoint layout: 0 = key absent (transitional, misses in a
    // separate file), 2 = misses embedded in this file.
    int cfg_format = 0;
};

static CkptState load_checkpoint() {
    CkptState ck;
    ifstream f(CKPT_FILE);
    if (!f.is_open()) return ck;

    // (14)/(18)/(23) Per-field presence tracking (cfg_threads is
    // informational and not tracked). A current-format checkpoint must carry
    // each required field EXACTLY once: a bare count of 18 cannot tell
    // "one field removed + another duplicated" from an intact file.
    enum Field { F_OS, F_OE, F_LV, F_TV, F_TT, F_QP, F_CC, F_CS, F_MR,
                 F_TQ, F_TM, F_UN, F_ES, F_VER, F_FMT, F_BB, F_PAL, F_SL,
                 F_QS, F_MPW, F_MPL, F_COUNT };
    static const char* FIELD_NAMES[F_COUNT] = {
        "original_start", "original_end", "last_verified", "total_verified",
        "total_total", "qhot_probes", "cold_candidates", "cold_survivors",
        "mr_calls", "total_qhot", "total_misses", "unresolved",
        "elapsed_sec", "cfg_version", "cfg_format", "cfg_block_bits",
        "cfg_p_anchor_limit", "cfg_sample_limit", "cfg_qhot_size",
        "cfg_minp_per_window", "cfg_minp_limit" };
    u32 seen_mask = 0;
    u32 dup_mask  = 0;
    auto mark = [&](int fbit) {
        u32 b = 1u << fbit;
        if (seen_mask & b) dup_mask |= b;
        seen_mask |= b;
    };
    constexpr u32 FULL_MASK = (1u << F_COUNT) - 1u;
    u64 parse_failures = 0;

    string line;
    while (getline(f, line)) {
        if (line.empty() || line[0] != '#') continue;
        auto eq = line.find('=');
        if (eq == string::npos) continue;
        size_t ks = 1;
        while (ks < eq && line[ks] == ' ') ++ks;
        string key = line.substr(ks, eq - ks);
        string val = line.substr(eq + 1);
        try {
            if      (key == "original_start") { ck.original_start = strict_u64(val); mark(F_OS); }
            else if (key == "original_end")   { ck.original_end   = strict_u64(val); mark(F_OE); }
            else if (key == "last_verified")  { ck.last_verified  = strict_u64(val); mark(F_LV); }
            else if (key == "total_verified")  { ck.total_verified = strict_u64(val); mark(F_TV); }
            else if (key == "total_total")     { ck.total_total    = strict_u64(val); mark(F_TT); }
            else if (key == "qhot_probes")     { ck.qhot_probes    = strict_u64(val); mark(F_QP); }
            else if (key == "cold_candidates") { ck.cold_cand      = strict_u64(val); mark(F_CC); }
            else if (key == "cold_survivors")  { ck.cold_surv      = strict_u64(val); mark(F_CS); }
            else if (key == "mr_calls")        { ck.mr_calls       = strict_u64(val); mark(F_MR); }
            else if (key == "total_qhot")      { ck.total_qhot     = strict_u64(val); mark(F_TQ); }
            else if (key == "total_misses")    { ck.total_misses   = strict_u64(val); mark(F_TM); }
            else if (key == "unresolved")      { ck.unresolved     = strict_u64(val); mark(F_UN); }
            else if (key == "elapsed_sec")     { ck.elapsed_sec    = strict_nonneg_double(val); mark(F_ES); }
            else if (key == "cfg_version")     { ck.has_version    = true; ck.cfg_version_str = val; mark(F_VER); }
            else if (key == "cfg_format")      { ck.cfg_format     = strict_int(val); mark(F_FMT); }
            else if (key == "cfg_block_bits")     { ck.cfg_block_bits     = strict_int(val); mark(F_BB); }
            else if (key == "cfg_p_anchor_limit") { ck.cfg_p_anchor_limit = strict_u64(val); mark(F_PAL); }
            else if (key == "cfg_sample_limit")   { ck.cfg_sample_limit   = strict_u64(val); mark(F_SL); }
            else if (key == "cfg_qhot_size")      { ck.cfg_qhot_size      = strict_int(val); mark(F_QS); }
            else if (key == "cfg_minp_per_window") { ck.cfg_minp_per_window = strict_u64(val); ck.has_minp_cfg = true; mark(F_MPW); }
            else if (key == "cfg_minp_limit")      { ck.cfg_minp_limit      = strict_u64(val); mark(F_MPL); }
            else if (key == "miss") {
                // (17) # miss=N,p,q,status -- authoritative miss record.
                size_t c1 = val.find(',');
                size_t c2 = (c1 == string::npos) ? string::npos : val.find(',', c1 + 1);
                size_t c3 = (c2 == string::npos) ? string::npos : val.find(',', c2 + 1);
                if (c3 == string::npos) { ++parse_failures; }
                else {
                    MissRecord m;
                    m.N = strict_u64(val.substr(0, c1));
                    m.p = strict_u64(val.substr(c1 + 1, c2 - c1 - 1));
                    m.q = strict_u64(val.substr(c2 + 1, c3 - c2 - 1));
                    const string status = val.substr(c3 + 1);
                    if      (status == "ANCHOR_EXHAUSTED") m.status = MISS_ANCHOR_EXHAUSTED;
                    else if (status == "UNRESOLVED")       m.status = MISS_UNRESOLVED;
                    else { ++parse_failures; continue; }
                    ck.misses.push_back(m);
                }
            }
        } catch (...) { ++parse_failures; }
    }

    // (14) Reject with a stated reason instead of the old rule
    // "last_verified > 0 means valid".
    auto reject = [&](const string& why) {
        cerr << "ERROR: checkpoint " << CKPT_FILE << " rejected: " << why << "\n";
        ck.valid = false;
        return ck;
    };
    if (parse_failures)
        cerr << "WARNING: " << parse_failures
             << " checkpoint field(s) failed to parse.\n";
    if (!(seen_mask & (1u << F_OS)) || !(seen_mask & (1u << F_OE)) ||
        !(seen_mask & (1u << F_LV)))
        return reject("missing required field(s): original_start,"
                      " original_end and/or last_verified");
    // (18)/(23) A checkpoint that declares the current format must carry
    // every field the current save writes, each exactly once, all parsed
    // cleanly. Missing and duplicated fields are named in the rejection.
    if (ck.has_version) {
        // (25) Transitional (format-1) files legitimately lack cfg_format;
        // everything else in the current field set is still required.
        // (31) Like cfg_format, the minp keys are transitional: absent in
        // checkpoints from older builds. When either is present, both are
        // required -- one without the other is field damage, not a format.
        u32 optional_absent = 0;
        if (!(seen_mask & (1u << F_FMT))) optional_absent |= (1u << F_FMT);
        if (!(seen_mask & (1u << F_MPW)) && !(seen_mask & (1u << F_MPL)))
            optional_absent |= (1u << F_MPW) | (1u << F_MPL);
        const u32 required_mask = FULL_MASK & ~optional_absent;
        if (dup_mask || seen_mask != required_mask) {
            ostringstream o;
            o << "current-format checkpoint field damage --";
            for (int i = 0; i < F_COUNT; ++i) {
                if (!(seen_mask & (1u << i))) o << " missing:" << FIELD_NAMES[i];
                if (dup_mask & (1u << i))    o << " duplicated:" << FIELD_NAMES[i];
            }
            return reject(o.str());
        }
        if (parse_failures)
            return reject("current-format checkpoint has unparseable field(s)");
        // (29) When the marker is present it must be exactly the format this
        // build writes; the transitional layout is represented by an ABSENT
        // key, so an explicit 1 (or any other value) is not a valid state.
        if ((seen_mask & (1u << F_FMT)) && ck.cfg_format != CKPT_FORMAT) {
            ostringstream o;
            o << "unsupported checkpoint format " << ck.cfg_format
              << " (this build writes " << CKPT_FORMAT
              << "; transitional files omit the key)";
            return reject(o.str());
        }
        // (25) An embedded-miss checkpoint must actually carry its misses:
        // counters claiming misses with no # miss= rows means the rows were
        // stripped or lost -- this must NOT fall back to the separate file
        // the way a genuine transitional checkpoint does.
        if (ck.cfg_format >= 2) {
            if (ck.total_misses > 0 && ck.misses.empty())
                return reject("embedded-miss format declares total_misses > 0"
                              " but carries no # miss= rows (rows stripped or"
                              " corrupted)");
            if ((u64)ck.misses.size() != ck.total_misses)
                return reject("embedded miss row count does not match"
                              " total_misses (single-artifact mismatch ="
                              " corruption, not loss)");
        }
    }
    if (ck.last_verified == 0)
        return reject("last_verified is 0");
    if (ck.original_start > ck.last_verified || ck.last_verified > ck.original_end)
        return reject("range invariant violated: need original_start <="
                      " last_verified <= original_end");
    // (18) Overflow-safe form: with corrupted values near 2^64, the additive
    // check total_verified + total_misses == total_total can wrap and pass.
    if (ck.total_verified > ck.total_total ||
        ck.total_misses != ck.total_total - ck.total_verified)
        return reject("counter invariant violated: total_verified +"
                      " total_misses != total_total");
    if (ck.unresolved > ck.total_misses)
        return reject("counter invariant violated: unresolved > total_misses");

    ck.valid = true;
    return ck;
}

// (8) Returns false if the checkpoint could not be written completely and
// durably. On any failure the previous checkpoint file is left untouched:
// nothing is renamed until the tmp file has been written, verified and
// fsynced. This matters on disk-full: an unchecked ofstream close would
// otherwise replace the only valid checkpoint with a truncated one.
static SaveResult save_checkpoint(const CkptState& s, const vector<string>& all_samples,
                            const vector<MissRecord>& all_misses,
                            int threads, int block_bits, u64 p_anchor_limit,
                            u64 sample_limit, int qhot_size,
                            u64 minp_per_window, u64 minp_limit) {
    string tmp = string(CKPT_FILE) + ".tmp";
    ofstream f(tmp, ios::out | ios::trunc);
    if (!f.is_open()) return SaveResult::FAILED_BEFORE_RENAME;

    f << "# Goldbach verification checkpoint + samples\n";
    f << "# cfg_version="   << PROG_VERSION      << "\n";
    f << "# cfg_format="    << CKPT_FORMAT       << "\n";
    f << "# original_start=" << s.original_start << "\n";
    f << "# original_end="   << s.original_end   << "\n";
    f << "# last_verified="  << s.last_verified  << "\n";
    f << "# total_verified=" << s.total_verified << "\n";
    f << "# total_total="    << s.total_total    << "\n";
    f << "# qhot_probes="    << s.qhot_probes    << "\n";
    f << "# cold_candidates="<< s.cold_cand      << "\n";
    f << "# cold_survivors=" << s.cold_surv      << "\n";
    f << "# mr_calls="       << s.mr_calls       << "\n";
    f << "# total_qhot="     << s.total_qhot     << "\n";
    f << "# total_misses="   << s.total_misses   << "\n";
    f << "# unresolved="     << s.unresolved     << "\n";
    f << fixed << setprecision(6) << "# elapsed_sec=" << s.elapsed_sec << "\n";
    // Run configuration, so the artifact records how it was produced.
    f << "# cfg_threads="    << threads          << "\n";
    f << "# cfg_block_bits=" << block_bits       << "\n";
    f << "# cfg_p_anchor_limit=" << p_anchor_limit << "\n";
    f << "# cfg_sample_limit="   << sample_limit   << "\n";
    f << "# cfg_qhot_size="  << qhot_size        << "\n";
    f << "# cfg_minp_per_window=" << minp_per_window << "\n";   // (31)
    f << "# cfg_minp_limit="      << minp_limit      << "\n";
    f << "# source tags: QHOT = accepted from the QHot ring;"
         " COLD = accepted on the ordered anchor scan (Miller-Rabin).\n";
    f << "# SIEVE is a legacy alias for COLD and is accepted on resume;"
         " no build ever used a segmented sieve to accept a witness.\n";
    // (31) Provenance: state the scan order explicitly in the artifact.
    // (34) The listed first anchors match this build's anchor set.
    f << "# witness scan order: the cold path tries anchor primes p in"
#ifdef SEXY_ANCHORS
         " ascending order over sexy starters (5, 7, 11, 13, 17, ...),"
#elif defined(SINGLE_GAP)
         " ascending order over g-starters, g = " GB_STR(GAP) ","
#else
         " ascending order over twin/cousin starters (3, 5, 7, 11, 13, ...),"
#endif
         " so a COLD witness p is the minimal anchor partner of its N.\n";
    // (17) Authoritative miss records. They live inside this file so that
    // checkpoint counters and miss details are committed by ONE atomic
    // rename: there is no crash window in which the counters claim misses
    // whose details exist only in a second, not-yet-written file.
    // goldbach_misses.csv is a derived convenience report.
    f << "# miss rows: authoritative; status ANCHOR_EXHAUSTED or UNRESOLVED;"
         " p,q are 0 when unresolved.\n";
    for (const auto& m : all_misses) {
        f << "# miss=" << m.N << ","
          << (m.status == MISS_ANCHOR_EXHAUSTED ? m.p : 0) << ","
          << (m.status == MISS_ANCHOR_EXHAUSTED ? m.q : 0) << ","
          << (m.status == MISS_ANCHOR_EXHAUSTED ? "ANCHOR_EXHAUSTED" : "UNRESOLVED")
          << "\n";
    }
    f << "N,p,q,source\n";
    for (const auto& line : all_samples) f << line << "\n";

    f.flush();
    if (!f.good()) return SaveResult::FAILED_BEFORE_RENAME;  // short write
    f.close();
    if (!f.good()) return SaveResult::FAILED_BEFORE_RENAME;  // close-time flush
    // ofstream::close flushes to the OS, not to disk. fsync the file before
    // the rename, and the directory after it, so the checkpoint survives a
    // power loss on either side of the swap.
    if (!fsync_path(tmp.c_str())) return SaveResult::FAILED_BEFORE_RENAME;
    if (rename(tmp.c_str(), CKPT_FILE) != 0) return SaveResult::FAILED_BEFORE_RENAME;
    // (30) The rename has happened: the new checkpoint IS the visible one.
    // A directory-sync failure only means its survival across an immediate
    // power loss is unconfirmed -- reporting "previous checkpoint kept"
    // here would be untrue.
    if (!fsync_path(".")) return SaveResult::RENAMED_NOT_DURABLE;
    return SaveResult::SUCCESS;
}

// ------------------ Argument parsing ------------------

// (11)/(29) CLI wrappers over the same strict parsers used for checkpoint
// fields -- one implementation, one set of rules -- adding the argument name
// to the error message.
static u64 parse_u64(const string& s, const char* what) {
    try { return strict_u64(s); }
    catch (const exception&) {
        throw invalid_argument(string(what) +
            " must be a non-negative integer, got '" + s + "'");
    }
}

static int parse_int(const string& s, const char* what) {
    try { return strict_int(s); }
    catch (const exception&) {
        throw invalid_argument(string(what) +
            " must be a non-negative integer in int range, got '" + s + "'");
    }
}

// ------------------ Usage ------------------

static void usage(const char* prog) {
    cout << "Usage: " << prog
         << " [--resume] <start> <end> <block_bits> <threads>"
            " [p_anchor_limit] [sample_limit] [qhot_size]\n\n";
    cout << "  --resume   resume from " << CKPT_FILE
         << " (accepted in any position)\n";
    cout << "  --allow-miss-loss   resume even if earlier miss details cannot"
            " be fully restored\n";
    cout << "  --minp-per-window=K   (31) log the minimal anchor partner of the"
            " first K even N\n"
            "                        of every " << SIEVE_SAMPLE_INTERVAL
         << " stretch, to " << MINP_FILE << "\n";
    cout << "  --minp-limit=M        cap on minimal-partner rows kept"
            " (default 1000000)\n\n";
    cout << "Examples:\n";
    cout << "  " << prog << " 4000000000000000000 4001000000000000000 24 44 100000 10000 1024\n";
    cout << "  " << prog << " --resume 0 0 24 44 100000 10000 1024\n";
}

// ------------------ Main ------------------

int main(int argc, char** argv) {
    signal(SIGINT,  signal_handler);
    signal(SIGTERM, signal_handler);

    // (5) Flags stripped wherever they appear; positional args read from what remains.
    bool resuming = false;
    bool allow_miss_loss = false;
    string minp_pw_arg, minp_lim_arg;   // (31)
    vector<string> pos_args;
    for (int i = 1; i < argc; ++i) {
        string a = argv[i];
        if (a == "--resume")      resuming = true;
        else if (a == "--allow-miss-loss") allow_miss_loss = true;
        else if (a.rfind("--minp-per-window=", 0) == 0) minp_pw_arg  = a.substr(18);
        else if (a.rfind("--minp-limit=", 0) == 0)      minp_lim_arg = a.substr(13);
        else if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
        else                      pos_args.push_back(a);
    }

    // (15) At most 7 positionals: a stray extra argument is a typo, not
    // something to ignore silently.
    if (pos_args.size() < 4 || pos_args.size() > 7) { usage(argv[0]); return 1; }

    u64 start, end;
    int block_bits, threads, qhot_size;
    u64 p_anchor_limit, sample_limit;
    u64 minp_per_window, minp_limit;   // (31)
    try {
        start          = parse_u64(pos_args[0], "start");
        end            = parse_u64(pos_args[1], "end");
        block_bits     = parse_int(pos_args[2], "block_bits");
        threads        = parse_int(pos_args[3], "threads");
        p_anchor_limit = pos_args.size() > 4 ? parse_u64(pos_args[4], "p_anchor_limit") : 1000000ULL;
        sample_limit   = pos_args.size() > 5 ? parse_u64(pos_args[5], "sample_limit")   : 10000ULL;
        qhot_size      = pos_args.size() > 6 ? parse_int(pos_args[6], "qhot_size")      : 1024;
        minp_per_window = minp_pw_arg.empty()  ? 0ULL       : parse_u64(minp_pw_arg,  "--minp-per-window");
        minp_limit      = minp_lim_arg.empty() ? 1000000ULL : parse_u64(minp_lim_arg, "--minp-limit");
    } catch (const exception& e) {
        cerr << "ERROR: could not parse arguments: " << e.what() << "\n";
        return 1;
    }

    if (block_bits < 1 || block_bits > 30) {
        cerr << "ERROR: block_bits must be in [1, 30]. Got " << block_bits << "\n";
        return 1;
    }
    if (threads < 1) {
        cerr << "ERROR: threads must be >= 1. Got " << threads << "\n";
        return 1;
    }
    if (qhot_size <= 0 || (qhot_size & (qhot_size - 1)) != 0) {
        cerr << "ERROR: qhot_size must be a positive power of two. Got " << qhot_size << "\n";
        return 1;
    }
    if (p_anchor_limit < 3 || p_anchor_limit > numeric_limits<u32>::max()) {
        cerr << "ERROR: p_anchor_limit must be in [3, 2^32). Got " << p_anchor_limit << "\n";
        return 1;
    }
    // (31) One window per stretch: 2K even-offsets must fit in the stretch.
    if (minp_per_window > SIEVE_SAMPLE_INTERVAL / 2ULL) {
        cerr << "ERROR: --minp-per-window must be <= " << (SIEVE_SAMPLE_INTERVAL / 2ULL)
             << " (one window per " << SIEVE_SAMPLE_INTERVAL << " stretch). Got "
             << minp_per_window << "\n";
        return 1;
    }

    CkptState st;
    vector<string> all_samples;
    vector<MissRecord> all_misses;
    vector<string> all_minp;   // (31) rows "N,p,q", ascending N

    if (resuming) {
        CkptState ck;
        try {
            ck = load_checkpoint();
        } catch (const bad_alloc&) {
            cerr << "ERROR: out of memory loading checkpoint.\n";
            return 1;
        } catch (const length_error&) {
            cerr << "ERROR: checkpoint is too large to load.\n";
            return 1;
        }
        if (!ck.valid) {
            cerr << "ERROR: --resume given but no valid checkpoint in " << CKPT_FILE << "\n";
            return 1;
        }

        // (33) The checkpoint must come from THIS build. Different builds
        // verify different claims (plain Goldbach, twin/cousin starters, sexy
        // starters); resuming one with another would silently mix claims in
        // a single artifact. Restricted builds always write cfg_version.
        if (!ck.has_version) {
            cerr << "ERROR: checkpoint " << CKPT_FILE << " has no cfg_version;"
                    " this build (" << PROG_VERSION << ") only resumes"
                    " checkpoints it wrote itself.\n";
            return 1;
        }
        if (ck.cfg_version_str != PROG_VERSION) {
            cerr << "ERROR: checkpoint " << CKPT_FILE << " was written by '"
                 << ck.cfg_version_str << "', not by this build ('"
                 << PROG_VERSION << "'). Refusing to resume: the two builds"
                    " verify different claims.\n";
            return 1;
        }

        // (9) The witness sample is a pure function of (start, end,
        // block_bits, p_anchor_limit, sample_limit, qhot_size). Resuming
        // with different values would silently mix two parameter sets into
        // one artifact, so refuse. threads is deliberately not checked --
        // guarantee (7) makes the sample independent of it.
        if (ck.cfg_block_bits != 0) {  // cfg keys present in the checkpoint
            bool mismatch = false;
            auto req = [&](const char* name, u64 ck_v, u64 cli_v) {
                if (ck_v != cli_v) {
                    cerr << "ERROR: --resume config mismatch: " << name
                         << " is " << ck_v << " in the checkpoint but "
                         << cli_v << " on the command line.\n";
                    mismatch = true;
                }
            };
            req("block_bits",     (u64)ck.cfg_block_bits, (u64)block_bits);
            req("p_anchor_limit", ck.cfg_p_anchor_limit,  p_anchor_limit);
            req("sample_limit",   ck.cfg_sample_limit,    sample_limit);
            req("qhot_size",      (u64)ck.cfg_qhot_size,  (u64)qhot_size);
            // (31) The minp artifact is deterministic in these too; refuse a
            // silent change mid-campaign. Older checkpoints without the keys
            // resume with a NOTE: coverage simply starts at the resume point.
            if (ck.has_minp_cfg) {
                req("minp_per_window", ck.cfg_minp_per_window, minp_per_window);
                req("minp_limit",      ck.cfg_minp_limit,      minp_limit);
            } else if (minp_per_window > 0) {
                cout << "NOTE: checkpoint predates minimal-partner logging;"
                        " MINP coverage starts at the resume point.\n";
            }
            if (mismatch) {
                cerr << "       Re-run with the checkpoint's values to keep the witness\n"
                        "       sample reproducible, or start a fresh run without --resume.\n";
                return 1;
            }
        } else {
            cout << "NOTE: checkpoint has no cfg record (older build); config"
                    " consistency cannot be verified on this resume.\n";
        }

        st    = ck;
        start = ck.last_verified + 2;
        end   = ck.original_end;

        try {
            // (20) Resume artifact loading allocates; fail with a
            // message instead of an uncaught bad_alloc.
            ifstream rf(CKPT_FILE);
            string line;
            while (getline(rf, line)) {
                if (line.empty() || line[0] == '#') continue;
                if (line.rfind("N,", 0) == 0) continue;
                all_samples.push_back(line);
            }

            // (31) Reload minimal-partner rows so the derived log survives
            // resume. Best-effort by design: a crash between the checkpoint
            // rename and the minp write leaves this file one chunk behind,
            // which shows as a gap, never as a duplicate.
            if (minp_per_window > 0) {
                ifstream mf(MINP_FILE);
                string ml;
                while (mf.is_open() && getline(mf, ml)) {
                    if (ml.empty() || ml[0] == '#') continue;
                    if (ml.rfind("N,", 0) == 0) continue;
                    if ((u64)all_minp.size() >= minp_limit) break;
                    all_minp.push_back(ml);
                }
            }

            // (17) Miss records come from the checkpoint itself -- same atomic
            // artifact as the counters, so no cross-file crash window exists.
            if (!ck.misses.empty()) {
                all_misses = std::move(ck.misses);
            } else if (st.total_misses > 0) {
                // Transitional: checkpoint written before misses moved into the
                // unified file. Fall back to the separate report.
                u64 malformed = 0;
                all_misses = load_miss_file(malformed);
                if (malformed)
                    cerr << "WARNING: " << malformed << " malformed row(s) in "
                         << MISS_FILE << " were skipped.\n";
                cout << "NOTE: miss records loaded from " << MISS_FILE
                     << " (checkpoint predates unified miss storage).\n";
            }

            // (17) Any disagreement between restored miss records and the
            // checkpoint counter means earlier miss details are already lost or
            // duplicated. Continuing would rewrite the record as if it were
            // complete, so this blocks unless the loss is explicitly accepted.
            if ((u64)all_misses.size() != st.total_misses) {
                if (!allow_miss_loss) {
                    cerr << "ERROR: checkpoint counts " << st.total_misses
                         << " miss record(s) but " << all_misses.size()
                         << " could be restored. Earlier miss details are"
                         " incomplete.\n"
                         "       Re-run with --allow-miss-loss to resume anyway"
                         " (counters are kept; lost details are NOT recovered"
                         " and the miss record stays permanently incomplete).\n";
                    return 1;
                }
                cerr << "WARNING: resuming with " << all_misses.size() << " of "
                     << st.total_misses << " miss record(s) (--allow-miss-loss)."
                     " The miss record is permanently incomplete.\n";
            }

            // (22) Content checks. Corruption (as opposed to loss) is never
            // overridable: a wrong record must not be carried into the artifact.
            {
                const bool accepted_loss =
                    (u64)all_misses.size() != st.total_misses;
                string why = validate_misses(all_misses, st.original_start,
                                             st.last_verified, st.unresolved,
                                             st.cfg_p_anchor_limit,
                                             accepted_loss);
                if (!why.empty()) {
                    cerr << "ERROR: restored miss records failed validation: "
                         << why << "\n";
                    return 1;
                }
            }
        } catch (const bad_alloc&) {
            cerr << "ERROR: out of memory loading resume artifacts.\n";
            return 1;
        }

        cout << "RESUMING from checkpoint\n"
             << "  Original range: [" << st.original_start << " .. " << st.original_end << "]\n"
             << "  Last verified:  " << ck.last_verified << "\n"
             << "  Resume at:      " << start << "\n"
             << "  Samples kept:   " << all_samples.size() << "\n"
             << "  Misses kept:    " << all_misses.size() << "\n"
             << "  Minp rows kept: " << all_minp.size() << "\n"
             << "  Cumulative:     " << st.total_verified << " verified, "
             << fixed << setprecision(1) << st.elapsed_sec << " s\n\n";
    } else {
        st.original_start = start;
        st.original_end   = end;
        // (19) No early deletion of a stale miss file here: the report is
        // now rewritten atomically (header-only if empty) at every
        // checkpoint, so a previous run's file survives argument typos and
        // setup failures, and is replaced only once this run actually saves.
    }

    const u64 original_start = st.original_start;
    const u64 original_end   = st.original_end;

    if (start & 1ULL) ++start;
    if (end & 1ULL)   --end;
    if (start < 6) {
        // (11) Be explicit instead of silently clamping: p = 2 is excluded
        // by design, so 4 = 2 + 2 is outside this verifier's claim.
        cout << "NOTE: verification starts at N = 6. N = 4 (= 2 + 2) is excluded"
                " by design (p = 2 is not an anchor prime).\n";
        start = 6;
    }
#ifdef SEXY_ANCHORS
    if (start < 8) {
        // (32) [SX] 6 = 5 + 1 and 6 = 3 + 3, but 1 is not prime and 3 is not
        // a sexy starter (9 is composite): N = 6 has no sexy-starter
        // partition. Excluded explicitly, like N = 4 above.
        cout << "NOTE: [SX] verification starts at N = 8. N = 6 is excluded:"
                " it has no sexy-starter partition (3 is not a sexy starter).\n";
        start = 8;
    }
#elif defined(SINGLE_GAP)
    if (start < 8) {
        // (37) For 6 | g, 3 + g is divisible by 3, so 3 is never a g-starter and
        // 6 = 3 + 3 is unavailable: N = 6 is always an exception. Other small
        // exceptions (e.g. N = 8 for g = 30) are genuine members of E(g) and are
        // reported as misses -- start windows well above them.
        cout << "NOTE: [GAP " << GAP << "] verification starts at N = 8. N = 6 is"
                " always an exception for 6 | g (3 is never a g-starter).\n";
        start = 8;
    }
#endif
    if (start > end) {
        if (resuming) {
            cout << "Checkpoint indicates the run is already complete.\n";
            // (21) Report the checkpointed outcome, not unconditional
            // success: automation must see the same exit status whether the
            // result came from a fresh run or a resume of its checkpoint.
            if (st.unresolved   != 0) return 4;
            if (st.total_misses != 0) return 2;
            return 0;
        }
        cerr << "ERROR: invalid range after even adjustment.\n";
        return 1;
    }

    const u64 BLOCK = 1ULL << block_bits;
    const u64 SPAN  = 2ULL * BLOCK;

    // (11) Within one span of 2^64 the block builder's b + SPAN - 2 would
    // wrap, creating an empty block whose range is silently never verified
    // while the checkpoint still advances past it. Refuse outright; nothing
    // this close to 2^64 is a sensible target for this tool anyway.
    if (end > numeric_limits<u64>::max() - SPAN) {
        cerr << "ERROR: end must be <= 2^64 - 2^(block_bits+1) ("
             << (numeric_limits<u64>::max() - SPAN) << ") to avoid block"
             " arithmetic overflow. Got " << end << "\n";
        return 1;
    }

    // (15) Cap the per-chunk block count. The per-block bookkeeping
    // (blocks, block_wit, block_misses) is allocated per chunk; block_bits=1
    // against a full checkpoint-sized chunk would try to build ~5e10 entries
    // and die in bad_alloc. 2^24 blocks ~= 1 GB of bookkeeping, a sane ceiling.
    {
        u64 evens_in_range = (end - start) / 2ULL + 1ULL;
        u64 chunk_evens    = std::min<u64>(CKPT_INTERVAL, evens_in_range);
        u64 blocks_per_chunk = (chunk_evens + BLOCK - 1ULL) / BLOCK;
        if (blocks_per_chunk > (1ULL << 24)) {
            cerr << "ERROR: block_bits=" << block_bits << " gives "
                 << blocks_per_chunk << " blocks per chunk (max 2^24)."
                 " Increase block_bits.\n";
            return 1;
        }
    }

    cout << PROG_VERSION << "\n";
    cout << "Range: [" << start << " .. " << end << "]\n";
    cout << "Block bits: " << block_bits << " | Threads: " << threads << "\n";
    cout << "P-anchor limit: " << p_anchor_limit << " | QHot size: " << qhot_size << "\n";
    cout << "Checkpoint interval: " << CKPT_INTERVAL << " even numbers\n";
    cout << "Sample limit: " << sample_limit << " | Escalation limit: " << ESCALATION_LIMIT << "\n";
    if (minp_per_window > 0)
        cout << "Minp per window: " << minp_per_window << " | Minp limit: "
             << minp_limit << " | Minp file: " << MINP_FILE << "\n";
    cout << "Building base primes up to " << p_anchor_limit << "..." << flush;

    auto setup0 = chrono::steady_clock::now();
    // (15) The sieve and the lookup table scale with p_anchor_limit (up to
    // ~4.5 GB combined at the 2^32 ceiling). Fail with a message, not an
    // uncaught bad_alloc, if the machine can't hold them.
    vector<u8> anchor_lookup;
    try {
        build_base_primes((u32)p_anchor_limit);
        anchor_lookup.assign((size_t)p_anchor_limit + 1, 0);
    } catch (const bad_alloc&) {
        cerr << "\nERROR: out of memory building the anchor tables for"
                " p_anchor_limit=" << p_anchor_limit
             << ". Lower p_anchor_limit.\n";
        return 1;
    }
    double setup_sec = chrono::duration<double>(chrono::steady_clock::now() - setup0).count();

#ifdef SEXY_ANCHORS
    cout << " done. Anchor primes (sexy starters only, p+6 prime): " << base_primes.size()
#elif defined(SINGLE_GAP)
    cout << " done. Anchor primes (g-starters only, g = " << GAP << "): " << base_primes.size()
#else
    cout << " done. Anchor primes (twin/cousin starters only): " << base_primes.size()
#endif
         << " | Setup: " << fixed << setprecision(3) << setup_sec << " s\n";

    for (u32 p : base_primes) {
        if ((u64)p <= p_anchor_limit) anchor_lookup[(size_t)p] = 1;
        else break;
    }

    // [B] Odd-indexed anchor bitset: bit k set iff p = 2k+1 is an anchor prime.
    // ~6 KB at p_anchor_limit = 1e5, plus two zero padding words.
    vector<u64> abits((size_t)((p_anchor_limit / 2ULL) / 64ULL) + 3, 0ULL);
    for (u32 p : base_primes) {
        if ((u64)p > p_anchor_limit) break;
        const u64 k = ((u64)p - 1ULL) >> 1;
        abits[(size_t)(k >> 6)] |= 1ULL << (k & 63ULL);
    }

    // Witness sampling: at most one accepted witness per N-distance bucket per
    // path. Both paths are sampled so the artifact certifies the QHot fast path
    // as well as the cold path.
    u64 nbuckets = (original_end >= original_start)
                 ? ((original_end - original_start) / SIEVE_SAMPLE_INTERVAL + 2ULL) : 1ULL;
    // (20) These scale with the range (~2 bytes per 10^10 of N) and can
    // reach gigabytes for extreme ranges; fail with a message, not an
    // uncaught bad_alloc.
    vector<u8> cold_bucket_seen, qhot_bucket_seen;
    try {
        cold_bucket_seen.assign((size_t)nbuckets, 0);
        qhot_bucket_seen.assign((size_t)nbuckets, 0);
    } catch (const bad_alloc&) {
        cerr << "ERROR: out of memory allocating " << nbuckets
             << " sampling buckets for this range.\n";
        return 1;
    }

    u64 legacy_sieve_rows = 0, unknown_tag_rows = 0;
    for (const string& s : all_samples) {
        size_t c1 = s.find(',');
        size_t c3 = s.rfind(',');
        if (c1 == string::npos || c3 == string::npos) continue;
        try {
            u64 sN = stoull(s.substr(0, c1));
            if (sN < original_start) continue;
            u64 b = (sN - original_start) / SIEVE_SAMPLE_INTERVAL;
            if (b >= nbuckets) continue;

            const string src = s.substr(c3 + 1);
            if (src == "QHOT") {
                qhot_bucket_seen[(size_t)b] = 1;
            } else if (src == "COLD") {
                cold_bucket_seen[(size_t)b] = 1;
            } else if (src == "SIEVE") {
                // Legacy alias. No build of this program ever used a segmented
                // sieve to accept a witness -- the cold path has always been
                // Miller-Rabin. Rows tagged SIEVE come from an older label and
                // are cold-path witnesses.
                cold_bucket_seen[(size_t)b] = 1;
                ++legacy_sieve_rows;
            } else {
                ++unknown_tag_rows;
            }
        } catch (...) {}
    }
    if (legacy_sieve_rows)
        cout << "NOTE: " << legacy_sieve_rows
             << " sample row(s) carry the legacy SIEVE tag; treated as COLD.\n";
    if (unknown_tag_rows)
        cerr << "WARNING: " << unknown_tag_rows
             << " sample row(s) have an unrecognised source tag and were ignored.\n";

    // ---- Mega-chunk loop ----
    u64 run_total = 0, run_verified = 0, run_qhot = 0, run_misses = 0, run_unres = 0;
    u64 run_probes = 0, run_cand = 0, run_surv = 0, run_mr = 0;
    // (38) Minimal-partner window statistics for this process (merged chunks only).
    struct MinpStat { u64 max_p = 0, N_at_max = 0, sum_p = 0, count = 0, none = 0; };
    MinpStat run_ms;
    double run_elapsed = 0.0;
    u64 chunk_start = start;
    int ckpt_count = 0;
    // (16) True whenever a valid checkpoint exists on disk -- loaded by
    // --resume or saved by this process. ckpt_count alone misreports a
    // resumed run interrupted before its first new save.
    bool have_checkpoint = resuming;
    // (28) Set only when THIS process wrote the derived report successfully;
    // the final summary must not point at a report that failed to write or
    // predates the last checkpoint.
    bool have_current_miss_report = false;

#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif

    while (chunk_start <= end) {
        try {   // (15) bad_alloc in per-chunk bookkeeping or merge => fatal,
                // chunk discarded, last checkpoint intact. The parallel
                // region itself cannot throw (10).
        u64 chunk_end = chunk_start + CKPT_INTERVAL * 2ULL - 2ULL;
        if (chunk_end > end || chunk_end < chunk_start) chunk_end = end;

        vector<pair<u64, u64>> blocks;
        for (u64 b = chunk_start; b <= chunk_end; ) {
            u64 blk_end = std::min<u64>(b + SPAN - 2ULL, chunk_end);
            blocks.push_back({b, blk_end});
            if (chunk_end - b < SPAN) break;
            b += SPAN;
        }

        u64 c_total = 0, c_verified = 0, c_qhot = 0, c_misses = 0;
        u64 c_probes = 0, c_cand = 0, c_surv = 0, c_mr = 0;
        vector<string> chunk_samples;
        vector<MissRecord> chunk_misses;
        vector<string> chunk_minp;   // (31)

        // (7) Per-block output. Each element is written by exactly one thread
        // (the owner of that block), so no synchronisation is needed and the
        // contents do not depend on scheduling.
        vector<BlockWitness> block_wit(blocks.size());
        vector<vector<MissRecord>> block_misses(blocks.size());
        // (31) Minimal-partner rows, written by the block's owning thread
        // only, merged serially in block order like misses and witnesses.
        vector<vector<string>> block_minp(blocks.size());
        vector<MinpStat> bstat(blocks.size());   // (38) per-block window statistics
        MinpStat chunk_ms;
#ifdef DUMP_ALL
        struct DumpRow { u64 N, p, q; u8 src; };   // src: 0 ring, 1 propagated, 2 cold
        vector<vector<DumpRow>> dump_rows(blocks.size());
#endif

        auto t0 = chrono::steady_clock::now();

#pragma omp parallel
        {
            QHot qhot;
#ifdef CARRY
            // [CARRY] Words carried between consecutive batches for the AVX-512 prefix.
            // When N0 advances by 128 and the ring is unchanged, each prefix entry's
            // partner grows by 128: its word index advances by exactly one and its
            // shift is unchanged, so this batch's low word is the previous batch's
            // high word. Valid only if the previous batch took the prefix at N0-128
            // with the same ring version.
            alignas(64) u64 carry_hi[KFIX];
#ifdef CARRYIDX
            alignas(64) u64 carry_idx[KFIX], carry_sr[KFIX], carry_sl[KFIX];   // [CARRYIDX] next hi word index, shifts
#endif
            u64 carry_N0 = 0, carry_ver = ~0ULL;
            bool carry_ok = false;
#endif
            // (10) init allocates; a bad_alloc here must not escape the
            // structured block. The thread still reaches the omp-for below
            // (all threads must, or the implicit barrier deadlocks) and
            // simply skips every block via the shutdown check.
            bool thread_ok = true;
            try {
                qhot.init(qhot_size);
            } catch (...) {
                thread_ok = false;
                g_fatal.store(1, memory_order_relaxed);
                g_shutdown.store(1, memory_order_relaxed);
            }
            // (7) Sampling is no longer restricted to thread 0. Every thread
            // records a candidate for every block it owns; the serial pass
            // below decides which candidates survive.
            const bool sampling = (sample_limit > 0);

            u64 l_total = 0, l_verified = 0, l_qhot = 0, l_misses = 0;
            u64 l_probes = 0, l_cand = 0, l_surv = 0, l_mr = 0;
#ifdef N3FIRST
            u64 l_n3_tests = 0, l_n3_hits = 0;   // (36)
#endif

#pragma omp for schedule(dynamic)
            for (size_t bi = 0; bi < blocks.size(); ++bi) {
                if (!thread_ok || g_shutdown.load(memory_order_relaxed)) continue;

                // (10) The only allocation in this body is the miss-record
                // push_back, but the catch is written to cover the whole
                // block so the claim holds by construction, not by audit.
                try {

                auto [blk, blk_end] = blocks[bi];

                // (7) Start every block from an empty ring. Without this, the
                // witness accepted for a given N depends on which blocks this
                // thread processed before it -- i.e. on the scheduler. Warm-up
                // costs a few thousand cold scans out of 2^block_bits numbers.
                qhot.reset();

                // Sampling is resolved per block, not per N: remember the first
                // witness of each kind, then map it to a bucket after the block.
                bool have_qs = false, have_cs = false;
                u64 qsN = 0, qsp = 0, qsq = 0, csN = 0, csp = 0, csq = 0;

                // ---- [B] Batched QHot: 64 consecutive even N per step ----
                // For N_i = N0 + 2i (i < 64) and a ring entry q, p_i = N0 - q + 2i
                // runs over consecutive odd numbers, so "which N_i does q certify?"
                // is one 64-bit slice of the odd-indexed anchor bitset (bit k <-> p =
                // 2k+1), starting at bit (N0 - q - 1)/2. The bitset (~6 KB at
                // p_anchor_limit = 1e5) stays in L1. Bits for p > p_anchor_limit are
                // zero by construction, so a slice never certifies a non-anchor p.
                for (u64 N0 = blk; N0 <= blk_end; N0 += 128ULL) {
                    const u64 left = (blk_end - N0) / 2ULL + 1ULL;
                    const unsigned cnt = (unsigned)std::min<u64>(64ULL, left);
                    const u64 full = (cnt == 64) ? ~0ULL : ((1ULL << cnt) - 1ULL);
                    l_total += cnt;

                    // (31) [MINP] Minimal-partner scan for windowed N. Runs
                    // before any ring walk or batch shortcut and consults
                    // neither, so the logged p is the true minimal anchor
                    // partner (same ascending order as the cold path) and a
                    // pure function of N -- reproducible across thread
                    // counts, schedules and builds. The batch-level overlap
                    // precheck keeps the non-window cost to one modulo and
                    // two compares per 64 N.
                    if (minp_per_window > 0) {
                        const u64 twoK = 2ULL * minp_per_window;
                        const u64 off0 = (N0 - original_start) % SIEVE_SAMPLE_INTERVAL;
                        if (off0 < twoK ||
                            off0 + 2ULL * (u64)(cnt - 1) >= SIEVE_SAMPLE_INTERVAL) {
                            u64 off = off0;
                            for (unsigned i = 0; i < cnt; ++i) {
                                if (off < twoK) {
                                    const u64 Nw = N0 + 2ULL * i;
                                    // [TC-INST] Three minima in one ascending
                                    // scan: the overall minimal starter
                                    // partner (p_min, as before), the minimal
                                    // TWIN-starter partner and the minimal
                                    // COUSIN-starter partner. A prime is only
                                    // primality-tested when it can still
                                    // improve an unfilled minimum, so the
                                    // extra cost over the old scan is small.
                                    u64 p_min = 0, q_min = 0, p_tw = 0, p_co = 0;
                                    u8  ty = 0;
                                    for (size_t pi = 0; pi < base_primes.size(); ++pi) {
                                        const u32 pr = base_primes[pi];
                                        if ((u64)pr > p_anchor_limit || (u64)pr >= Nw) break;
                                        const u8 t = starter_types[pi];
                                        if ((p_tw || !(t & 1)) && (p_co || !(t & 2)))
                                            continue;   // fills no open slot
                                        const u64 qw = Nw - pr;
                                        if (!survives_small_trial_division(qw)) continue;
                                        if (is_prime64(qw)) {
                                            if (!p_min) { p_min = pr; q_min = qw; ty = t; }
                                            if (!p_tw && (t & 1)) p_tw = pr;
                                            if (!p_co && (t & 2)) p_co = pr;
#ifdef SINGLE_GAP
                                            if (p_tw) break;   // [SX]/(37) single gear: first hit is final
#else
                                            if (p_tw && p_co) break;
#endif
                                        }
                                    }
                                    // (38) Window statistics, kept per block and
                                    // merged serially in block order, so they are
                                    // independent of threads and of --minp-limit.
                                    if (p_min) {
                                        if (p_min > bstat[bi].max_p) { bstat[bi].max_p = p_min; bstat[bi].N_at_max = Nw; }
                                        bstat[bi].sum_p += p_min; ++bstat[bi].count;
                                    } else {
                                        ++bstat[bi].none;
                                    }
                                    if (p_min && minp_limit > 0) {   // (38) no row strings when rows are not kept
                                        // type: T = twin-only, C = cousin-only,
                                        // B = both (p, p+2, p+4 all prime).
#ifdef SEXY_ANCHORS
                                        const char tc = 'S';   // [SX]
                                        (void)ty;              // (35) single type in SX
#elif defined(SINGLE_GAP)
                                        const char tc = 'G';   // (37)
                                        (void)ty;
#else
                                        const char tc = (ty == 1) ? 'T' : (ty == 2) ? 'C' : 'B';
#endif
                                        block_minp[bi].push_back(
                                            to_string(Nw) + "," + to_string(p_min) +
                                            "," + to_string(q_min) + "," + tc +
                                            "," + to_string(p_tw) + "," + to_string(p_co));
                                    }
                                    // No anchor partner: that N is a cold-path
                                    // miss and lives in the miss records; no
                                    // MINP row is written for it.
                                }
                                off += 2ULL;
                                if (off >= SIEVE_SAMPLE_INTERVAL) off -= SIEVE_SAMPLE_INTERVAL;
                            }
                        }
                    }

#ifdef N3FIRST
                    // (36) [N3FIRST] Lemma-first check, before the ring or any other
                    // probe: for N == 2 or 4 (mod 6), try N = 3 + (N - 3). 3 is both a
                    // twin starter (3, 5) and a cousin starter (3, 7), and by Lemmas
                    // 1-2 it is the only possible partner of the MISMATCHED gear in
                    // these classes. Hits are certified here and kept out of the
                    // QHot counters; the MR calls are counted in mr_calls. q is NOT
                    // recorded in the ring: that would bump the ring version on most
                    // batches and defeat the [CARRY] words.
                    u64 n3mask = 0;
                    {
                        unsigned r6 = (unsigned)(N0 % 6ULL);
                        for (unsigned i = 0; i < cnt; ++i) {
                            if (r6 != 0) {
                                const u64 N = N0 + 2ULL * i;
                                if (N >= 8) {
                                    ++l_n3_tests;
                                    const u64 q = N - 3ULL;
                                    if (survives_small_trial_division(q)) {
                                        ++l_mr;
                                        if (is_prime64(q)) n3mask |= 1ULL << i;
                                    }
                                }
                            }
                            r6 += 2; if (r6 >= 6) r6 -= 6;
                        }
                        l_n3_hits += (u64)__builtin_popcountll(n3mask);
                    }
                    u64 done = n3mask;
#else
                    // ---- Hot: walk the ring newest-first until the batch is covered ----
                    u64 done = 0;
#endif
#ifdef DUMP_ALL
                    u64 certq[64] = {0};   // [DUMP] q that first certified each bit
#endif
                    if (qhot.count > 0) {
                        int t0 = 0;
#if KFIX > 0
                        // [KFIX] Fixed-length, unrolled, branch-free prefix, after the
                        // tuned-I inner loop of Oliveira e Silva et al. (2014, Alg. 1.4):
                        // OR the first KFIX slices with no per-probe exit test, then test
                        // once. Taken only when all KFIX entries are usable: the ring is
                        // sorted (descending q), so buf[0] < N0 bounds every entry and
                        // entry KFIX-1 has the largest partner of the prefix.
                        if (qhot.count >= KFIX && qhot.buf[0] < N0 &&
                            N0 - qhot.buf[(size_t)(KFIX - 1)] <= p_anchor_limit) {
                            const u64* __restrict qb = qhot.buf.data();
#ifdef SIMD_CHECK
                            const u64 done_before = done;
#endif
#if defined(DUMP_ALL) || (!defined(SIMD_NEON) && !defined(KACC4) && !defined(SIMD_AVX512) && !defined(SIMD_AVX2))
#if defined(CARRY)
                            {   // [CARRY] scalar prefix: one load per entry when the previous batch's
                                // high words are still valid, two otherwise
                                const u64* ab = abits.data();
                                const bool use = carry_ok && carry_ver == qhot.ver && carry_N0 + 128ULL == N0;
                                u64 r = 0;
                                if (use) {
                                    #pragma GCC unroll 64
                                    for (int t = 0; t < KFIX; ++t) {
                                        const u64 b = (N0 - qb[t] - 1ULL) >> 1; const unsigned sh = (unsigned)(b & 63ULL);
                                        const u64 lo = carry_hi[t], hi = ab[(b >> 6) + 1];
                                        carry_hi[t] = hi;
                                        r |= sh ? ((lo >> sh) | (hi << (64U - sh))) : lo;
                                    }
                                } else {
                                    #pragma GCC unroll 64
                                    for (int t = 0; t < KFIX; ++t) {
                                        const u64 b = (N0 - qb[t] - 1ULL) >> 1; const unsigned sh = (unsigned)(b & 63ULL);
                                        const u64 lo = ab[b >> 6], hi = ab[(b >> 6) + 1];
                                        carry_hi[t] = hi;
                                        r |= sh ? ((lo >> sh) | (hi << (64U - sh))) : lo;
                                    }
                                }
                                carry_ok = true; carry_ver = qhot.ver; carry_N0 = N0;
                                done |= r;
                            }
#else
                            #pragma GCC unroll 64
                            for (int t = 0; t < KFIX; ++t) {
                                const u64 d = N0 - qb[t];
#ifdef DUMP_ALL
                                { u64 nb = anchor_slice(abits, (d - 1ULL) >> 1) & full & ~done;
                                  while (nb) { certq[__builtin_ctzll(nb)] = qb[t]; nb &= nb - 1ULL; } }
#endif
                                done |= anchor_slice(abits, (d - 1ULL) >> 1);
                            }
#endif
#elif defined(KACC4)
                            // [ACC4] four independent accumulators: breaks the serial OR
                            // chain on `done` so slice work from different entries overlaps.
                            {
                                u64 a0 = 0, a1 = 0, a2 = 0, a3 = 0;
                                const u64* ab = abits.data();
                                #pragma GCC unroll 16
                                for (int t = 0; t < KFIX; t += 4) {
                                    const u64 b0 = (N0 - qb[t]     - 1ULL) >> 1, b1 = (N0 - qb[t + 1] - 1ULL) >> 1;
                                    const u64 b2 = (N0 - qb[t + 2] - 1ULL) >> 1, b3 = (N0 - qb[t + 3] - 1ULL) >> 1;
                                    const unsigned s0 = b0 & 63, s1 = b1 & 63, s2 = b2 & 63, s3 = b3 & 63;
                                    const u64* p0 = ab + (b0 >> 6); const u64* p1 = ab + (b1 >> 6);
                                    const u64* p2 = ab + (b2 >> 6); const u64* p3 = ab + (b3 >> 6);
                                    // (x << 1) << (63 - s) == x << (64 - s), and is 0 when s == 0
                                    a0 |= (p0[0] >> s0) | ((p0[1] << 1) << (63 - s0));
                                    a1 |= (p1[0] >> s1) | ((p1[1] << 1) << (63 - s1));
                                    a2 |= (p2[0] >> s2) | ((p2[1] << 1) << (63 - s2));
                                    a3 |= (p3[0] >> s3) | ((p3[1] << 1) << (63 - s3));
                                }
                                done |= (a0 | a1) | (a2 | a3);
                            }
#elif defined(SIMD_AVX512)
                            // [AVX512] eight ring entries per 512-bit vector: two gathers
                            // (word w and w+1) and per-lane variable shifts. VPSRLVQ/VPSLLVQ
                            // give 0 for counts >= 64, so the s == 0 case needs no branch.
                            {
                                const long long* ab = (const long long*)abits.data();
                                const __m512i one = _mm512_set1_epi64(1), m63 = _mm512_set1_epi64(63);
                                const __m512i k64 = _mm512_set1_epi64(64), vN0 = _mm512_set1_epi64((long long)N0);
                                __m512i acc = _mm512_setzero_si512();
#ifdef CARRY
                                if (carry_ok && carry_ver == qhot.ver && carry_N0 + 128ULL == N0) {
                                    // [CARRY] one gather per vector: low words come from the previous batch
                                    #pragma GCC unroll 8
                                    for (int t = 0; t < KFIX; t += 8) {
#ifdef CARRYIDX
                                        // [CARRYIDX] index advances by one word per batch; shifts unchanged
                                        const __m512i idx = _mm512_add_epi64(_mm512_load_si512((const void*)(carry_idx + t)), one);
                                        _mm512_store_si512((void*)(carry_idx + t), idx);
                                        const __m512i sh = _mm512_load_si512((const void*)(carry_sr + t));
                                        const __m512i sl = _mm512_load_si512((const void*)(carry_sl + t));
                                        const __m512i lo = _mm512_load_si512((const void*)(carry_hi + t));
                                        const __m512i hi = _mm512_i64gather_epi64(idx, ab, 8);
                                        _mm512_store_si512((void*)(carry_hi + t), hi);
                                        acc = _mm512_or_si512(acc, _mm512_or_si512(_mm512_srlv_epi64(lo, sh), _mm512_sllv_epi64(hi, sl)));
#else
                                        const __m512i q = _mm512_loadu_si512((const void*)(qb + t));
                                        const __m512i b = _mm512_srli_epi64(_mm512_sub_epi64(_mm512_sub_epi64(vN0, q), one), 1);
                                        const __m512i w = _mm512_srli_epi64(b, 6);
                                        const __m512i sh = _mm512_and_si512(b, m63);
                                        const __m512i lo = _mm512_load_si512((const void*)(carry_hi + t));
                                        const __m512i hi = _mm512_i64gather_epi64(_mm512_add_epi64(w, one), ab, 8);
                                        _mm512_store_si512((void*)(carry_hi + t), hi);
                                        acc = _mm512_or_si512(acc, _mm512_or_si512(_mm512_srlv_epi64(lo, sh),
                                                                                   _mm512_sllv_epi64(hi, _mm512_sub_epi64(k64, sh))));
#endif
                                    }
                                } else
#endif
                                {
                                #pragma GCC unroll 8
                                for (int t = 0; t < KFIX; t += 8) {
                                    const __m512i q = _mm512_loadu_si512((const void*)(qb + t));
                                    const __m512i b = _mm512_srli_epi64(_mm512_sub_epi64(_mm512_sub_epi64(vN0, q), one), 1);
                                    const __m512i w = _mm512_srli_epi64(b, 6);
                                    const __m512i sh = _mm512_and_si512(b, m63);
                                    const __m512i lo = _mm512_i64gather_epi64(w, ab, 8);
                                    const __m512i hi = _mm512_i64gather_epi64(_mm512_add_epi64(w, one), ab, 8);
#ifdef CARRY
                                    _mm512_store_si512((void*)(carry_hi + t), hi);
#endif
#ifdef CARRYIDX
                                    _mm512_store_si512((void*)(carry_idx + t), _mm512_add_epi64(w, one));
                                    _mm512_store_si512((void*)(carry_sr + t), sh);
                                    _mm512_store_si512((void*)(carry_sl + t), _mm512_sub_epi64(k64, sh));
#endif
                                    acc = _mm512_or_si512(acc, _mm512_or_si512(_mm512_srlv_epi64(lo, sh),
                                                                               _mm512_sllv_epi64(hi, _mm512_sub_epi64(k64, sh))));
                                }
                                }
#ifdef CARRY
                                carry_ok = true; carry_ver = qhot.ver; carry_N0 = N0;
#endif
                                done |= (u64)_mm512_reduce_or_epi64(acc);
                            }
#elif defined(SIMD_AVX2)
                            // [AVX2] four ring entries per 256-bit vector (same scheme).
                            {
                                const long long* ab = (const long long*)abits.data();
                                const __m256i one = _mm256_set1_epi64x(1), m63 = _mm256_set1_epi64x(63);
                                const __m256i k64 = _mm256_set1_epi64x(64), vN0 = _mm256_set1_epi64x((long long)N0);
                                __m256i acc = _mm256_setzero_si256();
#ifdef CARRY
                                const bool use = carry_ok && carry_ver == qhot.ver && carry_N0 + 128ULL == N0;
#endif
#if defined(CARRY) && defined(CARRYIDX)
                                if (use) {
                                    // [CARRYIDX] AVX2: index advances one word per batch, shifts unchanged
                                    #pragma GCC unroll 16
                                    for (int t = 0; t < KFIX; t += 4) {
                                        const __m256i idx = _mm256_add_epi64(_mm256_load_si256((const __m256i*)(carry_idx + t)), one);
                                        _mm256_store_si256((__m256i*)(carry_idx + t), idx);
                                        const __m256i sh = _mm256_load_si256((const __m256i*)(carry_sr + t));
                                        const __m256i sl = _mm256_load_si256((const __m256i*)(carry_sl + t));
                                        const __m256i lo = _mm256_load_si256((const __m256i*)(carry_hi + t));
                                        const __m256i hi = _mm256_i64gather_epi64(ab, idx, 8);
                                        _mm256_store_si256((__m256i*)(carry_hi + t), hi);
                                        acc = _mm256_or_si256(acc, _mm256_or_si256(_mm256_srlv_epi64(lo, sh), _mm256_sllv_epi64(hi, sl)));
                                    }
                                } else
#endif
                                {
                                #pragma GCC unroll 8
                                for (int t = 0; t < KFIX; t += 4) {
                                    const __m256i q = _mm256_loadu_si256((const __m256i*)(qb + t));
                                    const __m256i b = _mm256_srli_epi64(_mm256_sub_epi64(_mm256_sub_epi64(vN0, q), one), 1);
                                    const __m256i w = _mm256_srli_epi64(b, 6);
                                    const __m256i sh = _mm256_and_si256(b, m63);
#ifdef CARRY
                                    const __m256i lo = use ? _mm256_load_si256((const __m256i*)(carry_hi + t))
                                                           : _mm256_i64gather_epi64(ab, w, 8);
#else
                                    const __m256i lo = _mm256_i64gather_epi64(ab, w, 8);
#endif
                                    const __m256i hi = _mm256_i64gather_epi64(ab, _mm256_add_epi64(w, one), 8);
#ifdef CARRY
                                    _mm256_store_si256((__m256i*)(carry_hi + t), hi);
#endif
#if defined(CARRY) && defined(CARRYIDX)
                                    _mm256_store_si256((__m256i*)(carry_idx + t), _mm256_add_epi64(w, one));
                                    _mm256_store_si256((__m256i*)(carry_sr + t), sh);
                                    _mm256_store_si256((__m256i*)(carry_sl + t), _mm256_sub_epi64(k64, sh));
#endif
                                    acc = _mm256_or_si256(acc, _mm256_or_si256(_mm256_srlv_epi64(lo, sh),
                                                                               _mm256_sllv_epi64(hi, _mm256_sub_epi64(k64, sh))));
                                }
                                }
#ifdef CARRY
                                carry_ok = true; carry_ver = qhot.ver; carry_N0 = N0;
#endif
                                const __m128i r = _mm_or_si128(_mm256_castsi256_si128(acc), _mm256_extracti128_si256(acc, 1));
                                done |= (u64)_mm_extract_epi64(r, 0) | (u64)_mm_extract_epi64(r, 1);
                            }
#else
                            // [NEON] two ring entries per 128-bit vector. USHL by 64 yields 0,
                            // so the s == 0 case needs no branch.
                            {
                                const u64* ab = abits.data();
                                uint64x2_t acc0 = vdupq_n_u64(0), acc1 = vdupq_n_u64(0);
#ifdef CARRY
                                const bool use = carry_ok && carry_ver == qhot.ver && carry_N0 + 128ULL == N0;
#endif
                                #pragma GCC unroll 16
                                for (int t = 0; t < KFIX; t += 4) {
                                    const u64 b0 = (N0 - qb[t]     - 1ULL) >> 1, b1 = (N0 - qb[t + 1] - 1ULL) >> 1;
                                    const u64 b2 = (N0 - qb[t + 2] - 1ULL) >> 1, b3 = (N0 - qb[t + 3] - 1ULL) >> 1;
                                    const u64* p0 = ab + (b0 >> 6); const u64* p1 = ab + (b1 >> 6);
                                    const u64* p2 = ab + (b2 >> 6); const u64* p3 = ab + (b3 >> 6);
#ifdef CARRY
                                    const uint64x2_t lo01 = use ? vld1q_u64(carry_hi + t)     : vcombine_u64(vcreate_u64(p0[0]), vcreate_u64(p1[0]));
                                    const uint64x2_t lo23 = use ? vld1q_u64(carry_hi + t + 2) : vcombine_u64(vcreate_u64(p2[0]), vcreate_u64(p3[0]));
#else
                                    const uint64x2_t lo01 = vcombine_u64(vcreate_u64(p0[0]), vcreate_u64(p1[0]));
                                    const uint64x2_t lo23 = vcombine_u64(vcreate_u64(p2[0]), vcreate_u64(p3[0]));
#endif
                                    const uint64x2_t hi01 = vcombine_u64(vcreate_u64(p0[1]), vcreate_u64(p1[1]));
                                    const uint64x2_t hi23 = vcombine_u64(vcreate_u64(p2[1]), vcreate_u64(p3[1]));
#ifdef CARRY
                                    vst1q_u64(carry_hi + t, hi01); vst1q_u64(carry_hi + t + 2, hi23);
#endif
                                    const int64x2_t s01 = vcombine_s64(vcreate_s64((int64_t)(b0 & 63)), vcreate_s64((int64_t)(b1 & 63)));
                                    const int64x2_t s23 = vcombine_s64(vcreate_s64((int64_t)(b2 & 63)), vcreate_s64((int64_t)(b3 & 63)));
                                    const int64x2_t k64 = vdupq_n_s64(64);
                                    acc0 = vorrq_u64(acc0, vorrq_u64(vshlq_u64(lo01, vnegq_s64(s01)), vshlq_u64(hi01, vsubq_s64(k64, s01))));
                                    acc1 = vorrq_u64(acc1, vorrq_u64(vshlq_u64(lo23, vnegq_s64(s23)), vshlq_u64(hi23, vsubq_s64(k64, s23))));
                                }
                                const uint64x2_t acc = vorrq_u64(acc0, acc1);
#ifdef CARRY
                                carry_ok = true; carry_ver = qhot.ver; carry_N0 = N0;
#endif
                                done |= vgetq_lane_u64(acc, 0) | vgetq_lane_u64(acc, 1);
                            }
#endif
#ifdef SIMD_CHECK
                            {   // [CHECK] recompute the prefix with the scalar slice and compare
                                u64 ref = done_before;
                                for (int t = 0; t < KFIX; ++t) ref |= anchor_slice(abits, (N0 - qb[t] - 1ULL) >> 1);
                                if (ref != done) { fprintf(stderr, "SIMD_CHECK MISMATCH at N0=%llu\n", (unsigned long long)N0); abort(); }
                                ++simd_checks;
                            }
#endif
                            l_probes += (u64)KFIX;
                            t0 = ((done & full) == full) ? qhot.count : KFIX;
                        }
#ifdef CARRY
                        else { carry_ok = false; }   // [CARRY] prefix skipped: carried words are stale
#endif
#endif
                        for (int t = t0; t < qhot.count; ++t) {
                            const u64 q = qhot.buf[(size_t)t];  // [SORT] descending q = ascending d
                            ++l_probes;
                            if (q >= N0) continue;              // cannot occur (q <= N0 - 5); guard
                            const u64 d = N0 - q;               // p for i = 0; odd
                            if (d > p_anchor_limit) break;      // [SORT] every later entry is staler
#ifdef DUMP_ALL
                            { u64 nb = anchor_slice(abits, (d - 1ULL) >> 1) & full & ~done;
                              while (nb) { certq[__builtin_ctzll(nb)] = q; nb &= nb - 1ULL; } }
#endif
                            done |= anchor_slice(abits, (d - 1ULL) >> 1);
                            if ((done & full) == full) break;
                        }
                        done &= full;
                    }
#ifdef N3FIRST
                    const u64 hot_from_ring = done & ~n3mask;   // (36) ring-only bits
#else
                    const u64 hot_from_ring = done;
#endif

                    // Sample: lowest ring-certified N, witness re-derived with the
                    // unchanged v4 per-N probe while the ring is still as walked.
                    if (hot_from_ring && !have_qs) {
                        const u64 Ns = N0 + 2ULL * (u64)__builtin_ctzll(hot_from_ring);
                        u64 sp = 0, sq = 0;
                        if (qhot.hit(Ns, p_anchor_limit, anchor_lookup, sp, sq) > 0) {
                            have_qs = true; qsN = Ns; qsp = sp; qsq = sq;
                        }
                    }

#ifndef DUMP_ALL
                    // [FULLSKIP] Batch fully covered by the ring: count it and skip the per-N
                    // loop, which would only count. Counters and witnesses are unchanged
                    // (sampling happens above). Suggested in external review, adopted in v6.1.
#ifdef N3FIRST
                    if (done == full) {   // (36) N-3 hits are not QHot hits
                        l_verified += cnt;
                        l_qhot += cnt - (u64)__builtin_popcountll(n3mask);
                        continue;
                    }
#else
                    if (done == full) { l_verified += cnt; l_qhot += cnt; continue; }
#endif
#endif
                    // ---- Cold: ascending N, identical to v4, plus propagation ----
                    // A q found for N_i is also OR-ed into later N_j (j > i) of this
                    // batch, mirroring v4 where that q enters the ring before N_i + 2.
                    for (unsigned i = 0; i < cnt; ++i) {
                        if ((done >> i) & 1ULL) {
#ifdef N3FIRST
                            ++l_verified;
                            if (!((n3mask >> i) & 1ULL)) ++l_qhot;   // (36)
#else
                            ++l_verified; ++l_qhot;
#endif
#ifdef DUMP_ALL
                            { const u64 Nd = N0 + 2ULL * i;
                              dump_rows[bi].push_back({Nd, Nd - certq[i], certq[i],
                                  (hot_from_ring >> i) & 1ULL ? (u8)0 : (u8)1}); }
#endif
                            continue;
                        }
                        const u64 N = N0 + 2ULL * i;

                        bool ok = false;
                        for (u32 pr : base_primes) {
                            if ((u64)pr > p_anchor_limit) break;
                            if ((u64)pr >= N) break;

                            u64 q = N - pr;
                            ++l_cand;
                            if (!survives_small_trial_division(q)) continue;
                            ++l_surv;

#ifdef COLDLEAN
                            // [COLDLEAN] q cannot be in the ring here: if it were, the ring walk
                            // would have tested exactly this pair (bit i of q's slice is "is p an
                            // anchor prime", true by construction) and N would be covered.
                            const bool known = false;
                            ++l_mr;
                            if (is_prime64_after_td(q)) {
#else
                            const bool known = qhot.contains_q(q);
                            if (!known) ++l_mr;
                            if (known || is_prime64(q)) {
#endif
                                ++l_verified;
                                if (!known) qhot.record(q);

                                if (!have_cs) { have_cs = true; csN = N; csp = pr; csq = q; }
                                // [B] propagate q to N_j, j > i: p_j = pr + 2(j - i)
                                if (i + 1 < cnt) {
#ifdef DUMP_ALL
                                    { u64 nb = (anchor_slice(abits, ((u64)pr - 1ULL) >> 1) << i) & full & ~done & ~(1ULL << i);
                                      while (nb) { certq[__builtin_ctzll(nb)] = q; nb &= nb - 1ULL; } }
#endif
                                    done |= (anchor_slice(abits, ((u64)pr - 1ULL) >> 1) << i) & full;
                                }
#ifdef DUMP_ALL
                                dump_rows[bi].push_back({N, (u64)pr, q, (u8)2});
#endif
                                ok = true;
                                break;
                            }
                        }

                        // ---- (6) Miss: identical to v4 ----
                        if (!ok) {
                            ++l_misses;
                            MissRecord m;
                            m.N = N;
                            u64 ep = 0, eq = 0;
                            u64 hi = std::min<u64>(ESCALATION_LIMIT, N / 2);
                            if (escalated_search(N, p_anchor_limit + 1, hi, ep, eq)) {
                                m.p = ep; m.q = eq; m.status = MISS_ANCHOR_EXHAUSTED;
                            } else {
                                m.status = MISS_UNRESOLVED;
                            }
                            block_misses[bi].push_back(m);
                        }
                    }
                }
                // (7) Record the candidates only. Bucket claiming and the
                // sample_limit cutoff are decided serially, in block order.
                if (sampling) {
                    BlockWitness& w = block_wit[bi];
                    if (have_qs) {
                        w.have_qhot = 1; w.qN = qsN; w.qp = qsp; w.qq = qsq;
                    }
                    if (have_cs) {
                        w.have_cold = 1; w.cN = csN; w.cp = csp; w.cq = csq;
                    }
                }

                } catch (...) {
                    // (10) Most plausibly bad_alloc from the miss buffer.
                    // Abort the run; the chunk is discarded below, so a
                    // partially processed block can never reach an artifact.
                    g_fatal.store(1, memory_order_relaxed);
                    g_shutdown.store(1, memory_order_relaxed);
                }
            }

#pragma omp atomic
            c_total    += l_total;
#pragma omp atomic
            c_verified += l_verified;
#pragma omp atomic
            c_qhot     += l_qhot;
#pragma omp atomic
            c_misses   += l_misses;
#pragma omp atomic
            c_probes   += l_probes;
#pragma omp atomic
            c_cand     += l_cand;
#pragma omp atomic
            c_surv     += l_surv;
#pragma omp atomic
            c_mr       += l_mr;
#ifdef N3FIRST
            g_n3_tests += l_n3_tests;   // (36)
            g_n3_hits  += l_n3_hits;
#endif
        }

        double chunk_sec = chrono::duration<double>(chrono::steady_clock::now() - t0).count();

        // (10) A fatal error inside the parallel region discards the whole
        // in-flight chunk, exactly like an interrupt: nothing from it is
        // resolved, merged or checkpointed.
        if (g_fatal.load(memory_order_relaxed)) {
            cerr << "\n*** FATAL: exception (most likely out of memory) inside the"
                    " parallel region. ***\n";
            cerr << (have_checkpoint
                     ? "The last saved checkpoint is intact; resume from it.\n"
                     : "No checkpoint was saved before the failure.\n");
            break;
        }

#ifdef DUMP_ALL
        {   // [DUMP] every certificate, in ascending N (debug builds only)
            FILE* df = fopen("goldbach_v62_dump.csv", "a");
            static const char* tag[] = {"RING", "PROP", "COLD"};
            for (size_t bi = 0; bi < blocks.size(); ++bi)
                for (const auto& r : dump_rows[bi])
                    fprintf(df, "%llu,%llu,%llu,%s\n", (unsigned long long)r.N,
                            (unsigned long long)r.p, (unsigned long long)r.q, tag[r.src]);
            fclose(df);
        }
#endif

        // ---- (7) Serial, block-ordered resolution ----
        // Everything that reaches an output file is decided here, single
        // threaded, walking blocks in ascending index (= ascending N). Bucket
        // ownership and the sample_limit cutoff therefore depend only on the
        // range and the block partition, never on thread count or scheduling.
        for (size_t bi = 0; bi < blocks.size(); ++bi) {
            for (const auto& m : block_misses[bi]) {
                chunk_misses.push_back(m);
                if (m.status == MISS_ANCHOR_EXHAUSTED)
                    cerr << "MISS(ANCHOR_EXHAUSTED) N=" << m.N
                         << " needs p=" << m.p << " q=" << m.q
                         << " (> p_anchor_limit=" << p_anchor_limit << ")\n" << flush;
                else
                    cerr << "MISS(UNRESOLVED) N=" << m.N
                         << " no starter partition found up to the escalation limit"   // (34)
                            " -- INVESTIGATE\n" << flush;
            }

            // (31) Minp rows merge in block order (= ascending N); the cap
            // is applied here, serially, so the kept set is deterministic.
            {   // (38) ascending block order + strict '>' keeps the smallest N at the max
                const MinpStat& b = bstat[bi];
                if (b.max_p > chunk_ms.max_p) { chunk_ms.max_p = b.max_p; chunk_ms.N_at_max = b.N_at_max; }
                chunk_ms.sum_p += b.sum_p; chunk_ms.count += b.count; chunk_ms.none += b.none;
            }
            for (auto& r : block_minp[bi]) {
                if ((u64)(all_minp.size() + chunk_minp.size()) >= minp_limit) break;
                chunk_minp.push_back(std::move(r));
            }

            if (sample_limit == 0) continue;
            const BlockWitness& w = block_wit[bi];

            if (w.have_qhot && w.qN >= original_start &&
                all_samples.size() + chunk_samples.size() < sample_limit) {
                u64 b = (w.qN - original_start) / SIEVE_SAMPLE_INTERVAL;
                if (b < nbuckets && !qhot_bucket_seen[(size_t)b]) {
                    qhot_bucket_seen[(size_t)b] = 1;
                    chunk_samples.push_back(to_string(w.qN) + "," + to_string(w.qp) +
                                            "," + to_string(w.qq) + ",QHOT");
                }
            }
            if (w.have_cold && w.cN >= original_start &&
                all_samples.size() + chunk_samples.size() < sample_limit) {
                u64 b = (w.cN - original_start) / SIEVE_SAMPLE_INTERVAL;
                if (b < nbuckets && !cold_bucket_seen[(size_t)b]) {
                    cold_bucket_seen[(size_t)b] = 1;
                    chunk_samples.push_back(to_string(w.cN) + "," + to_string(w.cp) +
                                            "," + to_string(w.cq) + ",COLD");
                }
            }
        }

        if (g_shutdown.load(memory_order_relaxed)) {
            cout << "\n*** Shutdown signal received. ***\n";
            cout << (have_checkpoint
                     ? "Last valid checkpoint covers up to the previous chunk.\n"
                     : "No checkpoint saved yet (interrupted during the first chunk).\n");
            cout << "Resume with: " << argv[0] << " --resume 0 0 " << block_bits << " "
                 << threads << " " << p_anchor_limit << " " << sample_limit << " "
                 << qhot_size << "\n";
            break;
        }

        for (auto& s : chunk_samples) all_samples.push_back(std::move(s));
        for (auto& r : chunk_minp)    all_minp.push_back(std::move(r));   // (31)
        if (chunk_ms.max_p > run_ms.max_p) { run_ms.max_p = chunk_ms.max_p; run_ms.N_at_max = chunk_ms.N_at_max; }   // (38)
        run_ms.sum_p += chunk_ms.sum_p; run_ms.count += chunk_ms.count; run_ms.none += chunk_ms.none;
        for (auto& m : chunk_misses) {
            all_misses.push_back(m);
            if (m.status == MISS_UNRESOLVED) ++run_unres;
        }

        run_total += c_total; run_verified += c_verified; run_qhot += c_qhot;
        run_misses += c_misses; run_probes += c_probes; run_cand += c_cand;
        run_surv += c_surv; run_mr += c_mr; run_elapsed += chunk_sec;

        CkptState s = st;
        s.last_verified  = chunk_end;
        s.total_total    = st.total_total    + run_total;
        s.total_verified = st.total_verified + run_verified;
        s.total_qhot     = st.total_qhot     + run_qhot;
        s.total_misses   = st.total_misses   + run_misses;
        s.unresolved     = st.unresolved     + run_unres;
        s.qhot_probes    = st.qhot_probes    + run_probes;
        s.cold_cand      = st.cold_cand      + run_cand;
        s.cold_surv      = st.cold_surv      + run_surv;
        s.mr_calls       = st.mr_calls       + run_mr;
        s.elapsed_sec    = st.elapsed_sec    + run_elapsed;

        // (8)/(30) Three outcomes, three truthful messages. Only a failure
        // BEFORE the rename leaves the previous checkpoint authoritative;
        // after the rename the new checkpoint is in place regardless of
        // whether its durability could be confirmed.
        const SaveResult sr = save_checkpoint(s, all_samples, all_misses,
                                              threads, block_bits,
                                              p_anchor_limit, sample_limit,
                                              qhot_size,
                                              minp_per_window, minp_limit);
        if (sr != SaveResult::FAILED_BEFORE_RENAME) {
            ++ckpt_count;
            have_checkpoint = true;
            if (sr == SaveResult::RENAMED_NOT_DURABLE)
                cerr << "WARNING: new checkpoint is in place, but the"
                        " directory sync failed -- its survival across an"
                        " immediate power loss is not confirmed.\n";
            // (24) The derived report follows the checkpoint that is now in
            // place, stamped with its last_verified.
            const SaveResult mr = write_miss_file(all_misses, original_start,
                                                  original_end,
                                                  s.last_verified,
                                                  s.total_misses);
            if (mr != SaveResult::FAILED_BEFORE_RENAME) {
                have_current_miss_report = true;
                if (mr == SaveResult::RENAMED_NOT_DURABLE)
                    cerr << "WARNING: miss report is in place but its"
                            " directory sync failed (durability"
                            " unconfirmed).\n";
            } else {
                have_current_miss_report = false;
                cerr << "WARNING: miss report save FAILED (disk full or I/O"
                        " error); " << CKPT_FILE << " remains authoritative.\n";
            }
            // (31) The minp log follows the checkpoint that is now in place.
            // Failure is tolerated: rows stay in memory and are retried at
            // the next checkpoint.
            if (minp_per_window > 0) {
                const SaveResult mnp = write_minp_file(all_minp, original_start,
                                                       original_end,
                                                       s.last_verified,
                                                       minp_per_window,
                                                       minp_limit,
                                                       p_anchor_limit);
                if (mnp == SaveResult::FAILED_BEFORE_RENAME)
                    cerr << "WARNING: minimal-partner log save FAILED; rows"
                            " kept in memory and retried at the next"
                            " checkpoint.\n";
                else if (mnp == SaveResult::RENAMED_NOT_DURABLE)
                    cerr << "WARNING: minimal-partner log is in place but its"
                            " directory sync failed (durability"
                            " unconfirmed).\n";
            }
        } else {
            cerr << "WARNING: checkpoint save FAILED before the swap (disk"
                    " full or I/O error). Previous checkpoint kept. Progress"
                    " since it is not persisted -- free space before stopping"
                    " this run.\n";
        }

        double tp  = s.elapsed_sec > 0.0 ? (double)s.total_total / s.elapsed_sec / 1e6 : 0.0;
        double pct = 100.0 * (double)(chunk_end - original_start + 2) /
                             (double)(original_end - original_start + 2);

        cout << fixed << setprecision(2)
             << "[CKPT " << ckpt_count << "] Up to " << chunk_end
             << " | " << pct << "% | " << s.total_verified << " verified | "
             << all_samples.size() << " samples | "
             << setprecision(1) << tp << " M/s | " << s.elapsed_sec << " s";
        if (c_misses > 0) cout << " | MISSES: " << c_misses;
        cout << "\n" << flush;

        if (chunk_end >= end) break;
        chunk_start = chunk_end + 2;
        } catch (const bad_alloc&) {
            g_fatal.store(1, memory_order_relaxed);
            cerr << "\n*** FATAL: out of memory in chunk setup or merge. ***\n";
            cerr << (have_checkpoint
                     ? "The last valid checkpoint is intact; resume from it.\n"
                     : "No checkpoint exists yet.\n");
            break;
        }
    }

    const u64 gt  = st.total_total    + run_total;
    const u64 gv  = st.total_verified + run_verified;
    const u64 gq  = st.total_qhot     + run_qhot;
    const u64 gm  = st.total_misses   + run_misses;
    const u64 gu  = st.unresolved     + run_unres;
    const u64 gpb = st.qhot_probes    + run_probes;
    const u64 gcc = st.cold_cand      + run_cand;
    const u64 gcs = st.cold_surv      + run_surv;
    const u64 gmr = st.mr_calls       + run_mr;
    const double ge = st.elapsed_sec  + run_elapsed;

    // (24) No final report rewrite: the report was last written alongside
    // the last SUCCESSFUL checkpoint and must not be rewritten from
    // in-memory state that may exceed what that checkpoint holds.

    const double dt = gt ? (double)gt : 1.0;
    cout << fixed << setprecision(6);
    cout << "\n===== FINAL SUMMARY =====\n";
    cout << "Full range:        [" << original_start << " .. " << original_end << "]\n";
    cout << "Coverage:          " << gv << " / " << gt << "\n";
#ifdef N3FIRST
    {   // (36) this process only: run_total, not the cumulative gt
        const double rt = run_total ? (double)run_total : 1.0;
        cout << "N-3 first (this run): " << g_n3_hits.load() << " hits / "
             << g_n3_tests.load() << " tests ("
             << 100.0 * (double)g_n3_hits.load() / rt << "% of all N)\n";
    }
#endif
    cout << "Misses:            " << gm
         << " (unresolved: " << gu
         << ", anchor-exhausted: " << (gm - gu) << ")\n";
    // (4) Four independent cost counters instead of one ambiguous tries/N.
    cout << "QHot hit rate:     " << 100.0 * (double)gq / dt << "%\n";
    cout << "QHot probes/N:     " << (double)gpb / dt << "\n";
    cout << "Cold cand/N:       " << (double)gcc / dt << "   (anchor primes considered)\n";
    cout << "Cold surv/N:       " << (double)gcs / dt << "   (survived trial division)\n";
    cout << "MR calls/N:        " << (double)gmr / dt << "   (Miller-Rabin invocations)\n";
    cout << "Total time:        " << ge << " s\n";
#ifdef SIMD_CHECK
    cout << "SIMD checks passed: " << simd_checks.load() << "\n";
#endif
    cout << "Throughput:        " << (ge > 0.0 ? (double)gt / ge / 1e6 : 0.0) << " M/s\n";
    cout << "Samples:           " << all_samples.size() << "\n";
    if (minp_per_window > 0)
        cout << "Minp rows:         " << all_minp.size()
             << " (" << MINP_FILE << ")\n";
    if (minp_per_window > 0) {
        // (38) Window summary: on screen always; appended to MINPSTATS_FILE only
        // for a run that finished its range (no interrupt, no fatal error).
        const double mean_p = run_ms.count ? (double)run_ms.sum_p / (double)run_ms.count : 0.0;
        cout << "Minp window stats: " << run_ms.count << " N with a partner, "
             << run_ms.none << " without; max minimal p = " << run_ms.max_p
             << " at N = " << run_ms.N_at_max << "; mean minimal p = "
             << setprecision(3) << mean_p << setprecision(6) << "\n";
        const bool finished = !g_shutdown.load(memory_order_relaxed) && !g_fatal.load(memory_order_relaxed);
        if (finished) {
            bool header = false;
            { ifstream probe(MINPSTATS_FILE); header = !probe.good(); }
            ofstream sf(MINPSTATS_FILE, ios::app);
            if (header)
                sf << "prog_version,range_start,range_end,resumed,minp_per_window,"
                      "p_anchor_limit,n_with_partner,n_without_partner,max_min_p,"
                      "N_at_max,mean_min_p\n";
            sf << PROG_VERSION << "," << original_start << "," << original_end << ","
               << (resuming ? 1 : 0) << "," << minp_per_window << "," << p_anchor_limit << ","
               << run_ms.count << "," << run_ms.none << "," << run_ms.max_p << ","
               << run_ms.N_at_max << "," << fixed << setprecision(4) << mean_p << "\n";
            sf.flush();
            if (sf.good()) cout << "Minp stats row:    appended to " << MINPSTATS_FILE << "\n";
            else cerr << "WARNING: could not append to " << MINPSTATS_FILE << "\n";
        }
    }
    cout << "Checkpoints:       " << ckpt_count << "\n";
    cout << "Output:            " << CKPT_FILE << "\n";
    if (have_current_miss_report)
        cout << "Miss report:       " << MISS_FILE << "\n";
    else if (have_checkpoint)
        cout << "Miss report:       unavailable or not rewritten this run;"
                " authoritative records are in " << CKPT_FILE << "\n";

    if (g_fatal.load(memory_order_relaxed)) return 5;  // exception in parallel region
    if (gu != 0) return 4;  // no partition found -- investigate
    if (gm != 0) return 2;  // anchor limit too small for some N
    if (g_shutdown.load(memory_order_relaxed)) return 3;
    return 0;
}
