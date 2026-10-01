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

#include <algorithm>
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

static const char* PROG_VERSION = "goldbach_v5_sorted";
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

static constexpr u64 CKPT_INTERVAL         = 100000000000ULL; // even numbers per checkpoint chunk
static constexpr u64 SIEVE_SAMPLE_INTERVAL = 10000000000ULL;  // one witness per this much N-distance
static constexpr u64 ESCALATION_LIMIT      = 100000000ULL;    // (6) how far past the anchor limit a miss is chased

// (30) Three-state save outcome. A plain bool cannot distinguish "the old
// artifact is still authoritative" from "the new artifact is already in
// place but its survival across a power loss is unconfirmed" -- and the
// caller's diagnostics must not claim the former when the latter happened.
enum class SaveResult {
    SUCCESS,              // written, renamed, directory sync confirmed
    FAILED_BEFORE_RENAME, // previous artifact untouched and authoritative
    RENAMED_NOT_DURABLE   // new artifact visible; crash durability unconfirmed
};

static const char* CKPT_FILE = "goldbach_v5sort_checkpoint.csv";   // [B] never the v4 campaign file
static const char* MISS_FILE = "goldbach_v5sort_misses.csv";        // [B]

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

// ------------------ Base primes (anchor table) ------------------

static vector<u32> base_primes;

static void build_base_primes(u32 limit) {
    vector<bool> is_prime((size_t)limit + 1, true);
    is_prime[0] = false;
    if (limit >= 1) is_prime[1] = false;

    for (u32 i = 2; (u64)i * i <= limit; ++i) {
        if (is_prime[i]) {
            for (u64 j = (u64)i * i; j <= limit; j += i)
                is_prime[(size_t)j] = false;
        }
    }

    base_primes.clear();
    // p = 2 is deliberately excluded: N is even, so q = N - 2 is even and > 2.
    // (15) u64 counter: with a u32 counter and limit = 2^32-1 (which the
    // validation range [3, 2^32) permits), i += 2 wraps and loops forever.
    for (u64 i = 3; i <= limit; i += 2)
        if (is_prime[(size_t)i]) base_primes.push_back((u32)i);
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
static bool escalated_search(u64 N, u64 p_from, u64 p_to, u64& p_out, u64& q_out) {
    u64 p = p_from | 1ULL;
    if (p < 3) p = 3;
    for (; p <= p_to && p < N; p += 2) {
        if (!survives_small_trial_division(p)) continue;
        if (!is_prime64(p)) continue;
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
    f << "# status ANCHOR_EXHAUSTED: partition exists but requires p > p_anchor_limit\n";
    f << "# status UNRESOLVED: no partition found up to the escalation limit "
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

    // (17) Miss records read back from # miss= rows in this file.
    vector<MissRecord> misses;
    // (18) True when the checkpoint declares cfg_version (current format);
    // such checkpoints are validated strictly for field completeness.
    bool has_version = false;
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
                 F_QS, F_COUNT };
    static const char* FIELD_NAMES[F_COUNT] = {
        "original_start", "original_end", "last_verified", "total_verified",
        "total_total", "qhot_probes", "cold_candidates", "cold_survivors",
        "mr_calls", "total_qhot", "total_misses", "unresolved",
        "elapsed_sec", "cfg_version", "cfg_format", "cfg_block_bits",
        "cfg_p_anchor_limit", "cfg_sample_limit", "cfg_qhot_size" };
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
            else if (key == "cfg_version")     { ck.has_version    = true;        mark(F_VER); }
            else if (key == "cfg_format")      { ck.cfg_format     = strict_int(val); mark(F_FMT); }
            else if (key == "cfg_block_bits")     { ck.cfg_block_bits     = strict_int(val); mark(F_BB); }
            else if (key == "cfg_p_anchor_limit") { ck.cfg_p_anchor_limit = strict_u64(val); mark(F_PAL); }
            else if (key == "cfg_sample_limit")   { ck.cfg_sample_limit   = strict_u64(val); mark(F_SL); }
            else if (key == "cfg_qhot_size")      { ck.cfg_qhot_size      = strict_int(val); mark(F_QS); }
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
        const u32 required_mask = (seen_mask & (1u << F_FMT))
                                ? FULL_MASK
                                : (FULL_MASK & ~(1u << F_FMT));
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
                            u64 sample_limit, int qhot_size) {
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
    f << "# source tags: QHOT = accepted from the QHot ring;"
         " COLD = accepted on the ordered anchor scan (Miller-Rabin).\n";
    f << "# SIEVE is a legacy alias for COLD and is accepted on resume;"
         " no build ever used a segmented sieve to accept a witness.\n";
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
            " be fully restored\n\n";
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
    vector<string> pos_args;
    for (int i = 1; i < argc; ++i) {
        string a = argv[i];
        if (a == "--resume")      resuming = true;
        else if (a == "--allow-miss-loss") allow_miss_loss = true;
        else if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
        else                      pos_args.push_back(a);
    }

    // (15) At most 7 positionals: a stray extra argument is a typo, not
    // something to ignore silently.
    if (pos_args.size() < 4 || pos_args.size() > 7) { usage(argv[0]); return 1; }

    u64 start, end;
    int block_bits, threads, qhot_size;
    u64 p_anchor_limit, sample_limit;
    try {
        start          = parse_u64(pos_args[0], "start");
        end            = parse_u64(pos_args[1], "end");
        block_bits     = parse_int(pos_args[2], "block_bits");
        threads        = parse_int(pos_args[3], "threads");
        p_anchor_limit = pos_args.size() > 4 ? parse_u64(pos_args[4], "p_anchor_limit") : 1000000ULL;
        sample_limit   = pos_args.size() > 5 ? parse_u64(pos_args[5], "sample_limit")   : 10000ULL;
        qhot_size      = pos_args.size() > 6 ? parse_int(pos_args[6], "qhot_size")      : 1024;
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

    CkptState st;
    vector<string> all_samples;
    vector<MissRecord> all_misses;

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

    cout << " done. Base primes: " << base_primes.size()
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

        // (7) Per-block output. Each element is written by exactly one thread
        // (the owner of that block), so no synchronisation is needed and the
        // contents do not depend on scheduling.
        vector<BlockWitness> block_wit(blocks.size());
        vector<vector<MissRecord>> block_misses(blocks.size());
#ifdef DUMP_ALL
        struct DumpRow { u64 N, p, q; u8 src; };   // src: 0 ring, 1 propagated, 2 cold
        vector<vector<DumpRow>> dump_rows(blocks.size());
#endif

        auto t0 = chrono::steady_clock::now();

#pragma omp parallel
        {
            QHot qhot;
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

                    // ---- Hot: walk the ring newest-first until the batch is covered ----
                    u64 done = 0;
#ifdef DUMP_ALL
                    u64 certq[64] = {0};   // [DUMP] q that first certified each bit
#endif
                    if (qhot.count > 0) {
                        for (int t = 0; t < qhot.count; ++t) {
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
                    const u64 hot_from_ring = done;

                    // Sample: lowest ring-certified N, witness re-derived with the
                    // unchanged v4 per-N probe while the ring is still as walked.
                    if (hot_from_ring && !have_qs) {
                        const u64 Ns = N0 + 2ULL * (u64)__builtin_ctzll(hot_from_ring);
                        u64 sp = 0, sq = 0;
                        if (qhot.hit(Ns, p_anchor_limit, anchor_lookup, sp, sq) > 0) {
                            have_qs = true; qsN = Ns; qsp = sp; qsq = sq;
                        }
                    }

                    // ---- Cold: ascending N, identical to v4, plus propagation ----
                    // A q found for N_i is also OR-ed into later N_j (j > i) of this
                    // batch, mirroring v4 where that q enters the ring before N_i + 2.
                    for (unsigned i = 0; i < cnt; ++i) {
                        if ((done >> i) & 1ULL) {
                            ++l_verified; ++l_qhot;
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

                            const bool known = qhot.contains_q(q);
                            if (!known) ++l_mr;
                            if (known || is_prime64(q)) {
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
            FILE* df = fopen("goldbach_v5sort_dump.csv", "a");
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
                         << " no partition found up to the escalation limit"
                            " -- INVESTIGATE\n" << flush;
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
                                              qhot_size);
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
    cout << "Throughput:        " << (ge > 0.0 ? (double)gt / ge / 1e6 : 0.0) << " M/s\n";
    cout << "Samples:           " << all_samples.size() << "\n";
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
