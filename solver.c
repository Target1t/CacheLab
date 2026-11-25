/*
 * solver.c - Cache Lab Bonus: reverse-engineer a mystery cache.
 *
 * solve() is called once per cache, in a fresh process. Only the cache_*
 * probes and report() are used. The four parameters fall out in this order:
 *
 *   b  the first power-of-two offset that misses ends the first block;
 *   E  blocks that all map to one set fill it up;
 *   s  the smallest power-of-two block stride whose blocks still collide;
 *   policy  touch a line, add one more line, see which line disappeared.
 *
 * Every test that depends on a random eviction is repeated: a deterministic
 * policy gives the same answer every time, RANDOM eventually gives another.
 */
#include "solver.h"

/* Repetitions for a probe whose answer RANDOM can get right by luck. */
#define TRIALS 256
/* Repetitions of the replacement-policy experiment. */
#define POLICY_TRIALS 64

/*
 * block_addr - address of the j-th block when consecutive blocks are
 *     2^shift blocks apart. A shift of 9 (512 blocks, more than the 512
 *     sets the largest cache has) makes every block land in set 0.
 */
static uint32_t block_addr(int b, int shift, int j)
{
    return (uint32_t)j << (b + shift);
}

/*
 * any_miss - Load n blocks of one set, then read the first one again.
 *     Returns 1 if the first block was evicted at least once, 0 if it was
 *     still there on every trial.
 */
static int any_miss(int b, int shift, int n, int trials)
{
    int i, j;

    for (i = 0; i < trials; i++) {
        cache_flush();
        for (j = 0; j < n; j++)
            cache_load(block_addr(b, shift, j));
        if (!cache_load(block_addr(b, shift, 0)))
            return 1;
    }
    return 0;
}

/*
 * find_b - The first load of a block misses, everything below its size hits.
 */
static int find_b(void)
{
    int k;

    cache_flush();
    cache_load(0);
    for (k = 1; k <= 7; k++) {
        if (!cache_load((uint32_t)1 << k))
            return k;
    }
    return 7;
}

/*
 * find_E - Blocks that all map to the same set fit E at a time, so the
 *     first n that pushes the oldest block out is E + 1.
 */
static int find_E(int b)
{
    int n, E = 1;

    for (n = 1; n <= 16; n++) {
        if (any_miss(b, 9, n, TRIALS))
            break;
        E = n;
    }
    return E;
}

/*
 * find_s - Blocks 2^(b+t) apart land in the same set once 2^t is a multiple
 *     of the number of sets, i.e. once t >= s. For t < s they spread over
 *     distinct sets and E + 1 of them never evict the first one.
 */
static int find_s(int b, int E)
{
    int t;

    for (t = 1; t <= 9; t++) {
        if (any_miss(b, t, E + 1, TRIALS))
            return t;
    }
    return 9;
}

/*
 * find_policy - Fill one set, touch its first line (LRU moves it, FIFO does
 *     not), then add one more line. LRU drops the second line, FIFO drops
 *     the first one, RANDOM drops either.
 */
static int find_policy(int b, int E)
{
    int i, j, missing;
    int lru_ok = 1, fifo_ok = 1;

    if (E == 1)
        return POLICY_LRU;   /* not graded for a direct mapped cache */

    for (i = 0; i < POLICY_TRIALS; i++) {
        cache_flush();
        for (j = 0; j < E; j++)
            cache_load(block_addr(b, 9, j));
        cache_load(block_addr(b, 9, 0));
        cache_load(block_addr(b, 9, E));
        missing = -1;
        for (j = 0; j < E; j++) {
            if (!cache_load(block_addr(b, 9, j))) {
                missing = j;
                break;
            }
        }
        if (missing != 1)
            lru_ok = 0;
        if (missing != 0)
            fifo_ok = 0;
        if (!lru_ok && !fifo_ok)
            break;
    }
    if (lru_ok)
        return POLICY_LRU;
    if (fifo_ok)
        return POLICY_FIFO;
    return POLICY_RANDOM;
}

/*
 * find_write_allocate - A store miss either brings the block in (the next
 *     load hits) or leaves the cache alone (the next load misses).
 */
static int find_write_allocate(void)
{
    cache_flush();
    cache_store(0);
    return cache_load(0);
}

void solve(void)
{
    int b, E, s, policy, write_allocate;

    b = find_b();
    E = find_E(b);
    s = find_s(b, E);
    policy = find_policy(b, E);
    write_allocate = find_write_allocate();

    report(b, s, E, policy, write_allocate);
}
