# NArrays

A collection of specialized array types and data structures for memory-efficient storage and manipulation of large datasets. Used as a dependency of NeoKCT.

## Data Structures

### DeltaArray

A delta-encoded, sorted array of unsigned integers. Instead of storing each value in full, it stores the difference from the previous element using a smaller type. Checkpoints at regular intervals make random access practical.

```julia
# Default: UInt64 values, UInt32 deltas, checkpoint every 256 elements
a = DeltaArray([1, 5, 6, 300, 301, 1000])

# Explicit types
a = DeltaArray{UInt64, UInt16}([1, 5, 6, 300], 128)

# Build incrementally, then encode from a sorted flat vector
a = DeltaArray()
encode!(a, sorted_vec)

# Input not sorted yet
a = DeltaArray(unsorted_vec; presorted=false)
```

**Memory layout:** Values are split into checkpoints (stored at full `C` width) and deltas (stored at `D` width). A checkpoint is inserted every `checkpoint_interval` elements, and also whenever a delta would overflow `typemax(D)`. The `regular_cp_idx` vector maps each chunk to its checkpoint position, allowing random access without scanning from the start.

**Complexity:**
- Indexing `a[i]`: O(`checkpoint_interval`), decodes from the nearest checkpoint
- Iteration: O(1) amortised, state carries the current checkpoint and value forward
- `searchfirst(a, target, lo, hi)`: linear scan in `[lo, hi]` starting from the checkpoint covering `lo`; exits early once the decoded value exceeds `target`

**API:**

| Function | Description |
|---|---|
| `DeltaArray(arr; presorted=false)` | Build from a vector (sorts if needed) |
| `DeltaArray{C,D}(interval)` | Empty array with explicit value/delta types |
| `encode!(a, arr)` | Re-encode in-place from a new sorted flat vector |
| `searchfirst(a, target, lo, hi)` | First index in range whose value equals `target`, or `0` |
| `a[i]` | Random access |
| `length(a)`, `size(a)` | Standard `AbstractVector` interface |

---

### PackedArray

A bit-packed array where multiple small values are packed into each unsigned-integer word. Each word is an independent container; `arr[i]` returns a `Vector{T}` of all values packed in the `i`-th word. A companion `BitVector` (the bitmap) tracks element boundaries within each word.

```julia
arr = PackedArray{UInt8, UInt64}()   # pack UInt8 values into UInt64 words

# Push a value into a new word
push!(arr, UInt8(42))

# Push a value into an existing word (spills to a new word if no space)
push!(arr, UInt8(7), word_id)

# Allocate an empty word and get its index
id = new_word!(arr)

# Read all values from word i
vals = arr[i]  # Vector{UInt8}

# Deduplicate words; returns (deduped_arr, perm, global_perm)
deduped, perm, global_perm = permdedup(arr)
```

**Memory layout:** Each word is a `W`-bit integer. The bitmap has one bit per bit-position per word (`length(words) * bitsizeof(W)` bits total). A `1` in the bitmap marks the last bit of a packed value, so the bitmap encodes both count and size of all values in a word without any per-element header.

**API:**

| Function | Description |
|---|---|
| `PackedArray{T,W}()` | Empty array, packing type `T` into words of type `W` |
| `push!(arr, elem)` | Pack `elem` into a new word |
| `push!(arr, elem, word_id)` | Pack `elem` into word `word_id` (spills if full) |
| `new_word!(arr)` | Allocate a new empty word; returns its index |
| `arr[i]` | All values in word `i` as `Vector{T}` |
| `prepack(elem, W)` | Compute packed representation and bit-width (helper) |
| `get_val(arr, word_id, elem_id)` | `elem_id`-th value in word `word_id` (as `W`) |
| `bitsizeof(e)` | Bit width of a value or type (`sizeof * 8`) |
| `permdedup(arr)` | Deduplicate words; returns `(deduped, perm, global_perm)` |
| `word_bitmap_slice(word_id, W)` | Bitmap index range for word `word_id` |
| `word_last_set(arr, word_id)` | Last occupied bit position in word `word_id` |
| `update_word_bitmap!(arr, word_id, bitsize)` | Mark a value boundary in the bitmap |

---

## Utilities

### Parallel Sort (`psort!`, `psortperm!`, `psortperm`)

A parallel merge sort using Julia's `Threads.@spawn`. Falls back to serial `MergeSort` for chunks smaller than 100 000 elements.

```julia
psort!(v)                  # sort v in-place
perm = psortperm(data)     # return a permutation that sorts data
psortperm!(perm, data)     # fill perm in-place
```

Adapted from the [Julia multithreading blog post](https://julialang.org/blog/2019/07/multithreading/).

---

## Installation

This package is not registered. Add it by path:

```julia
using Pkg
Pkg.add(path="/path/to/NArrays")
```

Or from a Julia environment that already depends on it, it is resolved automatically alongside NeoKCT.
