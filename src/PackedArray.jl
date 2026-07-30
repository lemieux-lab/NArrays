"""
    struct PackedArray{T, W <: Unsigned} <: AbstractVector{Vector{T}}
        words::Vector{W}
        bitmap::BitVector
    end

A bitpacked array of type T, using words of size W. Keeps a bitmap to remember element
boundaries within each word. Implements the AbstractVector interface, where each element
is a `Vector{T}` of all values packed into the corresponding word.
"""
struct PackedArray{T, W <: Unsigned} <: AbstractVector{Vector{T}}
    words::Vector{W}
    bitmap::BitVector
end

bitsizeof(e) = sizeof(e) * 8

PackedArray{T, W}() where {T, W<:Unsigned} = PackedArray{T, W}(W[0], BitVector(falses(bitsizeof(W))))
PackedArray{T}() where {T} = PackedArray{T, UInt64}(UInt64[0], BitVector(falses(64)))

## AbstractVector interface ##

Base.size(arr::PackedArray) = (length(arr.words),)
Base.IndexStyle(::Type{<:PackedArray}) = Base.IndexLinear()

# Return all values packed in word i as a Vector{T}
#
# Single pass over the bitmap window: get_val below re-derives (start, finish) from
# scratch for one element at a time, which made the naive "call get_val n times" loop
# O(n*W) per word (each call rescans up to W bits and recopies the bitmap slice).
# Walking the view once and extracting every value as its terminating bit is found
# drops that to O(W), independent of how many values are packed into the word.
function Base.getindex(arr::PackedArray{T, W}, word_id::Integer) where {T, W<:Unsigned}
    word = arr.words[word_id]
    bitmap_word = @view arr.bitmap[word_bitmap_slice(word_id, W)]
    assembly = T[]
    start = 0
    for (pos, bit) in enumerate(bitmap_word)
        bit || continue
        finish = pos
        push!(assembly, T(word << start >> (start + (bitsizeof(W) - finish))))
        start = finish
    end
    return assembly
end

## Internal helpers ##

# Computes the packed version of elem, returns (packed::W, bitsize::Int)
function prepack(elem::T, ::Type{W}=UInt64) where {T, W<:Unsigned}
    !isbits(elem) && throw("Type $T is not bits type and cannot be packed")
    pack = W(elem & typemax(W)) # /!\ Will round down to typemax, careful
    pack == 0 && return pack, 1
    return pack, bitsizeof(W) - leading_zeros(pack)
end

# Bitmap slice covering the word_id-th word
word_bitmap_slice(word_id::Integer, word_size::DataType=UInt64) =
    (word_id-1)*bitsizeof(word_size)+1 : (word_id-1)*bitsizeof(word_size)+bitsizeof(word_size)

# Last set bit position within the word_id-th word's bitmap slice (0 if none)
function word_last_set(arr::PackedArray{T, W}, word_id::Integer) where {T, W<:Unsigned}
    last_set = findlast(@view arr.bitmap[word_bitmap_slice(word_id, W)])
    return !isnothing(last_set) ? last_set : 0
end

function update_word_bitmap!(arr::PackedArray{T, W}, word_id::Integer, bitsize::Integer) where {T, W<:Unsigned}
    arr.bitmap[word_bitmap_slice(word_id, W)[1]+word_last_set(arr, word_id)+bitsize-1] = true
end

function update_word_bitmap!(arr::PackedArray{T, W}, word_id::Integer, bitsize::Integer, last_set::Integer) where {T, W<:Unsigned}
    arr.bitmap[word_bitmap_slice(word_id, W)[1]+last_set+bitsize-1] = true
end

# Return the elem_id-th value packed in word word_id (as W)
function get_val(arr::PackedArray{T, W}, word_id::Integer, elem_id::Integer) where {T, W<:Unsigned}
    word = arr.words[word_id]
    bitmap_word = @view arr.bitmap[word_bitmap_slice(word_id, W)]
    sum(bitmap_word) < elem_id && throw("Out of Bound")
    start = 0; finish = 0; counter = 0; ones = 0
    for bit in bitmap_word
        counter += 1
        ones == elem_id && break
        (start, finish, ones) = bit ? (finish, counter, ones+1) : (start, finish, ones)
    end
    return word << start >> (start + (bitsizeof(W)-finish))
end

## Push interface ##

# Pack elem into word at word_id, spilling to a new word if there is no space. Returns arr.
function Base.push!(arr::PackedArray{T1, W}, elem::T2, word_id::Integer) where {T1, T2, W<:Unsigned}
    elem, bitsize = prepack(elem, W)
    last_set = word_last_set(arr, word_id)
    if bitsizeof(W)-last_set >= bitsize
        arr.words[word_id] |= (elem<<(bitsizeof(W)-last_set-bitsize))
        update_word_bitmap!(arr, word_id, bitsize, last_set)
    else
        new_id = new_word!(arr)
        push!(arr, elem, new_id)
    end
    return arr
end

# Pack elem into a brand-new word. Returns arr.
function Base.push!(arr::PackedArray{T, W}, elem) where {T, W<:Unsigned}
    elem, _ = prepack(elem, W)
    new_id = new_word!(arr)
    push!(arr, elem, new_id)
    return arr
end

# Pack elem into word_id given its current fill pointer last_set (bits already used in
# that word), spilling to a new word if there is no space. Returns (word_id, new_last_set)
# for the word actually written to.
#
# Base.push!(arr, elem, word_id) re-derives last_set via word_last_set on every call, an
# O(W) findlast rescan of the word's bitmap. A caller building a word sequentially (as
# repack does, one word per source element/k-mer) already knows last_set from the previous
# call, so this skips the rescan entirely: O(1) amortized per push instead of O(W). At W=256
# and up, with billions of values pushed across a full repack, the rescan dominates runtime
# and lets the job run long enough for peak memory to balloon well past what the final
# structure needs.
function push_at!(arr::PackedArray{T1, W}, elem::T2, word_id::Integer, last_set::Integer) where {T1, T2, W<:Unsigned}
    elem, bitsize = prepack(elem, W)
    if bitsizeof(W) - last_set >= bitsize
        arr.words[word_id] |= (elem << (bitsizeof(W) - last_set - bitsize))
        update_word_bitmap!(arr, word_id, bitsize, last_set)
        return word_id, last_set + bitsize
    else
        new_id = new_word!(arr)
        return push_at!(arr, elem, new_id, 0)
    end
end

# Allocate a new empty word and return its index.
#
# resize! + fill! extends arr.bitmap's underlying chunk storage in place; the previous
# `append!(arr.bitmap, falses(bitsizeof(W)))` allocated a whole new W-bit BitVector just to
# copy it in, which at W=256/512 and one call per source element (hundreds of millions to
# billions in a full repack) was a large, avoidable source of transient garbage.
function new_word!(arr::PackedArray{T, W}) where {T, W<:Unsigned}
    push!(arr.words, W(0))
    old_len = length(arr.bitmap)
    resize!(arr.bitmap, old_len + bitsizeof(W))
    @view(arr.bitmap[old_len+1:end]) .= false
    return length(arr.words)
end

## Repacking ##

# Repacks arr into a PackedArray using word type W2, splitting any word whose
# values no longer fit into as many new words as needed (in order).
# Returns (repacked_arr, id_map), where id_map[i] lists the new word id(s)
# that word i's values were split across.
function repack(arr::PackedArray{T, W1}, ::Type{W2}) where {T, W1<:Unsigned, W2<:Unsigned}
    new_arr = PackedArray{T, W2}()
    id_map = Vector{Vector{UInt32}}(undef, length(arr.words))

    @showprogress "Repacking Words..." for i in eachindex(arr.words)
        wid = new_word!(new_arr)
        ids = UInt32[wid]
        last_set = 0
        for v in arr[i]
            new_wid, last_set = push_at!(new_arr, v, wid, last_set)
            if new_wid != wid
                wid = new_wid
                push!(ids, UInt32(wid))
            end
        end
        id_map[i] = ids
    end

    return new_arr, id_map
end

## Deduplication ##

# Deduplicates words; returns (deduped_arr, perm, global_perm).
# perm[i] = original word index of the i-th unique word.
# global_perm[i] = index of word i in the deduped array.
function permdedup(arr::PackedArray{T, W}) where {T, W<:Unsigned}
    uniq = Vector{W}()
    indexof = Dict{Pair{W, SubArray}, UInt32}()
    perm = Vector{UInt32}()
    global_perm = Vector{UInt32}(undef, length(arr.words))

    @inbounds @showprogress "Deduping Count Words..." for i in eachindex(arr.words)
        word = arr.words[i]
        bits = @view arr.bitmap[word_bitmap_slice(i, W)]
        j = get(indexof, word=>bits, 0)
        if j == 0
            push!(uniq, word)
            j = length(uniq)
            indexof[word=>bits] = UInt32(j)
            push!(perm, UInt32(i))
        end
        global_perm[i] = j
    end

    bitmap_slices = word_bitmap_slice.(perm, W)
    deduped_arr = PackedArray{T, W}(uniq, arr.bitmap[vcat(bitmap_slices...)])
    return deduped_arr, perm, global_perm
end
