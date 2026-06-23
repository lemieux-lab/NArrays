module NArrays

using ProgressMeter

include("DeltaArray.jl")
include("PackedArray.jl")

export DeltaArray, encode!, searchfirst
export PackedArray, bitsizeof, prepack, get_val, new_word!, permdedup, update_word_bitmap!, word_bitmap_slice, word_last_set
export psortperm, psort!

end # module NArrays
