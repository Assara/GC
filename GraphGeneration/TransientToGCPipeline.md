# Transient-to-GC conversion

`TransientToGCPipeline<V,E>` converts one new transient-stage file into one GC
file. It is single-threaded and separate from the transient generation and
legacy finalization pipelines.

```sh
make -f GraphGeneration/Makefile.pipeline gc-converter GC_VERTICES=4 GC_EDGES=6
./build/transient_to_gc_V4_E6 input/transient_V4_E6.gcg output/gc_V4_E6.gcg
make -f GraphGeneration/Makefile.pipeline test-gc-converter
```

```cpp
auto counts = GraphGeneration::TransientToGCPipeline<4, 6>{}.run(input, output);
// counts.total, counts.odd_gc, counts.even_gc
```

## Conversion rules

1. Read the mmap transient file and discard each graph containing a vertex of
   valence below 3. Both 1- and 2-valent transient vertices are excluded.
2. Canonicalize each retained graph with the existing final graph canonicalizer.
   Its automorphism parity check determines survival in both conventions:
   - **oddGC:** vertex-permutation sign multiplied by `(-1)^reversed_edges`
     (odd vertices, even edges, odd edge reversals; `C=1, D=1`). Permuting
     whole edges exchanges pairs of half-edges and contributes no extra sign.
   - **evenGC:** odd edges, even vertices (`C=0, D=1`).
3. Write the canonical graph and the two Boolean survival flags. No orientation
   sign or automorphism count is stored.

Input stage files are already deduplicated; conversion does not deduplicate
again. Graphs that vanish in both conventions remain in the output. `total`
counts all records after the valence filter; `odd_gc` and `even_gc` independently
count nonzero records in their respective complexes. Their sum need not equal
`total`.

The converter counts retained records in a first pass, then canonicalizes and
writes them directly into an exactly sized mmap file in a second pass. It does
not buffer the graph collection. An empty result still writes a valid header.
The input is unchanged, existing outputs are refused, and a temporary output
is renamed after writing and synchronizing the mapping.

## Byte format

This format has its own magic identifier and does not change old file formats.
The header is 64 bytes. All integer fields below are unsigned 64-bit values in
little-endian byte order.

| Offset | Field |
| --- | --- |
| 0 | Eight magic bytes: `GCGC2\r\n\0` |
| 8 | Format version: 2 |
| 16 | Vertex count |
| 24 | Edge count |
| 32 | Total record count |
| 40 | Nonzero oddGC count |
| 48 | Nonzero evenGC count |
| 56 | Record size: `2 * E + 1` |

Each record contains `2 * E` endpoint bytes followed by one flags byte: bit 0
means oddGC is nonzero; bit 1 means evenGC is nonzero. Other bits must be zero.
There is no padding between records. `MappedGCGraphReader<G>` exposes each
record and the persisted counts, and checks the format, dimensions, sizes,
count bounds, and flag bits.

Version 1 used vertex parity alone for oddGC and is rejected by the reader.
Regenerate GC files from their transient inputs. The legacy finalization
pipeline retains its original convention; the new converter explicitly enables
edge-reversal signs in the shared canonicalizer.

The pipeline has no explicit exception handling. Vertex-label preconditions
use assertions; an existing output aborts the process, including in release
builds. File operations retain the file layer's error handling.
