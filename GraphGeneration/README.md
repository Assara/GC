# Graph generation

The current pipeline generates triangle-covered seeds, then splits vertices at
a fixed loop number. It offers a shared-bucket splitter and a separate valence-group
splitter that keeps one target group in memory and reads its parents from disk.

See [TriangleSeeds.md](TriangleSeeds.md) for commands, filtering rules, output
formats, validation, and measured runtime and memory usage.

```sh
make -f GraphGeneration/Makefile.pipeline triangle-seeds TRIANGLE_MAX_LOOP=9
./build/triangle_seeds_unrestricted_L9 output/triangle_seeds_L9_new
make -f GraphGeneration/Makefile.pipeline triangle-splits-valence TRIANGLE_SPLIT_LOOP=9
OMP_NUM_THREADS=8 ./build/triangle_splits_valence_L9 \
  output/triangle_seeds_L9_new output/triangle_splits_valence_L9_new
make -f GraphGeneration/Makefile.pipeline test
```

Use fresh output directories. Generated graphs and benchmark files live under
`output/`; binaries live under `build/`.
