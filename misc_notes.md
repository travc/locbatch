# miscellaneous notes

sox +bash magic
------
Extrat a bit (first 20 sec here) of an audio file.
```bash
for n in $(seq -f "%02g" 1 14); do sox ../20120410_063059_${n}.wav foo_$n.wav trim 0 20 ; done
```

Segment/split those files into parts (each 1 sec long here).
```bash
for n in $(seq -f "%02g" 1 14); do sox foo_${n}.wav segmented/foo_${n}_.wav trim 0 1 : newfile : restart ; done
```
