sox (+bash) magic for extrating a bit (first 20 sec here) of an audio file..
------
for n in $(seq -f "%02g" 1 14); do sox ../20120410_063059_${n}.wav foo_$n.wav trim 0 20 ; done

and to segment (into 1 sec here) those files ....

for n in $(seq -f "%02g" 1 14); do sox foo_${n}.wav segmented/foo_${n}_.wav trim 0 1 : newfile : restart ; done
