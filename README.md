# Locbatch #
>Locbatch is a program to localize sounds recorded by arrays of microphones.  
> Specifically, it is being used (I hope) to localize bird calls recorded by arrays of microphones 
> distributed over a relatively wide area (10s to a few 100s of meters) which are all 
> plugged into a single multi-channel recorder.

## Table of Contents ##
* [Introduction](#introduction)
* [Installing](#installing)
* [Usage](#usage)
  * [Input Files](#input-files)
    * [Sound files](#sound-files)
    * [Microphone location file](#microphone-location-file)
    * [Annotation file](#annotation-file)
    * [Configuration file](#configuration-file)
  * [Running](#running)
  * [Interpreting Results](#interpreting-results)
* [Citation](#citation)
* [Acknowlegements](#acknowlegements)

## Introduction ##
Locbatch implements the "Correlation-Envelope-Sum" (CES) method which is described in detail in section 2.5 of [my PhD dissertation](http://taylor0.biology.ucla.edu/~travc/dissertation.pdf) and breifly in section F of a JASA paper ([journal](http://scitation.aip.org/content/asa/journal/jasa/128/1/10.1121/1.3425729), [preprint](http://taylor0.biology.ucla.edu/bibliography/pdf/TravisJASA10.pdf)).  
The primary benefit of CES over simple Correlation-Sum (aka. [Accumulated Correlaiton](http://www.ces.clemson.edu/~stb/research/acousticloc/)) is that CES produces a much smoother likelihood space so searching for the maximum is a lot simpler, faster, and robust (likely to find the global max instead of a local one).  In theory, it trades-off some precision, but I find that errors from surveying microphone positions and the variablity of sound propagation in air far exceed ths loss of precision of the method... In fact, I would assert that CES will frequently be more accurate because it doesn't try to be too precise, and is therefore less sensitive to these sorts of errors.

### Some terminology ###
Locbatch is based on code which was developed for a system where each node was actually a sub-array 
with multiple microphones, so some of the terminology may be a bit different than expected.

* **Node**: A node is a microphone.  Each node has an id (node_id) which is a *string* identifying it.  Be warned, `01` and `1` are different node_ids.  
* **Channel**: Channel is **not** a microphone. Channel should always be 1 here.
* **Event**: A sound which you want to localize.  At the most basic level, it is a start-time, duration (or end-time), low frequency, and high frequency box.  Within the code, it also includes lots of other information which might vary from one sound to the next, such as the temperature (used to calculate the speed of sound).

## Installing ##

Basic instlall instructions.  Power users can do a git clone and figure out the rest.

* Install python with pylab (ect.)  
I suggest just using the [Anaconda Python distribution](https://store.continuum.io/cshop/anaconda/) which includes everything you need.
* [Download locbatch.zip](https://github.com/travc/locbatch/archive/master.zip)
* Unzip it someplace
* (If using MS Windows) Make a shortcut to `locbatch.bat` and put that shortcut someplace easily accessible like the Desktop.
The program can be run by dragging a config file onto that shortcut (details below).

## Usage ##

### Input files ###

The input files locbatch needs are:
* [Sound files](#sound-files) (a collection of recordings and (optionally) a Syrinx-like metafile grouping consecutive (segmented) soundfiles)
* [Microphone location file](#microphone-location-file)
* [Annotation file](#annotation-file)
* [Configuration file](#configuration-file)

The input files should (probably) all be in the same directory.

#### Sound files ####
The actual sound files / recordings must be single channel files.
Wave (aka .wav) PCM files are supported natively and tested.  
Support for many different formats could be easily added by using [sox](http://sox.sourceforge.net/) (the code is already all there, but not enabled or tested).  
See the examples.

The sound files may be specified with a metafile OR a filename_pattern.

A metafile is the same as used by Syrinx, and lists **consecutive (no gaps)** recordings grouped by node.  The format (ignoring any blank lines, *italic* values vary) is something like:

<pre>
annotationfile=<i>ingored</i>
<i>node_01_file1.wav</i>
<i>node_01_file2.wav</i>
filegroup:
<i>node_02_file1.wav</i>
<i>node_02_file2.wav</i>
filegroup:
<i>node_03_file1.wav</i>
<i>node_03_file2.wav</i>
</pre>
Since Syrinx metafiles don't list node_ids, they are inferred from the order of the groups... `01`, `02` and `03` for this example.

Alternatively, a filename_pattern may be given if the sound files are named in a sensible way.  Specifically, the files must include the node_id and they must sort in the correct (by time) order.  
Again, the files matching the pattern (per node_id) must be consecutive (no gaps).

For example, the filename_pattern `20120410_*_{node_id}.wav` would work if your recordings are named like:

    20120410_063059_01.wav
    20120410_063059_02.wav
    20120410_063059_03.wav
    20120410_064400_01.wav
    20120410_064400_02.wav
    20120410_064400_03.wav
    20120410_065916_01.wav
    20120410_065916_02.wav
    20120410_065916_03.wav

#### Microphone location file ####
A file giving the position of each node/microphone is also required.  
The preferred format is one line per node listing: `node_id channel X Y Z` where `X`, `Y` and `Z` are position values in meters. NOTE: `channel` will always be `1` for wired arrays.  
Locbatch will also accept a file which just has one `X Y` pair per line,
from which it will infer the `node_id` from the line number (`01`, `02`, ect.) and set all the `Z` values to 0.  
See the examples.

#### Annotation file ####
The annotation file lists the events you want to localize.  
Locbatch currently supports two types of annotation files: Syrinx annotation files and Raven selection tables.  
I do intend to add support for other types of annotation files (at least Praat TextGrids and a locbatch/simple format).  
See the examples.

#### Configuration file ####
The confguration file provides locbatch all the information it needs to run.  
This includes specifying the events_filename (annotations), mics_filename (node/mic locations), ect.  
It also gives information needed for localization which may not be included in the annotations file under the `[Event defaults]` section.  
Finally, it specifies various options for how the localization is to be computed.

Two fiarly well commented example configuration files are included and should be used as a guide (one for Syrinx annotations, one for Raven):  
[examples/foo_syrinx.cfg](examples/foo_syrinx.cfg)  
[examples/foo_raven.cfg](examples/foo_raven.cfg)  
The [locbatch_code/default.cfg](locbatch_code/default.cfg) file gives all the possible options you can set.

Locbatch is intended to be run by creating a new config file (normally by copying, renaming an editing and existing one) for each new batch of localizations.  This way, the documentation of exactly what settings were used will be kept along with the results.  Any options not set in the config file will be set from the [locbatch_code/default.cfg](locbatch_code/default.cfg) file.

Results will automatically be output to a directory based on the config filename, so try to give it a meaningful name and be careful not overwrite existing files you want to keep.

Note: The config file is read by a python [SafeConfigParser](https://docs.python.org/2/library/configparser.html)
which means settings (known as options) are organized within sections.  It can also do some fancy stuff like interpolation where `%(var)s` will be replaced with the value of `var`, assuming `var` is set somewhere in the config file.

### Running ###

Under MS Windows, assuming you have installed the [Anaconda Python distribution](https://store.continuum.io/cshop/anaconda/), locbatch can be run by dragging a config file onto `locbatch.bat` (or a shortcut to `locbatch.bat`), assuming that the config file is in the same direcotry as the other files it needs.

Locbatch can also be run from the command line, though the exact command may depend on your python distribution.  It will be something like:  
`python something/locbatch/locbatch_code/locbatch.py -c your_config_file.loccfg`

## Interpreting results ##

[An example result figure](https://github.com/travc/locbatch/wiki)

### Things to watch out for ###
* **The error indicator/estimate value is very rough and not reliable.**  
Good vs bad localizations should normally be obvious from looking at the figure, but it is suprisingly tricky to quantify that into a metric.
* Spatial aliasing caused by some low frequency (< 200Hz) sounds, error in microphone positions, and/or error in the speed of sound can result in a 'bumpy' likelihood map.  
Try applying additional smoothing to the cross-correlation-evelopes (cenvs) by adding (in the `[Event defaults]` section of the config or per-event in the annotation file):  
`smooth_cenv_window: 75`  
This smooth the cenvs by convolving them with at 75 millisecond Hamming window.  Increase/decrease the `75` to appply more/less smoothing.
* ...


## Citation ##

For now, please cite my dissertation and or the JASA paper linked above.  
Hopefully I will get around to writing a little paper announcing this software officially.

## Acknowlegements ##

The example data was kindly provided by Gail L. Patricelli and Alan H. Krakauer at UC Davis ([lab website](http://www.eve.ucdavis.edu/gpatricelli/index.html)).


 
