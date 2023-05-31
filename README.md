
geftools : Tools for manipulating GEFs
==============

For a full documentation, see [geftools GitHub page](https://bgiresearch.github.io/geftools/).


## Installation
To install stereopy from source, you need:
- HDF5 1.12.1 or newer with development headers
- OpenCV 4.5.4 or newer with development headers
- A C compiler

You can specify build options for geftools as environment variables when you build it from source:
```shell
git clone https://github.com/STOmics/geftools.git
cd geftools
HDF5_ROOT=/path/to/hdf5 OpenCV_DIR=/path/to/opencv cmake .
make
make install
```

The supported build options are:
- HDF5_ROOT: To specify where to find HDF5.
- OpenCV_DIR: To specify where to find OpenCV.
- CMAKE_INSTALL_PREFIX: Install directory used by make install.
- GEFTOOLS_BUILD_DOC: Option to build documentation.


## LIST OF COMMANDS
```text
Command: bgef          Generate common bin GEF(.bgef) according to gem file or bin1 GEF
         cgef          Generate cell bin GEF(.cgef) according to common bin GEF and mask file
         view          View GEF
```


## COMMANDS AND OPTIONS
### geftools bgef [OPTION...]
```text
  -i, --input-file FILE   input gene expression matrix file(.gem/.gem.gz) or bin1 bGEF file [request]
  -o, --output-file FILE  output bin GEF file (.bgef) [request]
  -b, --bin-size STR      Set bin size by the comma-separated list [request] (default: 1,10,20,50,100,200,500)
  -r, --region STR        Restrict to a rectangular region. The region is represented by the comma-separated list of
                          two vertex coordinates (minX,maxX,minY,maxY) (default: "")
  -t, --threads INT       number of threads (default: 8)
  -s, --stat    BOOL      create stat group (default: true)
  -O, --omics   STR       input omics[request]
  -v, --verbose           Verbose output
      --help              Print help
```

Create bgef example:
```text
  ./geftools bgef -i xxx.bgef/xxx.bgem -b 1,10,20 -o xx.cgef -t 5 -O Transcriptomics
```


### geftools cgef [OPTION...]
Generate cell bin GEF (.cgef) according to common bin GEF (.bgef) file and mask file
```text
  -i, --input-file FILE   input bin GEF file [request]
  -m, --mask-file FILE    input mask file [request]
  -o, --output-file FILE  output cell bin GEF file (.cgef) [request]
  -b, --block FILE        Pre block size (default: 256,256)
  -t, --threads INT       number of threads
  -v, --verbose           Verbose output
      --help              Print help
```

Create cgef example:
```text
  ./geftools cgef -i xxx.bgef -m xxx.tif -o xx.cgef -t 5
```


### geftools view [OPTION...]
Show the contents of cell bin GEF
```text
  -i, --input-file FILE     Input bGEF/cGEF file [request]
  -o, --output-gem FILE     Output gem file (default: stdout)
  -d, --exp_data FILE       Input bGEF file to get exp_data.
  -b, --bin-size INT        Set bin size for bgef file, just support bGEF.
  -s, --serial-number STR   Input Serial number [request]
  -e, --exon INT            whether or not output exon (default: 1)
  --help                Print help
```

Create bgem example:
```text
  ./geftools view -i xxx.bgef -o xxx.gem -s SS200000135TL_D1 
```

Create cgem example:
```text
  ./geftools view -i xxx.cgef -d xxx.bgef -o xxx.gem -s SS200000135TL_D1 
```


## Coding Style Guide
See [here](docs/coding_style_guide.md).


