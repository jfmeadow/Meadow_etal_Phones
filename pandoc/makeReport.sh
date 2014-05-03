#!/bin/bash

export mainFile=phones
export outFile=Meadow_etal_Phones_Supplemental_R_Analysis
# export bibFile=surfaceDemo
export fromPath=../


# copy figures from Rmd folder
cp -r $fromPath/figure .

# produces default latex pdf with no citations
# pandoc $fromPath$mainFile.md -o $mainFile.pdf --highlight-style=tango

# produces formated latex with no citations.
pandoc $fromPath$mainFile.md -o $outFile.pdf -H margins.sty --highlight-style=tango