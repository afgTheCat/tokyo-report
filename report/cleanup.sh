#!/bin/sh
ls | grep "^.*\.\(aux\|log\|out\|bbl\|blg\|fdb_latexmk\|fls\|gz\|dvi\)" | xargs -I '{}' rm '{}'