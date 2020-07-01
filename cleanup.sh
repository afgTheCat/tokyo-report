#!/bin/sh
ls | grep "^.*\.\(aux\|log\|out\|bbl\|blg\|fdb_latexmk\|fls\|gz\)" | xargs -I '{}' rm '{}'